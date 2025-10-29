#!/usr/bin/env python3
"""
Test framework for TRASH.R - compares output file hashes between runs.
Uses pytest for testing infrastructure.

Usage:
    pytest test_trash.py -v                    # Run tests
    pytest test_trash.py -v --update-baseline  # Update baseline hashes
"""

import subprocess
import hashlib
import json
import shutil
from pathlib import Path
import pytest

# Configuration
TEST_FASTA = "testing_fastas/ath_Chr1_extraction_trc.fasta"
BASELINE_FILE = "tests/test_baseline.json"
OUTPUT_PATTERN = "ath_Chr1_extraction_trc.fasta_*"
TEST_OUTPUT_DIR = "test_output"  # Dedicated directory for test outputs
TRASH_CMD = ["Rscript", "src/TRASH.R", "-f", TEST_FASTA]
# Files to exclude from hash comparison (non-deterministic)
EXCLUDE_FILES = ["ath_Chr1_extraction_trc.fasta_run_time.csv"]


def compute_file_hash(filepath):
    """Compute MD5 hash of a file."""
    md5_hash = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def get_output_files(src_dir):
    """Get all output files matching the pattern."""
    return sorted(src_dir.glob(OUTPUT_PATTERN))


def compute_all_hashes(src_dir):
    """Compute hashes for all output files (excluding non-deterministic files)."""
    files = get_output_files(src_dir)
    hashes = {}
    for filepath in files:
        relative_path = filepath.name
        # Skip files that are expected to be non-deterministic
        if relative_path in EXCLUDE_FILES:
            continue
        hashes[relative_path] = compute_file_hash(filepath)
    return hashes


def clean_output_files(src_dir):
    """Remove previous output files."""
    for filepath in get_output_files(src_dir):
        filepath.unlink()


def run_trash(test_output_dir, project_root):
    """Run TRASH.R in the test output directory and return success status."""
    # Build command with absolute paths
    trash_script = project_root / "src" / "TRASH.R"
    fasta_file = project_root / TEST_FASTA
    cmd = ["Rscript", str(trash_script), "-f", str(fasta_file)]

    result = subprocess.run(
        cmd,
        cwd=test_output_dir,
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print(f"STDOUT: {result.stdout}")
        print(f"STDERR: {result.stderr}")
        raise RuntimeError(f"TRASH.R failed with exit code {result.returncode}")
    return result


@pytest.fixture(scope="session")
def project_root():
    """Get the project root directory."""
    return Path(__file__).parent.parent


@pytest.fixture(scope="session")
def baseline_path(project_root):
    """Get the baseline file path."""
    return project_root / BASELINE_FILE


@pytest.fixture(scope="session")
def test_output_dir(project_root):
    """Create and manage a dedicated test output directory."""
    output_dir = project_root / TEST_OUTPUT_DIR

    # Remove existing test output directory if it exists
    if output_dir.exists():
        shutil.rmtree(output_dir)

    # Create fresh test output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nCreated test output directory: {output_dir}")

    yield output_dir

    # Cleanup: remove test output directory after tests complete
    print(f"\nCleaning up test output directory: {output_dir}")
    if output_dir.exists():
        shutil.rmtree(output_dir)


@pytest.fixture(scope="session")
def current_hashes(test_output_dir, project_root):
    """Run TRASH and compute hashes of output files."""
    # Run TRASH.R in the isolated test directory
    print(f"\nRunning TRASH.R in {test_output_dir}")
    run_trash(test_output_dir, project_root)

    # Compute hashes from the test output directory
    hashes = compute_all_hashes(test_output_dir)

    return hashes


def test_trash_output_deterministic(current_hashes, baseline_path, request):
    """Test that TRASH.R produces deterministic output by comparing file hashes."""

    update_baseline = request.config.getoption("--update-baseline")

    if update_baseline:
        # Update baseline mode
        with open(baseline_path, 'w') as f:
            json.dump(current_hashes, f, indent=2, sort_keys=True)
        print(f"\n✓ Baseline updated: {baseline_path}")
        print(f"  Files tracked: {len(current_hashes)}")
        pytest.skip("Baseline updated successfully")

    # Normal test mode - compare against baseline
    if not baseline_path.exists():
        pytest.fail(
            f"No baseline found at {baseline_path}.\n"
            f"Run: pytest test_trash.py --update-baseline"
        )

    with open(baseline_path, 'r') as f:
        baseline_hashes = json.load(f)

    # Compare files
    all_files = set(current_hashes.keys()) | set(baseline_hashes.keys())

    differences = []
    for filename in sorted(all_files):
        current = current_hashes.get(filename)
        baseline = baseline_hashes.get(filename)

        if current is None:
            differences.append(f"  - Missing in current run: {filename}")
        elif baseline is None:
            differences.append(f"  + New in current run: {filename}")
        elif current != baseline:
            differences.append(
                f"  ✗ {filename}\n"
                f"    Expected: {baseline}\n"
                f"    Got:      {current}"
            )

    if differences:
        diff_msg = "\n".join(differences)
        pytest.fail(
            f"\nHash mismatch detected:\n{diff_msg}\n\n"
            f"To update baseline: pytest test_trash.py --update-baseline"
        )

    print(f"\n✓ All {len(current_hashes)} output files match baseline")


if __name__ == "__main__":
    # Allow running directly
    import sys
    sys.exit(pytest.main([__file__, "-v"] + sys.argv[1:]))
