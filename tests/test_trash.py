#!/usr/bin/env python3
"""
Test framework for TRASH.R - compares output file hashes between runs.
Uses pytest for testing infrastructure.
Tests are parameterized to run with multiple FASTA files (small and med).

Usage:
    pytest test_trash.py -v                    # Run all tests
    pytest test_trash.py -v -k small           # Run only small tests
    pytest test_trash.py -v -k med             # Run only med tests
    pytest test_trash.py -v --update-baseline  # Update all baseline hashes
"""

import subprocess
import hashlib
import json
import shutil
from pathlib import Path
import pytest

# Configuration
TEST_FASTAS = [
    {
        "name": "small",
        "fasta": "testing_fastas/ath_Chr1_extraction_trc.fasta",
        "baseline": "tests/test_baseline_small.json",
        "output_pattern": "ath_Chr1_extraction_trc.fasta_*",
        "exclude": ["ath_Chr1_extraction_trc.fasta_run_time.csv"]
    },
    {
        "name": "med",
        "fasta": "testing_fastas/Ath_Chr1_med.fasta",
        "baseline": "tests/test_baseline_med.json",
        "output_pattern": "Ath_Chr1_med.fasta_*",
        "exclude": ["Ath_Chr1_med.fasta_run_time.csv"]
    }
]
TEST_OUTPUT_DIR = "test_output"  # Dedicated directory for test outputs


def compute_file_hash(filepath):
    """Compute MD5 hash of a file."""
    md5_hash = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def get_output_files(src_dir, output_pattern):
    """Get all output files matching the pattern."""
    return sorted(src_dir.glob(output_pattern))


def compute_all_hashes(src_dir, output_pattern, exclude_files):
    """Compute hashes for all output files (excluding non-deterministic files)."""
    files = get_output_files(src_dir, output_pattern)
    hashes = {}
    for filepath in files:
        relative_path = filepath.name
        # Skip files that are expected to be non-deterministic
        if relative_path in exclude_files:
            continue
        hashes[relative_path] = compute_file_hash(filepath)
    return hashes


def clean_output_files(src_dir, output_pattern):
    """Remove previous output files."""
    for filepath in get_output_files(src_dir, output_pattern):
        filepath.unlink()


def run_trash(test_output_dir, project_root, fasta_path):
    """Run TRASH.R in the test output directory and return success status."""
    # Build command with absolute paths
    trash_script = project_root / "src" / "TRASH.R"
    fasta_file = project_root / fasta_path
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


@pytest.mark.parametrize("test_config", TEST_FASTAS, ids=[cfg["name"] for cfg in TEST_FASTAS])
def test_trash_output_deterministic(test_config, test_output_dir, project_root, request):
    """Test that TRASH.R produces deterministic output by comparing file hashes."""

    update_baseline = request.config.getoption("--update-baseline")

    # Run TRASH.R in the isolated test directory
    print(f"\nRunning TRASH.R for {test_config['name']} in {test_output_dir}")
    run_trash(test_output_dir, project_root, test_config['fasta'])

    # Compute hashes from the test output directory
    current_hashes = compute_all_hashes(
        test_output_dir,
        test_config['output_pattern'],
        test_config['exclude']
    )

    baseline_path = project_root / test_config['baseline']

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

    # Clean up output files for this test
    clean_output_files(test_output_dir, test_config['output_pattern'])


if __name__ == "__main__":
    # Allow running directly
    import sys
    sys.exit(pytest.main([__file__, "-v"] + sys.argv[1:]))
