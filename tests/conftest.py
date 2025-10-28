"""pytest configuration and hooks."""

def pytest_addoption(parser):
    """Add custom pytest command line options."""
    parser.addoption(
        "--update-baseline",
        action="store_true",
        default=False,
        help="Update the baseline hashes instead of comparing"
    )
