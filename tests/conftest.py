"""Pytest configuration: ensure the repo root is importable.

This file makes the tests work without requiring a project-level
`pytest.ini`, `pyproject.toml`, or installed package. Pytest auto-loads
`conftest.py` before collection, so the path manipulation here runs
before any test imports `metrics.*` or `lib.*`.

If you later add a `pyproject.toml` with a `[tool.pytest.ini_options]`
section (or a `pytest.ini`) and a proper package install, this file can
be reduced to just `import sys` or removed entirely.
"""
import sys
from pathlib import Path

# tests/conftest.py → repo root is one level up.
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
