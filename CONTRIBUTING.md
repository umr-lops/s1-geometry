# Contributing

1. Create a branch from `main` and open a pull request (never push to `main`).
2. Install the development environment:

   ```bash
   python -m pip install -e '.[parquet,test,docs]' pre-commit
   pre-commit install
   ```

3. Before pushing:

   ```bash
   pre-commit run --all-files                      # ruff + file hygiene
   python -m pytest --cov=s1_geometry              # tests
   sphinx-build -b html -W docs docs/_build/html   # docs build without warnings
   ```

4. Embedded model coefficients (`heading.py`, `angles.py`, `wv.py`,
   `orbit.py`) are validated deliverables: do not change them without
   refitting **and** revalidating, and document the new validation in the
   README and `docs/certification.rst`.
5. Add an entry to `CHANGELOG.md` (section *Unreleased*).

Releases: create a GitHub release with a `YYYY.MM.DD` tag; the version is
taken from the tag and the `Publish` workflow uploads to PyPI.
