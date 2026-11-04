# Equal-nine independent audit

Frozen publication allowlist: `REVIEW.md`, `README.md`, `audit.py`, `normal.json`, `optimized.json`, `SHA256SUMS`.

Run `python audit.py > normal.json`, `python -O audit.py > optimized.json`, then compare the outputs. Only the Python standard library is used. The script reads the named sibling source packets pinned in its receipts and replays `equal_nine_ambient/verify.py`; it does not import producer code.

PASS: every one of 2,373 retained coefficient values matches an independent physical-word/simplex-vertex reconstruction, including all 1,480 new order-eight/nine coefficients. The mathematical review separately checks the profile domain, strict scalar sign and actual physical arbitrary-ambient lift. The larger order-ten simplex's positive coefficient is only a method boundary.
