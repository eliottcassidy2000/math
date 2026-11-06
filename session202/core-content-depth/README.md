# Common-content core depth and a complete target-67 obstruction

**Independently reviewed. No eventual-coverage threshold is proved.**

- `CONTENT.md`: if every allowed rooted core has a common divisor c>=2 of both its root-deletion Hamilton count and two-bin count, then a word with m>=2 cores has c^(m-2) dividing its Hamilton count. Arbitrary O/Z singleton spacers are allowed. Prime outputs require at most two cores.
- `TRIANGLE_67.md`: a complete all-depth inverse for triangle-root cores and arbitrary O/Z spacers excludes exactly {7,55,67} among coprime-to-three odd targets through67.
- `SCOPE.md`: the infinite bank of one-triangle-plus-singleton deletion cores generates no additional language, by exact SCC factorization.
- `audit/AUDIT.md`: independent PASS review, including an alternative 229-state saturated inverse and physical profile checks.

Run from this directory:

    sha256sum -c MANIFEST.sha256
    PYTHONDONTWRITEBYTECODE=1 python3 verify_content.py > /tmp/content-normal.json
    PYTHONDONTWRITEBYTECODE=1 python3 -O verify_content.py > /tmp/content-optimized.json
    cmp /tmp/content-normal.json /tmp/content-optimized.json
    cmp /tmp/content-normal.json CONTENT-normal.json
    PYTHONDONTWRITEBYTECODE=1 python3 verify_triangle_67.py > /tmp/triangle67-normal.json
    PYTHONDONTWRITEBYTECODE=1 python3 -O verify_triangle_67.py > /tmp/triangle67-optimized.json
    cmp /tmp/triangle67-normal.json /tmp/triangle67-optimized.json
    cmp /tmp/triangle67-normal.json TRIANGLE67-normal.json
    PYTHONDONTWRITEBYTECODE=1 python3 audit/independent.py > /tmp/content-independent.json
    cmp /tmp/content-independent.json audit/independent-normal.json

All scripts use the Python standard library. The imported `verify.py` is included and pinned in the current manifest. Frozen input proofs and historical manifests are preserved under `originals/`; the current release manifest is the one to replay after status promotion.
