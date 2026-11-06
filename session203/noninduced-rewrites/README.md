# A physical 775 repair and a restricted source diagnostic

**Independently audited in the stated scopes. The eventual odd tournament Hamilton-path spectrum remains OPEN.**

- [Physical repair](FINITE_REPAIR.md): reverse one old arc of E(129) and add two vertices to obtain exactly 775 Hamilton paths. The same named recipe has a proved count at every shifted dyadic input and hits its intended target only at 129. This is not an all-input extractor or a minimum-edit theorem.
- [Source-ray diagnostic](SOURCE_RAY_BOUND.md): a finite exact test excludes one-state affine O/Z/J source digit grammars at bases 3,5,7,9. Several states, nonlinear encodings, and arbitrary graph rewrites are outside its scope.
- [Independent audit](audit/NONINDUCED_REPAIR_AUDIT.md): literal path enumeration, exact gap counting, and separate rational linear algebra.

Frozen dependencies and their repository commits are in dependencies/PROVENANCE.json. The original proof bytes are preserved in reviewed-originals/, and PROMOTION.json records the two status-only promotions. Exploratory searches are omitted.

## Reproduction

Python 3 standard-library independent checks:

    python3 audit/audit_repair775.py > /tmp/repair775.json
    cmp /tmp/repair775.json audit/REPAIR775_NORMAL.json
    python3 audit/audit_source_rays.py > /tmp/source_rays.json
    cmp /tmp/source_rays.json audit/SOURCE_RAYS_NORMAL.json

Author checks (the third uses SymPy):

    python3 verify775.py > /tmp/author775.json
    cmp /tmp/author775.json VERIFY775.json
    python3 check_repair_polynomial.py > /tmp/clone775.json
    cmp /tmp/clone775.json REPAIR_POLYNOMIAL.json
    python3 verify_source_rays.py > /tmp/author_rays.json
    cmp /tmp/author_rays.json SOURCE_RAYS_VERIFIED.json

All checks also pass under python3 -O. No external research files are read by these scripts.
