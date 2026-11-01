# Short fresh-block radix banks

Read `PROOF.md` for exact scope and proofs.

- **FINITE-EXACT:** among words of at most four operations in the authenticated 32-matrix fresh-transitive-block bank, no fixed output-normalized affine ray supports two distinct nonzero digits at any common integer radix >=2. Every integer eigenvalue is checked; signed affine offsets are allowed.
- **PROVED:** actual three-operation migration yields `F_k=(52*80^k+27)/79`, on `7k+2` vertices, starting with `1,53,4213,...`.
- The global spectrum and longer-word / changing-profile digit constructions remain **OPEN**.

Reproduce with Python 3, standard library only:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 verify.py
PYTHONDONTWRITEBYTECODE=1 python3 -O verify.py
```

`bank.json` is regenerated exactly from `migration.py`, a verbatim copy of the sibling constructor packet's audited compiler. That packet's proof is the inherited all-source operation theorem. `normal.json` and `optimized.json` are identical receipts. The complete frozen bank, spectral candidates, rays and actual first-iterate graph witness are retained.

Only files listed in `PUBLIC_FILES.txt` are publication inputs. Exploratory scripts are not dependencies.
