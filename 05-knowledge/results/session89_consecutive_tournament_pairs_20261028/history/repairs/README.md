# Minimal legacy square-codec repairs

The patches target exact files at `03408f84ac43fcdcfb0887e5d6bc1fc8792415cf`; both originals also match `b38ec1a1eaf7f4e31a7a38308e2764b5c7650d6d`.

- `tiling_square_s325.patch`: two executable coordinate changes, `hi-2` to `hi-1`, in the illustrative upper-grid loop and `build_square`.
- `new_sequences_s327.patch`: the same one-expression correction in the inline square constructor. No production function is introduced and no rank/determinant experiment is rerun.

Both add explicit provenance warnings: saved results from the old lossy maps were **not full-binary-matrix rank/determinant censuses**. Historical output is not rewritten. Other calculations, OEIS claims, and exploratory narrative are neither changed nor certified. The full corrected copies are staging artifacts; publication can apply the patches to the existing live files instead of duplicating them.

## Public regression checker

`audit_patched_codec.py` defaults to live repository-relative `04-computation/tiling_square_s325.py` and `04-computation/new_sequences_s327.py`, locating the repository from its own ancestors or the working directory. It never reads private snapshots. From a checkout:

```sh
python path/to/audit_patched_codec.py --repo-root .
python -O path/to/audit_patched_codec.py --repo-root .
```

For staged copies, pass both `--s325-source PATH` and `--s327-source PATH`; either source can also be overridden individually when the other is resolvable from `--repo-root`.

The checker AST-extracts S325's actual `build_square`; S327 has no named codec function, so it extracts the actual tile-list setup and inline matrix-filling block into a temporary callable. Structural guards fail closed if that source layout changes. It does **not** execute either script's top-level experiments, rank/determinant evaluations, or output logic.

For each codec, it checks all 2, 16, 512, and 65,536 input pairs for matrix sides 1–4 against an independent inverse-cell oracle, checks complete bijection, explicitly separates the old `000`/`101` collision, and checks zero, every coordinate basis vector, and all ones through side 12. Normal and optimized runs agree at **399,112 always-active checks** (including extraction guards). Finite tests supplement the all-size coordinate proof in `../REPAIRED_SQUARE_CODEC.md`; they do not prove arbitrary-size behavior by enumeration alone.

`regression.json` additionally records byte-for-byte patch replay, rejection of each original buggy codec, and successful repository-relative/default and explicit-source invocation. The old invalid-backslash escape in S325's unchanged docstring can emit a harmless `SyntaxWarning` during parsing; no unrelated cleanup is included.
