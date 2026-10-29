# Minimal publication manifest — 2026-10-29

Publish these **14 essential files**, plus this manifest and its JSON companion:

```
README.md
halving_certificate.json
binary_kernel_certificate.json
joint_candidate_3567.json
joint_candidate_3568.json
check_halving.py
check_kernel.py
check_joint_halves.py
phase2/README.md
phase2/extension_3567.json
phase2/extension_3568.json
phase2/inverse_transport_3567.json
phase2/inverse_transport_3568.json
phase2/verify_extension.py
```

The JSON manifest records exact sizes, SHA256 hashes, isolated test output
hashes, and output values. If the parent adds the planned dated phase2 pointer
to the first README before publication, refresh **that file's** size/hash in
this manifest; no mathematics or checker depends on the README hash.

## Runtime and isolation

Only `check_halving.py` requires **SymPy 1.14.0**. The other three checkers
use the Python standard library. All paths are resolved relative to `__file__`.
All four checkers were run normally and with `python -O` from `/tmp`, against
an isolated copy containing exactly these 14 files. All four output pairs
matched byte for byte.

No adjacent raw or compact rank database is required: the exact two source
curve models, original point lists, and binary character banks are embedded
in `halving_certificate.json`. In isolation `check_halving.py` reports 66,625
checks rather than 66,628, because the three optional adjacent source-file
SHA256 comparisons are skipped. All mathematical checks and hostile tests
still run. The original source SHA256 digests remain recorded in the JSON.

Point indices throughout are **zero-based** in the unchanged embedded source
order. Phase2 preserves all exact point transport and each elementary
index-two enlargement. The certificate supplies both forward rational and
inverse integral basis matrices.

## Exclude

The following are regenerators or transcripts, not checker dependencies:

- `produce_halving.py`;
- `inspect_kernel.py`;
- `factorization.log`;
- `phase2/extend_saturation.py`;
- `phase2/progress.log`;
- all `__pycache__` files.

The earlier README's mentions of these files describe research regeneration;
the minimal publication packet intentionally omits them.

## Independent audit and claim boundary

A third agent independently audited all 26 exact half identities, every
rank increase 17 through 30, the infinity/two-torsion character special cases,
both determinants and integral inverses, and both rootless good primes. Its
normal/-O audit passed; the parent holds that separate audit file and can
publish it alongside this packet.

The proved endpoint is **2-primary saturation of the two specified rank-30
subgroups**, with exact index 8192 and quotient `(Z/2Z)^13`. Thus the thirteen
original binary-kernel directions are genuine global **divisibility classes**,
not point linear relations. No full Mordell–Weil group, rank upper bound,
or odd-prime saturation is claimed.
