# Frozen phase7 independent release audit

**Status: PASS. Frozen 2026-10-29.** No producer modules are imported by the
four audit scripts below. All use only the Python standard library; normal
and optimized (`python -O`) outputs agree. Earlier frozen audit files are
unchanged. No raw database or research-paper PDF is required.

## Findings

### Isogeny separation and the unique finite prime bank

All66 compact models have independently recomputed good-reduction flags and
point-count traces. Every one of their2145 unordered pairs has an unequal
trace at a common good prime<=71, proving pairwise non-Q-isogeny. This does
not infer geometric non-isogeny, family independence, or a rank-density law.

The bank{19,29,31,37,43,47,53,59,61,67,71} covers every pair. Each of its eleven
primes has a pair for which it is the only witness anywhere in the recorded
window. These forced-pair certificates prove both minimality and uniqueness
within that window. No outside-window optimum is asserted.

### Individual versus joint divisibility

On3567 and3568, every individual source column has a nonzero local binary
character, excluding individual rational halves. The independently rebuilt
rank17 banks have13-dimensional kernels. Exhausting8192 words per kernel
recovers minimum supports6/four instances and5/one instance respectively.
Thus every smaller nonzero odd-coefficient support is excluded. The first
selected joint support on each curve is independently verified to have the
stored rational half, with sign retained. No inference from local vanishing
alone is used, and this auditor did not rerun the unnecessary characteristic-
zero quartic factorizations.

### Exact2-primary saturation

All26 replacement identities and their source-coordinate transports pass.
Every terminal bank has binary rank30; the derivative characters at reduced
two-torsion and neutral characters at infinity are handled explicitly.
Rootless good primes23 and3 exclude rational2-torsion. Hence each terminal
rank30 subgroup is2-saturated in E(Q), without assuming ambient rank30.

The transported matrix has determinant2^-13, entries with denominators at
most2, and an independently verified integral inverse. Explicit recovery
of each old point proves containment. Thus the enlarged subgroup has exact
index8192 and quotient C2^13, and is precisely the2-primary saturation of
the original subgroup. The initial thirteen-dimensional finite-bank kernel
therefore equals the actual global divisibility kernel for these inputs.
The final certificate also independently implies original point independence
through the invertible transport; priorF3 certificates are compatible but
not needed for this conclusion.

No full Mordell–Weil identification, exact ambient rank, or odd-prime
saturation is claimed.

### Tournament solutions and attachment memory

The proper-neighbourhood dependency-height induction is sound and gives
stabilization at every global iteration k>=h(T), with h(T)<=|T|-1. Independent
Tarjan/SCC and synchronous-iteration checks cover all1099 labelled tournaments
through order5, with all-vertex and maximum-outdegree seeds.

The five-vertex host is strong and rigid. Two globally equivariant rules
with identical singleton child-type rows produce source-winner sets of
sizes3 and5. All120 relabellings transport correctly. Rigidity makes the
rule extensions unambiguous, and child embeddings are the missing sidecar.
This is not TEQ nonuniqueness under one rule or an A000568 runtime claim.

## Exact release files and compact dependencies

All paths below are relative to the Session97 packet root.

| Audit script | Required mathematical inputs |
|---|---|
| `audit/check_isogeny_independent.py` | `ranks/phase2/highrank_inputs.json`; `isogeny/certificate.json`; `isogeny/minimum_prime_bank.json` |
| `audit/check_joint_halving_independent.py` | `halving/halving_certificate.json`; `halving/binary_kernel_certificate.json`; `halving/joint_candidate_3567.json`; `halving/joint_candidate_3568.json` |
| `audit/check_two_saturation_independent.py` | `halving/halving_certificate.json`; `halving/phase2/extension_3567.json`; `halving/phase2/extension_3568.json`; `halving/phase2/inverse_transport_3567.json`; `halving/phase2/inverse_transport_3568.json` |
| `audit/check_tournament_solution_independent.py` | None; the finite hostile and census are self-contained. |

Receipts: `audit/isogeny-normal.json`, `audit/halving-normal.json`,
`audit/saturation-normal.json`, `audit/solution-normal.json`.
Replay hashes and exact compact-input hashes are retained in
`audit/PHASE7_REPLAY_MANIFEST.json`.

Further proof detail is in `audit/ISOGENY_AND_TWO_SATURATION_AUDIT.md` and
`audit/TOURNAMENT_SOLUTION_AUDIT.md`. Exclude the exploratory
`audit/check_rank_data_independent.py` and the optional full-database integrity
checker from this release: both require a raw snapshot that is not distributed.

This completes the requested bounded audit. Remaining open infinite-failure,
generic-section and global-rank questions are not promoted by these checks.
