# LRC(14): THM-4009 Fourier prefix-order and owner audit

**Final status (2026-08-24): FINITE-EXACT bounded laboratory + PROVED loss
ledger.  LRC(14) remains OPEN.**  No theorem ID is reserved, no global Fourier
closure is claimed, and no weak boundary row is closed.

The inherited source is
[THM-4009](../01-canon/theorems/THM-4009-euclidean-covering-transference-short-relation-compression.md):
a hypothetical counterexample has a Graver relation of Euclidean square norm
at most `195`, `l1` norm at most `50`, and height at most `13`.  Its
support-two branch has 47 ratios, 3,666 unoriented ratio/support packets, and
7,332 oriented labelled assignments modulo global sign.  The target tested
here is a finite Hermitian danger matrix obtained from ordered prefixes of
cap-respecting scalar relations.

On the strict control `V_38={1,...,11,13,38}`, the Euclidean-shortest rows are
exactly thirty additive triples `x+y=z`, of square norm three.  For
`1<=k<=16`, every one of the sixty single-order prefix-bouquet matrices is
rigorously positive definite; every Arb `LDL^T` pivot is greater than `2/5`.
Combining both prefix orderings produces exactly two negative matrices.  For
`2+6=8`, a frozen integer polynomial of coefficient square norm `394` gives

```text
integral (N_V-1)|P|^2 = -5.9113725705681529025... < 0.
```

Thus shortness alone does not yield the certificate; cross-order phase is the
load-bearing new coordinate.  Centered Fejer orders one through seventeen
remain positive on the same row.

The equality hostiles AP and Goddyn--Wong remain positive definite on all
thirty two-order bouquets and have open danger cover almost everywhere.
Their valid `t=1/14` witness is atomic, with opposite-slope owners `(1,+)` and
`(13,-)`.  THM-4002's two rows with identical complete norm-three relation
fibre and opposite covariance also remain untouched by the new cap.

```text
source:      THM-4009 short Graver row
target:      finite Toeplitz/Hermitian danger matrix
map:         ordered closed-walk prefixes -> frequency bouquet
preserved:   sign, harmonic scale, prefix order, cross-order phase
destroyed:   endpoint owner, component incidence, arrival, parity
sidecar:     signed endpoint word / owner-address cross-correlation
```

The cheapest decisive next test is the exact SDP atlas on the 47 support-two
ratios crossed with the 17 THM-3910 types, retaining all prefix orderings and
the owner word, with AP, Goddyn--Wong, and both THM-4002 phase twins as hostile
controls.  Numerical SDP may discover vectors; only frozen integer vectors
with Arb sign certificates should carry claims.

Reproduction:

```bash
python3 -B 04-computation/lrc14_graver_bouquet_fejer_owner_audit_20260824.py
python3 -B -O 04-computation/lrc14_graver_bouquet_fejer_owner_audit_20260824.py
```

Both commands byte-match
`05-knowledge/results/lrc14_graver_bouquet_fejer_owner_audit_20260824.out`.
The detailed literature/guardrail audit is in
`05-knowledge/results/lrc14_thm4009_fourier_lp_bottleneck_audit_20260824.md`.
