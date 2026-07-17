---
id: THM-948
title: THE EXACT PACKET SUPPORT DECOMPOSITION AND SIGN AUDIT — exact rational Möbius masses and B5 identity on three packets; packet-level near-sharpness measured; universal positive association refuted by an infinite support-three family
status: The subset-Möbius decomposition and the direct B5 identity are exact rational equalities on all three audited packets. Their aggregate absolute-vs-signed slack is 1.000–1.001 on those packets only. The proposed universal law M(A) >= 0 for |A| >= 3 is refuted at supports 3, 4, and 5; support 3 has a complete infinite modulo-28 counterfamily.
source: kind-pasteur-2026-07-17-S128 cont.39; corrected by codex-S50/S52
depends_on:
  - THM-935 (the E_s identity this makes exactly computable)
  - THM-930 (depth-spectrum sweep, the direct side of the referee)
related:
  - THM-946 (corrected two-pole tail program; the first-pushed use of that number)
  - THM-944/945/947 (B5 scoreboard, moment wall, and arc wire)
  - THM-938/939/942/943 (the deviation-fee ledger)
script:
  - 04-computation/exact_support_decomposition_kps_S128c39.py -> 05-knowledge/results/exact_support_decomposition_kps_S128c39.out
  - 04-computation/exact_support_sign_counterexample_codex_S50.py -> 05-knowledge/results/exact_support_sign_counterexample_codex_S50.out
  - 04-computation/exact_support_triple_sawtooth_codex_S51.py -> 05-knowledge/results/exact_support_triple_sawtooth_codex_S51.out
---

# THM-948 — exact packet support decomposition and sign audit

## I. What is exact

For each support `A`, `2 <= |A| <= 5`, define

```text
M(A) = mu(intersection_(v in A) B_v) - p^|A|
       - sum_(S proper subset A, |S| >= 2) M(S) p^(|A|-|S|),
p = 1/7.
```

Every intersection measure is computed by an exact sweep on the common endpoint
grid `14*lcm(A)`. Summing `M(A)` over supports of each size gives the following
exact referee against the direct depth-spectrum computation of `B5`:

| packet | B5 direct | support identity | M2 | M3 | M4 | M5 |
|---|---|---|---|---|---|---|
| certified | +0.082053 | exact | +0.000311 | +0.006682 | +0.023545 | +0.030061 |
| opus-30Z | -0.601603 | exact | +0.000098 | +0.231821 | +0.365015 | +0.505866 |
| random | +0.088381 | exact | -0.000060 | +0.002007 | +0.019141 | +0.027255 |

Thus concrete packets can be autopsied support by support without estimating a
relation tail. This does not provide a universal tail theorem.

## II. What the three packets show—and do not show

The signed total versus its absolute envelope is respectively

```text
0.040039 vs 0.040083,   0.72370 vs 0.72371,   0.033711 vs 0.033711.
```

So the absolute support-size budget is nearly sharp on these three packets. That is
an observed packet property, not a sharpness law. In particular, these examples do
not imply that every individual support mass, or every aggregate support layer, has
a fixed sign.

## III. Exact refutation of positive association

The conjecture `M(A) >= 0` for every `|A| >= 3` is false. Exact witnesses occur at
each support size used by `B5`:

```text
M({1,2,15})       = -4/1715,
M({1,2,3,28})     = -5/4116,
M({1,2,3,4,32})   = -109/806736.
```

The support-three failure is an infinite quasiperiodic law. For every `N > 2`,

```text
M({1,2,N}) = k[N mod 28] / (686 N),
```

where

```text
k = (0,25,50,61,44,55,38,49,32,43,26,37,20,24,
     0,-24,-20,-37,-26,-43,-32,-49,-38,-55,-44,-61,-50,-25).
```

Hence the mass is positive for residues `1..13`, zero for `0,14`, and negative
for `15..27`. In particular,

```text
M({1,2,15+28m}) = -12 / (343(15+28m))   for every m >= 0.
```

The proof rescales the fixed `B1` and `B1 intersection B2` windows into integer-
and half-integer-centered tooth discrepancy functions. An endpoint sweep through
`N=2000` and an independent midpoint-cell engine through `N=300` referee the exact
closed formula.

## IV. Structural reading

The sign-bearing datum first appears on supports of size three. Pair mass is
symmetric, so orienting it by runner labels creates an arbitrary tournament and
destroys the Möbius sign. The faithful quotient is a signed support hypergraph
decorated, in the `{1,2,N}` family, by the cyclic 28-residue boundary walk. Its
antipodal involution `r -> 28-r` reverses the mass sign.

For the opus-30Z packet the exact autopsy remains useful: the leading fees come
from the Schur triple `{420,450,870}` and the AP clusters `{450,510,570}` and
`{450,570,690}`, not merely from the common gcd-30 dilate structure. The next
universal argument must retain signed support data; a blanket FKG/positive-
association reduction cannot be used.
