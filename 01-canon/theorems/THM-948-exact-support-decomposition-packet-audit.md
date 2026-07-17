---
id: THM-948
title: THE EXACT PACKET SUPPORT DECOMPOSITION AND SIGN AUDIT — exact rational Möbius masses and B5 identity on three packets; universal positive association refuted; the complete two-parameter doubling-triple Fourier law
status: The subset-Möbius decomposition and the direct B5 identity are exact rational equalities on all three audited packets. Their aggregate absolute-vs-signed slack is 1.000–1.001 on those packets only. The proposed universal law M(A) >= 0 for |A| >= 3 is refuted at supports 3, 4, and 5. Support three has both the complete {1,2,N} modulo-28 family and the exact {a,2a,b} two-residue Fourier law, including its full zero locus and sign; the complete finite Bernoulli coefficient law and its arbitrary-residue lifts are Lean-checked, while the analytic Fourier/integral identification remains external.
source: kind-pasteur-2026-07-17-S128 cont.39; corrected and strengthened by codex-S50/S52/S59
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
  - 04-computation/exact_support_a2a_b_fourier_codex_S59.py -> 05-knowledge/results/exact_support_a2a_b_fourier_codex_S59.out
  - 04-computation/lean/TournamentH7/TournamentH7/LRCExactDoublingTriple.lean
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCExactDoublingTriple.lean
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

## IV. Complete doubling-triple Fourier law

The one-parameter sawtooth is a slice of a two-parameter theorem.  Let `a,b`
be positive with `{a,2a,b}` distinct, and put

```text
d=gcd(a,b),                 A=a/d, B=b/d.
```

There is an explicit integer table `C` on
`(Z/14Z) x (Z/28Z)` such that

```text
M({a,2a,b})=C(A mod 14,B mod 28)/(686AB).               (1)
```

Its complete zero and sign laws are

```text
C=0  <=>  7|A or 14|B,
sign C = centeredSign_14(A) * centeredSign_28(B)         (C != 0).
```

Thus the experimentally suggested implication `7|A => M=0` is true, but it
has the additional and exact `14|B` branch; together they are the full
converse.

The proof is Fourier-exact.  For

```text
g(x)=1_{||x||<1/14}-1/7,       H(x)=g(x)g(2x),
```

the support mass is the third centered moment

```text
M({A,2A,B})=integral H(At)g(Bt)dt
            =sum_(k!=0) Hhat(Bk) ghat(Ak),
```

because `(A,B)=1`.  The coefficients are

```text
ghat(n)=sin(pi n/7)/(pi n),
Hhat(n)=[sin(pi n/14)-(1/7)sin(pi n/7)
         -(2/7)1_{2|n}sin(pi n/14)]/(pi n).
```

The absolutely convergent sine-pair series reduces by the periodic Bernoulli
`B2` identity to (1).  Exact rational evaluation of all `14*28` residue cells
proves the zero/sign statement.  Independent referees compare the closed form
with 5,028 endpoint-sweep triples, every two- and three-face midpoint cell in
267 further triples, and 54 common dilations.

`TournamentH7.LRCExactDoublingTriple` formalizes the denominator-cleared
Bernoulli side.  A single finite audit proves all `392` cells integral, their
zero/sign laws and both residue reflections; subsequent theorems lift these to
arbitrary natural residues, prove `C(A,B)=0` exactly when `7|A or 14|B`, and
recover the infinite negative `15 mod 28` family at the closed-form interface.
It is sorry-free, avoids `native_decide`, and audits to the standard
foundational trio.  Its scope deliberately begins after the analytic identity
between the concrete support mass and the Bernoulli expression.

For B5 this is a genuine relation-tail simplification: every support-three
term containing a doubling pair has an exactly known sign or vanishes on the
two divisor hyperplanes.  It is not a universal positivity theorem; outside
the zero locus half of the centered residue quadrants are negative.

## V. Structural reading

The sign-bearing datum first appears on supports of size three. Pair mass is
symmetric, so orienting it by runner labels creates an arbitrary tournament and
destroys the Möbius sign. The faithful quotient is a signed support hypergraph.
In the `{1,2,N}` family it carries the cyclic 28-residue boundary walk; on the
full doubling locus it carries the Fourier cell `(A mod 14,B mod 28)` and the
intersection of the frequency lattices `A*Z` and `B*Z`.  Its antipodal residue
involutions reverse the mass sign, while the Fourier zero divisors give the two
exact vanishing branches.

For the opus-30Z packet the exact autopsy remains useful: the leading fees come
from the Schur triple `{420,450,870}` and the AP clusters `{450,510,570}` and
`{450,570,690}`, not merely from the common gcd-30 dilate structure. The next
universal argument must retain signed support data; a blanket FKG/positive-
association reduction cannot be used.
