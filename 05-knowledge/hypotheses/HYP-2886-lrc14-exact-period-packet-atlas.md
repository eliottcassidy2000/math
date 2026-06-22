---
id: HYP-2886
status: EVIDENCE / proof-target; exact finite scout, no LRC14 proof claimed
source: codex-2026-06-22-S102
tags: [lrc14, exact-period, euler-totient, primitive-packets, residue-atlas, chi7, affine-lane, crt, tournament-analysis]
related:
  - HYP-2865
  - HYP-2876
  - HYP-2882
  - HYP-2883
  - HYP-2884
  - HYP-2885
  - HYP-2632
  - HYP-2636
  - HYP-2628
  - HYP-2864
  - OPEN-Q-108
results:
  - 04-computation/lrc14_exact_period_packet_atlas_codex_s102.py
  - 05-knowledge/results/lrc14_exact_period_packet_atlas_codex_s102.out
---

# HYP-2886: LRC14 exact-period packets are the residue atlas; fixed bases are only charts

The useful finite-denominator statement is not:

```text
there is a universal finite witness basis B.
```

HYP-2865/HYP-2876 already refute that.  The useful statement is:

```text
LRC14 witnesses should be tracked as exact-period unit packets a/D
until after the mod-7, affine-lane, and CRT-defect data are evaluated.
```

This is the LRC-side continuation of HYP-2882, HYP-2883, and HYP-2884.  Tournament
strong components tell us not to contract a component unless its `H` label is
retained.  LRC denominator packets say the same thing arithmetically: do not
project a rational witness grid to a scalar count `N(S,D)` until the unit
packet's denominator, residue, `chi_7` class, affine-pair class, and CRT defect
have done their work.

After the concurrent HYP-2885 additive-energy extremality note, the split is
cleaner.  HYP-2885 is the cap-bound/realizability branch: prove no integer
speed set out-covers the interval by routing `L_y` through Fejer/additive-energy
extremality.  HYP-2886 is the witness-packet branch: if a row survives that
cap analysis locally, keep exact-period packets rather than fixed denominator
bases, and prove the remaining scaled packets cannot all be covered.  These are
not competing closures; they are the two faces of the same OPEN-Q-108 gap.

## Exact packet predicate

For a speed row `S` and denominator `D`, define an exact-period packet as a unit
`a mod D`.  It is safe at LRC14 level exactly when

```text
14 * min(s*a mod D, D - (s*a mod D)) >= D
```

for every `s in S`.

The packet count

```text
N(S,D) = #{a in (Z/DZ)^* : a is safe for S}
```

depends only on `S mod D`, while the packet capacity is `phi(D)`.  This is the
same exact-period law as HYP-2628:

```text
sum_{d|q} phi(d) = q.
```

## Exact proved fragments

The first obstruction is completely elementary.

```text
If some speed s in S is divisible by D, then N(S,D)=0.
```

Indeed `s*a = 0 mod D` for every unit `a`.  Therefore every fixed finite
denominator basis `B` is killed by the primitive covering row

```text
{1,...,11,13,84*lcm(B)}.
```

The S102 scout makes this visible in a scaled form.  For the row

```text
divload_B90 = {1,...,11,13,84*lcm(1,...,90)}
```

all counts on the prompt-style basis vanish:

```text
N(21)=N(41)=N(53)=N(83)=N(89)=0,
```

and the first unit witness in the scan is at

```text
D=97, a=8, N=2.
```

Thus finite bases are not global closures.  They are local residue charts.

Covering rows also kill every reduced denominator `2..14`: for each such
denominator, the covering condition supplies a speed divisible by it.  This is
the exact packet form of the apex-7 obstruction.

## S102 finite evidence

Script:

- `04-computation/lrc14_exact_period_packet_atlas_codex_s102.py`
- output: `05-knowledge/results/lrc14_exact_period_packet_atlas_codex_s102.out`

The first exact-period unit witnesses in the tested covering/tower rows are:

```text
AP13_boundary: D=14, a=1,  N=6
cover_84:      D=41, a=17, N=2
tower_m6:      D=53, a=22, N=2
tower_m53:     D=55, a=23, N=2
AP12_182:      D=27, a=2,  N=2
floor_star:    D=23, a=4,  N=2
divload_B60:   D=67, a=28, N=2
divload_B90:   D=97, a=8,  N=2
```

The prompt-like basis counts on `(21,41,53,83,89)` are:

```text
cover_84:    (0,2,2,2,2)
tower_m6:    (0,0,2,2,2)
tower_m53:   (0,0,0,2,2)
divload_B90: (0,0,0,0,0)
```

This is exactly the corrected HYP-2865 picture.  Small fixed denominators work
as a sampled atlas for ordinary rows; divisor-loaded rows force the atlas to
move.

## The mod-7 packet signal

Even when the scalar count is only `N=2`, the safe packets are not random in
the mod-7 quotient.  Examples:

```text
cover_84, D=41:
  safe mod7 residues = {3:2}
  chi_7 = NQR only
  affine pair = 3/6 only

tower_m53, D=55:
  safe mod7 residues = {2:1, 4:1}
  chi_7 = QR only
  affine pairs = 0/2 and 4/5

floor_star, D=23:
  safe mod7 residues = {4:1, 5:1}
  affine pair = 4/5 twice
```

S102 ranks exact-period quotient lenses by aggregate variance explained on
mixed packet-safety cases:

```text
mod14        0.167119
mod7         0.080466
chi_x_affine 0.064834
affine_pair  0.041515
chi7         0.025836
parity       0.000000
```

This is the correct reading.  `mod14` wins because the witness band is a
literal `2*7` band.  But the LRC14 proof object is not "mod14 wins"; it is
that `mod7`, `chi_7`, and the affine-pair quotient retain nontrivial signal
before scalarization.  That is exactly the layer used by HYP-2632, HYP-2883's
locally balanced signed packet graph, and HYP-2884's first lifted-divergence
probe.

## Scaled CRT decomposition

At the natural scaled modulus `q=14*Vmax`, witnesses decompose into many exact
periods.  S102 reports:

```text
cover_84, q=1176:
  good a = 12
  top exact denominators = 1176:4, 392:4, 147:2, 98:2

tower_m6, q=7056:
  good a = 80
  top exact denominators = 7056:22, 3528:12, 2352:8, ...

tower_m53, q=62328:
  good a = 706
  top exact denominators = 62328:212, 31164:110, 20776:104, ...
```

So the scaled CRT route is not a bounded-denominator route.  It is a
many-packet survival route: the forbidden classes generated by the runners
must fail to cover all exact-period packets.

This makes HYP-2864 and the CRT route compatible.  If a finite sheet or
denominator argument fails, it should expose a bounded quotient/gcd or a
divisor-loaded packet, not a mysterious analytic obstruction.

## Multiplicativity defect as a strong-atom analogue

Tournament strong components satisfy the exact product law

```text
H(T) = product_i H(C_i).
```

Exact-period packet capacity satisfies the arithmetic product law

```text
phi(D1*D2) = phi(D1)*phi(D2),  gcd(D1,D2)=1.
```

But LRC packet safety is not multiplicative, because the safe band is
archimedean.  The defect should be kept as a labelled atom.

S102 examples:

```text
cover_84:
  rate(7*13)=1/36, product rate(7)*rate(13)=0, defect=1/36
  rate(7*83)=1/123, product=0, defect=1/123

AP12_182:
  rate(3*41)=1/20, product=0, defect=1/20
  rate(7*83)=5/246, product=0, defect=5/246
```

These defects are not bugs in the phi ledger.  They are exactly the LRC
analogue of strong-ear insertion profiles: multiplicative capacity is present,
but witness feasibility needs the cross-factor boundary label.

## Proof target

The next theorem should be stated at packet level.

```text
For a primitive covering 13-row S, after deleting denominators killed by
divisibility, the exact-period packet hypergraph generated by the runner
danger classes cannot cover all scaled packets.
```

More usefully for the current proof DAG:

```text
1. finite low denominators killed or certified by divisibility/resonance;
2. local HYP-2883/HYP-2884 signed-current conservation lifted on exact-period
   residue fibers;
3. CRT multiplicativity defects retained as labelled atoms;
4. high-denominator incoherent residue packets routed to the existing
   spectrum/L2/Part-A witness floor.
```

This would not bypass OPEN-Q-108.  It would give OPEN-Q-108 a sharper discrete
front end: before applying Beurling-Selberg, L2, or slow-fast Part A, the
residual packet object has already been stripped of false fixed-basis claims
and organized by exact-period residues.

## Assumption challenge

Candidate vertices considered:

- runners;
- denominators;
- exact-period units `a mod D`;
- reduced exact denominators in `q=14*Vmax`;
- mod-7 residues;
- `chi_7` classes;
- affine-pair classes under `r -> 2-r`;
- CRT factors;
- support-six packet classes;
- proof obligations.

Chosen quotient: exact-period unit packets and quotient lenses on those
packets.  This preserves the LRC predicate "a rational time is a level-1/14
witness" and the Euler phi packet law.  It destroys raw speed magnitude except
through residues mod `D`; this is intentional, because divisor loading is the
actual obstruction to fixed finite bases.

The challenged assumption is that a successful finite sampled basis should
become a universal certificate basis.  The correct upgrade is not "one finite
basis"; it is "an adaptive exact-period atlas plus a proof that the atlas
cannot be covered."
