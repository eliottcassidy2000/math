---
id: HYP-2167
status: SUPPORTED by S611 diagnostic + S612 bounded carry-fiber search + S654 high-carry apex/parity bridge; full carry theorem open
source: user-2026-06-03; codex-2026-06-03-S611; monad-compute-2026-06-03-S3; codex-2026-06-05-S654
related:
  - HYP-2175
  - HYP-2230
  - HYP-2166
  - HYP-2165
  - HYP-2164
  - HYP-2163
  - HYP-2162
  - HYP-2161
  - HYP-2101
  - HYP-2100
  - HYP-2099
  - HYP-2088
  - THM-401
  - THM-407
---

# HYP-2167: n=14 lift/CRT conservativity is a carry-fiber theorem

## Claim

The remaining lift/CRT gap after the HYP-2164 least-positive quotient,
HYP-2165 owner bridge, and HYP-2166 quotient tower is not just "try more
integer representatives."  It has a specific number-theoretic coordinate.

For any lifted speed over the `C=27` residue quotient, write

```text
v = r + 27 k,        1 <= r <= 26.
```

Since `27 == -1 (mod 14)`,

```text
v == r - k (mod 14).
```

Thus the carry vector `k` is the CRT glue between the `Res_27` shell coimage
and the `n=14` clock ledger.  A conservative n=14 proof object must remember
the fiber coordinate `k`, or an equivalent owner/certificate datum, not only
the least-positive residue section `r`.

## S611 Evidence

`04-computation/lrc_n14_carry_conservativity_s611.py` runs an exact diagnostic
on the known primitive floor rows AP and `V*`.

First it recovers the THM-407 shell fold:

```text
G=<2,-1> shell orbits modulo 27:
  gcd=1: (1, 2, 4, 5, 7, 8, 10, 11, 13)
  gcd=3: (3, 6, 12)
  gcd=9: (9,)
```

Then it compares actual unit scalar lifts with their least-positive
`Res_27` shadows:

```text
scalar probes: 36
actual scaled rows at floor: 36/36
least-positive shadow route histogram: {'floor shadow': 3, 'strict shadow': 33}
shadow below-floor rows: 0
```

The only floor shadows are exactly the HYP-2164 section-floor rows:

```text
AP:u=1
AP:u=2      = nonprimitive 2*AP
V*:u=1
```

All other AP/`V*` scalar lifts remain floor as actual integer rows by scaling
invariance, but their least-positive shadows become strict.  This proves the
least-positive section is not conservative for floor behavior by itself.

Finally S611 probes local carries inside the AP and `V*` residue fibers.  Add
one or two `+27` wraps to canonical AP or `V*` residues and compute exact
maximin:

```text
AP weight=1: 13/13 strict, min M=1/13
AP weight=2: 78/78 strict, min M=1/12
V* weight=1: 13/13 strict, min M=2/25
V* weight=2: 78/78 strict, min M=1/12
local floor rows: 0
local below-floor rows: 0
```

So isolated carry defects do not produce a new floor family near AP or `V*`.
Any new lifted floor row must use a globally coherent carry pattern.

## Structural Interpretation

The current n=14 improvements form a quotient tower:

```text
integer speed rows
  -> (Res_27 shadow, carry fiber k)
  -> least-positive Res_27 section
  -> G=<2,-1> shell strata {gcd 1, gcd 3, gcd 9}
  -> owner/certificate discharge.
```

HYP-2164 classifies the least-positive section: primitive floor is AP plus
`V*`.  HYP-2165 shows the fixed-boundary parity scaffold and the `C=27` owner
layer separate cleanly: parity lives in the 64 self-converse classes, while
cheap-pair/positive-measure certificates live in owner labels.

HYP-2167 adds the missing fiber coordinate.  In categorical language, the
least-positive `Res_27` representative is a section of the coimage map, not the
whole coimage object.  The carry vector is descent data in the fiber.  Yoneda
probes such as the `n`-clock and endpoint-owner ledgers detect that descent
data because `27 == -1 (mod 14)`.

In the anti-Poisson frame, this says the all-orders cancellation is rigid in
the fiber.  AP and `V*` have true floor cancellation.  Local carry perturbations
break it immediately.  The only way a lifted row can remain at the floor is if
the carry pattern is globally coherent enough to imitate a scalar AP/`V*`
carry, or if it routes to a known multiple/Cprime owner certificate.

## Proof Target

A plausible carry-conservativity theorem for n=14 is:

```text
Let V = R + 27K be a primitive lifted row over a classified Res_27 shadow R.
If M(V) <= 1/14 after the HYP-2163 no-multiple clock split, then either:

  1. (R,K) is cohomologous to the AP or V* scalar carry pattern, hence
     normalizes to AP or V*;
  2. the carry changes the n-clock ledger in a way that gives a clock witness;
  3. the owner layer gives a HYP-2165 cheap-pair or positive-measure exit; or
  4. the Cprime/multiple branch gives a bounded CRT contradiction.
```

The new object is the carry cocycle

```text
K_i = (v_i - r_i)/27.
```

The base quotient records `r_i`.  The n-clock sees `r_i - K_i`.  Pair-sum
pinches see both `r_i+r_j` and `27(K_i+K_j)`.  Owner certificates should be
phrased as constraints on this cocycle rather than as an arbitrary lift search.

## Tournament Analysis

S611 uses carry probes as tournament vertices.

```text
vertices: 36 scalar section probes + 182 local carry probes = 218
observable: route, exact margin above 1/14, carry span, clock-change count
switch: proof-burden order, with label tie Hamiltonian path
fingerprint: transitive ledger; 218 singleton SCCs; 0 directed 3-cycles
edge flips between proof-burden and carry-complexity gauges: 8708/23653
```

The edge flips matter.  Carry complexity is not itself the proof certificate;
the certificate is the interaction between carry data and the LRC floor probes.

## Assumption Challenge

Candidate tournament vertices considered:

```text
runners,
gaps,
fixed circle sections,
section boundaries,
wall-crossing events,
residues,
cover arcs,
Fourier modes,
matroid circuits,
proof obligations,
lift carries.
```

This hypothesis chooses lift carries because HYP-2164 already classifies the
least-positive base quotient and HYP-2165 already tests the owner bridge.  The
challenged assumption is that the base `Res_27` shadow is enough.  S611 shows
it is not: the same floor scalar lift can project to a strict section shadow.

## Honest Status

This is not a proof of LRC n=14.  It is a sharper formulation of the last
missing theorem.

What is proved by computation:

```text
1. the carry/CRT identity is exact;
2. scalar AP/V* floor lifts collapse to only three floor section shadows;
3. AP/V* local carry moves of Hamming weight one or two are all strict.
```

What remains open:

```text
prove global carry-fiber conservativity, especially for nonlocal carry patterns
and for lifted rows whose least-positive shadows are not AP or V*.
```

## S654 parity/apex carry update (codex-2026-06-05-S654)

S654 connects this carry theorem to the repo's even/odd pair-carrier work.  For
the same lifted speed

```text
v = r + 27k,
```

there are two exact projections:

```text
v mod 2  = r+k mod 2,
v mod 14 = r-k mod 14.
```

Thus `k` simultaneously toggles the parity word and decides the apex
zero-divisor condition:

```text
14 | v  iff  k == r mod 14.
```

Pair-sum denominators satisfy the paired version:

```text
14 | v_i+v_j  iff  k_i+k_j == r_i+r_j mod 14.
```

S654 tests the forced multiple-of-14 branch that any counterexample must enter.
For each AP and `V*` coordinate, it sets the minimal apex carry
`k_i == r_i mod 14` and computes exact maximin.  This reaches high-carry
residues `7..13` that S612's `L1<=6` search could not touch.

Every single-apex lift is strict:

```text
AP best single-apex tax: M=28/365, margin=27/5110, via r=13 -> k=13.
V* best single-apex tax: M=2/25,  margin=3/350,  via r=13 -> k=13.
```

Then S654 adds all one- and two-coordinate `+27` toggles around every minimal
apex lift:

```text
rows screened: 2054
approximate near-floor rows exact-checked: 0
exact floor rows in checked set: 0
exact below-floor rows in checked set: 0
```

The full Boolean parity lattice also has only the zero-carry base row at the
floor in the near-floor screen.  This does not prove global conservativity, but
it moves the residual: a new floor lift must use a larger globally coherent
carry cocycle or an owner/Cprime route, not the first parity/apex bridge.

See HYP-2230; `04-computation/lrc14_parity_carry_bridge_s654.py`;
`05-knowledge/results/lrc14_parity_carry_bridge_s654.out`.

## S612 compute update (monad-compute-2026-06-03-S3)

S611 only tested *isolated* wraps (local carry moves of Hamming weight 1 and 2,
each a single `+27`).  S612 closes that gap within a bounded radius by an
exhaustive carry-fiber search over both floor shadows, lifting `row = R + 27*k`
and computing the EXACT loneliness radius for every `k` in:

- the **L1-budget** fiber `sum(k) <= 6` (27132 vectors per shadow, all
  magnitudes — reaches the coherent / large-magnitude patterns S611 could not);
- the **full Boolean lattice** `{0,1}^13` (8192 vectors per shadow, every
  Hamming weight 0..13 at magnitude 1).

`04-computation/lrc_n14_carry_fiber_search_s612.py` /
`05-knowledge/results/lrc_n14_carry_fiber_search_s612.out`.  Method: a fast
float maximin screens all ~70k rows; every row whose float `M` lands within
`1e-6` of the floor is re-checked with the VERIFIED S611 exact oracle (the float
maximin was validated against exact on 305 rows, worst error `2e-15`).

Result (sanity gate AP, V*, 2*AP all at floor):

```text
                              AT floor (M=1/14)   below floor   min M over nonzero k
  L1<=6   over AP                    0                 0          1/13
  L1<=6   over V*                    0                 0          7/88
  {0,1}^13 over AP                   0                 0          1/13
  {0,1}^13 over V*                   0                 0          2/25
```

Every nonzero carry — isolated OR coherent, at every Hamming weight 1..13 — is
strictly above `1/14`, and none drops below the floor (no LRC(14) counterexample
in this radius).  So within this bounded fiber neighbourhood the least-positive
section AP / V* is the UNIQUE floor representative.  This upgrades the open
clause: the only surviving escape routes for a new floor lift are a *large*
coherent carry outside `L1<=6` / magnitude>1 spread over many coordinates, or a
multiple / owner certificate (HYP-2165) — exactly the residual proof targets.

## See

`04-computation/lrc_n14_carry_conservativity_s611.py`,
`05-knowledge/results/lrc_n14_carry_conservativity_s611.out`,
`07-reflections/lrc-n14-carry-fiber-conservativity-s611.md`,
`05-knowledge/hypotheses/HYP-2175-lrc-dimension-descent-salience.md`,
`05-knowledge/hypotheses/HYP-2166-lrc-n14-res27-quotient-tower-conservativity.md`,
`05-knowledge/hypotheses/HYP-2164-lrc-n14-res27-pinch-certificate.md`,
`05-knowledge/hypotheses/HYP-2165-lrc-n14-res27-fixed-class-bridge.md`,
`01-canon/theorems/THM-407-twisted-involution-shell-reduction-of-the-LRC-additive-residual.md`.
