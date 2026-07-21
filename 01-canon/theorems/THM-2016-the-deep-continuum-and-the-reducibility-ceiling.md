---
id: THM-2016
title: "THE DEEP CONTINUUM AND THE REDUCIBILITY CEILING — three worked threads on the continuum coordinates (THM-2013). (A) ABSOLUTE-|R| AND LOCAL-PROFILE AUDIT: in the n=7 c3=12 shell, (char_A,|R|,score) resolves 28/47 and 4&5-profiles raise this to 41/47, leaving six twin pairs; at c3=11 score raises 30→46 and profiles raise 46→50, leaving two pairs; at c3=13 profiles raise 13→14, leaving one pair. (B) REDUCIBILITY CEILING (PROVED + verified): max c₃ over reducible tournaments = c₃_max(n−1), by discrete convexity of c₃_max rather than the false bound of an SCC sum by its largest summand. (C) H IS A THERMOMETER: mean/min/max H increase monotonically with cyclic temperature, locating the H≥disc binding case at τ=1."
status: >
  VERIFIED (boxeph-2026-07-21-S200), exhaustive n≤7. (A) With the advertised absolute |R|, the n=7
  hot-shell resolution rows (base; +score; +4-profile; +4&5; +score+4&5) are
  c₃=12: (28,28,28,41,41), c₃=11: (30,46,46,50,50), and c₃=13: (13,13,13,14,14).
  Final unresolved groups are respectively six, two, and one twin pairs.  Thus local density adds
  beyond score in all three shells, while score also contributes at c₃=11. (B) max-reducible-c₃
  = c₃_max(n−1) exactly for n=4..7 (1,2,5,8 = c₃_max(3..6)); c₃_max(n)=1,2,5,8,14 (n=3..7,
  = n(n²−1)/24 odd, n(n²−4)/24 even); τ_red=c₃_max(n−1)/c₃_max(n)→1. Corrected proof: the increments
  c₃_max(n+1)−c₃_max(n)=T_floor(n/2) are nondecreasing, so concentrating an SCC-size partition can
  only increase its sum; hence Σ_SCC c₃_max(|SCC|)≤c₃_max(n−1). (C) H strictly increases in mean with τ across all 15 shells at
  n=7 (1,3,6.7,10.5,14.4,23.9,33.4,44,60.9,79.7,96,114.5,134.7,153.9,178.3); the spread widens toward
  mid-hot. Uses THM-2013 (temperature), THM-1926 (cycle spectrum), THM-1966 (|R|), THM-1862 (c₃ over
  SCCs), death-star-S84 (H≥disc quasirandom).
source: boxeph-2026-07-21-S200 (owner: work the S199 open threads + one or two picked up from other agents)
depends_on: [THM-462, THM-1862, THM-2000, THM-2005]
related:
  - THM-2013  # continuum coordinates (this resolves its L3 residue + condensation threshold)
  - THM-1979  # the spectrum (single point -> continuum)
  - THM-1966  # mac-mini |R| (the L2 coordinate that collapses in the deep continuum)
  - THM-1926  # my zeta / cycle spectrum (=char_A, the L1 coordinate)
  - THM-1862  # order-join: c₃ = Σ over strong components
  - THM-2000  # discrete-convex increment law for c₃_max; repairs the ceiling proof
  - THM-2005  # condensation product, harmonic-hazard word, and parity-shuffle tax
  - MISTAKE-216  # records the dropped-SCC-summands error and its repair
  - MISTAKE-217  # corrects signed-R computations that were labeled |R|
  - "07-reflections/coordinates-for-the-continuum-cyclic-temperature-and-the-cycle-spectrum-boxeph-S199.md"
script: 04-computation/continuum_threads_boxeph_S200.py (+ .out)
---

# THM-2016 — the deep continuum and the reducibility ceiling

## Thread A — local density is an L3 probe, after auditing the score sidecar

Extending the coordinate budget of THM-2013 (`L0` temperature, `L1` cycle spectrum = char_A, `L2`
`|R|`), we tested both the score sequence and **local structure** — the `k`-profile (census of induced
`k`-vertex subtournaments).  They must be compared explicitly: a profile cannot be credited for a
split already made by score.  In the n=7 hot shells:

| shell | (char_A,\|R\|) | +score | +4-profile | +4&5-profile | +score+4&5 |
|---|---:|---:|---:|---:|---:|
| τ=6/7 (c₃=12, 47 cls) | 28 | 28 | 28 | **41** | **41** |
| τ=11/14 (c₃=11, 52) | 30 | **46** | **46** | **50** | **50** |
| τ=13/14 (c₃=13, 15) | 13 | 13 | 13 | **14** | **14** |

Local density genuinely resolves beyond score in every tested shell: `28→41`
at `c₃=12`, `46→50` after score at `c₃=11`, and `13→14` at `c₃=13`.
The final keys leave six, two, and one **twin pairs**, respectively—not “3/47
classes.”  These are non-isomorphic tournaments agreeing on temperature, the
full cycle spectrum, the absolute Rédei magnitude, score, and all 4- and
5-vertex densities. The deep continuum is
**invariant-resistant**: no cheap global (spectral/Rédei-magnitude) or local (flag) coordinate separates it; it
needs essentially the whole object. This is the sharp form of "enumeration is forced at the hot core."

## Thread B — the reducibility ceiling (proved)

> **Reducibility ceiling.** For `n>=2`, the maximum `c₃` over **reducible**
> (non-strongly-connected) tournaments on `n` vertices is `c₃_max(n−1)`.

*Proof.* Put `f(m)=c₃_max(m)`, with `f(1)=f(2)=0`.  A 3-cycle lies inside one
strong component, so for SCC sizes `n₁,...,n_r`, where `r≥2`,

```text
c₃(T)=Σ_i c₃(S_i) ≤ Σ_i f(n_i).
```

THM-2000 proves the all-`m` increment law
`f(m+1)-f(m)=T_floor(m/2)`.  These increments are nondecreasing.  If
`a≥b≥2`, transferring one vertex from the smaller part to the larger therefore
does not decrease `f(a)+f(b)`.  Iterating over the SCC-size partition
concentrates it at `(n-r+1,1,...,1)` and gives

```text
Σ_i f(n_i) ≤ f(n-r+1) ≤ f(n-1).
```

Equality is attained by a max-cyclic `(n-1)`-tournament joined in order to a
singleton. ∎

Verified exactly: max-reducible-c₃ = `c₃_max(n−1)` = 1,2,5,8 for n=4..7. The attained
**reducible ceiling** (for `n>=3`, when `c₃_max(n)>0`)
```
        τ_red = c₃_max(n−1)/c₃_max(n)   (= 1/2, 2/5, 5/8, 4/7 for n=4..7)
```
is exact and is attained by reducible classes; for `τ > τ_red` every class is
strongly connected.  THM-462's all-`n` no-holes theorem realizes every integer
`c₃`-level from zero through `c₃_max(n)`.  Therefore the first all-strong shell
is

```text
        τ_all=(c₃_max(n−1)+1)/c₃_max(n)   (= 1, 3/5, 3/4, 9/14 for n=4..7).
```

Since `c₃_max(n) = n(n²−1)/24` (odd),
`n(n²−4)/24` (even), `τ_red → 1`: reducible tournaments persist to nearly the hot edge (as a vanishing
fraction, THM-1978), and only the sliver `(τ_red, 1]` is purely strong.

### The ceiling ratios form a parity-shuffled harmonic word

The same ceiling has an exact reciprocal-sequence interpretation. Put
`τ_red(n)=c₃_max(n−1)/c₃_max(n)` for `n>=4`. Then

```text
τ_red(n)=(n−1)/(n+2)  (n even),       τ_red(n)=(n−3)/n  (n odd),
3/(1−τ_red(n))=6,5,8,7,10,9,... .
```

The last sequence is the adjacent-pair reversal of `{5,6,7,...}`. Hence,
for `Re z>1`,

```text
sum_(n>=4)((1−τ_red(n))/3)^z = zeta(z)−sum_(m=1)^4 m^(−z).
```

The ratios also telescope:

```text
product_(j=4)^N τ_red(j)=1/c₃_max(N),
sum_(N>=3) product_(j=4)^N τ_red(j)=75/4−24log2.
```

Thus the maximum-`c₃` reciprocal mass is the partition sum of cumulative
reducibility hazards. Sorting the same hazard support gives partition sum `2`;
THM-2005 proves that the order-sensitive parity-shuffle tax is
`67/4−24log2>0`.

### The ceiling is a reindexing; the first strong shell is an analytic shift

Put `f(m)=c₃_max(m)`, `m>=3`.  The reducibility-ceiling denominator sequence
`c₃_max(n-1)`, `n>=4`, has exactly the same support and full Dirichlet profile
as `f`; THM-2000/2005 give its abscissa `1/3` and mass
`75/4-24log2`.  It is not a new harmonic support.

The first all-strong numerator sequence is genuinely new:

```text
g(m)=f(m)+1,                    m>=3.
```

It still has abscissa `1/3`, because `f(m)~m^3/24`.  More precisely,

```text
D_g(s)-D_f(s)=sum_(m>=3)[(f(m)+1)^(-s)-f(m)^(-s)]                  (RC)
```

extends holomorphically to `Re(s)>-2/3`, far across the common convergence
line.  Indeed, on compact subsets of that half-plane,

```text
|(a+1)^(-s)-a^(-s)|
 =|s integral_a^(a+1)x^(-s-1)dx| <= C a^(-Re(s)-1),
```

and `a=f(m)` is cubic.  Thus crossing the reducibility ceiling by one integer
shell changes the finite reciprocal mass but not the critical Dirichlet
singularity.  This is the analytic form of “one discrete level above the same
hot edge.”

## Thread C — H is a thermometer (locating the H≥disc binding case)

The Hamiltonian-path count H, the flagship #P-hard invariant, is **macroscopically a function of
cyclic temperature**: mean H over the n=7 shells rises monotonically
`1, 3, 6.7, 10.5, 14.4, 23.9, 33.4, 44, 60.9, 79.7, 96, 114.5, 134.7, 153.9, 178.3` as τ goes 0→1,
with the *spread* (the free-coordinate contribution) widening through the mid-hot shells. So `τ`
predicts H's scale and the cycle-spectrum/|R| coordinates predict its fine value — consistent with
"scores explain ~97% of H-variance." This locates death-star-S84's result: `H ≥ disc` is saturated at
the **hot center τ=1** (regular/quasirandom), where H is largest and the inequality is tightest —
exactly where enumeration fails and the coordinate description takes over.

## Net

The continuum coordinates (THM-2013) now have a sharp top and bottom: **H is the thermometer** (the
temperature is even readable off the flagship invariant), the **reducibility ceiling** `c₃_max(n−1)`
pins the cold edge exactly, and the **deep hot center is invariant-resistant** — the one place where
local densities, the spectrum, and `|R|` all fail together, and there is nothing left to do but hold
the whole tournament.
