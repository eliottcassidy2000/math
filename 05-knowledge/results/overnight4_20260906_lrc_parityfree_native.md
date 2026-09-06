# Independent native audit, physical coarea, and the sharp norm-four ceiling

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The complete arbitrary-parity ternary-unit theorem in the
[fourth-round proof](overnight4_20260906_lrc_parityfree_probe.md) passes the
root audit: `mu(F_w)<=min_i E_i(w)<=6/55`, with equality in either bound
against `6/55` exactly at the primitive triple `(1,10,11)`.
This report supplies the independent native **all-projection** finite head,
an exact physical coarea identity, and the sharper norm-four consequence.
The [separate referee](overnight4_20260906_lrc_parityfree_audit.md) reconstructs
every coefficient polygon, owner residue fiber, and literal physical head.
No arbitrary chart entry, body floor, synchronization, or LRC(14) follows.

## 1. Literal interval replay of every projection

The universe is exactly `1<=a<b<c<=63`, `gcd(a,b,c)=1`, `3` not dividing
any speed. There is no parity restriction, no inherited coefficient filter,
and no exclusion of additive or norm-four triples. In particular the loops
include the empty-defect row `(1,2,4)`, which an old odd loop starting at
five would omit. The raw producer's TSV is comparison data, not a source
of the independently enumerated universe or interval endpoints.

The [C++ source](../../04-computation/overnight4_20260906_lrc_parityfree_native.cpp)
retains the `Cursor` and `literal` implementation of the previously audited
`lrc14_universal_literal_empty_core_sep06.cpp`. It changes the universe and
consumer checks. On denominator `42abc`, it builds all six literal distinct
sheet assignments, scans their exact intervals, and computes every capped
pair intersection and physical triple intersection. Tied right endpoints
advance simultaneously; contacts of zero length contribute nothing.
It uses no carrier formula to compute these quantities.

The wrap argument survives removal of parity. If a zero-label half interval
of speed `d` meets a nonzero-label interval of speed `s`, then
`11/(42s)<1/(14d)`, hence `s>11d/3`. The latter interval is shorter than
the former half interval, so a wrap cut cannot lower the effective omitted
capacity. A nonwrapping partner lies on one side of the cut, preserving
pair intersection length as well. This is the same typed endpoint argument
as in the independent third-round audit, now with the parity hypothesis
explicitly removed. Fixed-sheet six-sparsity is also parity-free.

The [retained output](overnight4_20260906_lrc_parityfree_native.out) records:

| Exact check | Result |
|---|---:|
| Independently generated eligible triples | 10,074 |
| Empty native contact dictionaries | 1,028 |
| Native contacts, summed over all rows | 117,864 |
| Physical masses above the old odd-only `6/77` ceiling | 151 |
| Maximum physical mass and maximum selected projection | `6/55` |
| Equality triples for either maximum | only `(1,10,11)` |
| Optimization-live checks | 181,385 |

Every row's three projection numerators and physical numerator agree with
the raw producer on the same denominator. The native contact count is
exactly three times the full raw carrier count. The extremizer has
`E=(6/55,12/77,23/154)`, so the equality statement concerns the selected
projection; it does not say every projection obeys `6/55`.
At `(4,7,11)`, physical mass `215/2156` is strictly below the selected
projection `223/2156`. Both distinctions are explicit controls.
Eight inherited positive/hostile controls include heights beyond the head.

Compilation at `-O0` and `-O3 -DNDEBUG` gives byte-identical output. Reproduce
from the repository root, choosing an executable path outside tracked files:

```text
g++ -std=c++17 -O3 -DNDEBUG 04-computation/overnight4_20260906_lrc_parityfree_native.cpp -o C:/w/overnight4_parityfree/native_o3.exe
C:/w/overnight4_parityfree/native_o3.exe 05-knowledge/results/overnight4_20260906_lrc_parityfree_probe_head63.tsv
```

Raw LF-byte SHA-256 values:

```text
source 459f7080672d61e72df1e836c15a27304ad4d847858c89fe05c49982bc996eda
output 67a44d7881d90a14d52aed3a8d510841dc191c9784c8e13a1195cb1e111bb0f2
raw TSV c3d33fdd136245aafe512b04963a6eb6f1b5db6f1a572a3e8535ef59d01a09fa
```

## 2. The all-height inference is a separate audited argument

The head alone proves no unbounded statement. The proof retains every
primitive coefficient pattern, modulo signed coordinate symmetries,
including odd relation norms. The independent referee reconstructs all
747 patterns of maximum coefficient at most 18 by cube-edge geometry,
using an algorithm distinct from the producer's fractional knapsack.
Their exact slopes and intercepts agree record by record, with semantic
digest `cf808062354debbefc1d8ead8ad0d10e9da5427cb42b8f083b6af24d0059c87c`.

For a nonempty admissible defect list, the general sector gives
`N<15c/98+2S/7+4/3`. The count target is `14c/55`; their slope difference
is `547/5390`. The sufficient cutoff is therefore
`c>=1540S/547+21560/1641`. A primitive short relation has
`S<4sqrt(c/3)`. If `S<=18`, the cutoff is less than 64; if `S>=19`,
`3S^2/16-1540S/547-21560/1641` is positive at 19 and increasing.
The coefficient tail above 18 and the separate norm-five pattern have
the strict bounds stated in the proof. This covers the unbounded general
sector; the additive and norm-four physical profiles handle its low cases.

The empty-defect case is essential. Pattern `(0,1,2)` allows no integer
defect and has `N=F=B=0`; a blanket strict inequality `N<F+B` would assert
`0<0`. The repaired proof first discharges this case with exact emptiness.
This was a rejected provisional implication, not a retraction of the
correctly scoped odd theorem. The residue-coordinate gate and strict
support endpoints both survive the repair.

For the selected-network upgrade, the
[incoming additive theorem](lrc14_additive_parity_empty_core_sep06.md) is
a proved dependency. It controls the complete additive projection, not
just physical mass. Its normalized profiles have continuum selected bound
`39/392`, strict residue error `4/(7c)`, and network cutoff 60. The root
and separate referee reviewed these formulas, the monotone cutoff jumps,
and the crossing at `a/c=1/4`. On norm four, a single doubled-coefficient
roof is pointwise minimal on the entire line, so selection commutes with
summation for that specifically proved reason. It does not commute in
general. The original physical proof has its own sufficient additive tail
and does not depend on the incoming additive ceiling.

## 3. A saturated physical coarea identity

Here is the reusable exact mechanism behind the low-relation physical
integrals. Let `w` be any primitive positive integer triple and `v` any
primitive integer relation perpendicular to it. For any radius `r>0`,
let `f(t)` be the length of the set of real `y` for which

```text
t n_1-y w belongs to [-r,r]^3,
```

where `z` is an integer vector with `z.w=1` and `n_1=v cross z`.
Then `w cross n_1=v`, and the planar parameterization `(t,y)` has area
Jacobian `||v||_2`. In particular `(w,n_1)` is a saturated integer basis
of the plane lattice `v^perp`: its covolume equals that of a primitive
normal. Consequently

```text
integral_R f(t) dt = area([-r,r]^3 intersect v^perp)/||v||_2.       (A)
```

For the physical radius `r=3/14`, the interval-intersection identity says
`f(t)=max(0,L_w(tv))`. Formula (A) is independent of the speed ratios
inside this relation plane. Strict versus weak cube boundaries change no
area or interval integral. Coordinate permutations and sign changes
preserve the right side. Primitivity fixes the integer multiplier spacing;
discarding it would lose the finite residue-sampling rule.

This is related to the incoming
[variable-radius coarea theorem](variable_radius_empty_core_sep06.md),
which integrates the **carrier width over affine defects** and obtains
`4r^2 sum(w)`. Formula (A) instead integrates the **physical intersection
length along one carrier line**. The two integrals have different units
and must not be identified. No external priority is claimed.

Solving for a coordinate of maximum coefficient gives these exact areas:

| Absolute relation pattern | Physical line integral |
|---|---:|
| `(1,1,1)` | `3r^2` |
| `(1,1,2)` | `2r^2` |
| `(1,2,2)` | `7r^2/4` |

For `(1,1,1)`, the square in the other coordinates loses two triangles of
total area `r^2`. For `(1,1,2)`, the whole square is allowed and division
by the maximum coefficient two gives `2r^2`. For `(1,2,2)`, the square
loses two triangles of area `r^2/4` each, and the same division gives
`7r^2/4`. These identities concern physical line profiles. They imply
neither pointwise equality of different roofs nor a bound for an arbitrary
selected projection without an additional fixed-roof argument.

## 4. Sharp arbitrary-parity norm-four theorem

For primitive distinct positive ternary-unit triples admitting a relation
with absolute pattern `(1,1,2)`, the stronger exact theorem is

```text
mu(F_w)=min_i E_i(w)<=11/140,
equality iff w=(2,11,20).                                      (B)
```

The three sorted sectors are `c=2a+b`, `c=a+2b`, and `2b=a+c`.
The doubled coefficient is at `a` or `b`. Its complementary speeds include
`c`, and its capped roof is pointwise minimal throughout the live line;
the endpoint proof appears in Section 6 of the main proof. Since
`r||v||_1=6/7<1`, only defect zero is possible. The complete dictionary
is exactly the carriers `kv` with `3` not dividing `k` and `f(k)>0`.

The preceding coarea identity gives `I=9/98` exactly. For a nonnegative
even profile decreasing on the positive half-line, with `f(0)=q=3/(7c)`,
layer-cake integration of

```text
#{k in Z: 1<=k<T, 3 does not divide k} < 2T/3+2/3
```

gives the strict bound

```text
mu(F_w) < (2/3) I+(4/3)q = 3/49+4/(7c).                      (C)
```

For `c>=33`, (C) is strictly below `11/140`; its margin at 33 is
`1/32340`. The exact remaining universe comprises the **82** eligible
norm-four triples with `c<=32`. All belong to the independently matched
native head. Direct fixed-profile enumeration has maximum `11/140`
only at `(2,11,20)`. The
[independent norm-four audit](overnight4_20260906_lrc_norm4_audit.md)
checks the profile, strict quadrature, coarea constants, and finite base.
Thus (B) is an all-height theorem, not a height-32 extrapolation.

The old odd-domain maximum `6/77` at `(1,5,11)` remains a sharper result
on that smaller domain. Statement (B) settles the norm-four sector; the
incoming results now extend it to all nonadditive triples as follows.
For nonprimitive ternary-unit triples, divide by their common divisor,
which is prime to three and permutes sheet labels while preserving Haar
measure. Equality in (B) is then exactly those allowed positive dilates.

## 5. Combined nonadditive ceiling and exact old-threshold classification

**PROVED, relative to the named incoming theorems and independently audited
finite bases.** Continue to assume `1<=a<b<c`, `gcd(a,b,c)=1`, and `3` not
dividing `abc`. The complete physical mass and selected minimum satisfy:

```text
If a+b!=c: mu(F_w)<=min_i E_i(w)<=11/140,
with equality for either quantity iff w=(2,11,20).             (D)

For either Q(w)=mu(F_w) or Q(w)=min_i E_i(w):
Q(w)>6/77 iff [a+b=c and w!=(1,4,5)] or w=(2,11,20);
Q(w)=6/77 iff w=(1,5,11);
Q(w)<6/77 otherwise.                                         (E)
```

For (D), every nonadditive triple either has a signed `(1,1,2)` relation,
has a signed `(1,2,2)` relation, or has none of the three low patterns.
The first case is (B). Incoming
[THM-4441 / signed-122 sharp ray closure](../../01-canon/theorems/THM-4441-lrc14-signed-122-sharp-ray-closure.md)
gives `min E<=46/665<6/77` and physical mass `<=51/770<6/77` in the second.
Incoming [THM-4437 / all-parity reduction](../../01-canon/theorems/THM-4437-lrc14-all-parity-network-reduction-to-three-low-circuits.md)
gives `min E<6/77` in the third. These cases exhaust the domain even if
some relation classes overlap; neither smaller bound can introduce another
equality in (D). Both incoming results were promoted before `058a8ded9`.
This is a synthesis of proved components, not an extrapolation from a new
nonadditive census.

For (E), the generic and norm-five cases are already strict. On norm four,
(C) is strictly below `6/77` for `c>=35`; at 35 its upper bound is
`19/245`, with margin `1/2695`. The complete 88-row norm-four base `c<=34`
has exactly one triple above the threshold, `(2,11,20)`, and exactly one
equal, `(1,5,11)`, for both quantities.

On the additive line, the independent
[third-checkpoint discrepancy proof](overnight3_20260906_lrc_additive_audit.md)
gives `mu(F_w)>=9/98-6/(7c)>6/77` for `c>=62`; at 62 the lower bound is
`237/3038`, with margin `3/33418`. The complete 146-row additive base
`c<=61` has both quantities strictly above `6/77` except `(1,4,5)`, where
both equal `1/28`. Since `mu<=min E`, the physical lower bound also proves
the selected lower bound. These strict tails and exhaustive finite bases
establish both directions and the equality boundary in (E).

The [finite-base verifier](../../04-computation/overnight4_20260906_lrc_threshold_classification.py)
independently enumerates both coefficient-identity universes and reads the
hash-pinned full native-audited head. Its
[transcript](overnight4_20260906_lrc_threshold_classification.out) passes
475 always-active gates, identically in ordinary and optimized Python.
An independent norm-five reconstruction also checked these two subbases
and tail margins. Reproduce from the repository root:

```bash
python3 04-computation/overnight4_20260906_lrc_threshold_classification.py
python3 -O 04-computation/overnight4_20260906_lrc_threshold_classification.py
```

```text
SHA-256, LF bytes:
source da16bcf875e4b111fd31d3d7a341357c674e678b21e5d4ccffe3f52c34297e96
output 0629e2ab0cd320c08197219adf9377848976d66efb87f8b10bc077a4049eced9
```

For common ternary-unit dilates, apply (D)-(E) to the primitive reduction.
The map preserves Haar measure and permutes sheet labels. The classification
concerns a local tail statistic; it does not supply a body-safe anchor,
chart entry, or a proof of LRC(14).
