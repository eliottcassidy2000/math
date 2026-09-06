# Parity-free physical and selected-network 6/55 ceilings

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The proved theorem is:

```text
w=(a,b,c), 1<=a<b<c, gcd(a,b,c)=1, every w_i nonzero mod3:
mu(F_w)<=min_i E_i(w)<=6/55,
mu(F_w)=6/55 iff min_i E_i(w)=6/55 iff w=(1,10,11).           (1)
```

Here `F_w` is the literal scale-three union over the six distinct sheet
assignments, defined below. There is **no oddness hypothesis**. A complete
raw-address calculation verifies all 10,074 eligible triples with `c<=63`.
The all-height reduction has strict tails above this head: the general
coefficient sector starts at 64, the additional `(1,2,2)` sector at 34,
the additive sector at 50, and the norm-four sector at 18.

No claim about arbitrary entry or LRC(14) is made. Sections 1-5 establish
the physical proof independently of the concurrent additive
ceiling. Section 6 upgrades to the selected network projection using that
audited additive result and an explicit fixed-roof lemma for norm four.
It does not invoke the odd network ceiling after dropping parity.

## 1. Inheritance, hostiles, and the preserved physical quantity

The closest proved mechanism is
[THM-4434 / universal scale-three network projection bound](../../01-canon/theorems/THM-4434-lrc14-universal-scale-three-network-projection-bound.md),
verified as promoted in incoming `d2f64b809`. Its zonotope, scalar residue,
open-interval count, and projected-relation-lattice arguments do not use
oddness. Its coefficient universe and the parity step in the final cutoff
do use oddness; both must be replaced here.

The canonical parity hostile is `(2,5,7)`, whose physical mass is `22/245`.
Deleting additive triples does not restore the old `6/77` ceiling:
`(2,11,20)` has mass `11/140`. The current exact head through 60 had maximum
`6/55`, uniquely at `(1,10,11)`. The least-used useful sidecar is the
physical mass of the complete zero-defect line, whose integral retains
roof information that the carrier count discards.

The live concept board is: exact physical sheet mass; affine owner residues;
relation-line counts; coefficient-polytope convexity; shortest relation scale;
and central-section area. The new operation is to remove parity in the
coefficient compiler and separate its new short-relation exceptions before
applying the count consumer.

**Concurrent inheritance update.** Incoming `d7027ffb2` contains the
[audited all-height additive theorem](lrc14_additive_parity_empty_core_sep06.md),
which independently proves the full additive `6/55` ceiling, for both the
selected network projection and physical mass. Its sharper one-sided
layer-cake error gives physical cutoff 34. The independent physical branch
below keeps its sufficient cutoff 50, so this new result is compatible
additional evidence rather than a hidden premise. The incoming
[parity probe](lrc14_parity_empty_core_sep06.md) independently verifies all
20,648 triples through height 79. At that original `d7027ffb2` baseline,
THM-4437 was **RESERVED / UNPROVED EMPTY STUB**, so it was not used here.
By incoming `058a8ded9`, [THM-4437 / all-parity reduction](../../01-canon/theorems/THM-4437-lrc14-all-parity-network-reduction-to-three-low-circuits.md)
and [THM-4441 / signed-122 closure](../../01-canon/theorems/THM-4441-lrc14-signed-122-sharp-ray-closure.md)
are **PROVED**. They strengthen the generic and norm-five branches below;
the original independent proof remains valid. Together with our norm-four
refinement they give the [sharp nonadditive bound and exact old-threshold classification](overnight4_20260906_lrc_parityfree_native.md#5-combined-nonadditive-ceiling-and-exact-old-threshold-classification).

Explicitly define

```text
D_(w_i,s)={x mod1: ||w_i(x+s/3)||<1/14}, s in {0,1,2},
F_w=union_(pi in S_3) intersection_i D_(w_i,pi(i)).
```

The [parity-free physical address derivation](overnight3_20260906_lrc_additive_audit.md)
starts from these literal sheets. Put `y=3x`, `r=3/14`, and let `n` be the
unique nearest-integer vector with `|n_i-w_i y|<r`. Its owner labels are
`s_i=-w_i^(-1)n_i mod3`. The raw condition is one-periodic in `y`; its
three inverse x-branches cancel the length factor of one third. Thus
physical mass equals raw mass over one period.

The primitive cross-product address `C=w cross n` identifies `n mod Zw`
with the full rank-two relation lattice. The distinct-sheet condition is
exactly `C_i!=0 mod3` for all three coordinates. Write

```text
Lambda(w)={C in Z^3:C dot w=0, every C_i!=0 mod3,
                       14|C_i|<3(w_j+w_k) for each i},
q=3/(7c),
L_w(C)=min(q, [3(w_j+w_k)-14|C_i|]/(14w_jw_k): i=1,2,3).
```

The one-dimensional interval intersection identity gives, without parity,

```text
mu(F_w)=sum_(C in Lambda(w)) L_w(C) <= q N,
N=|Lambda(w)|.                                             (2)
```

Consequently `N<=14c/55` is a sufficient physical count gate. Retaining
all roofs rather than using this gate is essential in the two exceptional
shortest-relation sectors below.

## 2. The parity-free general count sector

Let `v dot w=0` be a primitive integer relation, `S=||v||_1`, and
`M=max_i|v_i|`. The exact integer defect is `delta=v dot n`. The affine
carrier line is characterized by `v cross C=delta w`; primitivity makes
its integer longitudinal step exactly `v`.

For completeness, the inherited owner decomposition has two types:

* If every `v_i` is a unit modulo three, only `delta=0 mod3` is live and
  two of the three longitudinal classes are live.
* Otherwise exactly one coefficient is zero modulo three, only
  `delta!=0 mod3` is live and one of three longitudinal classes is live.

Two zero residues would force the third to vanish because all speeds are
units, contradicting primitivity. The new standalone verifier checks all
192 finite residue fibers independently of any parity filter.

Let `f_v(delta)` be the real longitudinal width of the error-cube section.
As in THM-4434, set the width to zero at the strict support endpoints;
this matters for an actual zero coefficient. The inherited elementary
arguments first require an empty-list split. If the allowed integer defect
list is empty, there are no live carriers and `N=F_v=B_v=0`; discharge the
target directly. For a **nonempty** allowed defect list they give

```text
N < F_v(w)+B_v,
B_v <= 2S/7+4/3,
F_v(w)/c <= 6/49+4/(7M).                                (3)
```

These use the exact area `integral f_v=9(a+b+c)/49`, monotone section
quadrature, peak width `f_v(0)<=6c/(7M)`, and open-interval residue counts.
None uses `S` even. For `M>=19`, (3) gives

```text
F_v(w)/c <=142/931<15/98.                               (4)
```

For `M<=18`, remove the odd-theorem requirement that `S` be even. The
complete coefficient universe, up to coordinate permutation, is

```text
0<=p_1<=p_2<=p_3<=18,
at least two nonzero coordinates,
gcd(p_1,p_2,p_3)=1,
at most one coefficient zero modulo three,
p not in {(0,1,1),(1,1,1),(1,1,2),(1,2,2)}.             (5)
```

There are exactly **747** patterns, comprising 698 full-support and 49
actual-zero patterns. The first exclusion forces equal speeds. The next
two are the norm-three and norm-four physical branches. The final exclusion
is the new count exception handled immediately below.

The finite universe contains exactly one empty-defect pattern, `(0,1,2)`.
For instance `w=(1,2,4)` has relation `(0,-2,1)`; the owner type requires
`delta!=0 mod3`, while `|delta|<9/14` forces `delta=0`. Thus its dictionary
is empty. The independent referee caught the false blanket implication
`N<F_v+B_v`, which would read `0<0` here. The repaired split above retains
the strongest correct conclusion, `N=0`, and the verifier now checks both
the coefficient pattern and this speed hostile explicitly. The non-strict
slope and intercept bounds remain valid for the empty case.

For each pattern the standalone exact compiler retains all admissible
integer defects, each isolated-sign sector, and every vertex of
`{w in [0,1]^3:v dot w=0}`. Permutations are absorbed by permuting speed
coordinates; global negation exchanges signs and defects. There are no
other sign sectors with a positive speed solution.

For fixed error slice the width is a maximum of linear functions minus a
minimum, hence convex in `w`. Its sum has its maximum at a normalized speed
vertex. This compiler uses exact fractional knapsack on the cube slice,
reusing the independently audited width algorithm from our THM-4434 audit;
it imports no repository mathematics. All 747 maxima obey

```text
F_v(w)<=15c/98.                                        (6)
```

The sole normalized boundary maximum is pattern `(1,7,8)`. The semantic
digest of every exact `(pattern, defects, slope, intercept)` record is

```text
cf808062354debbefc1d8ead8ad0d10e9da5427cb42b8f083b6af24d0059c87c.
```

The excluded `(1,2,2)` pattern has just defect zero, exact maximal slope
`3/14`, and intercept `4/3`. Therefore

```text
N<3c/14+4/3 <14c/55 for c>=34.                         (7)
```

This is a count closure of that one coefficient pattern, not an assertion
that every norm-five relation is exceptional.

## 3. A short relation reduces the general head to 63

The projected relation lattice is

```text
{(x,y) in Z^2:ax+by=0 mod c}, with determinant c.
```

The projected `l1` ball has area

```text
[2c(ab+ac+bc)/((a+b)(a+c)(b+c))] L0^2 > (3/4)L0^2.
```

The elementary planar Minkowski argument therefore supplies a primitive
nonzero relation with

```text
S<4sqrt(c/3).                                         (8)
```

This step also requires no parity. It is important that we do **not** use
the even-norm alternative from the odd proof.

After discharging an empty defect list, outside the three separated patterns,
(2), (3), and (6) imply the desired
strict bound whenever

```text
c >= A S+B,
A=1540/547, B=21560/1641.                              (9)
```

Indeed `14/55-15/98=547/5390`. For every integer `c>=64`, (8) pays (9):

* If `S<=18`, then `A S+B<=104720/1641<64`.
* If `S>=19`, the function `g(S)=3S^2/16-A S-B` is already positive at 19,
  where `g(19)=27763/26256`, and has positive derivative thereafter.
  Equation (8) gives `c>3S^2/16>A S+B`.

Thus every triple whose chosen short relation belongs to the general
sector has `mu(F_w)<6/55` at `c>=64`. Equations (7) and the next section
cover every possible exceptional choice at these heights.

## 4. The two physical short-relation branches

For either primitive norm-three or norm-four relation, strict eligibility
implies `|v dot n|<r S<1`, so the defect is exactly zero. Therefore the
complete carrier dictionary is on `C=k v`, with `3` not dividing `k` and
all three strict roofs retained.

Put `f(t)=max(0,L_w(t v))`, extending it evenly. This is a nonnegative even
function decreasing on the positive half-line. The precise physical sum is

```text
mu(F_w)=sum_(k in Z) f(k)-sum_(k in Z) f(3k).
```

If `I=integral_R f`, rectangle comparison at steps one and three gives

```text
|mu(F_w)-(2/3)I|<=2f(0)=6/(7c).                        (10)
```

For `(1,1,1)`, positive sorted speeds force the additive identity `c=a+b`.
The complete roof calculation in the linked additive audit gives exactly
`I=27/196`. Consequently

```text
mu(F_w)<=9/98+6/(7c)<6/55 for c>=50.                   (11)
```

For `(1,1,2)`, retain the common cap `H=3/(7c)` and the roof whose relation
coefficient has absolute value two. If its complementary speeds are `p,q`,
that roof is `[3(p+q)-28t]/(14pq)`. Since `p,q<=c`, its intercept is at least
`H`. The cap-plus-one-roof envelope has half-integral

```text
H * ([3(p+q)/(14pq)]-H/2)/(2/(pq))
 =9[c(p+q)-pq]/(196c^2)
 <=9/196,
```

where the last inequality is exactly `(c-p)(c-q)>=0`. The actual physical
function lies below this envelope because it also retains the other roofs.
Thus `I<=9/98`, and (10) gives

```text
mu(F_w)<=3/49+6/(7c)<6/55 for c>=18.                   (12)
```

This bound was independently derived by the niche agent. Root also derived
the exact equality `I=9/98` from the saturated error-plane Jacobian; the
simpler one-roof upper bound above suffices for (12). Neither derivation
requires the older odd-only norm-four theorem. The known mixed-parity
norm-four hostile `(2,11,20)` is preserved rather than incorrectly removed.

## 5. Complete physical head and equality boundary

The remaining universe is precisely

```text
1<=a<b<c<=63,
every speed nonzero modulo three,
gcd(a,b,c)=1.
```

The standalone verifier enumerates every such speed triple. For each one,
it scans the first and third carrier coordinates inside their exact strict
integer roof bounds and solves the middle coordinate. It retains only the
complete owner residue and plane constraints. On denominator `42abc`, the
physical summand is the minimum of the three integer terms

```text
min(18ab, 3w_i[3(w_j+w_k)-14|C_i|]).
```

The output is:

```text
eligible triples                                  10,074
physical masses exceeding 6/55                         0
equality triples                           {(1,10,11)}
physical masses exceeding the old 6/77 ceiling          151
```

The next-largest head value is `383/3640` at `(5,8,13)`. Hostile/positive
controls include `(2,5,7):22/245`, `(1,7,8):31/392`, `(1,5,11):6/77`,
`(1,10,11):6/55`, and `(2,11,20):11/140`.

Every tail inequality above is strict. The audited analytic chain and this
exact finite head yield the unique equality in (1). The wildcard referee
independently reconstructed the complete physical head, all 747 coefficient
slopes by explicit cube-edge vertices, and all 192 owner residue fibers.
Its empty-defect correction is incorporated above. Root independently
generalized the literal native six-sheet engine to all parities and compared
all three projection columns, physical mass, and the factor-three contact
count on every one of the 10,074 rows: **PASS, 181,385 gates**. The native
head has 1,028 empty rows and 117,864 positive contacts. It reproduces the
sole equality and all 151 failures of the old ceiling. No huge height census
is a premise of this proof.

A further immediate corollary would remove primitivity: a common divisor
of ternary-unit speeds is also a ternary unit. Multiplication by that divisor
preserves circle measure and permutes the three sheet labels. Thus equality
would be exactly positive dilates of `(1,10,11)` by integers prime to three.
This follows from the actual physical definition, not from a gcd operation
that ignores sheet residues.

## 6. Upgrade to the selected network projection

Define `E_i` by summing the individual capped `i`-roof over the same complete
dictionary. The exact identity with the degree-zero sheet-network capacity
is parity-free: a fixed speed sheet has interval length `1/(7w_i)` and gaps
six times that length. Intersections preserve this sparsity, so the sharp
sparse-interval lemma in
[THM-4414 / six-separated contact capacity collapse](../../01-canon/theorems/THM-4414-lrc14-six-separated-contact-capacity-collapse.md)
identifies every max-flow with the sum of its edgewise cap minima. The raw
change of variables has three inverse x-branches with length factor one
third. Owner orientation classes remain separate. These arguments use the
ternary-unit sheet labels, not oddness. In particular `E_i<=qN` for each i.

The native cut at zero also preserves this roof formula without parity.
If the zero-label half interval of speed d meets a nonzero-label interval
of speed s, its latter left endpoint is at least `11/(42s)`, while the
half interval ends at `1/(14d)`. Thus `s>11d/3`, so the entire s-interval
is shorter than the half interval. When the wrapped interval is omitted,
the other pair therefore already lies below the half capacity; replacing
the full cap with the half cap changes nothing. When the wrapped interval
belongs to the pair, its nonwrapping partner lies on only one side of the
cut, so clipping changes no intersection length. The endpoint near one is
symmetric. The native three-pointer scan enumerates all six distinct label
assignments and advances every interval tied at the least right endpoint;
this enumerates every positive contact once, including mixed parity. The
all-parity speed loops must begin early enough to include `(1,2,4)`.

The general and `(1,2,2)` count branches therefore already bound all three
network projections strictly in their stated tails. The incoming
[audited additive-family theorem](lrc14_additive_parity_empty_core_sep06.md)
supplies `min_i E_i<=6/55` for the entire `(1,1,1)` branch, with equality
only at `(1,10,11)`. This is the one additional proved dependency of the
network upgrade; the independently derived physical theorem did not need it.

For norm four, there is a single fixed projection equal to physical mass.
Let j be the coordinate with coefficient magnitude two. In every sorted
positive sector j is a or b, so its complementary speeds are `p,c` for some
`p<c`. Write

```text
T(t)=max(0,min(H,[3(p+c)-28t]/(14pc))),
H=3/(7c), k0=3(c-p)/28, k1=3(c+p)/28.
```

The other two roofs lie above H on `[0,k0]` and above the doubled-coefficient
roof on `[k0,k1]`. This follows by checking the two affine endpoints:

| Sorted identity | p | k0 | k1 | Other roofs at k0 | Other roofs at k1 |
|---|---:|---:|---:|---|---|
| `c=2a+b` | b | `3a/14` | `3(a+b)/14` | both `3/(14a)>=H` | `H,0` |
| `c=a+2b` | a | `3b/14` | `3(a+b)/14` | both `3/(14b)>=H` | `H,0` |
| `2b=a+c` | a | `3(b-a)/14` | `3b/14` | `H,3/(7b)` | both `3/(14b)` |

All roofs decrease, so the k0 inequalities also pay the preceding plateau;
their differences are affine on the sloping segment. Beyond k1 the doubled
roof is nonpositive. Thus the **entire** physical profile is exactly T,
not merely a function with the same integral. Every live integer multiplier
has the same fixed roof j as its minimum, and consequently

```text
E_j(w)=mu(F_w)=min_i E_i(w)                             (13)
```

throughout the norm-four sector. The strict physical tail (12) therefore
also closes its selected network projection. This whole-line equality was
independently derived by the niche agent, and the endpoint table above was
checked directly again for this report.

The updated complete H63 head checks all three E values as well as mass.
Its selected maximum is also exactly `6/55`, uniquely at `(1,10,11)`.
The network upper bound in (1) follows. We retain the physical
distinction: at `(4,7,11)`, `mu=215/2156` while `min_i E_i=223/2156`, so a
general physical identity cannot replace the additive projection theorem.

## 7. Connection contract and verification

The source is the literal six-sheet union. The first map sends it to the
complete integer error address and preserves physical interval length and
owner residues. The second map slices by a primitive relation; width and
residue density bound the carrier count. That map loses individual roof
values, so the norm-three and norm-four sections retain the physical roofs
as a sidecar. The coefficient compiler removes parity honestly rather than
reusing the old filtered output. The shortest relation supplies the scale;
the enlarged physical threshold then collapses the exact head to 63.

The cheapest decisive hostile probes were the two omitted short patterns:
both `(1,1,1)` and `(1,1,2)` have count slope `2/7`, above the proposed
physical count gate `14/55`. Their physical central-section areas are small
enough to restore the ceiling. The extra `(1,2,2)` slope `3/14` needs its exact
short intercept, not a new physical classification.

Reproduction after filing in the repository:

```text
python -B 04-computation/overnight4_20260906_lrc_parityfree_probe.py
python -B -O 04-computation/overnight4_20260906_lrc_parityfree_probe.py
python -B 04-computation/overnight4_20260906_lrc_parityfree_probe.py --head-tsv 05-knowledge/results/overnight4_20260906_lrc_parityfree_probe_head63.tsv
```

The verifier executes **229,857** optimization-live gates. Normal and
optimized outputs are byte-identical. The
[complete per-row TSV](overnight4_20260906_lrc_parityfree_probe_head63.tsv)
is retained for an independent literal-sheet comparison.

```text
SHA-256, LF bytes:
source 1405ce106b74d533c3f8e330fc95855996594913e031e5f80ea323092445d866
output 9f0f60138b6b9716b545501ce297076cab916ea94e8ec8b52b9c2d867a97da69
head63 TSV c3d33fdd136245aafe512b04963a6eb6f1b5db6f1a572a3e8535ef59d01a09fa
```

The producer transcript retains its pre-audit candidate wording as historical
run provenance; the independent analytic, coefficient, physical-head, and
native-network audits above establish the current proved status. The niche
agent's independent whole-line norm-four roof proof and root's saturated
plane-area derivation supply separate checks of the short-relation mechanism.
The author of this lane did not reserve a theorem ID or mutate Git. Root's
checkpoint integrates the independently audited proof and the later incoming
promotions described in Section 1. THM-4437 is not a dependency of the
original proof; it is a proved input to the later combined corollary.
