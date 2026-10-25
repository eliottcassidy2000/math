# Higher clone cancellation: minimum cycle count versus invariant toggle rank

2026-10-25. **PROVED (elementary, independently audited):** the exact
cancellation order, equality criterion, Boolean operation sidecar, real
exterior-algebra interpretation, and an unbounded family of genuine tournament
examples. **FINITE-EXACT:** the explicitly listed actions and tournaments below.
No historical novelty, faster enumeration, black-path transport, or LRC claim.

## 0. Inheritance, concepts, and scope

The inherited result is
`05-knowledge/results/ports_regular_actions_hn_20261025/counting/CLONE_INVOLUTION.md`
(at the source checkpoint `bd128a`). It identifies a fixed tournament's
weighted inflation fibre with orbits of cut subsets in its canonical module
gaps. Its single invariant odd-sized gap toggle proves cancellation at `t=-1`.
The module-path input is **THM-4463**, file
`01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md`,
section 2. Both inputs were read; this note does not re-prove their full
inflation-fibre classification.

The corrected near miss is toggling one *labelled* gap, which need not descend
to the cut-set quotient. The canonical inherited hostile is
`C3[TT2,TT2,TT2]`; the new hostile is the regular `C3^2` gap action. The
least-used sidecar is the full permutation action on gaps, not just the number
of module orbits. The local portfolio is: Anchor, higher signed cancellation;
Niche, exact equivariant operations; Wildcard, an exterior-invariant model.

| Concept | Object / operation | Retained or missing coordinate |
|---|---|---|
| Inflation fibre | weighted quotient -> cuts modulo automorphisms | retain the fixed inflated tournament and literal gaps |
| Signed enumerator | subset orbit rank polynomial at `-1` | full one-element cycle profile controls root order |
| Boolean toggles | invariant symmetric differences | only common invariant subsets are legal translations |
| Exterior invariants | permutation representation -> exterior algebra | preserves graded dimensions, not subset multiplication |
| Tournament realization | observers encode coordinates | actual intrinsic arcs and all automorphisms must be verified |

The arguments first concern an arbitrary finite odd-order action. No abstract
action is silently declared to be a tournament action. Section 5 constructs
actual tournaments and proves their full automorphism groups.

A bounded source recovery search covered the inherited counting packet,
THM-4463, relevant THM-47xx files, and exact phrases about subset orbits,
minimum cycles, and vanishing order. It recovered the inherited cancellation
and standard exterior Molien framework, but no already-filed version of the
exact refinement in the inspected scope. This is not a completeness or novelty
claim. All refinement proofs are included below.

## 1. Exact root-order theorem for odd actions

Let a finite group `A` of odd order act on a finite set `G`. Write `n=|G|`,
`Omega=2^G/A`, and

```
P_A(t) = sum_[C] in Omega t^|C|,
c(g) = number of cycles of g on G,
d = min_(g in A) c(g),
r = number of A-orbits on G.
```

Empty-set conventions are `P=1`, `d=r=0`, empty products `1`. Faithfulness is
not required. Burnside applied separately at each subset size gives

```
P_A(t) = (1/|A|) sum_(g in A) product_(D a cycle of g) (1+t^|D|).      (1)
```

Every cycle has odd length. For odd `ell`,

```
1+(-1+x)^ell = ell*x + terms of degree >=2.
```

Thus the summand for `g` has exact order `c(g)` and positive leading
coefficient `product_D |D|`. The minimum-order summands cannot cancel:

```
ord_(t=-1) P_A(t) = d,                                             (2)
kappa := [x^d] P_A(-1+x)
       = (1/|A|) sum_(c(g)=d) product_D |D| > 0.                    (3)
```

In particular `kappa` is a positive integer, since `P_A` has integer
coefficients. This is a theorem for all finite odd-order actions, not an
extrapolation from finite checks.

Each `g`-cycle lies in an `A`-orbit. Hence

```
r <= d <= n,
d = r iff some single g in A is one cycle on every A-orbit of G.    (4)
```

The word *single* is essential: finding a suitable, possibly different,
element for every orbit separately does not suffice (section 6). All cycle
lengths and all `A`-orbit sizes are odd, so

```
d == r == n (mod 2).                                              (5)
```

Complement gives `P_A(t)=t^n P_A(1/t)`. Since `(1+t)^d` divides `P_A(t)`
over the integers, `2^d` divides `|Omega|=P_A(1)`. This stronger integer
divisibility does **not** construct `d` invariant symmetric-difference toggles.

For a nontrivial transitive action (`r=1`, `n>1`), Burnside's fixed-point
average is one. The identity has `n>1` fixed points, so some element has no
fixed points. All of its cycles have odd lengths at least three. Therefore

```
d <= n/3.                                                        (6)
```

The regular actions of `C3^k` attain (6), including the genuinely realized
tournament family in section 5. Thus their separation `d-r` is unbounded.
No simultaneous derangement on several point orbits is assumed in (6).

### Clone signed moments

For a fixed tournament `T`, let `q` be its number of canonical transitive
modules, `N=|V(T)|`, and `G` their gaps, so `n=N-q`. The inherited fibre map
has `M=q+|C|`, excess `E=n-|C|`, and sign `(-1)^E`.
Using (2)--(3) and (5),

```
sum_[C] (-1)^(n-|C|) binom(|C|,j) = 0       for 0<=j<d,
sum_[C] (-1)^(n-|C|) binom(|C|,d) = kappa.                        (7)
```

Indeed the `x^j` coefficient of `P(-1+x)` is
`sum_[C] binom(|C|,j)(-1)^(|C|-j)`; the conversion factor in (7) is
`(-1)^(n+j)`, which equals `+1` at `j=d`. Equation (7), rather than just
its zeroth moment, is the higher signed consequence in the actual fibre.

## 2. The maximal invariant translation sidecar

Let `O_1,...,O_r` be the `A`-orbits of gaps, with odd sizes `n_i`.
For each invariant subset `F` of `G`, define `tau_F(C)=C symmetric_difference F`.

- `tau_F` commutes with `A` iff `F` is invariant. Necessity follows by
  applying commutation to the empty set; sufficiency is immediate.
- More strongly, a fixed `tau_F` descends to `Omega` iff `F` is invariant,
  without assuming equivariance in advance. If it descends, the equivalent
  inputs `F` and `gF` must have equivalent outputs `empty` and
  `F symmetric_difference gF`. Only the empty set belongs to its own
  `A`-orbit, so `gF=F`. Sufficiency follows from equivariance.
- The invariant subsets are exactly unions of the `O_i`. Thus all such
  translations form a group `B=C2^r`; its coordinate generators reverse
  subset-cardinality parity.
- This action is free even on `Omega=2^G/A`, not merely on labelled subsets.
  If `[C symmetric_difference F]=[C]`, then `gC=C symmetric_difference F`
  for some `g`. Since `gF=F`, iterating through the odd order of `g` gives
  `C=C symmetric_difference F`, so `F` is empty.

This is the **exact rank of fixed symmetric-difference operations that
descend to the quotient**, equivalently of the equivariant ones.
It is not a classification of all possible nonlinear permutations of `Omega`.
The action preserves the literal inflated tournament, so in the clone setting
it operates on precisely the inherited weighted quotient fibre.

For a representative `C` of a `B`-orbit in `Omega`, put
`a_i=|C intersect O_i|`. The polynomial of that entire Boolean orbit is

```
product_(i=1..r) (t^a_i+t^(n_i-a_i)).                              (8)
```

Each factor has a simple zero at `-1`: its derivative there is
`(-1)^a_i (n_i-2a_i)`, which is nonzero since `n_i` is odd. Thus every
Boolean orbit separately has exact cancellation order `r`. Higher order
`d>r` necessarily arises by cancellation **between** these Boolean orbits.
That is why the signed multiplicity need not equal the legal toggle rank.

## 3. The exterior-algebra bridge and its exact boundary

**PROVED generalization, real representations; odd group order NOT needed.**
Let a finite group `A` act linearly on a finite-dimensional **real** vector
space `V`, and define

```
H_V(t) = sum_j dim_R((exterior^j V)^A) t^j.
```

The standard exterior Molien identity is

```
H_V(t) = (1/|A|) sum_(g in A) det(I+t g).                           (9)
```

For a primary-source reference to the standard identity, take `x=0` in
section 2 of Daryl R. DeFord and Peter G. Doyle,
[Cyclic groups with the same Hodge series, Rev. Un. Mat. Argentina 59 (2018), 241--254](https://inmabb.criba.edu.ar/revuma/pdf/v59n2/v59n2a02.pdf).
The exact root-order argument below is proved here and is not claimed as a
new historical theorem or attributed to that paper.

For completeness, (9) follows by expanding the determinant as the trace
polynomial of exterior powers and averaging the finite group action, which
is a projection onto invariants. Average an inner product to make the action
orthogonal. For each `g`, let `f(g)=dim V^g`. The factor `det(I+t g)` has
exactly `f(g)` factors `1+t`. Every other real eigenvalue is `-1`, contributing
`2` at `t=-1`; every nonreal pair `lambda,conjugate(lambda)` contributes
`|1-lambda|^2>0`. Consequently its remaining factor is positive at `-1`.
As in section 1,

```
ord_(t=-1) H_V(t) = min_g dim V^g.                                (10)
```

Write `s=dim V^A` and `W=(V^A)^perp`. As graded algebras,

```
(exterior V)^A = exterior(V^A) tensor (exterior W)^A,
H_V(t) = (1+t)^s H_W(t).                                         (11)
```

Thus `H_V/(1+t)^s` has nonnegative integer coefficients. The common fixed
space controls the natural exterior factor; the smallest *individual*
fixed space controls the root multiplicity. Those are distinct predicates.

**Bridge to odd subset actions.** Set `V=R^G` with its permutation action.
For odd-order `A`, every element of a subset stabilizer induces only odd
cycles on that subset, hence has sign `+1`. Averaging an ordered wedge for a
subset therefore supplies one nonzero invariant line per subset orbit.
Different subset orbits have disjoint wedge supports. Hence

```
P_A(t) = H_(R^G)(t),     dim V^A=r,     dim V^g=c(g).               (12)
```

The map preserves graded dimensions; it does not identify subset union or
complement with exterior multiplication. The subset stabilizer's orientation
character is the sidecar that (12) requires and that even actions may lose.
Equivalently, a permutation cycle of length `ell` contributes
`1-(-t)^ell` to (9), equal to `1+t^ell` exactly for odd `ell`.

Two sharp boundary controls:

1. A transposition on two points has subset-orbit polynomial `1+t+t^2`,
   root order zero, although its minimum cycle count is one. Its exterior
   invariant polynomial is instead `1+t`, whose root order is one as (10)
   requires. Do not extend the subset statement to arbitrary even actions.
2. The **real** hypothesis in (10) cannot simply be erased. Let `C3` act
   on `C^3` by the scalar primitive cube root `omega`. Nonidentity elements
   have no fixed vector, but the exterior-invariant polynomial is `1+t^3`,
   with a simple zero at `-1`. Its minimal fixed dimension is zero. The
   conjugate-pair positivity used above is absent.

## 4. A strict action example and the negative quotient

Let `A=C3^2` act regularly on its nine elements. The identity has nine
one-cycles; each of the eight nonidentity elements has three three-cycles.
Thus

```
P(t) = ((1+t)^9 + 8(1+t^3)^3)/9
     = 1+t+4t^2+12t^3+14t^4+14t^5+12t^6+4t^7+t^8+t^9.
r=1, d=3, kappa=24, |Omega|=64.                                 (13)
```

There is only one independent invariant gap toggle: complement all nine
gaps. Its 32 pairs have smaller cut counts `0,1,2,3,4` with multiplicities
`1,1,4,12,14`. Their first derivatives at `-1` sum to

```
9 - 7 + 4*5 - 12*3 + 14*1 = 0.
```

This exhibits cancellation between honest toggle pairs. The maximal-degree
root factor is not a positive exterior factor:

```
P(t)/(1+t)^3 = 1-2t+7t^2-4t^3+7t^4-2t^5+t^6.                    (14)
```

Negative coefficients rule out interpreting (14) as the Hilbert polynomial
of an additional nonnegatively graded multiplicity space. By contrast,
`P/(1+t)` has coefficients `1,0,4,8,6,8,4,0,1`, as (11) predicts.
No claim is made that some unrelated, noncanonical free `C2^3` action on a
64-element set is impossible; only the specified invariant translations and
positive degree-one exterior-factor interpretation are obstructed.

The cyclic regular action `C9` on the same nine points has `r=d=1`,
`kappa=6`, and 60 subset orbits. Orbit size and the group order alone do not
determine the higher cancellation order; element-cycle profiles are necessary.

## 5. Genuine tournament realization with unbounded separation

For each integer `k>=2`, let `X=F3^k`. Give `X` the lexicographic tournament:
for distinct `a,b`, look at the first coordinate `j` where they differ and
orient `x_a -> x_b` iff `b_j-a_j=1` in `F3`. For each `i=1,...,k`, adjoin
three observers `Y_i={y_(i,c):c in F3}`. Orient:

```
y_(i,c) -> y_(i,c+1),
Y_i -> Y_j for i<j,
x_a -> y_(i,c) iff a_i=c.                                       (15)
```

These rules specify exactly one direction for every distinct vertex pair.
This is an actual tournament `R_k`, not a tournament orientation imposed on
an otherwise symmetric quotient relation.

The outdegrees are

```
deg(x_a) = (3^k-1)/2+k,
deg(y_(i,c)) = 1+2*3^(k-1)+3(k-i).                              (16)
```

They distinguish all `Y_i` and distinguish them from `X`, since the smallest
observer degree minus the `X` degree is
`(3^(k-1)+3-2k)/2>0` for `k>=2` (it is one at `k=2` and strictly increases
thereafter). An automorphism therefore preserves every part. On `Y_i` it
is a rotation `c -> c+s_i`. The incidence in (15) forces its action on `X`
to be `a -> a+s`. Conversely every such translation preserves all the arcs,
including the first differing coordinate in the lexicographic tournament.
Therefore

```
Aut(R_k) = C3^k.                                                (17)
```

The core is pair-free:

- Distinct `X` vertices are separated by an observer of a differing coordinate.
- Two vertices in the same `Y_i` are separated by some vertex of `X`.
- Vertices from different observer parts have outdegrees differing by at
  least three; a homogeneous pair must have degrees differing by exactly one.
- An `X` vertex and a vertex in `Y_i` are separated by another observer part
  `Y_j`, available because `k>=2`: the `Y_i` vertex sees all of `Y_j` alike,
  while the `X` vertex has exactly one outneighbor there.

Inflate **every** vertex of `X` to `TT2`; leave all observers as singleton
blocks. Call the resulting tournament `T_k`. Every intended two-vertex block
is homogeneous. Any pair from distinct base blocks is separated by an
inflated vertex over a base separator, since `R_k` is pair-free. Thus the
canonical transitive modules are exactly these clone pairs and singleton
observers. Automorphisms of `T_k` must permute canonical modules and preserve
sizes; their induced maps are precisely weight-preserving automorphisms of
`R_k`, and each lifts uniquely through the transitive block orders. All
translations in (17) preserve the weights. It follows that

```
|V(T_k)| = 2*3^k+3k,
Aut(T_k)=C3^k,
G(T_k) is a regular 3^k-point C3^k-set.                           (18)
```

Consequently its *actual weighted clone fibre* has

```
P_(T_k)(t) = ((1+t)^(3^k)+(3^k-1)(1+t^3)^(3^(k-1)))/3^k,
r=1, d=3^(k-1),
kappa=(3^k-1)*3^(3^(k-1)-k).                                   (19)
```

This proves the unbounded separation within genuine tournaments and attains
the transitive bound (6). The case `k=2` is the 24-vertex realization of (13).
The restriction `k>=2` is load-bearing: the six-vertex `R_1` has three
unexpected homogeneous pairs and is not the pair-free core asserted above.

## 6. Genuine simultaneous-cycle hostile

For a second example, use `X=F3^2` but four ordered observer triples for
linear forms `a`, `b`, `a+b`, `a-b`, in that order. Keep the lexicographic
arcs on `X`, directed cycles within observers, and ordered arcs between
observer parts; replace the last rule in (15) by
`x_(a,b) -> y_(i,c)` iff the corresponding linear form equals `c`.

The `X` vertices have degree eight. The four observer parts have degrees
`16,13,10,7`, so all parts are distinguished. Rotations of the first two
observer parts force `X` translation by `(alpha,beta)`. Rotations of the
last two are then forced to be `alpha+beta` and `alpha-beta`; all translations
work. The full automorphism group is `C3^2`. The same separator arguments
prove this 21-vertex core pair-free.

Now inflate all twelve **observer** vertices to `TT2` and leave `X`
singleton. The resulting tournament has 33 vertices, full automorphism
group `C3^2`, and four orbits of three gaps. Every orbit separately admits
a full cycle. However each nonzero `(alpha,beta)` lies in exactly one of
the four distinct one-dimensional kernels, so it fixes all three points of
one orbit and is a three-cycle on each other orbit. Hence

```
r=4, d=6, kappa=24,
P(t)=((1+t)^12+8(1+t)^3(1+t^3)^3)/9,
P coefficient list = [1,4,10,28,63,96,108,96,63,28,10,4,1].       (20)
```

The first failed implication is replacing the single element in (4) by
independent choices on each orbit. The exact repair retains simultaneous
cycle profiles. This also blocks the naive sum of separate transitive
bounds: `sum_i |O_i|/3=4`, whereas the actual `d` is six.

## 7. Exact finite validation and limits

Run in any directory with **Python 3.10+**:

```bash
python3 check.py > normal.json
python3 -O check.py > optimized.json
cmp normal.json optimized.json
```

All gates are explicit exceptions, not `assert`; both modes passed and outputs
are byte-identical. The script imports only the standard library.

The finite universe is explicit: empty and trivial actions, regular `C3`,
regular `C9`, regular `C3^2`, the four-projection `C3^2` action, and the even
transposition hostile. For each odd action it checks group closure, literal
subset-orbit ranks against Burnside polynomials, shifted root coefficients,
the leading coefficient formula, parity, the simultaneous-cycle criterion,
integer divisibility, complement symmetry, signed moments, and the freeness
of every invariant Boolean translation. It also verifies (14), the even
exterior distinction, and the scalar-complex hostile.

For the constructed tournaments, a separate literal arc backtracking search
uses only outdegree and already-matched arc compatibility to enumerate **all**
automorphisms. It does not assume the proposed translation group. It also
enumerates all homogeneous pairs directly from their external incidences.

| Actual tournament | Vertices | Full automorphisms | Canonical gaps | Gap orbits | Minimum cycles | Full subset census? |
|---|---:|---:|---:|---:|---:|---|
| `T_2` | 24 | 9 | 9 | 1 | 3 | yes: 512 subsets, 64 orbits |
| four-observer example | 33 | 9 | 12 | 4 | 6 | yes: 4096 subsets, 512 orbits |
| `T_3` | 63 | 27 | 27 | 1 | 9 | **no**: structural/cycle check only |

The 63-vertex check does not enumerate its `2^27` subsets. Its asserted
root order follows from the all-action theorem and the exact gap action;
no unseen full census is reported. The code also confirms the three unwanted
pairs in `R_1`. The actual weighted-quotient/cut-orbit identification is the
inherited theorem, not a new factorial weighted-graph enumeration at order 33.

Scope stays within finite group actions, invariant exterior algebras, and the
specified weighted tournament inflation fibres. Neither extra cancellation
nor a dimension count supplies new local toggles, black moves, or transport.
