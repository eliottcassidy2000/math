# Session73: the guard-support equality fails, not the sharp-depth conjecture

**Status:** `REFUTED` — the support-rank equality in HYP-9419 clause 2, at
incoming commit `64a7985e7054821c5958dc85f6dd489c844111f2`.
`FINITE-EXACT` — the small rational ranks and witnesses below.
`PROVED` — the support sandwich and signed-cylinder replacement formula.
The separate Session72 sharp-depth conjecture remains **OPEN**, with order
seven its first untested order. No order-seven computation was attempted.

## 1. Exact target, scope, and provenance

The pinned hypothesis says:

> The support form proposed here: the rank of the span of `k`-fold products of
> named guards on `X_n` equals the dimension of the span of the monomials `x_S`
> over the free-pair sets `S` that are unions of at most `k` path supports of
> states, so that it is computable from Hamiltonian-path supports without
> elimination; it would explain `8, 44, 273`.

Source: `05-knowledge/hypotheses/HYP-9419-unifying-candidates-three-layer-composition-rank-laws-dissection-obstructions-pattern-clocks.md`,
clause 2, at the pin above. `SOURCE_PINS.json` records Git blob IDs and SHA-256
hashes of byte-exact snapshots, including the operating contract.

This is an audit of that new support form, **not** a refutation of the cited
Session72 assertion that the least full-recovery depth is
`ceil((n-2)/2)`. The closest inherited exact result already warns against the
proposed shortcut: Session72's
`05-knowledge/results/session72_even_graph_unifying_carriers_20261026/session72-response/REPORT.md`,
section 7, lines 458–468 at this pin, states that order five has 46 visible
linear-forest Walsh characters but rank 44, hence two further dependencies.
This packet credits that prior result and independently translates the
obstruction into ordinary Boolean monomial coordinates. The numbers 8 and
44 are inherited and independently reverified, not newly discovered here.
The order-six rank 273 is not recomputed or used as a dependency.

**Universe and gauge.** Vertices are `0,...,n-1`. The arcs
`0 -> 1 -> ... -> n-1` are fixed. The free edges are
`E_n={(i,j): 0<=i<j<n, j-i>=2}`, in lexicographic order; `x_ij=1` means
`i -> j`. Thus `X_n` is the entire cube `{0,1}^{M}`, with
`M=(n-1)(n-2)/2`. Every permutation `p` has guard `g_p` indicating that its
consecutive arcs are directed as required. A reversed forced edge makes the
**whole guard zero**; it is not merely an unconstrained adjacency to discard.
A possible path is one with nonzero guard, equivalently one occurring in at
least one state. Its full support is the set of its free consecutive pairs,
forgetting orientation.

All monomials are multilinear functions on this full cube, where distinct
`x_S` are linearly independent. The identity-path guard is 1, so products of
at most `k>=1` guards and exactly `k` guards span the same space.

## 2. Exact census and minimality

Let `U_k` be the distinct unions of at most `k` full possible-path supports,
and let `D_k={S:S subset T for some T in U_k}`. The proposed literal monomial
span has dimension `|U_k|`. The generous downward-closure reading has dimension
`|D_k|`. Enumeration includes every labelled state and every permutation:

| n | k | states | full unions `|U_k|` | downward closure `|D_k|` | actual Q guard-product rank |
|---|---|---|---|---|---|
| 3 | 1 | 2 | 2 | 2 | 2 |
| 3 | 2 | 2 | 2 | 2 | 2 |
| 4 | 1 | 8 | 6 | 8 | 8 |
| 4 | 2 | 8 | 6 | 8 | 8 |
| 5 | 1 | 64 | 31 | 46 | 44 |
| 5 | 2 | 64 | 53 | 64 | 64 |

At these tested depths every union has a compatible signed-path
representative; the general replacement below nevertheless retains
compatibility explicitly. At depth one the union of **actual nonzero
monomial coefficients** is exactly the downward closure in all three orders.
Thus interpreting “support” as polynomial-expansion support does not rescue
the order-five equality either.

There are respectively 3, 11, 53 possible permutations and 2, 9, 46 distinct
nonzero guards at orders 3, 4, 5. Including the zero guard gives 3, 10, 47
distinct guards; in particular “47 guards at five” is consistent with the
46 nonzero rows here and is not the 46-dimensional candidate monomial span.

**Minimal orders within `n>=3`.** The literal equality first fails at
`n=4,k=1`: 6 is not 8. The downward-closure repair first fails at
`n=5,k=1`: 46 is not 44. At order three the depth-one span and support family
are already full, so every larger depth remains full. At order four the
depth-one guard span and downward closure are full, so no larger depth can
produce a counterexample to the downward-closure reading there. These facts
establish the asserted minimal orders, not merely the first sampled orders.

**Filter hostile.** If impossible permutations are incorrectly included as
paths “of states,” their full support counts become 8 and 38 at orders four
and five. That error hides the order-four witness but still fails at order
five (38 versus 44); their downward closures remain 8 and 46. The checker
records this control rather than silently changing the path universe.

## 3. Small witnesses and the lost information

### Literal reading: a missing support, yet a recoverable monomial

At order four write `a=x_02`, `b=x_03`, `c=x_13`. Full support masks, with bit
order `(a,b,c)`, are `{0,2,3,5,6,7}`: neither `{a}` nor `{c}` is present.
Nevertheless

```
g_(2,0,1,3) = (1-a)c,
g_(3,1,2,0) = (1-c)(1-a),
a = 1 - g_(2,0,1,3) - g_(3,1,2,0).
```

Hence `a` lies in the guard span but not the literal proposed monomial span.
The signed factors `(1-x)` create lower-support terms that the full-support
list discards. The eight-dimensional rank is independently certified by an
8-by-8 integer evaluation minor of determinant -1.

### Downward reading: an occurring coefficient cannot be isolated

At order five use `(a,b,c,d,e,f)=(x_02,x_03,x_04,x_13,x_14,x_24)`.
The guard of `(0,1,3,4,2)` is `d(1-f)`, so the monomial `d` certainly appears
in an actual guard expansion and in the downward-closure span.

For a multilinear polynomial `q`, let `[m]q` denote its coefficient of `m`.
Define the linear functional

```
L(q) = [d]q + [ad]q + [abd]q + [cd]q + [acd]q
     + [df]q + [adf]q + [abdf]q + [cdf]q + [def]q + [adef]q.
```

Exact exhaustive expansion verifies `L(g_p)=0` for every one of the 53
possible paths, whereas `L(d)=1`. Thus `d` is **not** in the guard span,
even though it occurs in a guard. This separates the two spaces directly;
it is not only a mismatch of their dimensions.

`check.normal.json` gives a second independent annihilator in bitmask
coordinates and verifies both on all 53 rows. They imply the upper bound
`46-2=44`. A displayed 44-row/44-column minor has determinant 2, proving the
matching rational lower bound independently of the rank routines. Two
elimination routes (integer primitive rows and exact `Fraction` arithmetic
with reverse column order) also give 44, both in evaluations and coefficients.
The missing coordinate is **signed coefficient coupling**, not a missing
monomial support. Session72's Walsh-coordinate warning is the same obstruction
in a different invertible polynomial basis.

## 4. Strongest support-only survivor and exact replacement

The following elementary statement applies to every `n,k`, over **any field**
when using the `0/1` Boolean coordinates. It is not an all-order rank formula
without elimination.

For a possible guard write `A_p` for its positively required free coordinates
and `B_p` for its negatively required coordinates. Thus

```
g_p = product_(e in A_p) x_e * product_(e in B_p) (1-x_e).
```

For a tuple of at most `k` guards put `A=union A_p`, `B=union B_p`. If
`A intersect B` is nonempty the product is zero. Otherwise it is the signed
cylinder `x_A product_(e in B)(1-x_e)`, whose monomial coefficient vector is

```
C_(A,B),S = (-1)^|S\A|  if A subset S subset A union B,
           0            otherwise.
```

Let `U_k^compatible` consist of the distinct supports `A union B` of these
nonzero products, and let `D_k^compatible` be its downward closure. Then

```
|U_k^compatible| <= rank_F(span of guard products)
                 <= |D_k^compatible| <= |D_k|,
rank_F(span of guard products) = rank_F(C).
```

**Proof.** The displayed expansion gives both the upper inclusion into the
monomial span and the exact coefficient-matrix replacement; Boolean monomial
evaluation is invertible over every field. For the lower bound, select one
nonzero cylinder for each distinct support `T`. Its coefficient at `x_T` is
`+1` or `-1`, and every other monomial has support contained in `T`. Order the
selected support columns and rows by increasing support cardinality. Their
square submatrix is triangular with diagonal entries `+1` or `-1`, hence has
determinant `+1` or `-1` over the integers and is invertible over every field.
This proves the sandwich. The checker additionally verifies unimodular
support minors for the three depth-one cases. QED.

The order-four literal candidate mistakes this generally strict lower bound
for equality. The order-five downward candidate mistakes an upper bound for
equality. Signed orientations and compatible joint constraints must survive
in the replacement matrix. Its rank is not supplied by merely counting its
nonzero columns.

The inherited degree obstruction also survives: each guard has degree at
most `n-1`, so `k` guards have degree at most `k(n-1)`, and full recovery needs
`k>=ceil((n-2)/2)`. None of this addresses whether three restored tests suffice
at order seven. That is precisely the separate open sharp-depth question.

## 5. Reproduction, dependencies, and stop condition

Run in this directory:

```bash
python3 check.py --verify-manifest > /tmp/guard-support.replay.normal.json
python3 -O check.py --verify-manifest > /tmp/guard-support.replay.optimized.json
cmp /tmp/guard-support.replay.normal.json /tmp/guard-support.replay.optimized.json
cmp check.normal.json /tmp/guard-support.replay.normal.json
```

Only the Python standard library is used. All test failures raise explicitly;
no validation depends on `assert`. Guards are independently built by direct
arc evaluation and signed cylinders, and polynomial coefficients are checked
by Boolean Möbius inversion. Rational and primitive-integer ranks agree;
Bareiss determinants separately check the lower-rank certificates. Source
snapshots are immutable at the incoming pin; no shared repository file was
edited by this audit worker.

This bounded audit stops at the minimal witnesses, the exact repaired
statement, and replayable evidence. No unrequested order-seven computation,
new broad classification, or promotion of the sharp-depth conjecture is part
of this packet. The new question is whether coefficient dependencies admit a
structural description that is genuinely finer than the support family.
