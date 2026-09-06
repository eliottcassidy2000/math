# Complete odd-prime three-node three-Hasse-jet precision and Smith law

**PROVED BY EXACT SYMBOLIC IDENTITIES + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [independent proof and computation audit](overnight7_20260906_oddjets_independent_audit.md)
passes without a mathematical correction. This closes the odd-prime counterpart of incoming
`6f450d8cb7b776b243b1031a1042a168d5757114`'s proved
[dyadic symbolic-minor law](three-node-threejet-dyadic-smith-overnight-hexagon-sep05.md).
The [arbitrary-multiplicity inverse law in THM-4443](../../01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md)
is an inherited proved dependency. No external priority is claimed.

## 1. Statement and inheritance

Let p be an odd prime, and observe Hasse derivatives of orders 0,1,2
at three distinct integer nodes on polynomials of degree<9. Write their
pairwise p-adic depths as `(e,e,e+d)`, with integers e,d>=0. Let s_1,...,s_9
be the nondecreasing p-Smith exponents and
`D_j=s_1+...+s_j`, `D_0=0`. Put `delta=[p=3]`, `eta=[p=5]`.

**Close-pair case, d>=1.** The complete determinantal law is

```
D_0=D_1=D_2=D_3=0,
D_4=min(3e,e+delta),
D_5=min(6e,4e+delta,3e+d+2delta),
D_6=min(9e,7e+d+delta+eta),
D_7=min(13e+d+delta,12e+4d+2delta+eta),
D_8=19e+4d+delta,
D_9=27e+9d.                                        (1)
```

Thus s_j=D_j-D_(j-1), and the largest loss is
`L=s_9=8e+5d-delta`. For p!=3 the extra terms 3e in D_4 and 6e in D_5
are dominated; they make the p=3 shallow mask explicit in one formula.

**Equilateral case, d=0.** Put `epsilon=[p=3 and e>0]` and
`theta=[p=5 and e>0]`. The full list is

```
(s_1,...,s_9)
 =(0,0,0,e+epsilon,2e+epsilon,4e+theta,
             5e,7e-epsilon-theta,8e-epsilon).        (2)
```

In particular all factors are units when e=0. The largest loss may also
be written `L=max(6e,8e-[p=3])`. At e>=1 the exceptional lists are

```
p=3:  (0,0,0,e+1,2e+1,4e,5e,7e-1,8e-1),
p=5:  (0,0,0,e,2e,4e+1,5e,7e-1,8e),
p>=7: (0,0,0,e,2e,4e,5e,7e,8e).                     (3)
```

Consequently the **full odd-prime partition is metric-only**, not merely
its largest factor. Combining this result with the inherited dyadic law,
all unit dependence in this precise three-node uniform-three-jet problem
is confined to the dyadic final two factors. This does not extend to more
nodes, heterogeneous multiplicities, higher jet orders or full two-jet
partitions at four nodes.

The closest proved mechanism is the incoming residual 886-minor proof.
The canonical hostile is the dyadic isometric pair with losses18/19.
The corrected near miss is to transport that unit failure to every prime,
or infer intermediate ideals from the largest loss. The least-used sidecar
is the joint gcd of two equal-weight minors: one selected minor can gain
unit-dependent valuation while its companion attains the same target.
The older [THM-4429 three-node two-jet law](../../01-canon/theorems/THM-4429-arbitrary-three-node-two-jet-smith-form-and-metric-precision.md)
supplies a related small-prime comparison, not a higher-jet dependency.
Exact-statement and correction-ledger searches retained MISTAKE-547's
observer warning and found no prior odd-prime three-jet full law in the
local truth surfaces read for this task.

The concept board is: exact inverse precision; weighted minors; small
prime coefficient divisibility; residue-complete attainment; shallow-depth
masks; and lower convex envelopes.

## 2. Normalize locally without losing units

Translate one node to zero and multiply the variable by a p-adic unit.
This is a unimodular change of polynomial coefficients and complete Hasse
data over Z_p. In both cases the nodes become

```
0, p^e, p^e a.                                     (4)
```

For d>=1 choose zero as a closest-pair endpoint, so `v_p(a)=d` and a-1
is a unit. For d=0, both a and a-1 are units. The normalized a can be any
corresponding element of Z_p, not just an integer representative. All
polynomial lower bounds below apply to that complete p-adic domain.

The three Hasse rows at zero clear an identity block. The residual 6x6
matrix has columns j=3,...,8 and rows indexed by `(node,r)`, with node
in `{1,a}` and r=0,1,2:

```
R_((node,r),j)=p^((j-r)e) binom(j,r) node^(j-r).       (5)
```

Rows 0,1,2 are the jets at 1; rows 3,4,5 are the jets at a. Every rank-r
minor has the common dilation factor

```
p^(W e) P_(I,J)(a),
W=sum J-sum_(i in I) r_i,
P_(I,J) in Z[a].                                   (6)
```

The three cleared unit factors make its residual rank-r determinantal
valuation exactly D_(r+3). The code expands **every** residual minor of
ranks1,2,3,4: 36,225,400,225 minors, **886 in total**, retaining all
**2,623** nonzero monomials after exact coefficient cancellation.
The JSON manifest contains every row/column set, weight and polynomial.
This is a finite symbolic identity universe, not a sampled node bank.

## 3. Exact largest odd-prime loss

At a node with differences u,v to the others, reciprocal expansion gives

```
a0=1/(uv)^3,
a1=-3(u+v)/(uv)^4,
a2=3(2u^2+3uv+2v^2)/(uv)^5.                        (7)
```

THM-4443 proves that the largest p-Smith exponent is the maximum of the
negative valuations of all these coefficients; every proposed denominator
is attained in an inverse matrix entry.

In the close-pair case a closest endpoint has difference depths e and
e+d. Since p is odd and d>=1, the term of depth 2e in its quadratic
numerator has coefficient 2 and is uniquely least. Thus its order-two
candidate is `8e+5d-delta`. Its order-one candidate is `7e+4d-delta`,
and order zero is `6e+3d`; both are strictly smaller. The outsider's three
candidates are at most `6e,7e-delta,8e-delta`, so it cannot dominate.
This proves the close-case L in (1), including p=3,e=0.

In the equilateral case, after removing the common square scale, the
three quadratic numerators are

```
N0(a)=2+3a+2a^2,
N1(a)=7-7a+2a^2,
Na(a)=7a^2-7a+2.                                   (8)
```

At every odd prime at least one is a unit. For p!=5, a common zero of
N0,N1 would satisfy `5(2a-1)=0`, so a=1/2; but N0(1/2)=4, nonzero at
an odd prime. At p=5 these two quadratics agree modulo five, and their
discriminant is 3, a nonsquare in F5. Thus there is no common root even
there. This argument is over the residue field F_p of Z_p; no statement
about arbitrary residue-field extensions is being inferred.

The largest order-two candidate is therefore `8e-delta`. All order-zero
candidates are 6e, and every order-one candidate is at most `7e-delta`.
For e>=1 the order-two candidate dominates; for e=0 order zero gives 0.
This proves `L=max(6e,8e-delta)` and correctly retains the ternary
order-zero mask. No odd-prime unit class analogous to the dyadic Gamma
is needed for this target.

Finally the confluent Hasse determinant has valuation `27e+9d`.
Therefore D_8=D_9-L, giving the last two entries of (1) and, for d=0,
`D_8=19e+[p=3 and e>0]`. The largest-factor identity supplies rank five
of the residual; no intermediate-minor inference is made from it.

## 4. All-depth close-pair lower bounds and witnesses

For a nonzero monomial c a^b in (6), its valuation is
`W e+b d+v_p(c)`. With d=1+d', d'>=0, encode this affine cost as
`(W,b,b+v_p(c))`. At primes p>=7 the certificate uses the universal lower
bound `v_p(c)>=0`, not a representative sampled prime. The complete
coordinatewise-undominated sets of these triples are:

| residual rank | p>=7 coefficient floor | p=3 | p=5 |
|---|---|---|---|
| 1 | `(1,0,0)` | `(1,0,1),(3,0,0)` | `(1,0,0)` |
| 2 | `(3,1,1),(4,0,0)` | `(3,1,3),(4,0,1),(6,0,0)` | `(3,1,1),(4,0,0)` |
| 3 | `(7,1,1),(9,0,0)` | `(7,1,2),(9,0,0)` | `(7,1,2),(8,1,1),(9,0,0)` |
| 4 | `(12,4,4),(13,1,1)` | `(12,4,6),(13,1,2)` | `(12,4,5),(13,1,1)` |

For every monomial, the certificate checks coordinatewise domination by
one of its row's displayed triples. This proves an unbounded lower
envelope for all e,d'>=0; specialization cancellation can only raise a
minor valuation.

The extra quinary rank-three cost 8e+d is redundant even though it survives
coordinatewise Pareto pruning. Its shifted triple `(8,1,1)` is no smaller
than the average `(8,1/2,1)` of `(9,0,0)` and `(7,1,2)`. Equivalently,

```
8e+d >= (9e+(7e+d+1))/2 >= min(9e,7e+d+1)   for d>=1.             (9)
```

This holds on the full real cone e>=0,d>=1; integrality is not necessary.
A preliminary integer case split was sufficient but less informative.
Dropping d>=1 is invalid: at e=1/2,d=1/4 the extra cost becomes strictly
smaller than the other two. The proof retains the actual domain rather
than treating coordinatewise undominated costs as automatically necessary.

The following exact minors attain the required close-case costs. Constants
and factors are shown explicitly; W is the dilation weight from (6).

| rows I | degrees J | W | polynomial |
|---|---|---:|---|
| 0 | 3 | 3 | `1` |
| 2 | 3 | 1 | `3` |
| 0,1 | 3,4 | 6 | `1` |
| 1,2 | 3,4 | 4 | `6` |
| 2,5 | 3,4 | 3 | `18a(a-1)` |
| 0,1,2 | 3,4,5 | 9 | `1` |
| 1,2,5 | 3,4,5 | 7 | `30a(a-1)(2a-1)` |
| 0,2,5 | 3,4,5 | 8 | `6a(a-1)(5a-2)` |
| 0,1,2,5 | 3,4,5,6 | 13 | `3a(a-1)(5a^2-5a+1)` |
| 1,2,4,5 | 3,4,5,6 | 12 | `90a^4(a-1)^4` |

For v_p(a)>=1, the factors a-1, 2a-1 and 5a^2-5a+1 are units at all odd
primes. The extra W=8 witness is used at p=5, where 5a-2 is also a unit.
Thus every needed cost has an actual noncancelling minor, not merely a
cheap monomial. The constants give
`v_p(3)=delta`, `v_p(6)=delta`, `v_p(18)=2delta`,
`v_p(30)=delta+eta`, `v_p(90)=2delta+eta`.
These witnesses and the complete lower certificate prove D_4,...,D_7
in (1). Every displayed factorization is independently expanded and
compared with its determinant polynomial in the companion source.

## 5. Equilateral lower bounds and joint attainment

If e=0 the determinant is a unit, so all nine factors are units. Assume
e>=1. The least possible dilation weights for residual ranks1,2,3,4 are
respectively `(1,3,7,12)`, from the smallest possible column sums and
largest possible derivative-order sums in (6).

The complete coefficient certificate gives the following lower bounds:

```
p>=7: (D4,D5,D6,D7) >= (e,3e,7e,12e),
p=5:                  (e,3e,7e+1,12e+1),
p=3:                  (e+1,3e+2,7e+2,12e+2).          (10)
```

For p>=7 and p=5 use the original 2,623 monomials of all 886 minors.
For p=3, the entire admissible equilateral residue class is `a=2+3z`,
z in Z_3. Substitute this into every polynomial and retain all **8,051**
nonzero coefficients. Let s_r=(1,3,7,12) and let beta_r be respectively
`(0,0,0,0)`, `(0,0,1,1)` or `(1,2,2,2)`. Every resulting coefficient c
satisfies the exact inequalities

```
W>=s_r,    W+v_p(c)>=s_r+beta_r.                     (11)
```

For the generic regime v_p(c) is again replaced by its lower bound zero.
Writing e=1+E, E>=0, proves
`W e+v_p(c) >= s_r e+beta_r`. Powers of integral z cannot decrease the
valuation. Thus (11) is a finite polynomial proof of (10) for the entire
unbounded domain, not a check of finitely many a or e values.

Ranks1,2,4 attain (10) with the entries/minors
`3`, `18a(a-1)` and `90a^4(a-1)^4`, respectively. Since a and a-1 are
units, their valuations are fixed in each prime regime.

For rank three, use **both** equal-weight minors

```
I=(1,2,5), J=(3,4,5): 30a(a-1)(2a-1),
I=(2,4,5), J=(3,4,5): 30a^4(a-1)(a-2),
both W=7.                                          (12)
```

At p!=3 their two linear factors cannot both be nonunits, because
`(2a-1)-2(a-2)=3` is a unit. Hence their minimum valuation is v_p(30):
zero at p>=7 and one at p=5. At p=3, a=2 mod3 makes both linear factors
divisible by three, but the same identity prevents both being divisible
by nine. Their minimum valuation is exactly one, giving total valuation
two after the factor 30. Thus the pair of minors attains D_6 in every
allowed unit class, even when one selected minor vanishes or gains
arbitrarily many digits. This joint-attainment step is the essential
equilateral repair.

Together with Section 3's D_8,D_9, equality in (10) proves the full lists
(2)-(3).

## 6. Controls, implications and boundaries

The companion source checks 432 literal full 9x9 integer Hasse Smith
matrices: 180 close-pair cases at p=3,5,7,11; 248 equilateral cases with
complete small residue-lift ranges; and four signed unit-affine controls.
It separately retains the excluded dyadic 18/19 hostile and rejects p=2
in its odd-prime predictor. At equilateral e=1 the exact lists are

```
p3: (0,0,0,2,3,4,5,6,7),
p5: (0,0,0,1,2,5,5,6,8),
p7: (0,0,0,1,2,4,5,7,8).
```

The p=5 exception changes intermediate factors while its largest loss
remains the generic 8e. Thus proving only the precision law would miss
a real full-partition correction. At e=0 in the close case, the formula
gives six unit factors followed by `(d+1,3d,5d-1)` at p=3 and
`(d,3d,5d)` at other odd primes, matching the separated two-node block.

A shallow affine fit is still unsafe: at p=7,e=4,d=1, D_7=52, while the
single candidate 13e+d gives 53. The second rank-four slope is needed.
This is an odd-prime instance of the incoming dyadic failure mechanism,
not a new discovery of the general warning against finite-depth fits.

The source-target map normalizes the integer observer to its local
residual matrix and retains exact polynomial coefficients, dilation
weights and close-pair powers. Monomial valuation discards cancellation,
so actual factored witnesses and, in the equilateral rank-three case,
their joint residue obstruction restore attainment. The final metric
quotient is valid here because those witnesses cover every admissible
unit class. The determinant remains insufficient by itself, and the
metric quotient is not transported to a different observer or node count.

This joins the prior three-node two-jet metric law, the incoming dyadic
three-jet unit boundary, and the sixth-round order-aware jet envelope:
the general normalized jet carries units, but in this odd-prime uniform
three-node problem the complete minor families erase them at the target.
No terminal-pair-only cancellation budget at arbitrary node count is
inferred; the sixth-round unbounded local cancellation family remains a
separate valid obstruction.

## 7. Reproduction and frozen proof manifest

```
python -B 04-computation/overnight7_20260906_oddjets.py --certificate 05-knowledge/results/overnight7_20260906_oddjets.json
python -B -O 04-computation/overnight7_20260906_oddjets.py --certificate 05-knowledge/results/overnight7_20260906_oddjets.json
```

The script writes `overnight7_20260906_oddjets.json` beside itself; an
explicit `--certificate PATH` may redirect it. The manifest contains the
complete 886 polynomial minors, all close-pair frontiers and 11 exact
factorization witnesses. No repository code is imported. Normal and
optimized runs produce identical LF outputs and identical JSON, with
**18,102 optimization-live gates PASS**. The full symbolic universe and
attainment arguments prove the unbounded result; literal matrices are
independent controls, not a finite extrapolation premise.

```
source b894f4860c1bdd824d4b2ce5ef2b443bdff2b73a684a06d5ac699c3afafa45b1
output 1c2b6d20a5a58d8d5a2c7d7565c8a058c3bac54b90d9c0a407785bbc351f5e85
JSON   05fc22845ec10078a2c9d568dea59d00bfdbdc8588d9e314e4134fcbe13e7d11
```

The independent referee reconstructs every polynomial by Bareiss integer
determinants and degree-bounded interpolation, using no producer imports
or permutation expansion. All 886 polynomials and 2,623 coefficients
match the manifest. It separately verifies the 8,051 shifted ternary
coefficient floors, all frontiers and 11 witnesses, then compares the
complete formulas with 177 fresh odd-prime matrices using a standard-library
p-local pivot Smith reduction, plus both dyadic hostiles. Normal and
optimized audit runs pass **24,396 optimization-live gates** with identical
output. The referee also accepts the full normalization, all-depth lower
bounds, joint attainment, precision proof, and real-cone quinary redundancy.

```
audit source 6ba392a11651dfd36e0b336ff3763be87531ceb8bf12772fc66e54f0eb35ba77
audit output 19cfa766d56ebf7229d01306258d9e6e51aca1feb19f029cee77b737ce066a97
```

Root filed these audited artifacts in the seventh checkpoint. No theorem
ID or external-priority reservation was made by this lane.
