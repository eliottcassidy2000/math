# Independent second-lane audit: periodic lattice cap and exact scale dichotomy

**Verdict: PASS. No mathematical blocker.** The infinite claims below have been
checked by proof, not inferred from the finite tests. The finite tests are
**FINITE-EXACT**. This audit makes no historical-novelty claim and proves no
new bound for the full Euclidean plane.

Audited producer: `LATTICE_CAP_REPAIR.md`, SHA256
`2c6212ebc62d6ceb18c4d5e5b1b61a1b5cc2ec85a78453d57a378e4a087ae04f`,
plus its exact verification script and outputs. The two requested precision
changes are present: the field graph has vertices `K subset C`, not `K^2`,
and the odd polygon is divided by `n` before claiming unit-distance edges.

## 1. Claims accepted and precise scope

1. Every single positive-distance graph on every rank-two Euclidean lattice
   has a proper coloring with at most three colors, periodic under a
   finite-index subgroup of the lattice.
2. A form not proportional to a rational form has at most two unoriented
   directions at any attained positive norm, hence every such graph is
   bipartite and its maximum over attained layers is two.
3. For a rational-up-to-scale form, let its imaginary quadratic field be
   `K=Q(sqrt(-d))`, with `d>0` squarefree. Every norm layer has the explicit
   finite-period coefficient/top-digit coloring in the producer.
4. The maximum chromatic number over attained positive norm layers of that
   **fixed lattice** is two if the fundamental discriminant is even
   (`d=1,2 mod4`) and three if it is odd (`d=3 mod4`). The maximum statement
   does not say every odd-discriminant layer is three-chromatic.
5. For the graph on `K subset C` whose edges satisfy `N(z-w)=1`, the same
   even/odd fundamental-discriminant cases have chromatic number two/three.
6. For `Q=x^2+xy+6y^2`, the `D=36` layer has exactly the six displayed
   steps, chromatic number three and odd girth exactly 23.
7. The separate dense-flower control has 21 vertices, 57 norm-seven edges,
   25 triangles, degree histogram `5^18,8^3`, and chromatic number three.

The theorem concerns one distance and an additive rank-two **discrete**
lattice. It does not extend to arbitrary finite-rank additive subgroups,
multiple distances, arbitrary number fields, or the full plane. The
finite-period assertion is for lattices, not for the whole divisible
additive group `K`.

## 2. Reduction from real forms and finite periods in the nonarithmetic case

With rows `r(u)=(u1^2,u1*u2,u2^2)`, direct expansion gives

    det(r(u),r(v),r(w)) = det(u,v) det(u,w) det(v,w).

For equal nonzero norms, distinct unoriented directions are linearly
independent in pairs. Three pairs therefore force the normalized
coefficients `(a,b,c)/D` to be rational by Cramer's rule. This repairs the
old genericity argument and establishes the exact two-direction limit.

With one step pair `+/-u`, write `u=g*w`, where `g>0` and `w` is primitive,
then extend `w` to a unimodular basis `(w,z)`. The coloring
`floor(A/g) mod2` at `A*w+B*z` has periods `2g*w,z`, including for negative
`A`. With two pairs `+/-u,+/-v`, the edge-generated group is
`H=Zu+Zv`, of finite index in `Z^2`. In each `H` coset choose a fixed origin
and color by the parity of the unique `H` coordinates. Every element of
`2H` is a period. There is no silent assumption that `H=Z^2`.
An empty layer has the constant one-coloring.

For arithmetic forms, the producer's map

    z(x,y)=x+(B+sqrt(B^2-4AC))*y/(2A)

has norm `F(x,y)/A`; division by any attained edge `z(u)` gives norm
`F(x,y)/F(u)`. This is an injective additive Euclidean similarity into
`K`, not merely an abstract rational vector-space map.

## 3. Independent elementary proof of the local coefficient lemma

This section checks the essential ramification step without requiring any
number-field local-ring theorem. Let `Z_(p)` denote rational numbers with
p-coprime denominator, and let `v_p(0)=infinity`.

Choose `(p,j)` as follows: `(2,1)` if `d=1 mod4`; `(2,0)` if `d=2 mod4`;
and `(p,0)` for any odd prime divisor `p` of `d` if `d=3 mod4`.
Put `pi=j+sqrt(-d)`. The polynomial

    X^2-2jX+(j^2+d)

is Eisenstein at `p`: its constant has p-valuation exactly one and its
linear coefficient is p-divisible. This includes the easy-to-miss shifted
2-adic uniformizer for `d=1 mod4`; using `sqrt(-d)` there would be wrong.

Write a norm-one step as `w=A+B*pi`, with rational `A,B`. Its norm is

    1=A^2+2jAB+(j^2+d)B^2.                         (L)

If `j=0`, the two summand valuations are `2v_p(A)` and `1+2v_p(B)`;
they have different parity, so the smaller is the norm valuation zero.
The even valuation must be zero. Thus `v_p(A)=0` and `v_p(B)>=0`.

If `(p,j)=(2,1)`, put `a=v_2(A)`, `b=v_2(B)`. If `a<=b`, then `A^2`
is the unique term of smallest valuation in (L): the other two have
valuations at least `2a+1`. Hence `2a=0`, so `a=0,b>=0`.
If `a>b`, integrality of the valuation group gives `a>=b+1`; the last
term in (L) uniquely has smallest valuation `1+2b`, an odd integer that
cannot be zero. This case is impossible. Cases where a coefficient is
zero follow by omitting that term, or by the infinity convention.

Consequently **in every case** `A` is a p-unit and `B` is p-integral.
Reducing (L) gives `A^2=1 modp`, so `A=+/-1 modp`, with the same residue
one for both signs when `p=2`. This proves the producer's decisive lemma
independently. In particular it proves edge noncollapse, not just a
necessary norm congruence.

For the whole field, norm-one steps generate an additive group `Gamma`
contained in `Z_(p)+Z_(p)*pi`. On each coset choose an origin and reduce
the constant coefficient of the displacement. Every edge has residue
`+/-1`, giving a graph map to `C_p` for odd `p`, or `K2` at two.
Countable enumeration supplies coset origins without a choice axiom.
Arbitrary field vertices themselves need not be locally integral.

## 4. Why the finite-period formula works on every lattice vertex

Let the normalized lattice basis values have constant coefficients
`A_1,A_2` in the uniformizer basis. Choose `t>=0` making both `p^t A_i`
p-integral, set `m=p^(t+1)`, and define

    ell_i = p^t A_i mod m,
    r(x,y) = ell_1*x+ell_2*y mod m,  0<=r<m.

The `A_i` can have p-denominators; treating every vertex as locally
integral would fail. The scaling above explicitly fixes this. Only these
two constant coefficients need their denominators cleared.

For every complete norm-layer step `s`, the lemma gives
`A(s)=+/-1 modp`, hence `r(s)=+/-p^t modm`.
The target Cayley graph on `Z/mZ` with these two steps is precisely a
union of `p^t` cycles of length `p` (copies of `K2` when `p=2`). Each
component is indexed by the low `t` digits. The top-digit function

    k(r)=floor(r/p^t)

maps each component to the same `C_p`, and each edge changes `k` by
`+/-1 modp`, including the wrap edge. This is a graph map; **the top-digit
function is not asserted to be an additive group homomorphism**.
Alternating colors on `0,...,p-2` and assigning color two to `p-1` gives
three colors when `p` is odd. At two use the digit itself.

This validates the formula on all `Z^2`, not only the edge-generated
subgroup or a sampled connected patch. Periodicity under `mZ^2` is
immediate. The projection need not be surjective for pullback coloring.

## 5. Sharpness on every fixed odd-discriminant arithmetic lattice

The discriminant squareclass is invariant under a rational change of
basis and overall rational scaling, since those change discriminants by
squares. Thus the field in the theorem is well defined.

If `d=3 mod4`, set `n=(d+1)/4` and `theta=(1+sqrt(-d))/2`.
The norm in the basis `(1,theta)` is `x^2+xy+n*y^2`.
The vectors `U=(n,0)`, `V=(n-1,1)`, `W=(-n,1)` all have norm `n^2` and

    nV-nW-(2n-1)U=0.

The walk `V^n,(-W)^n,(-U)^(2n-1)` has length `4n-1=d`.
The first segment has points `((n-1)j,j)`; the second has points
`(n(n-1)+nk,n-k)`. Equating a positive-height point from the two gives
`k(2n-1)=0`, so their only common point is the join. The final segment
strictly decreases along the x-axis; the other two meet that axis only
at its endpoints. The polygon is therefore simple, including `n=1`.
Dividing its coordinates by `n` puts this odd cycle at unit distance in
`K`, proving the field lower bound three.

For a fixed arithmetic lattice `Lambda subset K`, its basis and
`(1,theta)` are rational bases of the same two-dimensional rational space.
If `T` has columns the lattice basis in order coordinates, choose a
positive integer `M` clearing every denominator of `T^(-1)`.
Then `M Z[theta] subset Lambda`; in particular the **entire** polygon,
not just its step vectors, embeds after multiplication by `M`.
Its edges all have one attained positive norm `M^2*n^2`. Pulling back
through the initial similarity gives an odd cycle at one norm in the
original fixed lattice. No ideal-class, maximal-order, shortest-vector,
or multiplicative-closure assumption on `Lambda` is used.

At even fundamental discriminant the prime-two upper bound is two,
and any attained positive norm has an edge, giving lower bound two.
The same edge lower bound applies to nonarithmetic lattices.
This proves all iff directions in the maximum claim.

## 6. The former residual and odd-girth certificate

Completing the square gives, for `x^2+xy+6y^2=36`,
`x^2<=864/23` and `y^2<=144/23`; thus `|x|<=6`, `|y|<=2`.
Exact enumeration gives only

    +/-(6,0), +/-(5,1), +/-(-6,1).

For the three positive generators `U,V,W` above, an integer relation
`aU+bV+cW=0` has `b+c=0` and `6a+11b=0`. Therefore it is an integer
multiple of `(-11,6,-6)`. In a closed walk, its length has the same
parity as `a+b+c`; an odd walk has odd multiplier and length at least
`|a|+|b|+|c|>=23`. The simple polygon for `n=6` attains length 23.
This is an analytic exact odd-girth proof, independent of a search.

At prime 23 and normalization by `(6,0)`, the projection is
`4x+2y mod23`, whose values on `U,V,W` are `1,-1,1`.
It gives the claimed three-color upper bound; the odd cycle gives the
lower bound. The independent program also propagates all exact reachable
lattice positions for walk lengths 1 through 23 and finds first odd
return at 23, a separate finite check with no torus identification.

## 7. Dense-flower cross-check and the information boundary

To avoid copying the producer coordinate arithmetic, the audit uses
`theta=1+omega`, satisfying `theta^2-theta+1=0` and
`N(a+b*theta)=a^2+ab+b^2`. The two multipliers become
`alpha=2+theta`, `beta=3-theta`; the three centers are `0,1,theta`.
Enumerating the six norm-one vectors and forming `alpha P+beta B`
gives 21 distinct points. Enumerating all 210 unordered pairs gives 57
norm-seven edges; all 1330 triples give 25 triangles. Degrees are
`5^18,8^3`. The color `a+2b mod3` is proper, which is the same color
`a+b mod3` in omega coordinates. The triangles prove exact chromatic
number three. Division by `alpha` makes every displayed edge unit length
while remaining inside `Q(sqrt(-3))`.

This verifies the density-versus-chromatic hostile. It does not prove
57 is globally extremal among every planar 21-point configuration. The
same norm-seven arithmetic can support a dense graph, a perfect-code
quotient and a seven-color tiling construction; these are distinct
objects with distinct graph predicates. A seven-chromatic quotient
cannot be reversed into a seven-chromatic source graph.

### Infinite-order relative rotation, with no lattice escape

The appended dense-flower commensurator argument and `SYNTHESIS.md`
Section 6 also pass. For nonzero `z`, commensurability of `zO` and `O`
implies that some positive integer multiple of `z` belongs to `O`, hence
`z in K`. Conversely, denominators of `z,z^(-1)` in the rational basis
`1,omega` can be cleared, giving mutual finite-index containments. For
unit-modulus z this is exactly the rotational commensurator `K intersect
S1`, not the finite rotational stabilizer `O^times` of six units.

Here `alpha+beta=5`, `alpha*beta=7`, and

    alpha/beta = alpha^2/7 = (11+5sqrt(-3))/14.

Its norm is one and its trace is `11/7`. A finite-order element would be
a root of unity, so the sum with its inverse would be an algebraic
integer. A rational algebraic integer is an ordinary integer; `11/7`
is not. Thus the rotation has infinite order (equivalently its angle
is not a rational multiple of pi). Yet the two ideal lattices are
commensurate of index seven, with the explicit ideal-sum certificate

    3(alpha+beta)-2alpha*beta = 1.

Therefore `alpha O+beta O=O`. Repeated powers of this rotation need not
preserve any fixed lattice, but all still lie in `K`; the field-plane
three-color cap remains in force. No claim that every irrational
rotation is harmless, no high-rank additive-group cap, and no global
extremal edge-count claim follows. The synthesis now explicitly includes isolated vertices in the zero-step
case alongside paths and grids; the correction was read back and checked.

## 8. Inherited errors not used as dependencies

The repaired proof does not inherit any of the following THM-419 claims:

- A nonarithmetic form need have only one unoriented direction. The
  irrational rhombus `x^2+y^2+c*sqrt(2)*xy`, for either `c=1` or `c=1/2`,
  has four norm-one vectors and is positive definite.
- An index-three mod-three-null sublattice must be similar to the
  original. For `x^2+6y^2`, the restricted form is `9u^2+6v^2` and its
  minimum is six, whereas an index-three similarity would give minimum
  three.
- The nonzero hyperbolic level `xy=1` over `F3` has four points. It has
  exactly two.
- The four inherited counts `10238,3825,362,61` sum to 16486. They sum to
  **14486**. The old census total is not reproduced or relied upon here.
- The saved output `05-knowledge/results/all_forms_chromatic_s708.out`
  supplies exact chromatic numbers. It prints `(1,0,1): chi-set=[2,3]`
  and assigns maximum three to all 448 forms. The square-lattice case is
  all-layer bipartite, so these labels cannot be interpreted as exact
  chromatic numbers. The generating program was not available at its
  historical path in this worktree; this audit does not speculate about
  which implementation detail produced the labels.
- The parity of an arbitrary **order** discriminant classifies the
  lattice maximum. It is the fundamental field discriminant that does:
  for instance `x^2+3y^2` has order discriminant `-12` but field
  discriminant `-3`, and the maximum is three.

The fresh census is explicitly a verification of coloring certificates,
not a census of exact chromatic numbers. No finite census proves the
all-lattice theorem. Off a rank-two lattice also does not mean off all
algebraic-coordinate configurations, and no rotation-rank or Dehn claim
is established.

## 9. Reproduction, finite universe and provenance

Run from this directory:

    python3 audit.py > normal.json
    python3 -O audit.py > optimized.json
    cmp normal.json optimized.json

Both modes pass and agree byte-for-byte: **715865 explicit gates**.
`audit.py` imports only standard-library modules and neither imports
producer code nor consumes producer output. Its complete universes are:

- Every squarefree `1<=d<=128`, every ramified prime of that field, and
  every reduced parameter `a/b` with `-12<=a<=12`, `1<=b<=12`, tested by
  the exact rational unit-circle parametrization: 28548 unit tests.
- All integer forms `1<=A<=5`, `-5<=B<=5`, `1<=C<=6`, `4AC>B^2` and all
  attained `1<=D<=50`: 280 forms, 4798 layers, 18848 directed steps.
  Steps are obtained independently by solving the quadratic in x for
  each bounded y; every step is used as a normalizer, not only one per
  layer. No primitive-vector filter is used.
- Five deep-denominator controls, with t values `12,10,6,3,2` at primes
  `2,2,3,7,23`; the largest modulus is 12167. Every scalar residue of
  every target modulus, for both signs and positive/negative integer
  representatives, is checked: 109688 directed target edges.
- The first 100 members of the general polygon family, with 300 exact
  embeddings into three nontrivial rational-basis lattices.
- Every reachable lattice position at lengths through 23 in the residual
  six-step graph; every pair and triple of the 21-point dense flower.

The producer script was separately replayed normal and optimized,
with byte-identical results matching its supplied output. Those receipts
are `producer-normal.json` and `producer-optimized.json`.
`SOURCE_PINS.json` records exact input hashes, including historical
sources as provenance only. `SHA256SUMS` closes the audit package.

The primary source [Madore, arXiv:1509.07023v1](https://arxiv.org/pdf/1509.07023)
was directly inspected on 2026-10-25, Proposition 3.2 and its proof on
printed pages 6–7. It supports valuation-ring integrality plus translation
of cosets before reduction. Its stated anisotropic sum-of-squares
hypotheses are not silently identified with the present ramified norm
specialization; the latter was independently proved in Section 3 above.
The source's 2015 introductory planar bound is historical, not a
current-status citation. No priority conclusion follows from this audit.

## 10. Publication correction-draft scope review

The current-theorem preface for `THM-419-2d-lattice-chromatic-cap.md`,
`MISTAKES_ENTRY.md`, and `hexagonal_forbidden_s116m.CORRECTION.md` in
`/tmp/session61-publication/correction-draft/` were reviewed. Their cap,
periodicity, graph-quotient direction, and correction/historical-status
scope pass. The publisher was additionally asked to include the explicit
square-lattice census-label contradiction from Section 8 above. The old
full proof may remain only beneath its explicit superseded-history banner.
