---
id: THM-4466
title: "Sharp tournament cubic bound for common-edge and substitution classes"
status: >
  PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY AUDITED. The sharp bound
  F <= 3C <= E^(3/2)/sqrt(3) holds when all positive cyclic triangles share
  one edge, and throughout the stated recursive constant-contact substitution
  grammar. In particular it holds at orders at most four. A separate finite
  scout is not a universal weighted-tournament theorem or a PDE conclusion.
source: tournament-continuation-20260908
depends_on:
  - THM-4465-weighted-tournament-shear-production-and-contact-kernel
script: 04-computation/tournament_cubic_norm_20260908.py
output: 05-knowledge/results/tournament_cubic_norm_20260908.out
script_sha256: fc90e74e0d17d7c32c97e342160fa380a8e7b4141353a169f9446421a64e2025
output_sha256: 887f835b728cdc2342262cf9d9f19fc99070ddd539de1c7da9ab0687c847f9eb
hash_basis: raw LF bytes
audit: >
  Independent tournament-repair agent checked the common-edge argument,
  six-vertex equality constants, and substitution scaling. Explicit exact
  checks survive Python -O; matrix traces and triple sums are independent
  computation paths. Finite scouting is kept separate from proved scope.
---

# THM-4466 -- Sharp tournament cubic bound for common-edge and substitution classes

**PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY AUDITED, within the
classes below.** A universal extension to arbitrary weighted tournaments
is **OPEN in this research pass / unproved here**; this wording does not
assert that the extension is an unsolved problem in the external literature.

## 1. Inheritance, types and exact statement

The closest mechanism is the weighted triple and substitution identity in
[THM-4465 / weighted-tournament-shear-production-and-contact-kernel](THM-4465-weighted-tournament-shear-production-and-contact-kernel.md).
Its hostile example has a reducible weighted tournament with positive
production, even with every edge weight between one half and one. Thus the
inherited unweighted reducibility ceiling cannot control the weighted case.
The least-used sidecar is the **squared amplitude on every edge**, grouped
by actual substitution blocks. The live concepts are cyclic triples,
transitive interface cost, common-edge support, energy sectors and block
substitution.

Let M be an n by n real matrix, n >= 1, with diagonal zero,
M_ij >= 0 and M_ij M_ji = 0. For each unordered pair, its positive entry
orients that edge and gives its weight w_ij. A zero pair has no orientation;
it may be completed either way without changing any expression here. The
vertices are coordinate axes; the pairwise relation is the actual support
orientation. Relabelling preserves all statements. Arbitrary changes of
basis need not preserve this nonnegative cone.

Define

```
E = sum_(i,j) M_ij^2 = ||M||_F^2,
C = sum over directed cyclic triples of the product of their three weights,
T = sum over transitive triples of the product of their three weights,
F = 3C-T = 4 tr(sym(M) skew(M)^2).
```

Zero products contribute nothing. The following are proved:

1. If every positive cyclic triple contains one fixed directed edge, or
   there is no positive cyclic triple, then

       F <= 3C <= E^(3/2)/sqrt(3).                         (1)

2. Every weighted tournament on at most four vertices satisfies (1).
3. The class satisfying the stronger bound on 3C is closed under
   constant-contact substitution, provided each child and the rescaled
   quotient specified below satisfies that bound. Consequently (1) holds
   for every recursively built tournament whose internal quotient is
   transitive or has the common-edge property in item 1. Quotients of
   order at most four are allowed in particular.

The constant is sharp already for a cyclic triangle with equal weights.
There are also equality examples with six vertices and several return paths.
The theorem includes zero amplitudes; the latter example lies on that
boundary and is not a strictly positive complete tournament.

## 2. Common-edge proof and the equality boundary

If C=0 the result is immediate. Otherwise let the common directed edge be
u -> v with weight a > 0. Each positive cyclic triple is its return
v -> x -> u. Write b_x and c_x for these two positive edge weights.
All these edges are distinct, so

```
C = a sum_x b_x c_x
  <= (a/2) sum_x (b_x^2+c_x^2)
  <= (a/2)(E-a^2).
```

The first inequality is the sum of (b_x-c_x)^2 >= 0. For fixed E>0,
a(E-a^2) on 0 <= a <= sqrt(E) is maximized at a^2=E/3, with value
2E^(3/2)/(3sqrt(3)). This proves (1), since T>=0.

Equality in the bound on 3C holds exactly when a^2=E/3, every return has
b_x=c_x, and there is no positive energy outside the common edge and
these returns. Equality in the bound on F also requires T=0; the displayed
support automatically has no positive transitive triple. E=0 is the
separate trivial equality case.

For a rational, non-triangle equality witness use six vertices u,v,x1,...,x4:

```
M_uv=2;
M_vxi=M_xiu=1 for i=1,...,4;
all remaining entries are zero.
```

Then E=4+8=12, C=4*2=8, F=24 and 3F^2=E^3=1728. Thus equality cannot be
classified as concentration on just one cyclic triangle.

At n=4 the number of cyclic triangles in any completed tournament is
4-sum_v binom(d_v,2). The four outdegrees sum to six, and this convex
integer sum is at least two (its balanced minimizer is 1,1,2,2).
There are at most two cyclic triangles. Any two distinct triangles on
four vertices share an edge, whose orientation is the same in both.
This proves item 2, including arbitrary zero patterns. The smaller orders
are immediate.

## 3. Constant-contact substitution and the stronger cyclic bound

Let the children B_i have sizes n_i >= 1. For every quotient edge i -> j,
give all n_i n_j cross-block entries the same magnitude lambda_ij >= 0.
Keep every internal entry of each child. Define the rescaled quotient
edge magnitude

```
a_ij = lambda_ij sqrt(n_i n_j).
```

Denote quantities for this rescaled quotient by E_Q,C_Q,F_Q. Then

```
E_total = sum_i E_i + E_Q,
C_total = sum_i C_i + C_Q.                                (2)
```

The energy identity is direct. A triple with two vertices in one block
and one in another is transitive whenever its product is positive. A
cyclic triple is therefore internal or meets three distinct blocks.
The latter contribution is correct because

```
a_ij a_jk a_ki = n_i n_j n_k lambda_ij lambda_jk lambda_ki.
```

If the stronger bound holds for every child and this rescaled quotient,
then (2) gives

```
3 C_total <= (sum_i E_i^(3/2) + E_Q^(3/2))/sqrt(3)
          <= (sum_i E_i + E_Q)^(3/2)/sqrt(3).
```

The last inequality follows from (x+y)^(3/2) >= x^(3/2)+y^(3/2) for
x,y>=0, iterated finitely. It is strict if two sectors have positive
energy. The induction begins with singletons. At each node, a transitive
quotient has C_Q=0; a common-edge quotient satisfies section 2 for every
choice of its magnitudes, including the required rescaling. This proves
the stated grammar at all orders and all recursion depths.

For comparison, THM-4465 retains the production information that (2) loses:

```
F_total = sum_i F_i + F_Q
          - sum_(i<j) lambda_ij^2 (n_j e_i+n_i e_j),
e_i = sum of the internal edge magnitudes of B_i.
```

The omitted term is nonnegative. It describes why internal activity can
decrease actual production even though cyclic weight adds exactly.
Using unscaled quotient weights in (2) would be incorrect. Allowing
nonuniform cross-block contacts also destroys this small-state law;
the boundary quadratic kernel in THM-4465 is then necessary.

For another all-order equality family take a cyclic quotient of three
edgeless blocks, choose positive cross magnitudes so that the three
rescaled a_ij are equal, and put no energy inside the blocks. There is
only one nonzero energy sector at the substitution step, and all
positive triple products are cyclic. Hence equality is attained.

## 4. Exact audit universe and limits

Reproduction from the repository root:

```
python3 04-computation/tournament_cubic_norm_20260908.py
python3 -O 04-computation/tournament_cubic_norm_20260908.py
```

The two modes agree with the frozen output. The **188,936 always-active
gates** use exceptions rather than Python assertions. They include:

- All 46,878 labelled orientation carriers at n=2,3,4 with magnitudes in
  {0,1,2}. Zero entries deliberately retain multiple orientation carriers;
  the count is not the number of distinct weighted matrices.
- An independent full matrix contraction for the cubic, alongside its
  directed-triple sum, and exact integer checks of 27C^2 <= E^3 and the
  production inequality (with F<=0 handled separately before squaring).
- Eighteen common-edge controls with one through six return paths;
  the six-vertex equality witness above; zero and transitive controls.
- 234 constant-contact recursive substitutions, including 216 height-two
  choices and 18 height-three choices. Both cyclic and production
  composition laws and energy decomposition are checked independently.

The source also labels a separate **FINITE-EXACT scout**: 6,000 seeded
draws, 1,000 at each order 5 through 10, edge magnitudes in {0,1,2,3,4}.
It found no counterexample to either candidate inequality. This is a
finite declared universe, not a proof for arbitrary weighted tournaments,
and no scouting outcome is used as a dependency of sections 1--3.

The underlying identity concerns instantaneous matrices in a specified
coordinate cone. Neither the norm estimate nor its equality case supplies
a pressure law, a spatially realizable trajectory, invariant cone,
self-similar Euler profile, or blowup/regularity conclusion.

## 5. Changed question and cheapest next test

The source-to-target map takes a weighted oriented matrix to its cyclic
triple hypergraph and energy allocation. It preserves C and E, but loses
transitive production cost and the locations of nonuniform contacts. The
needed sidecars are T for F and the weighted boundary kernel for such
contacts. The decisive local tests are exact three-/four-vertex controls,
the multi-return equality witness, and two same-state blocks whose next
contact has opposite signs (THM-4465).

The next candidate is the universal bound 3C <= E^(3/2)/sqrt(3).
The missing step is control of overlapping cyclic triples when there is
neither a common edge nor a constant-contact substitution decomposition.
Search at fixed E over a nontrivial cyclic support, retaining the full
edge amplitudes; any proposed concentration step must preserve or increase
C and must handle the multi-return equality family. A sign tournament or
strong-component label alone cannot supply that step. External priority
and the universal inequality remain unestablished in this pass.
