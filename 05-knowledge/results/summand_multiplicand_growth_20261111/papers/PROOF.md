# An explicit norm graph for growing sums and products

**PROVED HERE:** the graph identities and stabilization assertions below.
They are elementary finite-set statements. **CITED:** the supplied papers'
operator frameworks, precisely separated in [SOURCES.md](SOURCES.md).
**FINITE-EXACT:** the small literal controls in [graph_norm_check.py](graph_norm_check.py).
No claim of novelty is made for finite cyclic-set or necklace identities.

The live concepts are pair symmetry, sum/product evaluation, cyclic norms,
ghost observations, and compatible finite truncations. The closest inherited
mechanism is the Euler-product ghost atlas. The hostile example is a hidden
p-cycle that becomes visible after Frobenius. The useful corrected near
miss is that the two compositions FV and VF are different in this graph
model. The least-used sidecar is the actual cyclic orbit, rather than just
its numerical total.

## 1. A graph object on which addition and multiplication are exact

A finite permutation graph is a finite set X with a permutation σ, drawn
with one directed arc x→σ(x) at each vertex. Loops are permitted. It is a
disjoint union of directed cycles C_d; write c_d(X) for their multiplicities.
There is no forced tournament or arbitrary orientation of tied observations.

Addition is disjoint union. Multiplication is the Cartesian product with
diagonal permutation `(x,y)↦(σx,τy)`. In particular,

    C_a × C_b = gcd(a,b) copies of C_lcm(a,b).           (1)

Indeed each pair has orbit length lcm(a,b), and the product contains ab
vertices. This is a literal decomposition of the product graph.

Define the m-th ghost observation by

    g_m(X)=|Fix(σ^m)|=sum_(d|m) d c_d(X),    m≥1.       (2)

Then g_m(X+Y)=g_m(X)+g_m(Y) and g_m(X×Y)=g_m(X)g_m(Y).
All the ghost observations recover the cycle multiplicities by Möbius
inversion:

    c_n(X)=(1/n) sum_(d|n) μ(n/d) g_d(X).              (3)

Thus they determine this unlabelled permutation graph. They do not recover
extra vertex labels, pair supports, or an underlying tournament. Also,
nonnegative integral ghost values alone are not an existence certificate:
g_1=2,g_2=0 would require c_2=-1. The recovered cycle multiplicities must
be nonnegative integers.

The dynamic zeta series is

    Z_X(t)=product_d (1-t^d)^(-c_d(X)),
    t(d/dt)log Z_X(t)=sum_(m≥1) g_m(X)t^m.            (4)

These c_d are cycle/necklace multiplicities, not to be confused with the
standard Witt coordinates in a product of factors `(1-a_dt^d)^(-1)`.

## 2. Frobenius, induction and the norm are physical graph operations

For an integer r≥1, define:

- F_r(X,σ)=(X,σ^r): jump r steps along every old cycle.
- V_r(C_d)=C_rd, extended over disjoint unions. Concretely, replace each
  old vertex by r stages; advance through the stages and apply σ on wrap.
- N_r(X,σ)=(X^r,ν), where
  `ν(x_0,...,x_(r-1))=(σx_(r-1),x_0,...,x_(r-2))`.

The last operation is a twisted cyclic rotation, not ordinary Cartesian
multiplication. Its r-th power applies σ to every coordinate. Therefore

    F_r V_r(X)=rX,           F_r N_r(X)=X^r,          (5)
    N_r(X×Y)=N_r(X)×N_r(Y).                           (6)

Here rX means r disjoint copies. All equalities are isomorphisms of actual
permutation graphs. The projection formula is also exact:

    V_r(X×F_r(Y)) = V_r(X)×Y.                         (7)

One proof of (7) reads ghosts: for r∤m both sides vanish; for r|m they
both have ghost `r g_(m/r)(X) g_m(Y)`. Equation (3) then gives the graph
isomorphism. Alternatively use the r-stage induction construction directly.

For the norm, let d=gcd(r,m). Rotation by m partitions the r positions
into d orbits. Each orbit imposes `σ^(m/d)x=x`, giving

    g_m(N_r X)=g_(m/d)(X)^d,       d=gcd(r,m).          (8)

For prime p this becomes

    g_m(N_pX)=g_m(X)                 if p∤m,
              g_(m/p)(X)^p          if p|m.            (9)

These are the explicit finite graph counterparts of the norm/Frobenius
ghost relations highlighted by Sulyma. They are proved here without any
cohomological input.

**Two necessary boundaries.** First, `V_p F_p=p` is false in this model:
on one fixed vertex it gives C_p, whereas p times that vertex is pC_1.
Thus the model shares (5) and (7), not the whole Raynaud algebra in Zhang.
Second, `N_pX≡X mod V_p` is an identity after group completion, not always
a positive graph difference. For X=C_2, N_2X=C_4, so

    N_2(C_2)-C_2 = V_2(C_2-C_1).                     (10)

The correction on the right is virtual. The positive growth rule in the
next section avoids this loss of positivity by comparing X+Y with X and Y
together.

## 3. The additive failure of the norm is the positive carry graph

For a prime p, every nonconstant two-color word of length p has a free
cyclic color orbit. Split `(X+Y)^p` by the number a of X-colored positions.
There are `binom(p,a)/p` mixed color orbits of each type. After p rotations,
the monodromy is σ on each X coordinate and τ on each Y coordinate.
Consequently

    N_p(X+Y)=N_pX+N_pY
      +sum_(a=1,...,p-1) [binom(p,a)/p]
                          V_p(X^a×Y^(p-a)).          (11)

This proves the integrality of the displayed coefficients combinatorially.
More importantly, each correction term is an actual positive graph.
At p=2,

    N_2(X+Y)=N_2X+N_2Y+V_2(X×Y).                    (12)

Let A_n be the n fixed vertices labelled 1,...,n. Adding one integer gives

    N_2(A_n+1)=N_2(A_n)+1+V_2(A_n).                 (13)

Thus the next alphabet size adds one loop and n directed two-cycles.
The full graph has n fixed diagonal pairs and `n(n-1)/2` two-cycles.
Forgetting the permutation gives n² vertices; applying F_2 makes all of
them fixed. Equation (13) then reads `(n+1)²-n²=2n+1`, with the one
diagonal vertex and the n off-diagonal pairs still distinguished.

For prime p, N_p(A_n) has n loops and `(n^p-n)/p` p-cycles. For composite
r there are genuine intermediate stabilizer strata: for every d|r the
number of d-cycles is

    (1/d) sum_(e|d) μ(e)n^(d/e).                     (14)

For example N_6(A_2) contains 2 loops, one two-cycle, two three-cycles and
nine six-cycles, totalling 64 vertices. The prime/composite distinction is
therefore the presence of intermediate orbit periods, not a decorative
numerical label.

## 4. Doubles and squares are two evaluations of the same fixed diagonal

On N_2(A_n), the permutation swaps the ordered pair (a,b). Two natural
observations are

    s(a,b)=a+b,           t(a,b)=ab.                  (15)

Both are constant on swap orbits. The fixed points are exactly a=b; their
images are 2a under s and a² under t. This proves the common origin of the
double/square dichotomy.

The separate observations lose information: 1+4=2+3, while 1·6=2·3.
The joint observation loses only the already-quotiented order:

    {a,b} ↦ (s,t)=(a+b,ab)                           (16)

is injective. Its inverse is the pair of roots of `z²-sz+t`. Precisely,
positive integer unordered pairs correspond to integers s,t with t>0,
`Δ=s²-4t≥0` a square, and `s>sqrt(Δ)` with matching parity. The two roots
are `(s±sqrt(Δ))/2`. The fixed diagonal is exactly Δ=0.

In the joint chart, the birth of label n contributes

    (s,t)=(n+a,na),   1≤a≤n,
    t=ns-n²,         n+1≤s≤2n.                      (17)

It is a finite line segment ending at `(2n,n²)` on the parabola `t=s²/4`.
The endpoint is the new fixed pair; every other point is a two-cycle.
This is a literal integer growth graph with two value projections.

The bivariate inventory retains all these data. Put

    B_n(u,v)=sum_(1≤a<b≤n) u^(a+b)v^(ab),
    D_n(u,v)=sum_(1≤a≤n) u^(2a)v^(a²),
    P_n(u,v)=2B_n(u,v)+D_n(u,v).                     (18)

Then

    B_n-B_(n-1)=sum_(a=1,...,n-1) u^(a+n)v^(an),
    D_n-D_(n-1)=u^(2n)v^(n²).                        (19)

The coefficient two remembers the free swap orbit. Erasing it changes
ordered representation counts into orbit counts.

There are two growth filtrations, and they must be distinguished. Alphabet
birth is `max(a,b)`; output activation is a+b or ab. The pair (n,n) is born
at alphabet size n but becomes visible at output 2n or n². The joint
label `(max(a,b),a+b,ab)` records both clocks.

At a fixed output m, the sum fibre has m−1 ordered positive pairs and
`1_(2|m)` fixed pairs; the product fibre has τ(m) ordered pairs and
`1_(m is a square)` fixed pairs. Their numbers of unordered orbits are
respectively `floor(m/2)` and `(τ(m)+1_square)/2`.

## 5. Exact stabilization, not an inferred infinite-limit miracle

For a fixed output bound M≥2, every positive pair with a+b≤M already lies
in A_(M-1)². Every positive pair with ab≤M already lies in A_M². Therefore
the corresponding output-restricted labelled graphs stop changing exactly
by these respective alphabet bounds. For sums the pair (1,M−1) shows the
bound is sharp; for products the pair (1,M) does.

This gives a strict transition-map and stabilization certificate for the
growing representation graphs. It is the elementary object that the
inverse-system viewpoint from Shao asks us to retain. No theorem about
almost modules is required here.

The positivity hypothesis matters: allowing all integers gives infinitely
many representations 0=k+(-k); allowing zero gives infinitely many
representations 0=0·k. A fixed output can then receive infinitely many new
alphabet births, so the displayed stabilization bounds fail.

For a more general tower of finite partial-witness sets W_n with actual
restriction maps, all nonempty levels imply a compatible infinite witness
by the finitely branching tree argument. Finiteness and the maps are the
sidecars. They cannot be replaced by a sequence of unrelated counts.
For example the constant two-element sets with identity transitions and
with constant-zero transitions have the same level sizes but different
inverse limits; the value 1 survives forever only in the first. Without
finiteness, even nonemptiness at every level is insufficient: the nested
sets `{n,n+1,...}` have empty inverse limit.

## 6. The exact finite ghost interface and its look-ahead requirement

Let π_N retain the cycle counts c_1,...,c_N, equivalently ghosts g_1,...,g_N.
The operations +, ×, V_r and N_r descend to this finite observation:
(2), (8) and (3) show that output ghosts up to N use only input ghosts up
to N. No bound on the total number of vertices is being imposed.

F_r is different. Since `g_m(F_rX)=g_rm(X)`, computing π_NF_rX requires
π_(rN)X. This look-ahead is necessary: the empty graph and C_(rN) have
the same π_N when r>1, but their Frobenius images differ at cycle length N.
In particular a hidden p-cycle is invisible at level one and becomes p
fixed points after F_p.

This supplies an explicit finite interface for a recursive program:
store nonnegative cycle multiplicities through the required look-ahead,
use (11) for positive additive carries, and lower the observation cutoff
when applying Frobenius. Closing an F_r operation at the old cutoff by
discarding the extra data is not justified.

## 7. What transports, and what still needs a separate constructor

| source | map | preserved predicate | information to retain |
|---|---|---|---|
| ordered labelled pairs | swap quotient | exact unordered representations | fixed versus free orbit, labels |
| unordered pairs | joint sum/product map | existence and full pair identity | both coordinates, discriminant |
| permutation graph | all ghosts | cycle graph up to isomorphism | positivity of Möbius coordinates |
| finite ghost prefix | norm/induction | prefix output counts | Frobenius look-ahead |
| growing positive alphabet | bounded-output restriction | eventual exact representation graph | alphabet/output clocks |

If X is a set of Hamilton paths of a tournament, the graph N_2(X) is a
graph on ordered pairs of paths. Its cardinality H² does not by itself
construct a tournament with that path count, nor does its positive carry
construct a tournament sum. Physical tournament operations need their own
endpoint/support constraints. The old OCF/Witt defect work gives a concrete
reason: closed-walk ghost data can forget simple-cycle overlap information.

The contribution of the papers is the choice of a faithful operator object
and its sidecars. The graph construction and the joint chart make that
choice concrete, while leaving any additional tournament-realization claim
as a separate mathematical obligation.

Reproduce the finite controls with:

    python3 05-knowledge/results/summand_multiplicand_growth_20261111/papers/graph_norm_check.py

The receipt states its finite universe and includes the hostile controls;
the proofs above are valid for every finite permutation graph and every
positive integer alphabet size.
