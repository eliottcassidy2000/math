---
id: THM-2615
title: "Endpoint-pair parabolic transvection and translation-gauge boundary"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE; INDEPENDENT HOSTILE AUDIT IN
  FLIGHT.  On the THM-2334 mod-thirteen target quotient, retain the left and
  right endpoint vectors L,R rather than only their difference q=L-R.  Their
  determinant Delta is uniformly erased by common endpoint translation: for
  q nonzero, every Delta occurs thirteen times in the 169-point translation
  orbit.  The 169 difference characters are only the anti-diagonal slice of
  the 28,561-character joint endpoint transform, and an exact two-point
  hostile has all difference characters equal while Delta changes from zero
  to one.  For fixed q nonzero and Delta nonzero, the thirteen endpoint pairs
  form exactly the graph of the order-thirteen projective transvection
  x -> x + det(q,x)q/Delta on P^1(F_13) minus its unique fixed point [q].
  Fixed q and Delta make projective adjacency equivalent to exact endpoint
  gluing.  Endpoint reversal and target swap generate an intrinsic V4 and
  act by inversion/conjugation on the transvection.  This constructs an exact
  relation-address-side transition grammar, not a positive current or an
  intertwiner with the physical LRC root chart; seven fixed transvection
  edges do not close, while thirteen do.  LRC(14) remains open.
source: mac-mini-2026-07-28-endpoint-transition-lift
depends_on: []
related:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2315-marked-target-gain-corolla-and-pairwise-composition-boundary
  - THM-2321-prescribed-root-character-bispectrum-slice-positivity
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2602-commutative-vertex-insertion-and-ordered-transition-curvature-no-go
  - THM-2603-hurwitz-projective-root-owner-atlas-and-nonabelian-seven-edge-trivialization
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2610-chronological-paired-slice-marked-triangle-graft-and-action-axis-boundary
script: 04-computation/lrc14_endpoint_pair_parabolic_transvection_thm2615.py
output: 05-knowledge/results/lrc14_endpoint_pair_parabolic_transvection_thm2615.out
script_sha256: 04c8f2f7bb05e05a0ea0b179a72ca93b0994b6d041c0e960ce3e385d61f3937f
output_sha256: 21114b0dfafda1cd3d6ebc3a32302d247ab5102882bced223f3737afc416e266
hash_basis: working-tree bytes (LF)
---

# THM-2615 -- the target difference hides a parabolic endpoint transition

**RESERVED / PROVISIONAL PROOF CANDIDATE; INDEPENDENT HOSTILE AUDIT IN
FLIGHT.**  Nothing in this file is available as a proved dependency until an
independent audit checks the exact companion, the endpoint typing, and every
scope boundary, and the banner is explicitly promoted.

The post-THM-2602 obstruction is not a shortage of vertex labels.  It is the
absence of an ordered pair

```text
incoming root  ->  outgoing root.
```

THM-2334 already begins with an ordered pair, but then subtracts it.  Its
left and right endpoint modes determine a relation address; the target
quotient retains their difference and discards the common endpoint origin.
The discarded coordinate is exactly where a projective transition lives.

This theorem isolates that finite object.  It yields all three structures
requested repeatedly by the recent frontier:

```text
two target roles:       a two-dimensional endpoint plane;
two endpoint sides:     an intrinsic row/column V4;
thirteen root states:   one parabolic projective orbit;
fourteenth state:       the target direction [q], fixed at the boundary.
```

The construction is exact and intrinsic on relation-address data.  It is not
yet physical LRC chronology: its projective vertices are endpoint target
directions, not predecessor-sheet roots, and no current is proved nonzero on
one fixed transition fibre.

## 1. The endpoint matrix refines the THM-2334 target quotient

Put

```text
k=F_13,                       V=k^2.                       (1)
```

In the specialization of THM-2334 §6, both endpoint frequency vectors lie in
the mod-thirteen relation hyperplane because `13|X,Y`; the delayed-word
dilation is target-neutral because `13|R`.  Apply the canonical quotient map

```text
pi: K_13 -> K_13/L_13 = V
```

separately to the left and right endpoints before taking their difference.
Write the resulting row vectors as

```text
L in V,                       R in V,

             [ L ]
E(L,R)  =   [ R ].                                      (2)
```

Linearity gives the target-vector label already retained by THM-2334:

```text
q=L-R in V.                                               (3)
```

Define the second endpoint coordinate

```text
Delta=det(L,R) in k.                                      (4)
```

Equation (2) is a termwise refinement.  It does not assert that the sum of
current terms in any prescribed `(L,R)` cell is nonzero.  That distinction is
load-bearing because THM-2331 proves term occurrence, while THM-2334 §10 gives
an exact cancellation hostile for grouped relation addresses.

## 2. The difference quotient erases the determinant uniformly

The common endpoint translation group `V` acts freely by

```text
t:E(L,R) -> E(L+t,R+t).                                  (5)
```

The target difference is invariant, but the determinant obeys

```text
det(L+t,R+t)=Delta+det(q,t).                              (6)
```

If `q!=0`, the linear functional

```text
t -> det(q,t)
```

is surjective, with kernel `span(q)`.  Hence, among the `169` translates of
one endpoint pair, every value of `Delta` occurs exactly thirteen times:

```text
#{t in V:det(L+t,R+t)=d}=13             for every d in k. (7)
```

Thus `Delta` is not a function of the target vector, its projectivization, or
any common-translation-invariant statistic.  The minimal witness is

```text
q=(1,0),

(L,R)=((1,0),(0,0))             gives Delta=0,

(L,R)=((1,1),(0,1))             gives Delta=1.            (8)
```

Both pairs have the same difference `q`.

This is an information-loss statement, not a symmetry assertion about the
analytic coefficients.  The current weights need not be invariant under
(5).  What is exact is that once the endpoint pair is pushed forward only by
(3), no later operation on that marginal can recover which translate carried
the mass.

## 3. The 169 target twists are a diagonal slice of a 28,561 transform

Let `J(L,R)` be any complex joint endpoint array and let

```text
A(q)=sum_(L-R=q) J(L,R).                                  (9)
```

For a character `alpha in V^`, finite Fourier transformation gives

```text
Ahat(alpha)
 =sum_(L,R) e_13(alpha.(L-R))J(L,R)
 =Jhat(alpha,-alpha).                                    (10)
```

Therefore the `169` difference characters see only the anti-diagonal

```text
{(alpha,-alpha):alpha in V^}
```

inside the `13^4=28,561` joint endpoint characters.  At every Abel parameter
`rho<1`, THM-2334's absolute convergence permits this regrouping for its
regularized current.  No boundary-limit claim for all `28,561` independently
twisted cells is needed here.

The two point masses in (8) have identical values at all `169` characters in
(10), because their differences agree.  The off-diagonal character

```text
(alpha,beta)=((0,1),(0,0))                               (11)
```

separates them.  Thus the diagonal target-twist bank cannot detect even the
binary question `Delta=0` versus `Delta!=0`.  This is stronger and more typed
than a cardinality warning: it gives an exact kernel witness for the existing
Fourier observable.

As a linear map, (9) is surjective from a `28,561`-dimensional joint space to
a `169`-dimensional difference space.  Its kernel has dimension

```text
28,561-169=28,392.                                       (12)
```

The missing chronology is therefore off-diagonal endpoint information, not
another transform of the same 169 target marginal.

## 4. Every nondegenerate endpoint fibre is a parabolic graph

Fix

```text
q!=0,                         Delta!=0.                   (13)
```

Writing `L=R+q`, the determinant condition becomes

```text
det(L,R)=det(q,R)=Delta.                                  (14)
```

The affine line

```text
H_(q,Delta)={R in V:det(q,R)=Delta}                       (15)
```

has exactly thirteen points.  Neither `R` nor `R+q` lies on `span(q)`, so
projectivization maps (15) bijectively onto

```text
A_q=P(V) minus {[q]},                  |A_q|=13.           (16)
```

Indeed, every projective line other than `[q]` has a unique representative
whose determinant with `q` is `Delta`.

Define the covector and transvection

```text
phi_(q,Delta)(x)=det(q,x)/Delta,

T_(q,Delta)(x)=x+phi_(q,Delta)(x)q.                       (17)
```

Because `phi(q)=0`, the rank-one operator `N=q phi` satisfies `N^2=0`.
Consequently, in characteristic thirteen,

```text
T^m=I+mN,

T^13=I,                   T^m!=I for 1<=m<=12.            (18)
```

The only projective fixed point is `[q]`: a fixed line is an eigendirection,
and the unique eigenspace of this nontrivial unipotent is `ker(phi)=span(q)`.
The remaining thirteen points therefore form one orbit of length thirteen.

For every pair in the fixed fibre, (14) gives `phi(R)=1`, hence

```text
T_(q,Delta)(R)=R+q=L.                                    (19)
```

Thus projectivization identifies the complete fixed-`(q,Delta)` endpoint
fibre with the directed graph

```text
[R] -> T_(q,Delta)[R],                 [R] in A_q.         (20)
```

The target direction `[q]` is exactly the omitted fourteenth point and the
unique fixed boundary.  This is the first exact relation-address-side reason
for the recurrent `13+1` projective shape.

The construction is homogeneous in the endpoint matrix.  Scaling both rows
by `c!=0` sends

```text
(q,Delta) -> (cq,c^2 Delta)                              (21)
```

and leaves `T_(q,Delta)` unchanged.  There are exactly

```text
168=14*12                                                (22)
```

nonidentity transvections: twelve with each fixed projective point.  The
`2,016` parameter pairs `(q,Delta)` with both entries nonzero have scalar
fibres of size twelve under (21).

## 5. Fixed determinant turns projective adjacency into exact gluing

Orient an endpoint matrix from its right row to its left row:

```text
R -> L=R+q.                                               (23)
```

Consider consecutive matrices with the same nonzero `q,Delta`.  Exact
chronological gluing is

```text
R_(i+1)=L_i.                                              (24)
```

Because both sides of (24) have determinant `Delta` with `q`, each is the
unique normalized representative of its projective line.  Therefore

```text
[R_(i+1)]=[L_i]       iff       R_(i+1)=L_i.              (25)
```

This is the gain over the bare projective target line: the determinant is a
normalization sidecar that makes projective composition literal.  A glued
chain satisfies

```text
R_i=R_0+i q,

[R_i]=T_(q,Delta)^i[R_0].                                (26)
```

The same calculation gives the sharp hostile:

```text
seven fixed edges:       R_7=R_0+7q != R_0,

thirteen fixed edges:    R_13=R_0.                       (27)
```

Thus the parabolic endpoint grammar does not itself cancel the seven-clock
holonomy.  It realizes a genuine ordered transition and simultaneously shows
why a clock-dependent conjugation or inverse action is still necessary.

## 6. Endpoint reversal and target swap form an intrinsic `V4`

There are two independent involutions on (2):

```text
rho: swap the endpoint rows L,R;

kappa: swap the two target columns a,b.                   (28)
```

They commute and generate `C_2 x C_2`.  Let `S` denote the coordinate swap on
`V`, whose determinant is `-1`.  Their effects are

| operation | target vector | determinant | transvection |
|---|---|---|---|
| identity | `q` | `Delta` | `T_(q,Delta)` |
| endpoint reversal `rho` | `-q` | `-Delta` | `T_(q,Delta)^(-1)` |
| target swap `kappa` | `Sq` | `-Delta` | `S T_(q,Delta) S^(-1)` |
| double swap | `-Sq` | `Delta` | `S T_(q,Delta)^(-1) S^(-1)` |

The identities in the last column are the explicit formulas

```text
T_(-q,-Delta)(x)=T_(q,Delta)^(-1)(x),

T_(Sq,-Delta)(Sx)=S T_(q,Delta)(x).                      (29)
```

Therefore `Delta` is the sign character that changes under either basic
involution and survives their product.  This is a genuine `V4` because both
binary relations are intrinsic: endpoint order and target role.  It is not
the quartic matching `V4` of THM-2598/2606, though the same origin-loss pattern
is present.  Nor does (20) define a tournament: it supplies one directed
neighbour at each vertex, not a pairwise observable on every unordered pair.

## 7. Degenerate boundary and exact census

The boundary cases are completely sharp.

If `q=0`, then `L=R` and `Delta=0`; there are `13^2=169` such matrices.  If
`q!=0` but `Delta=0`, equations (14)--(15) force

```text
R=cq,                       L=(c+1)q,                    (30)
```

with thirteen choices of `c`.  Two choices contain a zero row, and the other
eleven projectivize both endpoints to `[q]`.  No nontrivial parabolic edge is
defined there.

Over all `13^4=28,561` endpoint matrices, the exact partition is

```text
q=0:                                      169,
q!=0, Delta=0:                          2,184,
q!=0, Delta!=0:                        26,208.             (31)
```

For each fixed `q!=0`, every determinant fibre has thirteen matrices, in
agreement with both (7) and (15).

## 8. What this changes at the LRC frontier

The theorem separates four projective lines that had the same cardinality but
different types:

```text
THM-2315: target-gain directions;
THM-2321: root-bispectrum shapes mapped abstractly to those gains;
THM-2603: physical-root projective atlas in PSL_2(F_13);
here:     projectivized endpoint residues of one relation term.            (32)
```

Equal cardinality is not a map.  The new exact statement is narrower and more
useful: on the fourth object, the target direction `[q]` is the fixed boundary
and `(q,Delta)` determines an ordered parabolic transition on the other
thirteen points.

This gives a typed interface with THM-2602.  A successful current refinement
must retain, before marginalization,

```text
(clock ell, common ancestry, L, R, q=L-R, Delta),         (33)
```

prove survival in one `q!=0,Delta!=0` sector, and glue the right endpoint of
one clock to the left endpoint of the next.  Then (25) turns its projective
edge table into an honest ordered kernel.  The existing `A(q)` and its 169
twists cannot do this by Section 3.

Two further maps are still absent:

1. an **action axis** from the two-dimensional target vector `q` to the
   scalar root correction required by THM-2602/2610;
2. an **ancestry-preserving projective intertwiner** from the endpoint line in
   (20) to the physical predecessor-root line of THM-2603.

Such an intertwiner must send `[q]` to the chosen physical cusp, conjugate the
parabolic transition to the physical root action, preserve endpoint order and
clock, and carry a nonzero current rather than only labels.  THM-2603's
nonabelian seven-edge norm suggests how clock-dependent conjugation can close;
(27) proves that a fixed endpoint transvection cannot.

The determinant is also an exact origin invoice.  Choosing a representative
inside one common-translation fibre chooses `Delta`; forgetting the origin
makes every determinant equally possible by (7).  This parallels the
affine-origin sidecar of THM-2606, but over `F_13^2` rather than `V4=F_2^2`.

## 9. Exact companion and scope

Run

```bash
python3 04-computation/lrc14_endpoint_pair_parabolic_transvection_thm2615.py
python3 -O 04-computation/lrc14_endpoint_pair_parabolic_transvection_thm2615.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_endpoint_pair_parabolic_transvection_thm2615.out.
```

The dependency-free companion enumerates all `28,561` endpoint matrices and
checks `605,997` exact assertions, including:

1. the complete census (31) and every translation identity (6)--(7);
2. the 169-character diagonal hostile (8)--(11);
3. every fixed-`(q,Delta)` fibre and all `26,208` nondegenerate edges;
4. all projective fixed points, thirteen-cycles, and exact gluing identities;
5. all 168 transvections and their twelve-fold homogeneous parameter fibres;
6. every row/column `V4` sign law and transvection conjugacy; and
7. the seven-step nonclosure and thirteen-step closure in (27).

No floating point, random fixture, external package, or literature claim is
used.

The theorem proves a finite algebraic transition object and a sharp loss
theorem.  It does not prove that a prescribed joint endpoint current survives
Abel cancellation, that a positive Boolean carrier selects `Delta!=0`, that
successive endpoint terms share ancestry, that endpoint directions are
physical roots, that an action axis exists, that seven clocks close, that any
scalar row is impossible, or that LRC(14) holds.

QED, conditional on independent promotion audit.
