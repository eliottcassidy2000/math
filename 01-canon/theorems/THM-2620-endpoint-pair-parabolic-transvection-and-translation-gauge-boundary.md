---
id: THM-2620
title: "Endpoint-pair parabolic transvection and translation-gauge boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  After declaring a
  deep-leg allocation
  gauge on the THM-2334 mod-thirteen target quotient, retain endpoint vectors
  L,R rather than only q=L-R.  Their
  determinant Delta is uniformly erased by common endpoint translation: for
  q nonzero, every Delta occurs thirteen times in the 169-point translation
  orbit.  The 169 difference characters are only the anti-diagonal slice of
  the 28,561-character joint endpoint transform, and an exact two-point
  hostile has all difference characters equal while Delta changes from zero
  to one.  For fixed q nonzero and Delta nonzero, the thirteen endpoint pairs
  form exactly the graph of the order-thirteen projective transvection
  x -> x + det(q,x)q/Delta on P^1(F_13) minus its unique fixed point [q].
  Modulo common scalar, the pointed fibres are PGL2(F13); the determinant-square
  half is PSL2(F13), splitting the transvections into two 84-element classes.
  Left owner C7 and right endpoint C13 actions split PSL2 into twelve free
  91-point torsors without an internal order-91 element.  Endpoint reversal
  and target swap generate a V4 action, not globally a torsor.  The canonical
  THM-2603 descending trace-two SL lift closes at -I, exposing a central C2
  scale invoice.
  This is relation-address algebra, not a positive current or physical-root
  intertwiner; seven fixed transvection edges do not close.  LRC(14) remains
  open.
source: mac-mini-2026-07-28-endpoint-transition-lift
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2603-hurwitz-projective-root-owner-atlas-and-nonabelian-seven-edge-trivialization
related:
  - THM-2315-marked-target-gain-corolla-and-pairwise-composition-boundary
  - THM-2321-prescribed-root-character-bispectrum-slice-positivity
  - THM-2602-commutative-vertex-insertion-and-ordered-transition-curvature-no-go
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2610-chronological-paired-slice-marked-triangle-graft-and-action-axis-boundary
  - THM-2615-physical-diagonal-toric-kernel-and-dipole-radon-invoice
  - THM-2616-cross-time-target-future-diagonal-and-principal-action-no-go
  - THM-2619-psl2f13-seven-edge-norm-minimal-projective-kernel-and-retained-target-obstruction
script: 04-computation/lrc14_endpoint_pair_parabolic_transvection_thm2620.py
output: 05-knowledge/results/lrc14_endpoint_pair_parabolic_transvection_thm2620.out
script_sha256: 7edadae89f99f190bc313c11759eb0fa0e252d2b1efb6baa226965fca520dabd
output_sha256: 85f1e62c8dc43b079696051517dd69405db0f69a138f79b1d9d662991aadf8f4
hash_basis: working-tree bytes (LF)
---

# THM-2620 -- the target difference hides a parabolic endpoint transition

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The post-THM-2602 obstruction is not a shortage of vertex labels.  It is the
absence of an ordered pair

```text
incoming root  ->  outgoing root.
```

THM-2334's atomic address can be allocated to an ordered endpoint pair before
the target quotient subtracts it.  The quotient retains their difference and
discards the common endpoint origin.  The discarded coordinate is exactly
where a projective transition lives.

This theorem isolates that finite object.  It yields all three structures
requested repeatedly by the recent frontier:

```text
two target roles:       a two-dimensional endpoint plane;
two endpoint sides:     an intrinsic row/column V4;
thirteen root states:   one parabolic projective orbit;
fourteenth state:       the target direction [q], fixed at the boundary.
```

The construction is exact after fixing the allocation gauge.  It is not
yet physical LRC chronology: its projective vertices are endpoint target
directions, not predecessor-sheet roots, and no current is proved nonzero on
one fixed transition fibre.

## 1. The endpoint matrix refines the THM-2334 target quotient

Put

```text
k=F_13,                       V=k^2.                       (1)
```

In the specialization of THM-2334 §6, rename the delayed-word dilation `D`
to avoid confusing it with the right endpoint, and write one atomic relation
address as

```text
r=ell_0+d-r_0,

ell_0=u+D beta,             d=m e_(c_3),             r_0=v.   (2a)
```

Here `13|D`, while the deep leg `d` is target-visible on the live branch.
Apply the canonical quotient map

```text
pi: K_13 -> K_13/L_13 = V
```

to the complete address.  To express its target label as an endpoint
difference one must first declare where the deep leg lives.  This theorem
uses the **left-deep allocation gauge**

```text
L=pi(ell_0+d),                 R=pi(r_0),

             [ L ]
E(L,R)  =   [ R ].                                     (2b)
```

Then, and only with this declared allocation, linearity gives the
THM-2334 target-vector label

```text
q=pi(r)=L-R in V.                                         (3)
```

The natural unshifted coefficient endpoints instead differ by
`q-pi(d)`, not by `q`; since `pi(d)!=0` on the live branch, silently calling
them `(L,R)` would be false.  More generally every allocation of the same
difference has the form

```text
L_t=pi(ell_0+d)+t,               R_t=pi(r_0)+t,          (3a)
```

and THM-2334's decomposition gauge realizes every `t in V`.  The balanced
choice `pi(ell_0)+pi(d)/2, pi(r_0)-pi(d)/2` is
reversal-equivariant in the quotient, but need not be represented by two
literal integer frequency modes.  The allocation is therefore a gauge
choice, not an intrinsic splitting of the atomic address.

Define the second endpoint coordinate

```text
Delta=det(L,R) in k.                                      (4)
```

Equation (2b) is a termwise refinement.  It does not assert that the sum of
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

This homogeneity is only a classification of the projective transition.
The target vectors `q` and `cq` are distinct THM-2334 target fibres, and the
analytic current weights need not descend under common scaling.

### 4.1 The pointed fibres are `PGL_2`, with a `PSL_2` square-class half

Write columns rather than endpoint rows and put

```text
B=[q R],                   det(B)=Delta,

U=[[1,1],[0,1]].                                           (22a)
```

Equation (17) is exactly

```text
T_(q,Delta)=B U B^(-1).                                   (22b)
```

Modulo common scalar, `B` is a point of `PGL_2(F_13)`.  The map

```text
[B] -> (B U B^(-1), [B e_2])

PGL_2(F_13)
  -> {(T,x): T nonidentity transvection,
             x notin Fix(T)}                              (22c)
```

is a bijection.  The first component fixes `[q]=[Be_1]`, while the marked
point is `[R]=[Be_2]`; conversely a transvection and one nonfixed point recover
the unique projective frame.  Hence

```text
2,184=168*13.                                             (22d)
```

Projection to `T` is a principal right `C_13=<U>` bundle: right multiplication
`B -> B U^s` preserves `T` and advances the marked point to `T^s[R]`.

The determinant square class sharpens this to `PSL_2`.  Set

```text
C=[q,Delta^(-1)R] in SL_2(F_13).
```

Then

```text
C^(-1) T_(q,Delta) C=U^(Delta^(-1)).                     (22e)
```

The `168` transvections therefore split into two `PSL_2` conjugacy classes of
size `84`, distinguished exactly by the quadratic character of `Delta`.
Square `Delta` gives the class of `U`; nonsquare `Delta` gives the class of
`U^2`.  Each half has

```text
1,092=84*13,                                              (22f)
```

The square half is `PSL_2` itself as a pointed endpoint bundle; the nonsquare
half is the other `PGL_2/PSL_2` coset, a `PSL_2` torsor rather than the group
with a chosen identity.  In the square half, THM-2603's abstract homogeneous
carrier `Omega=PSL_2/<C>` becomes, after conjugating `<C>` to `<U>`, the
**unpointed** transvection class; the right `C_13` fibre restores the missing
position on its thirteen-cycle.

Finally let `A` be THM-2603's order-seven owner.  On `PSL_2`, left
multiplication by `<A>` and right multiplication by `<U>` commute as actions
and act jointly freely.  Thus

```text
1,092=12*7*13,                                            (22g)
```

giving twelve `C_7 x C_13` torsor orbits.  This is a shared carrier for the
two cyclic grammars, not an element or internal subgroup of order `91`.

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

## 6. Endpoint reversal and target swap give a `V4` action, not a torsor

There are two independent involutions on the **allocated matrix** (2b):

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
involution and survives their product.  This is a genuine `V4` action, but it
is not free.  In the operation order `(1,rho,kappa,rho kappa)`, exact Burnside
censuses are

| carrier | fixed counts | orbit-size census |
|---|---:|---:|
| all `28,561` endpoint matrices | `(28561,169,169,169)` | `1:13, 2:234, 4:7020` |
| `26,208` nondegenerate matrices | `(26208,0,0,144)` | `2:72, 4:6516` |
| `2,016` nonzero `(q,Delta)` pairs | `(2016,0,0,144)` | `2:72, 4:468` |
| `168` transvections | `(168,0,0,24)` | `2:12, 4:36` |

The `144` nondegenerate double-swap fixed matrices have the form

```text
[[a,b],[b,a]],                 (a-b)(a+b)!=0.             (29a)
```

Thus `V4` is a symmetry action with genuine size-two orbits, not a torsor.
It is not the quartic matching `V4` of THM-2598/2606, though the same
origin-loss pattern is present.

The row-swap formula applies literally to the retained left-deep allocated
matrix.  Reversing the analytic Fourier term and then reapplying the left-deep
convention is a different operation: the deep leg moves, so an additional
common translation may change `Delta`.  Likewise target-column swap is
covariant only when the deep vector is swapped with its target coordinates.
No analytic current symmetry is asserted without those sidecars.  Nor does
(20) define a tournament: it supplies one directed neighbour at each vertex,
not a pairwise observable on every unordered pair.

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

### 7.1 Projective closure retains a central `C_2` scale invoice

After determinant normalization, exact vector frames require an `SL_2` lift
rather than only a `PSL_2` class.  In THM-2603's affine chart put

```text
U=[[1,1],[0,1]],

G=[[7,5],[10,11]],                  G^7=-I.              (31a)
```

For the canonical trace-`+2` endpoint transvections, the descending ordered
norm is

```text
product_(i=6,...,0) G^i U G^(-i)=-I.                    (31b)
```

The actual THM-2603 lift satisfies `P C P^(-1)=-U`; its seven negative
factors multiply instead to `+I`.  Thus the same projective norm can carry
opposite exact central signs.

The complete signed refinement of the normalized six-state atlas is

| state | `+I` exponents | `-I` exponents |
|---|---|---|
| `3F` | `{10,11}` | `{6,8,9}` |
| `3R` | `{2,3}` | `{4,5,7}` |
| `5F` | `{2,12}` | `{8,10,11}` |
| `5R` | `{1,11}` | `{2,3,5}` |
| `6F` | `{1,3}` | `{9,11,12}` |
| `6R` | `{10,12}` | `{1,2,4}` |

Here `F` is the algebraic order `A_0...A_6` and `R` is `A_6...A_0`, with
`A_i=g_t^i U^a g_t^(-i)`.  The moving trace is `t-a` in the forward order
and `t+a` in the reverse order.  Cayley--Hamilton gives seventh power `-I`
for signed traces `{3,5,6}` and `+I` for their negatives; the telescope then
gives the table.

Consequently any intertwiner that lifts the projective atlas to literal
endpoint vectors needs a central `C_2` scale/orientation sidecar.  Projective
identity is not exact vector identity.

## 8. What this changes at the LRC frontier

The theorem separates four projective lines that had the same cardinality but
different types:

```text
THM-2315: target-gain directions;
THM-2321: root-bispectrum shapes mapped abstractly to those gains;
THM-2603: abstract Hurwitz root-owner projective atlas in PSL_2(F_13);
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
   (20) to THM-2603's abstract Hurwitz chart, followed by a separate map from
   that chart to a physical predecessor-root carrier.

Such an intertwiner must first send `[q]` to an abstract Hurwitz cusp and only
then, through a separate proved map, to a physical cusp.  It must conjugate
the parabolic transition to the physical root action, preserve endpoint order
and clock, retain the central sign in Section 7.1, and carry a nonzero current
rather than only labels.  THM-2603's nonabelian seven-edge norm suggests how
clock-dependent conjugation can close; (27) proves that a fixed endpoint
transvection cannot.

The determinant is also an exact origin invoice.  Choosing a representative
inside one common-translation fibre chooses `Delta`; forgetting the origin
makes every determinant equally possible by (7).  This parallels the
affine-origin sidecar of THM-2606, but over `F_13^2` rather than `V4=F_2^2`.

THM-2615's lawful `13 by 13` present/bare-endpoint square is a one-dipole
projection of the joint `V x V` endpoint array in Section 3.  Its diagonal
Radon sums recover one scalar target colour, but that projection necessarily
forgets the two-axis determinant `Delta`.  Thus the two new objects are
successive resolutions, not competitors:

```text
full allocated endpoint array J(L,R)
   -> one-dipole 13 by 13 square
   -> diagonal Radon target colour.                       (33a)
```

A constructive proof must populate the full array on one positive
chronological carrier before using either the transvection grammar or the
Radon consumer.  THM-2616's cross-time diagonal no-go rules out inferring
that array from separately positive source/future marginals.

## 9. Exact companion and scope

Run

```bash
python3 04-computation/lrc14_endpoint_pair_parabolic_transvection_thm2620.py
python3 -O 04-computation/lrc14_endpoint_pair_parabolic_transvection_thm2620.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_endpoint_pair_parabolic_transvection_thm2620.out.
```

The dependency-free companion enumerates all `28,561` endpoint matrices and
checks `672,925` exact assertions, including:

1. the complete census (31) and every translation identity (6)--(7);
2. the 169-character diagonal hostile (8)--(11);
3. every fixed-`(q,Delta)` fibre and all `26,208` nondegenerate edges;
4. all projective fixed points, thirteen-cycles, and exact gluing identities;
5. all 168 transvections, their `84+84` determinant-square split, and the
   pointed `PGL_2/PSL_2` bundle bijections;
6. the twelve free left-`C_7`/right-`C_13` orbits on `PSL_2`;
7. every row/column `V4` law, all four Burnside/orbit censuses, and every
   transvection conjugacy;
8. the six signed `SL_2` closure rows and the canonical central `C_2` defect;
   and
9. the seven-step nonclosure and thirteen-step closure in (27).

No floating point, random fixture, external package, or literature claim is
used.

The theorem proves a finite algebraic transition object and a sharp loss
theorem.  It does not prove that a prescribed joint endpoint current survives
Abel cancellation, that a positive Boolean carrier selects `Delta!=0`, that
successive endpoint terms share ancestry, that endpoint directions are
physical roots, that an action axis exists, that the declared deep-leg
allocation is selected analytically, that the central sign is transported,
that seven clocks close, that any scalar row is impossible, or that LRC(14)
holds.

The final independent hostile audit checked the deep-leg allocation, every
PGL/PSL and square-class assertion, the twelve left/right torsor orbits, all
four Burnside tables, the six signed closure rows, and every scope boundary.
Normal and optimized runs independently byte-match the stored transcript and
the recorded hashes.

QED.
