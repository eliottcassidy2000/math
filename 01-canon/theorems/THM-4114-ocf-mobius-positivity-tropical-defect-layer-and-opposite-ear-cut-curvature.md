---
id: THM-4114
title: "OCF Mobius positivity, tropical defect layer, and opposite ear-cut curvature"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. OCF normalization makes
  every induced-subset Mobius atom and mixed Boolean difference nonnegative.
  A fixed base word has a matching-controlled first repair layer. By
  contrast, every one-vertex ear response is a submodular directed-cut
  quadratic. Full transfer or interaction rank does not force an interval,
  while stability and Lorentzianity already fail on minimal faces.
source: codex-frontier-synthesis-creative-20260825g
depends_on:
  - THM-002-ocf
  - THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension
  - THM-4099-squarefree-gap-transfer-and-mixed-insertion-boundary
related:
  - THM-1975-the-path-cover-polynomial-is-the-refined-compositional-invariant
  - THM-2602-commutative-vertex-insertion-and-ordered-transition-curvature-no-go
  - THM-4111-uniform-ear-average-and-recursive-selected-bank-growth
  - THM-4115-uniform-ear-cut-walsh-variance-and-sharp-growth-refinement
  - THM-4118-ear-response-lattice-and-stateful-unit-component-intervals
script: 04-computation/tournament_gap_polynomial_ocf_curvature_candidate.py
output: 05-knowledge/results/tournament_gap_polynomial_ocf_curvature_candidate.out
independent_audit_script: 04-computation/tournament_ocf_mobius_ear_curvature_thm4114_independent_audit.py
independent_audit_output: 05-knowledge/results/tournament_ocf_mobius_ear_curvature_thm4114_independent_audit.out
script_sha256: 81d2ce6dd8f1ac1a3cd3b56a457bc0fef21a463377c95200d0c1033d1a80958f
output_sha256: 9c28f557b99966837ca0226dc96829f3300375282050c2b5247c88a8443aebc9
semantic_sha256: 0eeda799bfc4764678b3aa8c5067b04b7a74944e5b0b6dbf519f427b09dd56b5
independent_audit_script_sha256: 4fc9763c589c60979d53793b071eaf4dd6393fd932e5b040158d0ae3b40a05ef
independent_audit_output_sha256: 2ecbaaadce2303ee9d33daddca4d8f9f50e9ad01e30c9f3931dd9e2994f720e7
independent_semantic_sha256: 3917956b00891bdfb012076298867a484f4b0a41950e0e18326a782067131c5e
hash_basis: raw LF bytes
primary_audit: >
  PASS. A self-contained implementation checks OCF data on all 33,867
  labelled tournaments through order six and 2,131,018 induced masks. The
  Mobius census separately covers 220,387 atoms and 811,255 mixed differences
  through order five. It also checks tensor controls through k=6, the strict
  degree-delay and algebraic hostiles, and every ear cut through order five.
  Normal, optimized, and frozen outputs byte-match.
independent_audit: >
  ACCEPT. A clean-room implementation imports no primary code and independently
  verifies OCF Mobius positivity and parity, deletion/promotion, zeta Smith
  forms and characteristic-dependent ranks, the matching defect layer, both
  stability hostiles, the full-cube ear cut, and the gapped full-rank code-8
  witness. It reproduces every stated finite census with the Mobius scope
  explicitly bounded through order five. Normal, optimized, and frozen
  outputs byte-match; the smallest theorem failure is none.
---

# THM-4114 -- OCF positivity and the opposite ear-cut cube

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4099 makes finite vertex insertion compositional by retaining a
squarefree gap polynomial. This theorem identifies two different Boolean
cubes hidden in that construction.

- In the **presence cube**, a coordinate decides whether a fixed inserted
  vertex is present. OCF makes all additive mixed differences nonnegative.
- In the **ear cube**, a coordinate reverses the incident arc between one new
  vertex and one old vertex. THM-4097 makes the response a directed cut, so
  every pair difference is nonpositive and every higher difference vanishes.

The sign reversal is not a contradiction: the coordinates perform different
operations. It is the precise obstruction to moving positivity from
THM-4099 directly into the selected-ear intervals of THM-4102/4104.

## 1. OCF normalization of the presence cube

Let a finite tournament `W` have a fixed decomposition

```text
V(W)=B disjoint_union I,             B nonempty.             (1)
```

For `S subseteq I`, put

```text
h(S)=H(W[B union S]),
Z_(W/B)(X)=sum_(S subseteq I) h(S)X_S                       (2)
```

in `R_I=Z[X_i:i in I]/(X_i^2:i in I)`. By THM-4099, this is the
relative gap-transfer polynomial.

Let `Gamma` range over collections of pairwise vertex-disjoint
directed odd cycles in `W`, including the empty collection. Write

```text
A_I(Gamma)=V(Gamma) intersect I,
c_A=sum_(Gamma:A_I(Gamma)=A) 2^|Gamma|,
Phi_(W/B)=sum_(A subseteq I)c_A X_A,
P_I=product_(i in I)(1+X_i).                                (3)
```

> **Theorem 1 (OCF-normalized face).** In `R_I`,
>
> ```text
> Z_(W/B)=P_I Phi_(W/B)
>        =sum_(A subseteq I)c_A X_A
>           product_(i in I minus A)(1+X_i).                 (4)
> ```

### Proof

THM-002 gives

```text
h(S)=sum_(Gamma in W[B union S])2^|Gamma|.                  (5)
```

A packing of `W` occurs in `W[B union S]` exactly when
its inserted support lies in `S`. Therefore

```text
h(S)=sum_(A subseteq S)c_A.                                 (6)
```

The coefficient of `X_S` on either right side of `(4)` is
the same zeta sum. **QED.**

The symbol `Phi` deliberately does not reuse THM-002's `Omega(W)`, which is
the odd-cycle **conflict graph**, not a polynomial. Boolean Mobius inversion
gives

```text
c_A=sum_(T subseteq A)(-1)^(|A|-|T|)h(T)>=0.                (7)
```

For disjoint `A,S subseteq I`, define

```text
Delta_A h(S)
 =sum_(T subseteq A)(-1)^(|A|-|T|)h(S union T).             (8)
```

Substitution of `(6)` gives

```text
boxed: Delta_A h(S)=sum_(D subseteq S)c_(A union D)>=0.      (9)
```

Thus `h` is a totally monotone Boolean capacity: every additive mixed
difference is nonnegative. If `A` is nonempty, every term on the right of
`(9)` is even, so every nonempty mixed difference is even as well. The parity
of the atoms is exact:

```text
c_emptyset=H(W[B]) is positive and odd,
c_A is even for every nonempty A.                           (10)
```

The first identity counts OCF packings lying entirely in `B`. Every packing
counted by nonempty `c_A` has at least one cycle and therefore even weight.
The atom is a weighted packing sum; it is not generally twice an unweighted
cycle count.

## 2. Deletion, promotion, and flattening rank

Fix `i in I`. Let `Z_del` be the relative polynomial
after deleting `i`, and let `Z_promote` be the relative
polynomial after moving `i` into the base. Then

```text
Z_(W/B)=Z_del+X_i Z_promote,                                (11)
Z_(W/B)|_(X_i=0)=Z_del,
Coeff_(X_i) Z_(W/B)=Z_promote.                              (12)
```

Here `Coeff_(X_i)` is the normal-form linear operation
`f_0+X_i f_1 -> f_1`. Ordinary differentiation is not being asserted to be a
derivation of the squarefree quotient. For the normalized polynomials,
cancelling `P_(I minus {i})` in `(11)` gives

```text
Phi_B
 =Phi_del+X_i(Phi_promote-Phi_del).                         (13)
```

The coefficient of `X_J` in the difference is
`c_(J union {i})`, counting packings that use `i`.
Hence the difference is coefficientwise nonnegative.

Now split `I=U disjoint_union V` and form

```text
M_(A,D)=h(A union D),             A subseteq U, D subseteq V,
C_(J,K)=c_(J union K),            J subseteq U, K subseteq V. (14)
```

If `zeta_U[A,J]=1[J subseteq A]`, then `(6)` is

```text
boxed: M=zeta_U C zeta_V^T.                                 (15)
```

Boolean zeta matrices are integral unimodular matrices. Thus `M` and `C`
have the same Smith normal form over the integers and the same rank over each
specified field. The numerical rank may still depend on the characteristic;
no Boolean, tropical, or nonnegative-rank claim is made. This preserves
linear separability, not arithmetic spacing.

## 3. The first tropical layer of one base word

Fix one permutation `P` of the base in THM-4099, and let
`Z_P` be its product of local gap factors before summing over base
words. Let

```text
D(P)={internal gaps g=(a,b) of P with b->a},
delta(P)=|D(P)|,
L_g(X)=sum_(i in I)1[a->i->b]X_i.                           (16)
```

Write `MA` for deletion of monomials containing a square. Then

```text
boxed:
[Z_P]_(degree delta(P))
 =MA(product_(g in D(P))L_g).                               (17)
```

Every bad gap has zero constant coefficient, so minimum degree assigns one
distinct inserted vertex to every bad gap and leaves all other gaps empty.
For `|S|=delta(P)`,

```text
[X_S]Z_P=per(A[D(P),S]),       A[g,i]=1[a->i->b].           (18)
```

Thus the degree-`delta(P)` layer is live exactly when the
compatibility graph has a matching saturating `D(P)`.

The matching is only the singleton shadow. The exact first live degree is

```text
r(P)=min {
  sum_(g in D(P)) |S_g| :
  S_g nonempty, kappa_g(S_g)>0,
  and the S_g are pairwise disjoint },                      (19)
```

with value infinity if no allocation exists. This is a disjoint-hyperedge
repair problem, and necessarily `r(P)>=delta(P)`.

The smallest strict delay occurs at order four. In lexicographic tournament
code `42`, take

```text
P=(0,1), I={2,3},
arcs: 1->0, 0->2, 3->0, 1->2, 3->1, 2->3.                  (20)
```

Its bad gap has local coefficients

```text
(empty,{2},{3},{2,3})=(0,0,0,1).                           (21)
```

More explicitly, its left, internal, and right factors are

```text
1+V+UV,             UV,             1+U+UV,                (21a)
```

whose squarefree product is `UV`. Hence `delta(P)=1` but `r(P)=2`, repaired
by the block `0->2->3->1`. Polynomial degree here measures inserted-footprint
size, not cycle count, total cycle length, or the number of bad gaps.

There is a parallel 2-adic boundary. Every coefficient `h(S)` is
odd, so direct coefficient valuations of `Z` vanish. For a nonzero
OCF atom put

```text
m_A=min {|Gamma|:A_I(Gamma)=A}.                             (22)
```

Then `nu_2(c_A)>=m_A`, with equality exactly when the number of
minimum-cardinality packings with support `A` is odd. Tropical data
becomes informative only after OCF normalization.

## 4. Full transfer rank can have a sparse value image

For `1<=i<=k`, take directed triangles

```text
b_i -> u_i -> v_i -> b_i                                   (23)
```

and form their **ordinal sum**, orienting every inter-gadget arc from the
earlier gadget to the later one. Use base
`B={b_i}` and inserted banks `U={u_i},V={v_i}`. Every
Hamiltonian path concatenates gadget paths, so the face flattening is

```text
M=[[1,1],[1,3]] tensor_power k.                             (24)
```

Consequently

```text
rank_Q(M)=2^k,
|det M|=2^(k 2^(k-1)),                                     (25)
```

with the same full rank over fields of odd characteristic. Its Smith factors
are `2^j` with multiplicity `binom(k,j)`: the base matrix has Smith form
`diag(1,2)`, and unimodular tensor powers preserve the diagonal tensor form.
In characteristic two its rank is therefore `1`, not `2^k`. Yet the entry
formula and value set are only

```text
M_(A,D)=3^|A intersect D|,
{3^j:0<=j<=k}.                                              (26)
```

Thus maximal separator rank can coexist with an exponentially sparse image.
Rank loses spacing, gcds, and the overlap graph needed for interval filling.

## 5. Positivity does not imply stability or Lorentzianity

For a directed triangle with singleton base,

```text
Z=1+U+V+3UV.                                                (27)
```

This is not real stable: `U=sqrt(-1)` and
`V=(-2+sqrt(-1))/5` both lie in the open upper half-plane and make
`Z=0`.

The OCF-normalized footprint polynomial is

```text
Phi=1+2UV,                                                  (27a)
```

which also is not real stable: `U=V=i/sqrt(2)` is an upper-half-plane zero.

The ordinary homogenization

```text
T^2+TU+TV+3UV                                               (28)
```

has Hessian eigenvalues `-3,1,4`, hence two positive eigenvalues
and fails the degree-two Lorentzian Hessian-signature criterion.

Even a homogeneous slice can fail. For the order-five tournament of lex code
`10` with base `{0}`, the degree-two face has Hessian

```text
[0 3 1 3]
[3 0 1 1]
[1 1 0 3]
[3 1 3 0],                                                 (29)
```

with characteristic polynomial

```text
(t^2-4t-13)(t^2+4t-1).                                     (30)
```

Both factors have one positive root. Nonnegative OCF atoms and total Boolean
monotonicity therefore do not imply stability or Lorentzianity.

## 6. Every ear cube is a directed cut

Let `T` be a tournament on `V`. For
`S subseteq V`, adjoin `x_S` with
`x->v` exactly when `v in S`. Use THM-4097's
nonnegative data `Start(v),End(v),Q(a,b)`.

Build a directed network on `{s,t} union V` with capacities

```text
s->b: Start(b),
a->t: End(a),
a->b: Q(a,b).                                               (31)
```

For source side `{s} union (V minus S)`, the three classes of
crossing arcs are exactly the three sums in THM-4097. Hence

```text
boxed: H(T+x_S)=capacity of that directed s-t cut.           (32)
```

For distinct `u,v` disjoint from the context `S`,

```text
Delta_u Delta_v H(T+x_S)
 =-[Q(u,v)+Q(v,u)]<=0,                                     (33)

Delta_A H(T+x_S)=0,              |A|>=3.                   (34)
```

Thus every ear response is a submodular quadratic set function on the **full**
Boolean cube. The nonconstant signatures that preserve a strong base form a
punctured domain, not a lattice, so submodularity there must always mean the
restriction of the full-cube identity.

On a directed triangle, the presence-face pair difference is `+2`, whereas
all three ear-cube pair differences are `-2`. Higher differences separate as
well: order-four code `42` with base `{0}` has footprint atoms

```text
(1,0,0,0,0,0,2,2),                                        (34a)
```

so its presence-cube third difference is `2`, while every ear-cube third
difference is zero. Presence and ear cubes use different coordinates—vertex
selection versus incident-arc reversal—so `(9)` cannot be imported into
ear-interval propagation.

## 7. Full ear interaction rank still misses an interval

The smallest strong hostile is the order-five lex-code-`8`
tournament. It has

```text
H(T)=9,
{H(T+x_S):emptyset proper_subset S proper_subset V}
 ={15,17,19,23,25,27,29,33,37,41}.                         (35)
```

The internal odd values `21,31,35,39` are missing. Nevertheless its symmetric
interaction matrix

```text
K_(u,v)=Q(u,v)+Q(v,u),   K_(u,u)=0,

[0  6  4  4  8]
[6  0 10 10  4]
[4 10  0 10  4]
[4 10 10  0  6]
[8  4  4  6  0]                                           (36)
```

has

```text
det K=-24320,        rank_Q K=5.                            (37)
```

The directed `Q`-matrix itself also has rational determinant `5520`. A
`{0,2}` versus `{1,3,4}` response flattening has Smith form `(1,2,2,8)`, so
it has full row rank over every odd characteristic while the response image
remains gapped.

Hence strongness, full interaction rank, and exact cut curvature do not make
one parent's ear image solid. Exhaustively:

| base order | labelled strong bases | nonsolid individual ear images |
|---:|---:|---:|
| `3` | `2` | `0` |
| `4` | `24` | `0` |
| `5` | `544` | `544` |

The missing coordinate is arithmetic overlap—atom sizes, gcds, or a
unit-step adjacency graph. THM-4111's mean forces growth of selected maxima;
neither that mean nor rank supplies the overlap sidecar.

## 8. Exact variance companion

THM-4115 is a **PROVED** companion to the ear-cube curvature here.
Writing `K=Q+Q^T=2w`, its Walsh calculation gives

```text
hat X({i,j})=-K_ij/4,
Var(X)=1/4 sum_i h_i^2+||K||_F^2/32.                        (38)
```

The hostile in `(35)` is code `1015` in THM-4115's LSB-first upper-pair
encoding after reversing the five vertex labels. That proved audit gives

```text
H=9, F_1=30, W=33, Var=305/4,
variance floor=994/33, oddceil=31, actual maximum=41.        (39)
```

Its image still misses `21`. Thus adding exact variance to the full-rank
curvature control does not supply interval overlap. The exact
response-sufficient quotient is `(Start,End,Q)->(H,h,w)`; passing further to
the scalar variance loses the labelled coefficient arrangement.

## 9. Exact verification and scope

The primary dependency-free referee declares

```text
OCF: 33,867 labelled tournaments through order six,
     2,131,018 induced masks, 314,690 directed odd cycles;
Mobius: 220,387 atoms and 811,255 mixed differences through order five;
ear cut: 1,099 tournaments through order five,
         33,866 cut values, 83,506 pair differences,
         41,480 higher differences.                         (40)
```

It also freezes tensor powers through `k=6`, literal gadget
tournaments through `k=3`, the degree-delay witness, both Hessian
hostiles, and the strong-ear census.

Run

```text
python3 -B 04-computation/tournament_gap_polynomial_ocf_curvature_candidate.py
python3 -B -O 04-computation/tournament_gap_polynomial_ocf_curvature_candidate.py
PYTHONHASHSEED=0 python3 -B 04-computation/tournament_gap_polynomial_ocf_curvature_candidate.py
python3 -B 04-computation/tournament_ocf_mobius_ear_curvature_thm4114_independent_audit.py
python3 -B -O 04-computation/tournament_ocf_mobius_ear_curvature_thm4114_independent_audit.py
PYTHONHASHSEED=0 python3 -B 04-computation/tournament_ocf_mobius_ear_curvature_thm4114_independent_audit.py
```

This theorem proves structural identities and sharp no-go boundaries. It
does not prove a solid interval at a new order, compress an arbitrary gap
interface below its full squarefree state, or prove the global
Hamiltonian-path spectrum conjecture. THM-4115 separately proves the exact
Walsh variance and a sharper growth theorem. **QED.**
