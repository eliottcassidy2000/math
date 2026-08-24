---
id: THM-4015
title: "First-kind character-sensitive transference via Foster resistance"
status: >
  PROVED ALGEBRA + VERIFIED-EXACT, with a CITED low-dimensional consequence.
  If a rank-d Euclidean lattice admits an obtuse superbase and c is a
  nonzero half-lattice point modulo the lattice, then the distance to the
  lattice times the shortest dual vector odd at c is at most sqrt(d)/2. The
  constant is sharp for Z^d with the all-half centre. Consequently the bound
  holds for every Euclidean lattice in dimensions at most three. THM-4019
  refutes the same sharp constant for arbitrary lattices in every dimension
  at least seven; the weaker arbitrary-lattice d/2 bound remains OPEN.
source: root + transference_independent_audit / LRC14 frontier session, 2026-08-24
audit: >
  PASS. The proof identifies parity representatives with cuts of the conorm
  graph, crossing-edge dual vectors with integral odd characters, their
  squared norms with effective resistances, and the bound with weighted
  Foster trace. An independent Fraction-arithmetic companion checks all 771
  connected unweighted labelled graphs through five vertices, 28 nonuniform
  weighted controls through eight vertices, 799 Foster identities, 12,187
  cuts, dual-Gram resistance, and the sharp orthogonal-star family through
  dimension twelve. Normal and optimized runs match the frozen LF output
  after platform-newline normalization.
  A separate deterministic exact search found no violation of the stronger
  first-kind constant in D4, A4, D5 or 5,728 certified random
  characteristic instances in ranks four and five; this is finite evidence,
  not a general theorem. THM-4019 subsequently gives the exact E7 hostile
  A=6, B=3/2 and an integral counterexample family in every rank at least
  seven; the Foster proof itself and its first-kind scope are unchanged.
depends_on: []
related:
  - THM-4009-euclidean-covering-transference-short-relation-compression
  - THM-4014-lrc14-diagonal-polar-ellipsoid-fastest-coordinate-relation-compression
  - THM-4019-e7-character-transference-sharp-constant-counterexample
audit_script: 04-computation/character_sensitive_transference_first_kind_foster_audit_20260824.py
audit_output: 05-knowledge/results/character_sensitive_transference_first_kind_foster_audit_20260824.out
audit_report: 07-reflections/character-sensitive-transference-first-kind-foster-boundary-20260824.md
search_script: 04-computation/character_sensitive_transference_random_exact_search_20260824.py
search_output: 05-knowledge/results/character_sensitive_transference_random_exact_search_20260824.out
audit_script_sha256: 69f8b927590f4fd13152e36bc0989c86d3c9a266884385f4d253af9be83e07db
audit_output_sha256: 81f7494dd66a9b5ce789901487671c71764d82371b5ed807646f0c4ef247f81d
audit_report_sha256: 5b0c541bbb52793e68b4815bdbe850f009acb1cb1b31e5ec45f4ecaedc45f1f4
search_script_sha256: 15d686582a859f59cacaf5d0550ad9a4472b530a0c97eac50aa1da6c3abde7b1
search_output_sha256: 4f2f10fadabb032530ad1a1da539df5fd77b0df6c3192324ad45cc1d4ffab856
hash_basis: raw LF bytes
---

# THM-4015 -- a sharp odd-character theorem for first-kind lattices

**PROVED ALGEBRA + VERIFIED-EXACT.** Let `L` be a rank-`d` Euclidean lattice
admitting an obtuse superbase, and let

```text
c in (1/2)L \ L,
delta=dist(c,L),
lambda_odd=min{||y||:y in L*, <y,c> in 1/2+Z}.         (1)
```

Then

```text
delta lambda_odd<=sqrt(d)/2.                           (2)
```

The constant is sharp for `L=Z^d` and
`c=(1/2)(1,...,1)`. No extension to arbitrary lattices is claimed;
THM-4019 refutes that sharp extension in every rank at least seven.

## 1. Parity normalization

Write `2c=v in L`, choose `ell in L` closest to `c`, and set

```text
p=v-2ell=2(c-ell).                                     (3)
```

Then `p` is a shortest vector in the nonzero parity class `v+2L`,

```text
||p||=2delta,                                          (4)
```

and the odd vectors in (1) are exactly

```text
{y in L*:<y,p> is odd}.                                (5)
```

## 2. The conorm graph

Let

```text
b_0+b_1+...+b_d=0                                     (6)
```

be an obtuse superbase: `b_1,...,b_d` is a basis of `L` and
`b_i dot b_j<=0` for `i!=j`. Define nonnegative conductances

```text
w_ij=-b_i dot b_j.                                     (7)
```

They form a connected graph on vertices `0,...,d`. Its Laplacian is the
rank-`d` Gram matrix `Q=(b_i dot b_j)`.

Every nonzero class of `L/2L` has a cut representative. For a nonempty
proper subset `S` of the vertices, put

```text
v_S=sum_(i in S)b_i.                                   (8)
```

Complementary subsets give opposite vectors and hence the same parity class.
Choose `S` so that `v_S=p mod 2L`. Selling's identity follows directly from
(6)--(7):

```text
||v_S||^2=sum_(i in S,j notin S)w_ij=:C(S).            (9)
```

Since `p` is shortest in its parity class,

```text
||p||^2<=C(S).                                         (10)
```

## 3. Crossing edges are integral odd characters

For an oriented edge `e=(i,j)`, let `a_e=e_i-e_j` in the zero-sum subspace
of `R^(d+1)`. The map

```text
y -> (<y,b_0>,...,<y,b_d>)                             (11)
```

is an isomorphism onto that subspace. Let `y_e` map to `a_e`. Its pairings
with the basis `b_1,...,b_d` are all integral, so `y_e in L*`. If `e` crosses
the cut `S`, then

```text
<y_e,v_S>=+/-1.                                        (12)
```

Because `p=v_S mod 2L`, this makes `y_e` an odd vector in (5).

The squared norm of `y_e` is the effective resistance of the edge:

```text
||y_e||^2=a_e^T Q^+ a_e=:R_eff(e),                     (13)
```

where `Q^+` is the inverse on the zero-sum subspace. Equivalently, ground
vertex zero; the grounded Laplacian is the Gram matrix of
`b_1,...,b_d`, and (13) is the corresponding incidence quadratic in the
inverse dual Gram matrix.

## 4. Foster trace and sharpness

Weighted Foster is the exact trace identity

```text
sum_e w_e R_eff(e)
 =sum_e w_e a_e^T Q^+ a_e
 =tr(QQ^+)=rank(Q)=d.                                  (14)
```

Restricting to edges crossing `S`, whose conductances sum to `C(S)`, gives
at least one crossing edge with

```text
R_eff(e)<=d/C(S).                                      (15)
```

Equations (10), (12), and (15) imply

```text
||p||^2 lambda_odd^2<=d.                               (16)
```

Using (4) proves (2).

For `Z^d`, take `b_i=e_i` for `1<=i<=d` and
`b_0=-(1,...,1)`. The conorm graph is the unit star. The all-leaves cut has
capacity `d`, every crossing resistance is one, and the all-half centre has

```text
delta=sqrt(d)/2,              lambda_odd=1.            (17)
```

Thus equality holds for every `d`.

## 5. Cited low-dimensional consequence and LRC boundary

Conway and Sloane give a new proof of Voronoi's theorem that every Euclidean
lattice of dimension at most three is of the first kind:
[J. H. Conway and N. J. A. Sloane, *Low-dimensional lattices VI: Voronoi
reduction of three-dimensional lattices*, Proc. R. Soc. Lond. A 436 (1992),
55--68](https://doi.org/10.1098/rspa.1992.0004). Hence (2) holds for every
lattice in dimensions `1,2,3`.

For the LRC lattice `Lambda=P(Z^13)`, the projected coordinate vectors
`g_i=P(e_i)` satisfy

```text
g_i dot g_j=-n_i n_j/||n||^2<0,
sum_i n_i g_i=0.                                      (18)
```

This is an obtuse **weighted circuit**, not an unweighted unimodular obtuse
superbase. Put

```text
b_i=n_i g_i,                 P_n=product_i n_i.        (19)
```

Then `sum b_i=0` and every off-diagonal product is negative. Any twelve
members form a basis only of the sublattice `L_0=<b_i>`, and the quotient is

```text
Lambda/L_0 isomorphic to direct_sum_i Z/n_i Z,
[Lambda:L_0]=P_n.                                     (20)
```

Indeed this is the quotient of `Z^13/DZ^13`, where
`D=diag(n_1,...,n_13)` and `Zn` lies inside `DZ^13`. Thus natural crossing-
edge functionals need not lie in `Lambda*` or retain the odd centre
character.

The determinant defect is visible without coordinates. Put `N=||n||^2` and
write an alleged full superbase as `P(u_i)`. Its unscaled conorm numerator is

```text
q_ij=(u_i dot n)(u_j dot n)-N(u_i dot u_j)>=0.
```

The matrix-tree polynomial of a full obtuse superbase would have to satisfy

```text
tau(q)=N^11.                                           (21)
```

Indeed its conorms are `q_ij/N`, each spanning tree has twelve edges, and
the squared covolume of `Lambda` is `1/N`; Kirchhoff therefore reads
`tau(q)/N^12=1/N`.

For the natural circuit, `q_ij=n_i^2 n_j^2`, and the weighted complete-graph
tree formula gives

```text
tau(q)=P_n^2 N^11.                                     (22)
```

The factor `P_n^2` is exactly the squared index loss.

This does not prove that the full lattice is not of the first kind. For
`n=(1,2,3)`, write `x=P(e_2)` and `y=P(e_3)`. The three vectors

```text
-x-2y,                 y,                 x+y          (23)
```

form a full obtuse superbase, with pair products `-2/7,-1/7,-1/14`, even
though the natural sublattice has index six. In general, `Lambda` admits an
obtuse superbase exactly when there are `u_1,...,u_12 in Z^13` such that

```text
det[u_1 ... u_12 n]=+/-1,       u_13=-sum_(i=1)^12 u_i,
||n||^2(u_i dot u_j)<=(u_i dot n)(u_j dot n)           (24)
```

for all distinct `i,j` in `{1,...,13}`. This is an exact finite certificate
once the `u_i` are supplied, but its search is unbounded. The Smith-index
loss in (20)--(22) is the first failed implication, not a nonexistence proof.

Conditionally, if the rank-twelve LRC lattice did admit an obtuse superbase,
then THM-4009's counterexample inequality `delta>3/7` and (2) would give an
odd-sum integer relation `a` with

```text
||a||_2<7sqrt(3)/3,       sum a_i^2<=16,
||a||_1<=14.                                           (25)
```

A shortest odd row can also be chosen Graver: in a conformal decomposition
of an odd-sum row, one summand is odd and strictly shorter. Statement (25) is
**CONDITIONAL**; no obtuse superbase for a hypothetical LRC counterexample is
known.

The proposed arbitrary-lattice inequality

```text
dist(c,L) lambda_odd<=d/2                              (26)
```

remains **OPEN** here. Ordinary Banaszczyk transference may select the even
index-two dual sublattice and does not prove (26).

The sharper arbitrary-lattice candidate was the same constant as (2). In
basis coordinates, for positive-definite `G` and nonzero `u in F_2^d`, it is

```text
min_(x=u mod 2) x^T G x
  * min_(u dot z=1 mod 2) z^T G^(-1)z <=d.             (27)
```

THM-4019 **REFUTES** (27). For the standard E7 Cartan form, the unique
nonzero mod-two Gram radical has exact primal and odd-dual minima

```text
A=6,                    B=3/2,                    AB=9>7.   (28)
```

The integral family `diag((3E7)^r,2I_k)`, with the repeated radical
characteristic and `d=7r+k`, has

```text
A=18r+2k,               B=1/2,                    AB=d+2r>d. (29)
```

Thus (27) is false in every rank at least seven. It remains OPEN in ranks
four through six; the earlier rank-four/five searches and THM-4019's complete
A6,D6,E6 atlas are **FINITE-EXACT** evidence only. The weaker inequality (26),
equivalently the product bound `AB<=d^2`, remains OPEN.

## 6. Reproduction

The exact companion audits the graph algebra, effective resistance, cut
bound, and sharp boundary. Run

```text
python3 -B 04-computation/character_sensitive_transference_first_kind_foster_audit_20260824.py
python3 -B -O 04-computation/character_sensitive_transference_first_kind_foster_audit_20260824.py
python3 -B 04-computation/character_sensitive_transference_random_exact_search_20260824.py
python3 -B -O 04-computation/character_sensitive_transference_random_exact_search_20260824.py
```

Both normal/optimized pairs match their frozen LF outputs after
platform-newline normalization. **QED.**
