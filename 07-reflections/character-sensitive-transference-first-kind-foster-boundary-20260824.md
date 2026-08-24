# Character-sensitive transference: a sharp Foster theorem and the first non-graphic wall

## Status and verdict

The proposed general statement remains **OPEN** in this audit:

```text
L a rank-d Euclidean lattice,
c in (1/2)L \ L,
delta=dist(c,L),
lambda_odd=min{||y||: y in L*, <y,c> in 1/2+Z}

question: delta*lambda_odd <= d/2 ?
```

THM-4019 finds no counterexample to this weaker `d/2` bound.  Ordinary
Banaszczyk transference does not prove it: that theorem controls the shortest
unphased vector of `L*`, which may lie in the even-character index-two
sublattice.

There is, however, a decisive positive result on a large and structurally
natural class.

> **PROVED ALGEBRA (sharp).** If `L` is a lattice of Voronoi's first kind,
> i.e. it has an obtuse superbase, then
>
> ```text
> delta*lambda_odd <= sqrt(d)/2.                       (1)
> ```
>
> Equality holds for `L=Z^d` and
> `c=(1/2)(1,...,1)`.  Thus the constant in `(1)` is sharp on this class.

Since `sqrt(d)/2<=d/2`, this proves the proposed inequality on every lattice
of Voronoi's first kind.  By the classical low-dimensional theorem of
Voronoi/Selling, reproved by Conway--Sloane, every Euclidean lattice of rank
at most three is of the first kind.  Consequently the proposed inequality is
true in dimensions `1,2,3`, with the stronger sharp bound `(1)`.  The cited
source for the low-dimensional classification is:

- J. H. Conway and N. J. A. Sloane, *Low-dimensional lattices. VI. Voronoi
  reduction of three-dimensional lattices*, Proc. R. Soc. Lond. A 436
  (1992), 55--68, [DOI 10.1098/rspa.1992.0004](https://doi.org/10.1098/rspa.1992.0004).

The Foster/effective-resistance argument below is proved here and is not
imported from that citation.

## 1. Exact parity-coset formulation

Write `2c=v in L`, and choose `ell in L` closest to `c`.  Put

```text
p=v-2ell=2(c-ell).
```

Then `p` is the shortest vector in the nonzero class
`alpha=v+2L in L/2L`, and

```text
||p||=2 delta.                                         (2)
```

Moreover, the odd dual vectors are exactly

```text
{y in L*: <y,p> is odd}.                               (3)
```

Indeed `p` and `v` differ by an element of `2L`, and
`<y,c>=<y,v>/2` modulo integers.

In basis coordinates `L=A Z^d`, with Gram matrix `G=A^T A` and a nonzero
parity vector `u in F_2^d`, define

```text
m_G(u)=min{x^T G x: x=u mod 2},
m_G^odd(u)=min{z^T G^(-1) z: u dot z=1 mod 2}.          (4)
```

Then

```text
delta*lambda_odd=(1/2)sqrt(m_G(u)m_G^odd(u)).           (5)
```

Thus the proposed general bound is the quadratic-form inequality

```text
m_G(u)m_G^odd(u) <= d^2.                               (6)
```

The first-kind theorem proves the much stronger right side `d`.  A
provisional extrapolation of that constant to arbitrary lattices is now
**REFUTED** by THM-4019:

> **REFUTED SHARP EXTENSION.** For the standard E7 Cartan matrix and its
> unique nonzero mod-two Gram radical `u`,
>
> ```text
> m_G(u)=6,                 m_G^odd(u)=3/2,
> m_G(u)m_G^odd(u)=9>7.
> ```

The exact witness amplifies integrally to every rank `d>=7`: for
`d=7r+k`, the Gram form `diag((3E7)^r,2I_k)` with the repeated radical
characteristic has product `d+2r`.  This closes the arbitrary-lattice
`sqrt(d)/2` route, but it does not approach the still-open `d^2` threshold
in equation `(6)`.  Ranks four through six for the sharp extension also
remain open.

## 2. Obtuse superbase as a conductance graph

Let

```text
b_0+b_1+...+b_d=0
```

be an obtuse superbase for `L`: `b_1,...,b_d` is a lattice basis and
`b_i dot b_j<=0` for `i!=j`.  Put

```text
w_ij=-b_i dot b_j >=0.                                 (7)
```

These are conductances on a graph with vertices `0,...,d`.  The graph is
connected because its Laplacian is the rank-`d` Gram matrix

```text
Q=(b_i dot b_j)_(0<=i,j<=d).                           (8)
```

Every nonzero parity class has a canonical cut representative.  Namely,
after using `b_1,...,b_d` as basis, choose a nonempty proper subset
`S subset {0,...,d}` and set

```text
v_S=sum_(i in S)b_i.                                   (9)
```

Complementary subsets give opposite vectors, and the other subsets give all
nonzero classes of `L/2L`.  Selling's identity is immediate from `(7)`:

```text
||v_S||^2=sum_(i in S,j notin S)w_ij=:C(S),            (10)
```

the capacity of the cut.  If `v_S` represents the class of the shortest
vector `p` in `(2)`, then simply

```text
||p||^2<=C(S).                                         (11)
```

Notice that the proof needs only a canonical representative, not a separate
claim that it is the unique shortest representative.

## 3. A crossing edge is an odd dual vector

For an oriented edge `e=(i,j)`, let `a_e=e_i-e_j`, regarded as a vector in
the zero-sum subspace of `R^(d+1)`.  The map

```text
y -> (<y,b_0>,...,<y,b_d>)
```

is an isomorphism from the ambient Euclidean space to that zero-sum
subspace.  Let `y_e` be the unique vector mapped to `a_e`.  Every pairing
with a superbase vector is `0`, `1`, or `-1`; hence

```text
y_e in L*.                                             (12)
```

If `e` crosses the cut `S`, then

```text
<y_e,v_S>=sum_(r in S)(a_e)_r=+/-1.                    (13)
```

Since `p=v_S mod 2L`, `(13)` says that `y_e` is an odd-character vector for
the original half-lattice centre.

The squared norm of `y_e` is exactly the effective resistance of the edge:

```text
||y_e||^2=a_e^T Q^+ a_e=:R_eff(e),                     (14)
```

where `Q^+` is the Moore--Penrose inverse on the zero-sum subspace.  This is
also obtained without pseudoinverses by grounding vertex zero: the grounded
Laplacian is the Gram matrix of `b_1,...,b_d`, and the corresponding incidence
quadratic form in its inverse is `(14)`.

## 4. Foster's trace identity closes the cut

The weighted Foster identity is one line of exact linear algebra:

```text
sum_e w_e R_eff(e)
 =sum_e w_e a_e^T Q^+ a_e
 =tr(Q Q^+)
 =rank(Q)
 =d.                                                    (15)
```

Restrict `(15)` to the crossing edges of `S`.  Their conductances sum to
`C(S)`, so at least one crossing edge satisfies

```text
R_eff(e)<=d/C(S).                                      (16)
```

Equations `(11)`--`(16)` give

```text
||p||^2 lambda_odd^2
 <= C(S) min_(e crosses S)R_eff(e)
 <= d.
```

Using `(2)` proves `(1)`.

## 5. Equality and exact controls

For `L=Z^d`, take the obtuse superbase

```text
b_i=e_i (1<=i<=d),       b_0=-(1,...,1).
```

Its conorm graph is the unit star centred at zero.  For
`c=(1/2)(1,...,1)`, the leaf-versus-centre cut has capacity `d`, every
crossing effective resistance is `1`, and

```text
delta=sqrt(d)/2,        lambda_odd=1.                  (17)
```

Thus equality holds in `(1)` for every `d`.

The independent standard-library companion is:

```text
04-computation/character_sensitive_transference_first_kind_foster_audit_20260824.py
SHA-256 69f8b927590f4fd13152e36bc0989c86d3c9a266884385f4d253af9be83e07db
```

with frozen output:

```text
05-knowledge/results/character_sensitive_transference_first_kind_foster_audit_20260824.out
SHA-256 81f7494dd66a9b5ce789901487671c71764d82371b5ed807646f0c4ef247f81d
```

Normal and optimized runs agree byte-for-byte.  Fraction arithmetic checks
all `771` connected unweighted labelled graphs through five vertices, `28`
nonuniform weighted graph controls through eight vertices, all `799` Foster
identities, all `12,187` noncomplementary cuts, the incidence/dual-Gram
effective-resistance identity, and the sharp star boundary through `d=12`.

Reproduction:

```bash
python3 -B 04-computation/character_sensitive_transference_first_kind_foster_audit_20260824.py
python3 -B -O 04-computation/character_sensitive_transference_first_kind_foster_audit_20260824.py
```

The computation audits the algebra; the proof is Sections 1--4.

## 6. Exact stopping boundary for LRC(14)

For the LRC projected lattice

```text
Lambda=pi(Z^13) subset n^perp,
```

the projected coordinate vectors `g_i=pi(e_i)` satisfy

```text
g_i dot g_j=-n_i n_j/||n||^2<0       (i!=j),
sum_i n_i g_i=0.                                      (18)
```

Thus LRC supplies an enticing **obtuse weighted circuit**.  But an obtuse
superbase requires an unweighted relation `sum b_i=0` in which any twelve
members form a unimodular lattice basis.  Equation `(18)` does not provide
that: the speed coefficients `n_i` carry arithmetic index, and the edge
incidence vectors from the Foster proof need not lift to integral vectors of
`Lambda*` with odd centre character.

This is the first failed implication:

```text
obtuse weighted frame
  DOES NOT YET GIVE
unimodular obtuse superbase / integral crossing-edge characters.          (19)
```

No claim is made that the LRC lattice is or is not of Voronoi's first kind;
that must be tested or proved separately.  A valid extension would need one
of:

1. a character-sensitive transference theorem for arbitrary index-two
   lattice pairs at the surviving `d/2` scale (the first-kind
   `sqrt(d)/2` scale is false by THM-4019);
2. a weighted-circuit Foster theorem retaining the Smith-index/integrality
   sidecar in `(18)`; or
3. a proof that the particular projected lattice admits some different
   unimodular obtuse superbase.

Conditionally, if the first-kind theorem did apply in dimension `12`, then
THM-4009's `delta>3/7` would give

```text
lambda_odd < (sqrt(12)/2)/(3/7)=7sqrt(3)/3,
sum y_i^2<=16,
||y||_1<=14.                                           (20)
```

That would produce an odd-character relation dramatically shorter than the
unphased square-norm-`195` row.  Equation `(20)` is **CONDITIONAL**, not a
current LRC reduction.

The research target is now sharply typed: preserve the finite Smith-index
payload while converting the weighted obtuse circuit `(18)` into an
effective-resistance odd character.  Reusing ordinary unphased transference
would simply drop the coordinate that this question was designed to retain.
