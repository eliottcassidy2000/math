# Every planar lattice norm layer is periodically three-colourable

**Status: PROVED by the self-contained arguments below; exact controls are
FINITE-EXACT.** This is a repair/completion of the repository's THM-419 route,
using the classical ramified-prime colouring mechanism already recovered in
THM-4558. It is not a claim of historical novelty. It says nothing new about
six- or seven-colourability of the full Euclidean plane.

## 1. Inherited mechanism and precise correction target

- Closest proved mechanism: THM-418, `01-canon/theorems/THM-418-lattice-chromatic-dichotomy-of-the-plane.md`, for square/triangular lattices; THM-4558, `01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md`, Lemma R, for conjugation-fixed primes.
- Corrected near miss: THM-419, `01-canon/theorems/THM-419-2d-lattice-chromatic-cap.md`, lists the all-lattice claim as partly conjectural but its consequence prose states it unconditionally. Its nonarithmetic one-direction assertion and its degenerate mod-three similarity assertion are false. The cap itself is true, for the different proof below.
- Canonical hostile: the irrational rhombus has two, not one, unoriented equal-norm directions; the arithmetic form `x^2+6y^2` has a mod-three kernel sublattice that is not similar to the original lattice.
- Least-used sidecar: a single attained norm vector normalizes the whole layer into one imaginary-quadratic **field** plane. It is not necessary for a ramified-prime sublattice to be similar to the original lattice.

Concept board: equal-norm directions / rational Gram field / ramified valuation /
periodic quotient / odd relations / hourglass root-normal carrier.

## 2. The all-lattice cap, including arbitrary real Gram matrices

Let `Q(x,y)=a x^2+b xy+c y^2` be a positive-definite real quadratic form.
For `D>0`, let `U(Q,D)` be the graph on `Z^2` joining `z,z+s` exactly when
`Q(s)=D`. If no such `s` exists, its chromatic number is one. Otherwise:

**Theorem A.** `chi U(Q,D)<=3`, and there is a proper colouring periodic
under a finite-index subgroup of `Z^2`. If `Q` is not proportional to a
rational form, every norm layer with edges is bipartite.

### 2.1 Three direction pairs force arithmeticity

For `u=(u1,u2)`, put `r(u)=(u1^2,u1*u2,u2^2)`. Direct expansion gives

```
det(r(u),r(v),r(w)) = det(u,v) det(u,w) det(v,w).       (1)
```

Here rows are ordered as displayed. Equal nonzero norm vectors on the same
real line are opposites, so three distinct **unoriented** norm directions
make the right side nonzero. The equations

```
[r(u);r(v);r(w)] (a,b,c)^t = D (1,1,1)^t
```

then imply `(a,b,c)/D in Q^3`. This is exact, not a generic-position claim.
Consequently a genuinely nonarithmetic form has zero, one or two direction
pairs at every achieved norm.

One direction pair gives components isomorphic to the doubly-infinite
path. Two pairs `+/-u,+/-v` have real-linearly independent `u,v`; hence
`m u+n v=0` with integers implies `m=n=0`. Each component is a translate
of `Zu+Zv` and is coloured by `(m+n) mod2`. No assertion that `u,v` form a
basis of the ambient lattice is made.

These colourings can be periodic. With two directions, choose representatives
modulo `H=Zu+Zv` and use the same coefficient parity in each coset; `2H` is
a finite-index period lattice. With one direction, write `u=g w` with `w`
primitive, extend `w` to a unimodular basis `(w,z)`, and colour `A w+B z`
by `floor(A/g) mod2`. Its period lattice is `2g Zw+Zz`.

### 2.2 Arithmetic layers embed in imaginary-quadratic norm-one graphs

If `Q` is proportional to a rational form, rescale it to
`F(x,y)=A x^2+B xy+C y^2` with integers `A,B,C`, `A>0`, and
`Delta=B^2-4AC<0`. An achieved squared norm rescales to the integer
`h=F(u)>0`. Put

```
K=Q(sqrt(Delta)),
z(x,y)=x+(B+sqrt(Delta))*y/(2A),
alpha=z(u),
Phi(x,y)=z(x,y)/alpha.
```

These are injective additive maps into `K`, and

```
N(z(x,y))=F(x,y)/A,
N(Phi(x,y))=F(x,y)/h.                                (2)
```

Thus every edge difference maps to norm one. Dividing by the attained
`alpha` is indispensable: unnormalized field reduction need not preserve
the target unit norm.

### 2.3 Ramification gives a three-colour field plane

For squarefree `d>0`, let `K=Q(sqrt(-d)) subset C`. The graph `G(K)` has
vertex set `K` and edges `z~w` exactly when `N(z-w)=1`; it is not the
algebraic Euclidean graph on `K^2`. Choose a ramified prime as follows:

| d | p | uniformizer pi |
|---|---|---|
| `d=1 mod4` | 2 | `1+sqrt(-d)` |
| `d=2 mod4` | 2 | `sqrt(-d)` |
| `d=3 mod4` | any odd prime dividing d | `sqrt(-d)` |

In each case the minimal polynomial of `pi` is Eisenstein at `p`. The
unique prime `P` over `p` is fixed by conjugation, its residue field is
`F_p`, and conjugation acts trivially on that residue field. Normalize the
valuation by `v_P(pi)=1`, so `v_P(p)=2`.

Write `w=a+b*pi`, `a,b in Q`. The two summands have valuations
`2 v_p(a)` and `2 v_p(b)+1`; when nonzero these have opposite parity.
Therefore there is no valuation cancellation, and

```
v_P(w)=min(2 v_p(a),2 v_p(b)+1).
```

If `N(w)=1`, conjugation invariance gives `2v_P(w)=0`; hence `a,b` are
`p`-integral. Reduction is `w -> a modp`, and `a^2=1 modp`. Every norm-one
step therefore reduces to `+/-1` (to `1` when `p=2`). This proves both
integrality and noncollapse of edges, not only a formal residue equation.

For the whole field graph, let `Gamma` be the additive group generated by
all norm-one vectors. Then `Gamma` is contained in the local ring. On each
coset `x0+Gamma`, reduce `x-x0`; this is a graph map to
`Cay(F_p,{+/-1})`, which is the cycle `C_p` for odd `p`, and `K2` for
`p=2`. Coset representatives may be chosen as least elements of a fixed
countable enumeration, so no choice principle is needed here. Odd cycles
admit three colours. Pullback along (2) proves the arithmetic cap.

This is the classical local-colouring mechanism; Madore's valuation-ring
formulation is discussed in the primary-source record below.

## 3. Explicit finite-period certificates, without choosing components

The same proof gives a direct finite formula for every arithmetic layer.
Write the normalized basis values in the uniformizer basis as

```
Phi(e_i)=A_i+B_i*pi,  A_i,B_i in Q.
```

Choose `t>=0` with `p^t A_1,p^t A_2` p-integral, put `m=p^(t+1)`, and let
`ell_i` be the residues of `p^t A_i` modulo `m`. Rational residues here
use the inverse of their denominator, which is coprime to `p`. Define

```
r(x,y)=(ell_1*x+ell_2*y) mod m,   0<=r<m,
k(x,y)=floor(r(x,y)/p^t),        0<=k<p.
```

For an edge step `s`, Section 2.3 proves that its constant coefficient is
p-integral and `+/-1 modp`. Thus `r(s)=+/-p^t modm`. Adding a step preserves
the lower `t` digits and changes the top digit by `+/-1 modp`, including
wraparound. Hence pull back any proper `C_p` colouring through `k`.
An explicit choice is

```
p=2: c(k)=k;
p odd: c(k)=k mod2 for 0<=k<=p-2, and c(p-1)=2.
```

This is periodic under `mZ^2`. Only the constant coefficients need to be
cleared for this formula; local integrality of **both** unit-vector
coefficients was already proved. The quotient bound has the right
direction: a colouring of the finite quotient pulls back; a large
chromatic number of some other quotient would not give a lower bound on
the lattice.

## 4. The cap is sharp: exact maximum over all scales of a fixed lattice

**Theorem B.** For an arithmetic planar lattice, form the imaginary
quadratic field `K=Q(sqrt(-d))` from its rational Gram discriminant, with
`d` squarefree. The maximum chromatic number over its achieved positive
norm layers is

```
2, if d=1 or 2 mod4 (even fundamental discriminant);
3, if d=3 mod4 (odd fundamental discriminant).
```

For a nonarithmetic lattice this maximum is two. This does **not** say
every layer of an odd-discriminant lattice is three-chromatic.

The upper bounds follow above. For the odd case put `n=(d+1)/4` and
`theta=(1+sqrt(-d))/2`; its norm form is

```
Q_n(x,y)=N(x+theta*y)=x^2+xy+n*y^2.
```

The three vectors

```
U=(n,0), V=(n-1,1), W=(-n,1)
```

all have norm `n^2`, and

```
n V - n W - (2n-1) U=0.                              (3)
```

Following `n` copies of `V`, then `n` copies of `-W`, then `2n-1` copies
of `-U` gives a simple polygon of odd length `4n-1=d`. Simplicity follows
directly: the first ascending segment has points `((n-1)j,j)`, the second
has points `(n(n-1)+nk,n-k)`, and the last lies on the x-axis; the two
positive-height segments meet only at their join.

The original arithmetic lattice, after the similarity used in (2), is a
rank-two Z-lattice `Lambda` in `K` spanning `K` over `Q`. Since
`(1,theta)` and any basis of `Lambda` are two rational bases, clearing
matrix denominators supplies an integer `M>=1` with
`M Z[theta] subset Lambda`. Multiplying the whole polygon by `M` therefore
puts it in this **same** lattice at one attained norm. It forces three
colours, proving the maximum claim. No shortest-vector condition is used.

Dividing every polygon vertex by `n` makes each edge have norm one inside
`K`. This normalized polygon and the local upper bound also give the exact
field-plane classification

```
chi G(Q(sqrt(-d)))=2 for d=1,2 mod4;
chi G(Q(sqrt(-d)))=3 for d=3 mod4.
```

This is recorded as an elementary consequence of the inherited mechanism,
not a priority claim.

## 5. A complete former residual: odd girth 23

Take `F=x^2+xy+6y^2`, `h=36`. This is the former hyperbolic-mod-three,
`9|h` branch: `Delta=-23=1 mod3`. Completing the square in either
coordinate gives `|x|<=6`, `|y|<=2`. Exhausting this tiny complete box gives
exactly

```
S={+/-U,+/-V,+/-W},
U=(6,0), V=(5,1), W=(-6,1).
```

The integer relation lattice on these three generators is generated by
`(-11,6,-6)`. Indeed the y equation says `c=-b`; the x equation becomes
`6a+11b=0`, so `(a,b,c)=t(-11,6,-6)` with integer `t`.
For any closed walk, its signed net generator counts have this form.
Walk-length parity equals `a+b+c mod2`, so an odd walk has odd `t` and
length at least `|a|+|b|+|c|>=23`. Equation (3), at `n=6`, realizes a simple
23-cycle. Hence the **odd girth is exactly 23**, not just at most 23.

At the ramified prime 23, normalization by `alpha=6` gives

```
ell(x,y)=4x+2y mod23;
ell(U)=1, ell(V)=-1, ell(W)=1.
```

Composing with the displayed colouring of `C23` gives an explicit proper
three-colouring, verified on the entire `23x23` quotient. Thus this norm
layer is exactly three-chromatic and triangle-free. It illustrates why
neither the triangular-lattice triangle nor a small mod-three linear
colouring is a complete account of the cap.

## 6. Hostile controls and correction anatomy

1. **False one-direction claim.** `Q=x^2+y^2+(sqrt(2)/2)xy` is positive
   definite and nonarithmetic. `Q=1` has exactly `+/-e1,+/-e2`: irrational
   and rational coefficients separate. Its graph is a square grid, not a
   union of paths. Surviving claim: bipartiteness, by the two-direction proof.
2. **False scaled-similarity step.** For `Q=x^2+6y^2`, mod-three null
   vectors lie in `x=0 mod3`. The restricted form is `9u^2+6v^2`, with
   minimum six. An index-three similar sublattice would scale squared
   lengths by three and have minimum three. It does not. Surviving claim:
   graph components lie in cosets of that sublattice; similarity does not
   follow. The field-prime argument repairs the cap without it.
3. **Incorrect level-set count.** In `F3`, `xy=1` has two points, one
   projective line, not four points/two lines. The avoided-kernel argument
   still works, but the old uniform count is false.
4. **Irrational angles do not prove lattice escape.** The existing
   2026-10-22 correction gives square-lattice vectors `(1,2),(2,1)` of
   equal norm with cosine `4/5`, hence angle an irrational multiple of pi.
   Only an infinite rotation-invariance requirement is excluded by a
   lattice's finite orthogonal symmetry group. Classical Dehn invariants
   need a polyhedral carrier, absent here.
5. **Arithmeticity has two meanings.** A graph of chromatic number at least
   four cannot lie in a single planar rank-two lattice. It can nevertheless
   have entirely algebraic coordinates in a number field: the Moser field
   is already an inherited counterexample to an 'off all arithmetic'
   conclusion. No rotation-rank formula or class-number ladder follows.

## 7. Actual bridge to Session60 and surviving limits

The hourglass root normals `+/-(1,0),+/-(0,1),+/-(1,1)` form a single
rank-two Eisenstein lattice for the norm `a^2-ab+b^2`. The repaired theorem
explains why retaining only those root-normal additive directions cannot
supply chromatic forcing above three, at **any** chosen norm scale.
This is a real obstruction on the target graph, not a match of numerals.

Connection contract: source = root-normal lattice; target = its fixed-norm
Cayley graph; map = the imaginary-quadratic embedding and normalized
ramified-prime quotient; preserved predicate = unit-edge incidence and
pullback of proper colourings; destroyed data = branch incidence, individual
square-root characters, cover singularities, polarizations and periods;
needed sidecar = the labelled root-cover equations from Session60. The
K3/Campedelli geometry does not turn this low-chromatic carrier into a
six- or seven-chromatic unit-distance configuration.

No Hadwiger–Nelson resolution, new planar bound, magic-square integer-point
claim, or rational descent of the CM/Kummer identification is asserted.

## 8. Verification and sources

Run:

```
python3 verify_lattice_cap.py > normal.json
python3 -O verify_lattice_cap.py > optimized.json
cmp normal.json optimized.json
```

The script uses only the Python standard library, exact integers and
fractions, and explicit checks that survive optimization. The exact
universe and all counts are printed in JSON: all positive-definite integer
forms `1<=A<=5,-5<=B<=5,1<=C<=6`, every attained `1<=h<=50`, and complete
connection sets via proven coordinate bounds. It checks norm normalization,
Eisenstein uniformizers, local coefficient integrality, unit residues,
finite cyclic step certificates, periods, all small Veronese determinant
controls, and the three negative controls. This census tests the formula;
the infinite theorem is the proof in Sections 2–4.

Primary source inspected 2026-10-25:
David A. Madore, *The Hadwiger–Nelson problem over certain fields*,
arXiv:1509.07023v1 (23 September 2015), Proposition 3.2, printed pp.6–7,
https://arxiv.org/pdf/1509.07023 . Its valuation-ring reduction principle is
classical; the present norm-one ramified specialization is proved above
rather than silently attributing every statement here to that proposition.
Madore's 2015 planar lower-bound introductory status is historical, not a
current-status citation.

Repository correction dependencies: MISTAKE-576/579 refute the old
Heegner/class-number chromatic roadmap and record the classical provenance
of the local-colouring lemma; the 2026-10-22 irrational-angle warning
corrects the Dehn/lattice-escape gloss. THM-4628's quotient-chromatic
obstructions must not be reversed into real-plane lower bounds.

A bounded literature search did not locate a primary source stating exactly
Theorems A/B in this lattice language. That is **not** novelty evidence.
The mathematical contribution to this session is a repaired route, exact
constructive certificates, and the explicit former-residual odd-girth-23
control. The local-colouring mechanism is inherited/classical.
