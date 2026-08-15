---
id: THM-3406
title: "Affine-modification power jets and principal-part transgression"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY DERIVATION-AUDITED.  For any field
  K, nonzero g in K[x], and P=f(x)+g(x)z, write B=K[P,x]=K[x,t],
  t=P-f=gz, I=(g,t), R=K[x,z]=B[t/g], S=R[g^-1], and Pi=S/R.  Then
  {a in B:a/g^q in R}=I^q, so Q_q=(g^-q B+R)/R is B/I^q, Pi is their
  injective direct limit under multiplication by g, and
  Q_q/Q_(q-1)=B/(g,t^q).  Factorization of g gives exact CRT triangular
  jet arms.  In characteristic zero, under (P_x,g)=R, the canonical
  localized Bezout class eta lies in ker(Dbar); it does not map to the
  divergence class mu by Dbar.  The connecting homomorphism of the
  two-term localization sequence is injective and sends eta to mu, giving
  K[P]eta ~= K[P]mu ~= K[T]/(chi).  For nonconstant g, eta has exact
  filtration depth two.  If psi is the minimal polynomial of f modulo
  rad(g), then K[P]eta intersect Q_1=psi K[P]eta and there is an exact,
  generally nonsplit support/thickness sequence
  0 -> K[T]/(chi/psi) -> K[T]/chi -> K[T]/psi -> 0.  Characteristic p,
  mixed-generator, radical, and duplicate-factor hostiles are sharp.  This
  is an integral response/principal-parts theorem and proves no new JC(2)
  case.
source: root-2608-jc-affine-modification-2026-08-15
audit: q8-multiblocker derivation independently reconstructed; exact arbitrary-q coefficient proof; two-term-complex type audit; CRT, characteristic, collision, and filtration hostile replay
depends_on:
  - THM-3386-linear-z-canonical-divergence-minimal-polynomial-collision-law
related:
  - THM-3404-factorized-danielewski-principal-parts-and-finite-cover-obstruction
  - THM-3348-linear-z-generic-puncture-response-and-one-root-valuation
  - THM-3326-linear-in-z-unit-response-trichotomy-and-jet-torsion
  - THM-2063-one-fiber-linear-planar-keller-pairs
script: 04-computation/jc_affine_modification_transgression_thm3406.py
output: 05-knowledge/results/jc_affine_modification_transgression_thm3406.out
script_sha256: ead69cdc09f73f0f4d897f9f381c21ee9fa19a908455f2884018f8416fd98723
output_sha256: ce508e7ab787e9ee385c9c335b789652ae8415fc45f8b6e30d81f75763720049
semantic_sha256: 065f584f197ff000d9f42175ffa858de223ee42551b19d996a54da37cfcb88f7
hash_basis: LF-normalized bytes
---

# THM-3406 -- affine-modification power jets and principal-part transgression

**PROVED + VERIFIED-EXACT + INDEPENDENTLY DERIVATION-AUDITED.**

## 1. Inheritance, corrected map, and connection contract

[THM-3386](THM-3386-linear-z-canonical-divergence-minimal-polynomial-collision-law.md)
computes the exact annihilator of the canonical divergence class for a
linear-in-`z` Hamiltonian.  Its mechanism is a localized rational primitive;
the canonical hostile is [MISTAKE-374](../MISTAKES.md), where generic
exactness fails to produce an integral polynomial mate.

[THM-3404](THM-3404-factorized-danielewski-principal-parts-and-finite-cover-obstruction.md)
shows that the correct sidecar for a localized boundary is a full
principal-part module, not only a divisor class or pole order.  Here the same
theme appears in a different operation: the polynomial plane is an affine
modification of the `(P,x)`-plane.  Its principal parts form a filtered tower,
not THM-3404's Laurent-degree direct sum.

The first tempting map is incorrectly typed.  If a localized primitive `h`
satisfies `D(h)=m in R`, then its quotient class `eta=[h]` satisfies
`Dbar(eta)=0`; it is not sent to `[m]` by `Dbar`.  The correct map is the
connecting homomorphism of a short exact sequence of two-term complexes.

| item | exact content |
|---|---|
| source | the affine modification `R=B[I/g]` and the localized complex `[R -> R]` |
| target | the filtered principal parts `Pi=S/R` and integral response cokernel `C_P` |
| map | denominator filtration `Q_q`, then the connecting map `delta([h])=[D(h)]` |
| preserved | denominator depth, factor arm, truncated coefficient jet, `K[P]`-annihilator, and collision maximum |
| lost by the associated graded | extension glue between reduced root-value support and multiplicity thickness |
| required sidecar | the filtered module `Pi`, not only `gr(Pi)`, together with `delta` |
| cheapest decisive tests | `I^2=(g^2,gt,t^2)`, `g=x^2`, and `P=x+x^p z` over `F_p` |

## 2. Exact affine-modification filtration over every field

Let `K` be any field, put

```text
R=K[x,z],              P=f(x)+g(x)z,              0!=g in K[x],       (1)
t=P-f(x)=g(x)z,        B=K[P,x]=K[x,t],           I=(g,t) in B.       (2)
```

The equality `B=K[x,t]` follows because `P,x` are algebraically independent:
the highest `P`-degree of a proposed relation has a nonzero highest
`z`-coefficient divisible by the corresponding power of `g`.  Inside
`B[g^(-1)]`,

```text
R=B[t/g]=B[I/g],
S=R[g^(-1)]=B[g^(-1)],
Pi=S/R.                                                               (3)
```

For every `q>=0`, define

```text
N_q={a in B:a/g^q belongs to R},
Q_q=(g^(-q)B+R)/R  subset Pi.                                        (4)
```

Write a unique polynomial

```text
a=sum_(j>=0) a_j(x)t^j.                                               (5)
```

After substituting `t=gz`, the coefficient of `z^j` in `a/g^q` is
`a_j g^(j-q)`.  It is polynomial exactly when

```text
g^(q-j) divides a_j       for every j<q.                              (6)
```

But `(6)` is also the coefficient criterion for membership in

```text
I^q=(g,t)^q=(g^q,g^(q-1)t,...,gt^(q-1),t^q).                          (7)
```

Therefore, for every `q`,

```text
N_q=I^q,
B/I^q  -> Q_q,             a mod I^q |-> [a/g^q]                     (8)
```

is an isomorphism.  In particular `Q_0=0`.  Replacing denominator `g^q`
by `g^(q+1)` gives the injective transition

```text
B/I^q -> B/I^(q+1),                  [a] |-> [ga],                    (9)
```

because the same coefficient test proves

```text
(I^(q+1):g)=I^q.                                                      (10)
```

Every element of `S` acquires a numerator in `B` after increasing its
`g`-denominator past its `z`-degree.  Hence

```text
Pi=union_(q>=0) Q_q
  =colim_q (B/I^q, multiplication by g).                              (11)
```

Under `(8)`, the image of `Q_(q-1)` in `Q_q` is `(gB+I^q)/I^q`.
Since

```text
I^q+(g)=(g,t^q),                                                      (12)
```

one obtains, for every `q>=1`,

```text
Q_q/Q_(q-1) = B/(g,t^q).                                              (13)
```

Thus the associated graded remembers the reduced denominator support at
each height, while the transition maps in `(9)` remember how the heights
are glued.

If `d=deg(g)`, the coefficient normal form gives

```text
dim_K Q_q=d q(q+1)/2,
dim_K(Q_q/Q_(q-1))=dq.                                                (14)
```

Equivalently, the filtration and its associated graded have Hilbert series

```text
sum_(q>=1) dim_K(Q_q)y^q=d y/(1-y)^3,
sum_(q>=1) dim_K(Q_q/Q_(q-1))y^q=d y/(1-y)^2.                         (14a)
```

The formulas include a constant unit `g`: then `I=B`, `S=R`, and every
module in `(8)--(14)` is zero.

## 3. CRT triangular jet arms

Aggregate the factorization over `K` as

```text
g=c product_(i=1)^r p_i^(e_i),
c in K^*,       p_i monic pairwise distinct irreducibles.             (15)
```

Put `J_i=(p_i^(e_i),t)`.  These ideals are pairwise comaximal and

```text
I=intersection_i J_i,
I^q=intersection_i J_i^q.                                            (16)
```

Indeed `B/I=K[x]/(g)` is the CRT product of the `B/J_i`, and powers of
pairwise comaximal ideals remain comaximal.  Consequently

```text
B/I^q = direct_sum_i B/J_i^q,                                        (17)
```

and each factor arm has the triangular `K[x]`-module normal form

```text
B/J_i^q
 = direct_sum_(j=0)^(q-1)
   t^j K[x]/(p_i^(e_i(q-j))).                                        (18)
```

On arm `i`, transition `(9)` is multiplication by `p_i^(e_i)` times a
unit.  Formula `(13)` similarly splits as

```text
B/(g,t^q)
 =direct_sum_i K[x,t]/(p_i^(e_i),t^q).                               (19)
```

The arms are indexed by aggregated irreducible factors.  Repeating the same
factor does not create a second arm.

## 4. The connecting homomorphism, not the quotient derivation

Assume from now on that `char(K)=0` and that the gradient row is unimodular:

```text
(P_x,g)=R.                                                           (20)
```

Let

```text
D=D_P=P_x partial_z-g partial_x,                C_P=R/D(R).           (21)
```

Choose a polynomial Bezout row

```text
A P_x+Cg=1                                                           (22)
```

and define

```text
m=A_x+C_z,       mu=[m] in C_P,
h=-A/g in S,     eta=[h] in Pi.                                      (23)
```

Any second row differs by `(rg,-rP_x)`.  It changes `h` by `-r in R` and
changes `m` by `-D(r)`.  Thus both `eta` and `mu` are canonical.
Differentiating `(22)` in `z` gives

```text
A_z P_x+A g'+C_z g=0,                                                (24)
```

and direct substitution yields

```text
D(h)=m.                                                              (25)
```

In particular, the induced derivation on the quotient has

```text
Dbar(eta)=0,                                                         (26)
```

not `mu`.  Consider instead the short exact sequence of two-term
`K[P]`-linear complexes

```text
0 -> [R --D--> R] -> [S --D--> S] -> [Pi --Dbar--> Pi] -> 0.          (27)
```

The connecting homomorphism is explicitly

```text
delta:ker(Dbar) -> C_P,
delta([s])=[D(s)].                                                    (28)
```

It is well-defined because a different lift changes `D(s)` by an element of
`D(R)`.  In characteristic zero, holding `P` fixed in
`S=K[P,x,g^(-1)]` gives

```text
D=-g partial_x,
ker(D:S->S)=K[P].                                                     (29)
```

To see the last equality, first take constants in `K(P)(x)`, then use the
UFD identity

```text
K(P) intersect K[P,x,g^(-1)]=K[P].                                  (30)
```

The kernel in `R` is the same `K[P]`.  The first kernel map in the long
exact sequence of `(27)` is therefore an isomorphism, so

```text
delta is injective,
im(delta)=ker(C_P -> S/D(S)),
delta(eta)=mu.                                                       (31)
```

For `F in K[T]`, the Bezout identity makes `A` a unit modulo `g`, and
`P=f mod g`.  Therefore

```text
F(P)eta=0
iff F(P)h belongs to R
iff g divides F(P)
iff g divides F(f).                                                   (32)
```

Let `chi` be the monic generator of the evaluation kernel

```text
ker(K[T] -> K[x]/(g)),                 F |-> F(f mod g).              (33)
```

For constant `g`, set `chi=1`.  Equations `(31)--(33)` give canonical
isomorphisms

```text
K[P]eta  --delta (isomorphism)--> K[P]mu,
K[P]eta ~= K[P]mu ~= K[T]/(chi).                                    (34)
```

This recovers THM-3386's annihilator while locating its rational primitive
inside the exact localization sequence.  The quotient differential only
certifies that `eta` is a cycle; `delta` performs the integral transgression.

## 5. Exact depth two and the support/thickness sequence

Assume `g` is nonconstant.  In characteristic zero, `(20)` is equivalent to

```text
g=c product_i p_i^(e_i) with every e_i>=2,
gcd(f',g)=1.                                                          (35)
```

Indeed `f'+g'z` is a unit in `(K[x]/(g))[z]` exactly when `f'` is a unit
and `g'` is nilpotent modulo `g`.  Under `(35)`, `g` divides `(g')^2`.

Choose `u in K[x]` with

```text
u f'=1 mod g                                                         (36)
```

and take the first Bezout coordinate

```text
A=u-u^2 g'z.                                                         (37)
```

Expanding `A(f'+g'z)` and using `(35)--(36)` proves it is `1 mod g`, so
some polynomial `C` completes `(22)`.  In the affine-modification coordinate
`t=gz`,

```text
h=-A/g=a/g^2,
a=-ug+u^2g't.                                                        (38)
```

Thus `eta in Q_2`.  Its symbol under `(13)` is

```text
sigma_2(eta)=u^2g't in B/(g,t^2).                                   (39)
```

It is nonzero because `u` is a unit modulo `g` and a nonconstant
characteristic-zero polynomial never divides its derivative.  Hence

```text
eta belongs to Q_2 but not Q_1.                                     (40)
```

Let

```text
s=rad(g)=product_i p_i,
psi=monic generator of ker(K[T] -> K[x]/(s)),
                              F |-> F(f mod s).                       (41)
```

Since `P=f+t`, multiplying `(39)` by `F(P)` gives the top symbol

```text
u^2 g' F(f)t        in B/(g,t^2).                                   (42)
```

For every factor `p_i`,

```text
ord_(p_i)(g')=e_i-1.                                                 (43)
```

Therefore

```text
F(P)eta belongs to Q_1
iff g divides g'F(f)
iff s divides F(f)
iff psi divides F.                                                    (44)
```

Since `(chi)` is contained in `(psi)`, `psi` divides `chi`.  Combining
`(34)` and `(44)` gives

```text
K[P]eta intersect Q_1=psi(P)K[P]eta                                 (45)
```

and the exact sequence

```text
0 -> K[T]/(chi/psi) --times psi--> K[T]/(chi)
  -> K[T]/(psi) -> 0.                                                (46)
```

The injection in `(46)` follows from cancellation in the PID `K[T]`.
In fact `(46)` never splits for nonconstant `g` under `(20)`.  Indeed
`rad(chi)=psi`, and every irreducible factor `rho` of `psi` occurs in `chi`
with exponent `M>=2`.  Localizing a hypothetical splitting at `(rho)` would
split the uniserial module `K[T]_(rho)/(rho^M)` into its length `M-1` and
length-one pieces, impossible because its `rho`-socle has dimension one while
the proposed direct sum has socle dimension two.  For `f=x,g=x^2`, this is
the standard nonsplit sequence with middle term `K[T]/(T^2)`.

Over an algebraic closure, if `alpha` runs through roots of `g`, put

```text
M_b=max{ord_alpha(g):f(alpha)=b}.                                    (47)
```

Then THM-3386's collision law and `(41)` become

```text
chi(T)=product_b (T-b)^(M_b),
psi(T)=product_b (T-b),
chi(T)/psi(T)=product_b (T-b)^(M_b-1).                               (48)
```

Thus the top quotient in `(46)` is reduced root-value support and the
submodule is multiplicity thickness.  Coincident factor arms combine by
least common multiple, hence by maximum multiplicity, while the nonsplit
middle term retains their infinitesimal glue.

## 6. Sharp hostiles and failure boundaries

### 6.1 The mixed generator at `q=2`

For every nonunit `g`,

```text
I^2=(g^2,gt,t^2).                                                     (49)
```

The term `gt` cannot be discarded:

```text
gt/g^2=t/g=z belongs to R,
gt does not belong to (g^2,t^2).                                    (50)
```

Nor may `I^2` be replaced by `I`, since

```text
g/g^2=1/g does not belong to R.                                     (51)
```

For `g=x^2`, replacing `g` by its radical is also false:
`x/g=1/x` is not polynomial.

### 6.2 Aggregation before CRT

For `g=x^2`, the local algebra `K[x]/(x^2)` has only the idempotents `0,1`.
Writing `x^2=x*x` therefore does not make two arms.  CRT splitting in
`(16)--(19)` requires distinct irreducible factors with their multiplicities
already aggregated.

### 6.3 Characteristic `p` kills injectivity

Over `K=F_p`, take

```text
P=x+x^p z,             g=x^p,             A=1, C=0.                  (52)
```

Then `(P_x,g)=(1)`, but

```text
h=-x^(-p) notin R,       D(h)=0,       m=mu=0.                       (53)
```

The class `eta` is nonzero with annihilator `(T^p)`, yet `delta(eta)=0`.
Here `ker(D:S->S)` strictly contains `K[P]` because Frobenius creates the
extra constant `x^(-p)`.  This is the sharp failure of `(29)--(34)` in
positive characteristic.  The affine-modification identities
`(3)--(19)` themselves remain valid.

### 6.4 Scope at the JC boundary

The theorem computes an integral sidecar for the already-tame
linear-in-`z` stratum; THM-2063 supplies that tameness.  It does not show
that an arbitrary planar Keller pair admits presentation `(1)--(3)`, that a
terminal filtration is an ordinary power filtration rather than a saturation
or integral closure, or that the response class controls an inverse cover.

To transfer the mechanism to a nonmonomial terminal branch one must still
prove an actual affine-modification chart, identify its physical denominator,
retain all boundary and infinity arms, and show that the connecting class is
the relevant Keller obstruction.  No new `JC(2)`, `DC(2)`, `A4/S4`, or
nonmonomial terminal case follows here.

## 7. Exact companion

The standard-library companion uses rational sparse polynomial arithmetic
and no floating point or assertion-dependent truth gate.  It checks `4,026`
denominator/ideal comparisons and transition colons, `30` length cells,
`150` CRT arm cells, five split/nonsplit/collision response geometries, every
small ternary annihilator/top-symbol polynomial in the displayed ranges, and
four Frobenius characteristics.  It freezes `(49)--(53)`, the radical error,
and the duplicate-arm idempotent hostile.  Normal and optimized outputs are
byte-identical.

Reproduce with

```text
python3 04-computation/jc_affine_modification_transgression_thm3406.py
python3 -O 04-computation/jc_affine_modification_transgression_thm3406.py
```

Artifact and semantic hashes are pinned in the frontmatter.

**QED.**
