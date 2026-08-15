---
id: THM-3449
title: "Noncommuting smooth Hensel pairs: abelian cutoff and Heisenberg orbit law"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Two equal-depth
  congruence-to-identity automorphisms
  commute universally through level 2c; their first commutator carry is the
  Lie bracket, and through level 3c the exact replacement is a class-two
  Heisenberg lattice.  The dyadic depth-one square carries and the first
  higher-bracket boundary are explicit.
source: root-noncommuting-hensel-heisenberg-2026-08-15
audit: >
  independent pullback/right-action sign and filtration proof; mixed-modulus
  matrix reconstruction of the presentation, normal forms, order, and centre;
  minimum-effective-depth, bank, dependence, dyadic restricted-square,
  A3/A4, and G_m hostile audits; normal/optimized/stored and clean-room exact
  replay; dependency, hash, AST/security, ID, routing, and documentation
  gates clean
depends_on:
  - THM-3442-smooth-hensel-fibre-vector-field-orbit-law
  - THM-3444-commuting-smooth-hensel-vector-field-lattice-action
related:
  - THM-3446-weighted-depth-commuting-hensel-lattice-action
script: 04-computation/noncommuting_smooth_hensel_heisenberg_orbits_thm3449.py
output: 05-knowledge/results/noncommuting_smooth_hensel_heisenberg_orbits_thm3449.out
script_sha256: b1ab908fc42dbd95824b0535ef3cbbaf5d55d68faf77a959a38ab6a348df790f
output_sha256: dfe2528430981840feac984a47a45ca726ae3ec05a52966ba263230ae6356005
semantic_sha256: fa4dbd16a9f56fb350388fde53c8b3ef39cdfbc6f3092cd31aed83c506bc8154
hash_basis: LF-normalized bytes
---

# THM-3449 -- noncommuting smooth Hensel pairs: abelian cutoff and Heisenberg orbit law

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The proof and exact companion passed an independent derivation, hostile, and
replay audit.

## 1. Exact equal-depth statement

Let `p` be prime and `c>=1`.  Let `X` be a smooth affine finite-type
`Z_p`-scheme of pure relative dimension `d>=1`, and let

```text
g,h in Aut_(Z_p)(X)                                    (1)
```

be scheme-theoretically equal to the identity modulo `p^c`.  On the
coordinate ring put

```text
G=g^*=I+p^c D,       H=h^*=I+p^c E.                    (2)
```

Their first-carry vector fields are

```text
delta_g=D mod p,       delta_h=E mod p.                 (3)
```

Use the operator commutator convention

```text
K=[G,H]=G H G^(-1) H^(-1),
B=[delta_g,delta_h]=delta_g delta_h-delta_h delta_g.    (4)
```

For a point `x` and a pullback operator, write

```text
x dot G=g(x),       x dot H=h(x).                       (5)
```

This is a **right** action:

```text
x dot (G H)=(x dot G) dot H.                            (6)
```

Consequently the geometric automorphism corresponding to `K` is

```text
k=h^(-1) g^(-1) h g.                                   (7)
```

This convention is load-bearing for the sign and word order below.

### A. Universal commutator carry and abelian cutoff

For every prime and every `c>=1`, with no dyadic exclusion,

```text
K=I+p^(2c)(D E-E D) mod p^(3c).                        (8)
```

Thus `K=id mod p^(2c)`, and its depth-`2c` first carry is exactly `B`.
Moreover, for every `a>=c`,

```text
G^(p^(a-c))=H^(p^(a-c))=I mod p^a.                     (9)
```

For `c<=a<=2c`, put `N=a-c`.  On every depth-`c` fibre the action
therefore factors through the abelian exponent group

```text
A_N=(Z/p^N Z)^2.                                       (10)
```

Fix `xbar in X(F_p)` and a lift `x_c`.  If

```text
delta_g(xbar), delta_h(xbar) are independent,           (11)
```

then `(10)` acts freely on the fibre above `x_c`, and

```text
fibre size       =p^(dN),
every orbit size =p^(2N),
number of orbits =p^((d-2)N).                          (12)
```

Condition `(11)` is exact at the first lift.  If it holds at every
`F_p`-point, `(12)` holds on every nonempty depth-`c` fibre throughout the
abelian range.

The endpoint `2c` is the sharp **universal** cutoff.  If

```text
B(xbar)!=0,                                             (13)
```

then at level `2c+1` the geometric commutator `(7)` translates the lifts
above every depth-`2c` point over `xbar` by `B(xbar)`.  Hence `g` and `h`
actually fail to commute on that fibre, so the naive abelian exponent action
does not factor there.  If `B(xbar)=0`, this first obstruction vanishes and
the particular fibre may remain abelian farther; no unconditional
pair-specific failure is asserted.

### B. The class-two Heisenberg replacement

Assume now

```text
p is odd, or c>=2.                                      (14)
```

For

```text
2c<a<=3c,       N=a-c,       S=a-2c,                   (15)
```

define

```text
H_(N,S)(p)=<x,y,z |
  x^(p^N)=y^(p^N)=1, z^(p^S)=1,
  z=x y x^(-1) y^(-1), [x,z]=[y,z]=1>.                 (16)
```

Every element has a unique abstract normal form

```text
x^i y^j z^r,
i,j mod p^N,       r mod p^S.                          (17)
```

Equivalently, triples multiply by

```text
(i,j,r)(i',j',r')
 =(i+i',j+j',r+r'-j i'),                               (18)
```

with horizontal and central coordinates reduced modulo `p^N` and `p^S`,
respectively.  In particular

```text
|H_(N,S)(p)|=p^(2N+S).                                 (19)
```

Since `K` has depth `2c`, the filtration commutator law gives

```text
[G,K],[H,K]=I mod p^(3c).                              (20)
```

Thus

```text
x |-> G,       y |-> H,       z |-> K                  (21)
```

defines a pullback representation of `(16)` and hence a right action on the
point fibre.

If at `xbar` the three vectors

```text
delta_g(xbar), delta_h(xbar), B(xbar)                  (22)
```

are linearly independent, the action `(21)` is free on the fibre above
`x_c`.  Consequently

```text
fibre size       =p^(dN),
every orbit size =p^(2N+S),
number of orbits =p^((d-2)N-S).                        (23)
```

Triple independence forces `d>=3`.  When `d=3`, the exponent in the last
line is

```text
N-S=c,                                                   (24)
```

so exactly `p^c` orbit banks persist.  The delayed central digit cannot make
the action transitive, even at full visible tangent rank.

Condition `(22)` is exact at the first nonabelian level.  At `a=2c+1`, a
dependence

```text
alpha delta_g+beta delta_h+gamma B=0                   (25)
```

at `xbar` produces the nonzero abstract element

```text
x^(alpha p^c) y^(beta p^c) z^gamma                     (26)
```

which fixes the whole fibre over `x_c`.  Therefore `(22)` at every special
point is equivalent to freeness of `(21)` on every depth-`c` fibre at level
`2c+1`, and implies freeness throughout `(15)`.

### C. The dyadic depth-one repair

At `p=2,c=1`, assertions `(8)--(13)` and the structural factor-through
`H_(2,1)(2)` at level three remain valid.  What fails is preservation of an
ordinary first carry under squaring.  Define the intrinsic depth-two square
carries

```text
sigma_g=((G^2-I)/4) mod 2=delta_g+delta_g^[2],
sigma_h=((H^2-I)/4) mod 2=delta_h+delta_h^[2],          (27)
```

where `delta^[2]=delta compose delta` is the restricted square of a
characteristic-two derivation.

On a fibre above `xbar`, the `H_(2,1)(2)` action is free at level three if

```text
delta_g,delta_h are independent, and
sigma_g,sigma_h,B are independent                      (28)
```

at `xbar`.  These two conditions are together equivalent to freeness at
both levels two and three.  Indeed, a normal form with an odd horizontal
exponent is detected at depth one.  If `i=2u` and `j=2v`, its depth-two carry
is exactly

```text
u sigma_g+v sigma_h+r B.                               (29)
```

Under `(28)`, a level-three fibre has `2^(2d)` points, free orbit size `32`,
and `2^(2d-5)` orbit banks.

The same local statements hold for smooth `p`-adic formal schemes locally
topologically of finite presentation.  Global statements require the
displayed independence at every residue point.

## 2. Proof

### 2.1 Filtration algebra

Modulo `p^(3c)`, the inverse series are finite:

```text
G^(-1)=I-p^cD+p^(2c)D^2,
H^(-1)=I-p^cE+p^(2c)E^2.                               (30)
```

Multiplying the four factors in `(4)` cancels the depth-`c` terms and leaves

```text
I+p^(2c)(D E-E D),                                     (31)
```

proving `(8)`.  Reducing `(31)` modulo `p` after division by `p^(2c)` gives
the Lie bracket of the derivations `(3)`.  More generally, automorphisms of
depths `u,v` have commutator depth at least `u+v`; this proves `(20)`.

For `(9)`, use

```text
(I+p^cD)^(p^N)=sum_k binom(p^N,k)p^(ck)D^k.             (32)
```

The elementary bound

```text
v_p(binom(p^N,k))>=N-v_p(k)                            (33)
```

and `c(k-1)>=v_p(k)` put every nonconstant term in depth at least `c+N`.
This proves exponent factor-through for every prime, including `p=2,c=1`.

Under `(14)`, raising a depth-`s` operator with nonzero first carry to its
`p`th power increases the depth by exactly one and retains that carry.  For
odd `p`, all higher binomial terms are one layer later.  For `p=2`, the
quadratic term is
one layer later once `s>=2`.  Thus for `t>=0` and a unit `q mod p`,

```text
G^(p^t q)=I+q p^(c+t)D_t mod p^(c+t+1),
D_t mod p=delta_g,                                     (34)
```

and similarly for `H`.  The depth-`2c` operator `K` obeys the same rule.

### 2.2 Abelian freeness

For a nonzero `(i,j) in A_N`, let

```text
L=min(c+v_p(i),c+v_p(j))<a.                            (35)
```

Formula `(34)` says that the first visible carry of `G^iH^j` at depth `L`
is the corresponding nonzero linear combination of `delta_g,delta_h`.
Condition `(11)` keeps it nonzero.  A fixed point would force that tangent
vector to vanish after evaluation and cancellation of `p^L`, a
contradiction.  At `p=2,c=1`, the only nontrivial abelian level is the first
lift `a=2`, where the same statement is the tangent-torsor translation law
of THM-3442.  Smoothness supplies the fibre count in `(12)`.

Finally, `(8)` makes `K` the first-lift translation by `B(xbar)` from depth
`2c` to `2c+1`.  This proves the fibrewise failure under `(13)` and the sharp
universal abelian boundary.

### 2.3 Presentation and normal-form uniqueness

The bilinear cocycle `-j i'` in `(18)` is well-defined modulo `p^S` because
`N>=S`, and its cocycle identity proves associativity.  Direct calculation
gives `z=xyx^(-1)y^(-1)` and the presentation `(16)`.  Hence every abstract
element has exactly one triple `(17)` before any geometric representation is
considered.

For a nonzero normal form in the representation, put `v_p(0)=infinity` and

```text
L=min(c+v_p(i), c+v_p(j), 2c+v_p(r))<a.                (36)
```

Using `(34)` for `G,H,K`, its depth-`L` carry is a nonzero coefficient
combination of the tied members of `(22)`.  Every product term involving two
positive-depth displacements lies in depth at least `2L>=L+1`.  In
particular, the central word-area correction in `(18)` cannot enter or cancel
the first-carry tie.  Independence `(22)` therefore makes the represented
normal forms distinct and excludes every stabilizer.  This proves
`(19)--(23)`.

If `(25)` holds, the three factors in `(26)` all have effective depth `2c`;
their first carries sum to zero.  The first-lift translation law then fixes
every lift from depth `2c` to `2c+1` above `xbar`.  This proves the converse.

### 2.4 Dyadic normal forms

When `p=2,c=1`, direct squaring gives

```text
(I+2D)^2=I+4(D+D^2),                                   (37)
```

which proves `(27)`.  At level three, a normal form with at least one odd
horizontal exponent has the ordinary depth-one carry.  With both horizontal
exponents even it has `(29)`; all cross terms are already zero modulo eight.
This proves `(28)` and its exact two-level converse.

## 3. Inheritance and connection ledger

THM-3442 supplies the one-generator first-carry translation law.  THM-3444
shows that pairwise commutation plus independent carries gives a free
lattice action.  The present theorem identifies exactly how long
commutation is automatic and restores the first lost coordinate rather than
discarding it.  THM-3446 is the proved unequal-depth commuting analogue; no
unequal-depth noncommuting theorem is imported here.

The canonical hostile to the near miss “independent first carries are enough”
is a pair with independent `delta_g,delta_h` but nonzero bracket.  The
corrected near miss also respects MISTAKE-249: in a noncommutative setting an
ordered word is not an abelian or two-sided exponent congruence.  The
least-used sidecar is the delayed central word-area digit, followed by the
remaining orbit-bank label.

| field | exact content |
|---|---|
| source | two equal-depth smooth Hensel automorphisms |
| target | `A_N` through depth `2c`, then `H_(N,S)(p)` through `3c` |
| map | the pullback word `x^i y^j z^r -> G^i H^j K^r`, inducing a right point action |
| preserved | base point, horizontal exponents, and oriented commutator area |
| destroyed by abelianization | word order and the nontrivial central action |
| required sidecars | `r mod p^S` and one of the orbit-bank labels in `(23)` |
| cheapest positive | the affine Heisenberg pair in Section 4.1 |
| cheapest failure test | evaluate `B(xbar)`, then the triple carry matrix |

## 4. Equality and sharp boundaries

### 4.1 Exact affine Heisenberg equality

On `A^3`, with `q=p^c`, take

```text
g(x,y,z)=(x+q,y,z),
h(x,y,z)=(x,y+q,z+q x).                                (38)
```

Then

```text
delta_g=partial_x,
delta_h=partial_y+x partial_z,
B=partial_z,                                            (39)
```

and the geometric commutator is exactly `z |-> z+q^2`, central at every
level.  All normal forms are distinct by their `x,y,z` displacements.  The
model attains every invoice in `(12)` and `(23)`; for `a>2c` it has exactly
`p^c` orbit banks even beyond the theorem's universal class-two range.

### 4.2 The universal `3c` class-two range is sharp

On `A^4`, take the denominator-free triangular pair

```text
g(x,y,z,w)=(x+q,y,z,w),
h(x,y,z,w)=(x,y+q,z+q x,w+q z).                        (40)
```

Its first fields are

```text
X=partial_x,
Y=partial_y+x partial_z+z partial_w,
[X,Y]=partial_z,
[Y,[X,Y]]=-partial_w.                                  (41)
```

The first three are independent everywhere.  Exact composition gives

```text
k=h^(-1)g^(-1)hg:(x,y,z,w)
  |->(x,y,z+q^2,w-q^3).                                (42)
```

This `k` commutes with `g` exactly and with `h` modulo `p^a` precisely
through `a<=3c`; at level `3c+1` their `w`-coordinate defect is exactly
`q^3`.  Thus a higher bracket genuinely replaces the Heisenberg law one
level beyond the claimed range, for every prime.

### 4.3 Dyadic ordinary-carry hostile

On `G_m x A^2` over `Z_2`, take

```text
g(x,y,z)=(-x,y,z),
h(x,y,z)=(x,y+2,z+2x).                                 (43)
```

On the special fibre,

```text
delta_g=x partial_x,
delta_h=partial_y+x partial_z,
B=x partial_z                                          (44)
```

are independent, but `g^2=1`, so `sigma_g=0`.  On the depth-one fibre over
`(1,0,0)`, level two has eight points and two free orbits of size four.
Level three has `64` points, but the nominal order-`32` Heisenberg group has
kernel `x^2`; the image has four orbits of size `16`.  This makes the square
carry in `(27)--(29)` indispensable.

### 4.4 Scope boundary

The unequal-depth noncommuting thresholds suggested by the same filtration
are `c_1+c_2` and `c_1+c_2+min(c_1,c_2)`.  They are a roadmap, not part of
this equal-depth theorem and not a consequence asserted from THM-3446.

A Hensel word-area digit is not an LRC owner, physical time, boundary
current, or Keller response.  No LRC(14), `JC(2)`, or boundary-response
consequence follows.

## 5. Exact companion

The standard-library companion checks the abstract presentation, exact
orders, centre sizes, and complete bilinear cocycle identities for three
small `H_(N,S)(p)` groups.  Independent normal-form images and permutation
orbits verify seven affine equality rows, including `65,536` exhaustive
right-action/product checks.  It exhausts the level-two and
level-three `A^4` orbit controls at `p=2,3`, checks exact `3c+1` defects
`8,27,125,64`, and exhausts the dyadic `G_m` hostile.  It prints orbit and
bank counts—the theorem's consequence objects—not only injected formulas.

Reproduce with

```text
python3 -B 04-computation/noncommuting_smooth_hensel_heisenberg_orbits_thm3449.py
python3 -B -O 04-computation/noncommuting_smooth_hensel_heisenberg_orbits_thm3449.py
```

QED.
