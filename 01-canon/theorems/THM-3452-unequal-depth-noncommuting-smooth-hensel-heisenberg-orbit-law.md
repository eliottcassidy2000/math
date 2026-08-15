---
id: THM-3452
title: "Unequal-depth noncommuting smooth Hensel Heisenberg orbit law"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For two
  congruence-to-identity automorphisms at depths c1,c2, the theorem proves
  the weighted abelian range through c1+c2 and the mixed-modulus class-two
  range through c1+c2+min(c1,c2), with exact orbit-bank tariffs and a sharp
  higher-bracket boundary.  The depth-one dyadic channel has a separate
  square-carry iff criterion.
source: root-unequal-depth-noncommuting-hensel-2026-08-15
audit: >
  independent pullback/right-action sign, filtration, mixed-cocycle,
  normal-form, freeness, dependence, bank, dyadic across-range iff, A3/A4,
  and G_m hostile derivations; independent exact-control reconstruction;
  normal/optimized/stored replay; dependency, hash, semantic, AST/security,
  ID, scope, routing, diff, and documentation gates clean after the
  MISTAKE-400 provisional-status repair
depends_on:
  - THM-3442-smooth-hensel-fibre-vector-field-orbit-law
related:
  - THM-3446-weighted-depth-commuting-hensel-lattice-action
  - THM-3449-noncommuting-smooth-hensel-heisenberg-orbit-law
script: 04-computation/unequal_depth_noncommuting_hensel_heisenberg_orbits_thm3452.py
output: 05-knowledge/results/unequal_depth_noncommuting_hensel_heisenberg_orbits_thm3452.out
script_sha256: 9918d0edff5f16eccd64a3b81cf9685bccbc3917aefb3fe9f1ba8d3c0c0bfb50
output_sha256: 020803b42eee75331f70a5fe02be5251d85533e3a1c4dd6343a513ff285d9e70
semantic_sha256: c0c00380198dca87af2940d7986199885633b9a042aeae25b557656abadb6d4b
hash_basis: LF-normalized bytes
---

# THM-3452 -- unequal-depth noncommuting smooth Hensel Heisenberg orbit law

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The proof and exact companion passed an independent immutable-file derivation,
hostile, replay, dependency, security, and routing audit.

## 1. Setup, convention, and universal filtration

Let `p` be prime.  Let `X` be a smooth affine finite-type `Z_p`-scheme of
pure relative dimension `d>=1`.  Let

```text
g,h in Aut_(Z_p)(X)                                    (1)
```

be scheme-theoretically equal to the identity modulo `p^c1` and `p^c2`,
respectively, where `c1,c2>=1`.  On the coordinate ring write

```text
G=g^*=I+p^c1 D,       H=h^*=I+p^c2 E.                  (2)
```

Put

```text
C=min(c1,c2),       M=max(c1,c2),       T=c1+c2,       (3)

delta_1=D mod p,    delta_2=E mod p,
B=[delta_1,delta_2]=delta_1 delta_2-delta_2 delta_1.   (4)
```

The declared depths need not be exact.  If an automorphism is actually
deeper, its displayed first carry is zero and the independence hypotheses
below cannot hold.

Use the pullback commutator

```text
K=[G,H]=G H G^(-1) H^(-1).                             (5)
```

For a point `x`, write `x dot G=g(x)` and `x dot H=h(x)`.  Pullback
multiplication induces the right action

```text
x dot (G H)=(x dot G) dot H.                            (6)
```

Consequently the geometric automorphism represented by `(5)` is

```text
k=h^(-1) g^(-1) h g.                                   (7)
```

This convention fixes both the commutator sign and the cocycle below.

For every prime, including the dyadic depth-one case,

```text
K=I+p^T(D E-E D) mod p^(T+C),                          (8)
[G,K]=I mod p^(T+c1),       [H,K]=I mod p^(T+c2).      (9)
```

Thus the first possible commutator carry is `B` at depth `T`, and every
triple commutator vanishes through depth `T+C`.  These are sharp universal
bounds: a particular pair can have a deeper commutator when the displayed
carry vanishes.

## 2. Weighted abelian range

For every

```text
M<=a<=T,                                                (10)
```

the action on a depth-`C` fibre factors through

```text
A_a=Z/p^(a-c1)Z x Z/p^(a-c2)Z,                         (11)
```

where a cyclic factor of order one is understood as trivial.  This
factor-through assertion holds for every prime.

Assume now

```text
p is odd, or C>=2.                                     (12)
```

Fix `xbar in X(F_p)`.  If

```text
delta_1(xbar), delta_2(xbar) are independent,          (13)
```

then `(11)` acts freely on every depth-`C` fibre above `xbar`.  Such a fibre
has

```text
fibre size       =p^(d(a-C)),
every orbit size =p^(2a-T),
number of orbits =p^((d-2)(a-C)+M-C).                  (14)
```

In particular, full visible tangent rank `d=2` still leaves the persistent
tariff `p^(M-C)` when the depths differ.  If `(13)` holds at every residue
point, `(14)` holds on every nonempty depth-`C` fibre throughout `(10)`.

## 3. Mixed-modulus class-two range

For

```text
T<a<=T+C,                                              (15)
```

put

```text
N_1=a-c1,       N_2=a-c2,       S=a-T.                 (16)
```

Define the finite class-two group

```text
H_a(c1,c2;p)=<x,y,z |
  x^(p^N_1)=1, y^(p^N_2)=1, z^(p^S)=1,
  z=x y x^(-1) y^(-1), [x,z]=[y,z]=1>.                (17)
```

Every abstract element has a unique normal form

```text
x^i y^j z^r,
i mod p^N_1,       j mod p^N_2,       r mod p^S,       (18)
```

and normal-form triples multiply by

```text
(i,j,r)(i',j',r')
 =(i+i',j+j',r+r'-j i').                              (19)
```

The coordinates in `(19)` are reduced modulo their respective moduli.
Because `S<=N_1,N_2`, the cocycle is well-defined.  In particular,

```text
|H_a(c1,c2;p)|=p^(N_1+N_2+S)=p^(3a-2T).               (20)
```

Equations `(8)--(9)` and the exponent bounds give the pullback
representation

```text
x |-> G,       y |-> H,       z |-> K,                 (21)
```

and hence the corresponding right point action.  This structural
factor-through assertion again holds for every prime.

Under `(12)`, if at `xbar`

```text
delta_1(xbar), delta_2(xbar), B(xbar)
are linearly independent,                              (22)
```

then `(21)` is free on every depth-`C` fibre above `xbar`.  Consequently

```text
fibre size       =p^(d(a-C)),
every orbit size =p^(3a-2T),
number of orbits =p^((d-3)(a-C)+2M-C).                 (23)
```

Triple independence forces `d>=3`.  At full rank `d=3`, exactly

```text
p^(2M-C)                                               (24)
```

orbit banks persist.  The delayed central digit does not make the weighted
action transitive.

The first nonabelian level also gives the exact dependence test.  At
`a=T+1`, a relation

```text
alpha delta_1+beta delta_2+gamma B=0                   (25)
```

at `xbar`, with coefficients in `F_p` not all zero, produces the nonzero
abstract word

```text
x^(alpha p^c2) y^(beta p^c1) z^gamma.                 (26)
```

It fixes every lift from depth `T` to depth `T+1` over `xbar`.  Therefore
triple independence at `xbar` is equivalent to freeness on that fibre at
the first nonabelian level, and it implies freeness throughout `(15)`.

## 4. Dyadic shallow-depth repair

The remaining unequal-depth case is

```text
p=2,       C=1<M.                                      (27)
```

Relabel the generators so that `g_s` has depth one and `g_d` has depth `M`.
Write

```text
U=g_s^*=I+2A,       V=g_d^*=I+2^M E,                  (28)

delta_s=A mod 2,       delta_d=E mod 2,
B=[delta_s,delta_d].                                  (29)
```

The intrinsic shallow square carry is

```text
sigma_s=((U^2-I)/4) mod 2
       =delta_s+delta_s^[2],                           (30)
```

where `delta^[2]=delta compose delta` is the restricted square of a
characteristic-two derivation.

On a fixed depth-one fibre above `xbar`, consider every active exponent
action beginning at level two; before level `M`, only the shallow cyclic
channel is active.  Then the following are exact iff statements.

1. Every active action through level `M+1` is free iff

   ```text
   delta_s(xbar)!=0, and
   sigma_s(xbar), delta_d(xbar) are independent.       (31)
   ```

2. Every active action through level `M+2` is free iff

   ```text
   delta_s(xbar)!=0, and
   sigma_s(xbar), delta_d(xbar), B(xbar)
   are independent.                                   (32)
   ```

At the two-generator levels `M<=a<=M+1`, a free action has the abelian
counts `(14)` with `C=1`.  At `a=M+2`, a free action has the class-two
counts `(23)` with `C=1` and `T=M+1`.

The converse words make the iff directions explicit.  A dependence

```text
alpha sigma_s+beta delta_d=0                           (33)
```

at level `M+1` gives

```text
x^(alpha 2^(M-1)) y^beta.                              (34)
```

At level `M+2`, a dependence

```text
alpha sigma_s+beta delta_d+gamma B=0                   (35)
```

gives the nonzero word

```text
x^(alpha 2^M) y^(2 beta) z^gamma.                      (36)
```

The separate condition `delta_s!=0` is detected at level two by an odd
shallow exponent.  Thus `(31)--(32)` are across-level criteria, not a claim
that the repaired square alone remembers the first lift.

When `C=M=1`, both horizontal generators require square carries.  That is
the distinct THM-3449 criterion involving ordinary horizontal independence
and the triple `sigma_g,sigma_h,B`; it is not obtained by setting `M=1` in
`(31)--(32)`.

## 5. Proof candidate

### 5.1 Commutator and exponent filtration

Expand the inverse series only to the first layer after `T`.  Every term
beyond `p^T(DE-ED)` contains one extra factor of depth at least `C`, so

```text
G H G^(-1) H^(-1)
 =I+p^T(DE-ED) mod p^(T+C),                            (37)
```

proving `(8)`.  The general filtration rule

```text
[I+p^u R,I+p^v S]=I mod p^(u+v)                        (38)
```

applied again proves `(9)`.

For any depth-`c` operator `I+p^cR` and `n>=0`,

```text
(I+p^cR)^(p^n)=I mod p^(c+n).                          (39)
```

Indeed

```text
v_p(binomial(p^n,k))>=n-v_p(k),
c(k-1)>=v_p(k).                                        (40)
```

This proves all exponent relations, even for `p=2,c=1`.  Equations
`(8)--(9)` then prove the abelian and class-two factor-through ranges.

Under `(12)`, taking a `p`th power raises the depth by exactly one and keeps
the same nonzero first carry.  Thus, for a unit `q mod p`,

```text
G^(q p^n)=I+q p^(c1+n)D_n mod p^(c1+n+1),
D_n mod p=delta_1,                                    (41)
```

and similarly for `H` and for the depth-`T` operator `K`.

### 5.2 Mixed normal forms and freeness

The bilinear cocycle `-j i'` in `(19)` satisfies the cocycle identity, and
direct calculation gives `z=xyx^(-1)y^(-1)`.  This proves the abstract
normal forms and order `(18)--(20)` independently of the geometric
representation.

For a nonzero represented normal form, put `v_p(0)=infinity` and

```text
L=min(c1+v_p(i), c2+v_p(j), T+v_p(r))<a,               (42)
```

omitting the central term in the abelian range.  Formula `(41)` says that
the first visible carry is the nonzero coefficient combination of the
members of `(13)` or `(22)` tied at depth `L`.  Every cross term has depth
at least `2L>=L+1`; in particular, the word-area correction cannot cancel
the first-carry tie.  Independence therefore excludes a fixed point for
every nonidentity normal form.  Smoothness supplies the fibre cardinality,
and division by `(11)` or `(20)` gives `(14)` and `(23)`.

At `a=T+1`, all three factors of `(26)` have effective depth `T`.  Relation
`(25)` kills their first translation on the tangent torsor from depth `T`
to `T+1`, proving the dependence direction.

### 5.3 Dyadic normal forms

Direct squaring gives

```text
(I+2A)^2=I+4(A+A^2),                                   (43)
```

which proves `(30)`.  An odd shallow exponent is detected by `delta_s` at
depth one.  For an even exponent `i=2u`, the first possible carry comes from
`(U^2)^u` and is the unit coefficient of `sigma_s`.  At level `M+1`, this
can tie only with `delta_d`; at level `M+2`, the new tied coordinate is the
commutator carry `B`.  Conditions `(31)--(32)` therefore exclude every
nontrivial normal-form stabilizer.  Conversely, `(33)` and `(35)` align at
depths `M` and `M+1` and give `(34)` and `(36)`.  Level two separately
detects `delta_s`.  This proves both dyadic iff claims.

## 6. Depth-vector and central-digit ledger

For the safe range, the normal form exposes every scale explicitly.

| coordinate | modulus at level `a` | first visible depth | retained carry |
|---|---:|---:|---|
| `i` | `p^(a-c1)` | `c1+v_p(i)` | `delta_1` |
| `j` | `p^(a-c2)` | `c2+v_p(j)` | `delta_2` |
| `r` | `p^(a-T)` | `T+v_p(r)` | `B` |

In the dyadic depth-one lane, an odd `i` retains `delta_s`; an even `i=2u`
has visible depth `2+v_2(u)` and retains `sigma_s`.  This is the missing
coordinate that a flat ordinary-carry ledger destroys.

| field | exact content |
|---|---|
| source | two smooth Hensel automorphisms with ordered depth vector `(c1,c2)` |
| target | weighted abelian group `(11)`, then mixed Heisenberg group `(17)` |
| map | `x^i y^j z^r -> G^i H^j K^r`, inducing the right point action `(6)` |
| preserved | both horizontal exponent scales and oriented commutator area |
| destroyed by abelianization | word order and the central digit `r mod p^(a-T)` |
| destroyed by depth flattening | the tariffs `M-C` and `2M-C` |
| required sidecars | full depth vector, central digit, and an orbit-bank label |
| cheapest positive | unequal-depth affine `A^3` model in Section 7.1 |
| cheapest hostile | evaluate the bracket, then the shallow higher bracket or dyadic square carry |

## 7. Equality and sharp failure boundaries

### 7.1 Unequal-depth affine equality

Let `q=p^C`, `Q=p^M`, and on `A^3` take the shallow and deep maps

```text
g_s(x,y,z)=(x,y+q,z+q x),
g_d(x,y,z)=(x+Q,y,z).                                  (44)
```

With `G=g_s^*`, `H=g_d^*`,

```text
delta_s=partial_y+x partial_z,
delta_d=partial_x,
B=-partial_z,                                          (45)

k=g_d^(-1)g_s^(-1)g_d g_s:
  (x,y,z) |-> (x,y,z-qQ).                              (46)
```

The three carries are independent everywhere, and `k` is exactly central.
This model attains every orbit and bank invoice in `(14)` and `(23)`.

### 7.2 Sharp first higher-bracket boundary

On `A^4`, replace the shallow map by

```text
g_s(x,y,z,w)=(x,y+q,z+q x,w+q z),
g_d(x,y,z,w)=(x+Q,y,z,w).                              (47)
```

Exact composition gives

```text
k=g_d^(-1)g_s^(-1)g_d g_s:
 (x,y,z,w) |-> (x,y,z-qQ,w+q^2 Q).                    (48)
```

This commutator commutes with `g_d` exactly and with `g_s` modulo `p^a`
precisely through `a<=T+C`.  At level `T+C+1`, the `w`-defect is the unit
multiple

```text
q^2 Q=p^(T+C).                                         (49)
```

Thus the universal class-two endpoint is sharp for every prime.

### 7.3 Global dyadic positive control

For `M>1`, on `A^3/Z_2` take

```text
g_s(x,y,z)=(x,y+2,z+2y),
g_d(x,y,z)=(x+2^M y,y+2^M,z).                         (50)
```

On the special fibre,

```text
delta_s=(0,1,y),       sigma_s=(0,1,y+1),
delta_d=(y,1,0),       B=(1,0,1).                     (51)
```

The ordinary triple `delta_s,delta_d,B` has rank two at both values of `y`,
whereas `sigma_s,delta_d,B` has rank three and the repaired pair has rank
two.  Thus `(32)`, not ordinary triple independence, is the exact positive
criterion.  At level `M+2`, the companion verifies

```text
M=2:  512=8*64,       M=3:  4096=32*128               (52)
```

on every one of the eight depth-one fibres.

### 7.4 Dyadic square-collapse hostile

On `G_m x A^2/Z_2`, take

```text
g_s(x,y,z)=(-x,y,z),
g_d(x,y,z)=(x,y+2^M,z+2^M x).                         (53)
```

At every special point,

```text
delta_s=x partial_x,
delta_d=partial_y+x partial_z,
B=x partial_z                                         (54)
```

are independent, but `g_s^2=1` and `sigma_s=0`.  At level `M+2`, the
nominal group order is `2^(M+4)`, while the image orbit has size `16` and
the kernel has size `2^M`.  The exact rows are

```text
M=2:  512=32*16,       M=3:  4096=256*16.              (55)
```

This makes the repaired square carry indispensable.

## 8. Scope and non-consequences

The closest proved mechanisms are THM-3446 (commuting unequal depths) and
THM-3449 (noncommuting equal depths).  The canonical hostile is the
denominator-free `A^4` higher bracket `(47)`, the corrected near miss is an
ordinary dyadic triple without `sigma_s`, and the least-used sidecars are
the depth vector and delayed central digit.

The result concerns local Hensel fibres and finite nilpotent exponent
groups.  It supplies no LRC owner, physical time, boundary current, Keller
response, or transport functor.  No LRC(14), `JC(2)`, or other physical
consequence follows.

The same local proof applies to smooth `p`-adic formal schemes locally
topologically of finite presentation.  Global statements require the
displayed independence conditions at every residue point.

## 9. Exact companion

The deterministic standard-library companion checks three complete
mixed-modulus presentations, their orders, centres, inverses, commutators,
and cocycle identities.  It verifies seven unequal-depth affine `A^3` orbit
rows, including the abelian cutoff and class-two invoices, and freezes
`65,536` mixed right-action/product checks.  Four swapped `A^4` controls
have exact first higher-bracket defects `16,81,128,3125`.

The dyadic positive model is exhausted on all eight special fibres through
the claimed ranges for `M=2,3`; the ordinary and repaired ranks are computed
independently.  The `G_m` hostile exhausts `512` and `4,096` points and
recovers image size `16` and kernels `4,8`.  The script rejects `assert`,
dynamic evaluation, and dependency drift, and prints a frozen semantic
digest.

Reproduce with

```text
python3 -B 04-computation/unequal_depth_noncommuting_hensel_heisenberg_orbits_thm3452.py
python3 -B -O 04-computation/unequal_depth_noncommuting_hensel_heisenberg_orbits_thm3452.py
```

QED.
