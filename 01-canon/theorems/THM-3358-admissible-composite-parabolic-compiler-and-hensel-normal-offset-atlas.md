---
id: THM-3358
title: "Admissible composite parabolic compiler and Hensel-normal offset atlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every admissible
  odd composite grade has a centered Hensel-normal atlas, and every exact-grade
  normal lift has an explicit fixed-prefix/unary/fixed-suffix Berggren
  compiler.  The jet gauge gives the intrinsic first-normal unit bundle;
  unitary divisors compile the oriented allocation star, while determinant
  charge and parity sharply obstruct a universal affine-shell unary compiler.
source: codex-2026-08-14-composite-parabolic-compiler
audit: >
  Independent arithmetic and word-order audits reconstructed the centered
  Hensel lift, normal-jet gauge, full-alpha compiler, Gaussian content pair,
  unitary-divisor boundary, affine-shell extension, and determinant-charge
  no-go.  They required the primitive-allocation hypothesis, an explicit
  coarse exact-grade root in the affine no-go, the alpha/jet gauge split, and
  ramified, partial-power, charge, and parity hostiles.  Ordinary, optimized,
  and stored exact transcripts are byte-identical after those repairs.
depends_on:
  - THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation
  - THM-3346-u-spine-prime-toggle-root-atlas-and-conjugation-monodromy
  - THM-3353-split-prime-parabolic-branch-transplant-and-unary-transducer-compiler
  - THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar
related:
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
  - THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors
  - THM-3345-berggren-gaussian-source-groupoid-prime-toggle-circuit-and-ancestry-cost
  - THM-3347-u-spine-signed-prime-clock-gram-and-projective-fold-boundary
  - THM-3356-primitive-affine-determinant-shells-parabolic-orbits-and-prime-clock-resultants
script: 04-computation/admissible_composite_parabolic_compiler_thm3358.py
output: 05-knowledge/results/admissible_composite_parabolic_compiler_thm3358.out
script_sha256: 9d56f18fa79502363228f824aa50e5406d6d72825807d73b3a2fd756e7f57e7b
output_sha256: bd7f7c765267a2fed591b682350b257e0b3e0b5eb2b4ea01b2c469df36fff99c
hash_basis: LF-normalized bytes
---

# THM-3358 -- admissible composite parabolic compiler and Hensel-normal offsets

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The composite, prime-power, word-order, and exact-grade claims below have
passed independent type, proof, and hostile audits.  The theorem generalizes
the proved prime compiler while retaining its Gaussian allocation, normal
coordinate, and ancestry-loss sidecars.

## 1. Inheritance and connection contract

[THM-3353](THM-3353-split-prime-parabolic-branch-transplant-and-unary-transducer-compiler.md)
turns either root of `C_t=0 mod p`, for a split prime `p`, into one explicit
`p^2`-spaced Berggren/U-spine branch.  It stops at two selected non-Hensel
lifts and leaves the other `p-1` valuation-one lifts open.
[THM-3346](THM-3346-u-spine-prime-toggle-root-atlas-and-conjugation-monodromy.md)
proves that the roots of

```text
C_t=2t^2+2t+1                                           (1)
```

at an admissible grade form the CRT prime-toggle cube.
[THM-3336](THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation.md)
supplies the primitive Gaussian content invoice.

The corrected near miss is the raw offset in THM-3353: its two gates are not
antipodes modulo `p^2`.  The least-used sidecar is the normal coordinate of a
root modulo `N^2`.  The proved connection is

| field | contract |
|---|---|
| source | a primitive Gaussian representation of an admissible odd grade `N`, a gate `X`, and a normal unit `alpha mod N` |
| target | one fixed-prefix / unary / fixed-suffix Berggren branch |
| map | center the gate offset to a root `ell mod N^2`, then use `t=N^2s+ell+alpha N` |
| preserved | exact `N`-adic grade, Gaussian allocation toggle, hypotenuse, and primitivity of the compiled target |
| destroyed | raw Gaussian unit/conjugation gauge, absolute height under reduction, and global ancestry uniformity |
| sidecar | centered Hensel root, normal unit `alpha`, representation/gate label, and fixed word data |
| hostile | for `N=5`, raw offsets `{1,8}` are not antipodes, while centered roots `{21,3}` are |

The new claim is not merely that composite roots exist.  It identifies every
exact-grade normal lift and gives each one an explicit Berggren transducer.

## 2. Every oriented Gaussian allocation has an odd cusp

Let

```text
N=product_(j=1)^k p_j^e_j >1,       p_j=1 mod 4,        (2)
```

where the `p_j` are distinct odd primes and `e_j>=1`.  Call such `N`
**admissible**.  For odd `N`, this is equivalent to the existence of a
primitive representation as a sum of two squares: an inert prime divisor
would divide both coordinates, while choosing one Gaussian prime above each
split factor constructs the converse.  Choose such a representation

```text
N=a^2+b^2,        a>b>0,        gcd(a,b)=1,             (3)
```

whose coordinates necessarily have opposite parity.  Put

```text
d=a^2-b^2,       e=2ab,       kappa=sign(d-e),
H=[d e; -kappa e kappa d],
h=H(1,1)=(d+e,|d-e|).                                  (4)
```

Then

```text
H^T H=N^2 I,             det(H)=kappa N^2.              (5)
```

The two entries of `h` are positive, odd, coprime, and unequal.  Hence the
odd-cusp descent of THM-3353 gives a unique word `R_g` such that

```text
B(R_g)(3,1)=h.                                          (6)
```

Here `g=a+ib` records the chosen Gaussian representation, and the parameter
matrices and root-to-child word convention are

```text
B_U=[2 -1;1 0],       B_A=[2 1;1 0],       B_D=[1 2;0 1],
B(w_1...w_l)=B_(w_l)...B_(w_1).                         (7)
```

There are `2^(k-1)` choices (3): for every split prime power one allocates
`pi_j^e_j` or its conjugate, modulo simultaneous conjugation and units.

Use the cusp collision `B_A(1,1)=B_D(1,1)=(3,1)` and define

```text
S_X=X R_g,       M_X=B(S_X),       X in {A,D}.           (8)
```

Then

```text
M_X(1,1)=h,       det(M_X)=eta_X in {+1,-1},
eta_A=-eta_D.                                             (9)
```

Let `P=[0 1;1 0]`.  Choose the unique

```text
H_X in {H,HP}       with det(H_X)=eta_X N^2,             (10)
```

and let `c_X` be the second column of `M_X`.

## 3. Centering reveals the square-modulus Hensel root

Define the centered intercept `ell_X` by

```text
(ell_X+1,ell_X)=-H_X^T c_X.                             (11)
```

The right side really has consecutive coordinates.  Write the columns of
`H_X` as `k_1,k_2` and let `J(x,y)=(-y,x)`.  Equations (5), (9), and (10)
give

```text
k_1-k_2=-eta_X J h,       det(h,c_X)=eta_X.             (12)
```

Taking the coordinate difference of `-H_X^T c_X` in (12) gives one.
Taking norms in (11) gives the hidden exact law

```text
C_(ell_X)=N^2 ||c_X||^2.                                (13)
```

Thus `ell_X mod N^2` belongs to

```text
R_(N^2)={ell mod N^2:C_ell=0 mod N^2}.                  (14)
```

Equivalently `x_X=2ell_X+1` satisfies `x_X^2=-1 mod N^2`.
The Gaussian gcd of `N` with `x_X+i` recovers the corresponding norm-`N`
factor allocation up to its fixed unit/conjugation convention.  Thus the
classical square root of `-1` attached to a primitive sum of two squares is
literally the centered compiler offset.

The map

```text
(g,X) |--> ell_X mod N^2                                (15)
```

is a bijection.  Indeed the `2^(k-1)` primitive representations, together
with the two determinant-matched gates, give all `2^k` oriented Gaussian
factor allocations.  In complex notation, the ungauged matrix `H` acts by

```text
z |-> conjugate(g)^2 z       if kappa=+1,
z |-> g^2 conjugate(z)       if kappa=-1,               (15a)
```

while `P` acts by `z |-> i conjugate(z)`.  Thus the two gates supply the two
oriented lifts of the allocation represented by `g`.  Reduction of (11)
modulo every `p_j^e_j` identifies these with the two local square roots of
`-1`; the determinant-one relation `det(h,c_X)=eta_X` prevents a deeper or
vanishing local choice.  More explicitly, (5) and (11) give

```text
H_X(ell_X+1,ell_X)=-N^2 c_X.                            (15b)
```

The primitive Gaussian ideal `(N,(ell_X+1)+i ell_X)` has norm `N` and is the
oriented allocation encoded by `H_X`, up to the fixed unit convention.
Different allocations give different ideals and hence different roots.
Hensel uniqueness lifts each allocation modulo `p_j^(2e_j)`.  CRT gives
exactly `2^k` roots in (14), so this injection and the count prove exhaustion.

Equivalently, for a coarse root `tau in R_N`, its unique lift in (14) is

```text
ell=tau+Nk,
k=-(C_tau/N) [2(2tau+1)]^(-1) mod N.                   (16)
```

The derivative is a unit because `(2tau+1)^2=-1 mod N`.

## 4. The full normal atlas and the exact-grade shell

Every residue `t mod N^2` satisfying `N|C_t` has a unique expression

```text
t=ell+alpha N,       ell in R_(N^2),       alpha in Z/NZ. (17)
```

This is Hensel projection followed by the normal displacement from its unique
lift.  Taylor expansion at `ell` gives

```text
C_(ell+alpha N)/N = 2alpha(2ell+1) mod N.               (18)
```

Consequently

```text
gcd(N,C_t/N)=gcd(N,alpha),                              (19)
```

and, locally for `p^e || N`,

```text
min(v_p(C_t),2e)=e+min(v_p(alpha),e).                   (20)
```

In particular,

```text
v_(p_j)(C_t)=e_j for every j
  iff alpha in (Z/NZ)^x.                                (21)
```

Call (21) the **exact-grade normal shell**.  It contains

```text
|R_(N^2)| phi(N)=2^k phi(N)                             (22)
```

residue rays.  Equations (17)--(21) also prove the equality boundary:
nonunit `alpha` adds precisely its truncated prime-power valuation, and
`alpha=0` returns to the deeper square-modulus root.

There is a more intrinsic first-normal coordinate.  On the exact-grade shell
define

```text
nu_N(t)=C_t/N mod N in (Z/NZ)^x.                        (22a)
```

For `t=ell+alpha N`, equation (18) says

```text
nu_N(t)=2alpha(2ell+1) mod N.                           (22b)
```

The centered root can be recovered intrinsically by the Newton retraction

```text
Ret_N(t)=t-C_t [C'_t]^(-1) mod N^2,                    (22b1)
```

where `C'_t=2(2t+1)` is a unit.  On the exact-grade shell the map

```text
t |-> (Ret_N(t),nu_N(t))                               (22b2)
```

is a bijection onto `R_(N^2) x (Z/NZ)^x`, with inverse
`(ell,nu) |-> ell+N nu [C'_ell]^(-1)`.  Thus, after the standard CRT
prime-toggle gauges are chosen, the full exact-grade shell is a torsor for

```text
F_2^k x (Z/NZ)^x.                                      (22c)
```

The first factor toggles Hensel-root branches and the second rescales the
nonzero normal coordinate.  This is the unit-frame bundle of the first
infinitesimal normal line to the simple-root locus in the parameter line, not
another copy of the cube or a Euclidean normal bundle.  Arithmetically,
`nu_N(t)` is the complementary
hypotenuse factor `C_t/N` reduced modulo `N`; it records the transverse
coprimality that the root allocation alone forgets.

The gauge matters.  For a subset `S` of prime-power coordinates, let
`chi_S in (Z/NZ)^x` be the CRT sign that is `-1` on `S` and `+1` off `S`, and
let `T_S` be the corresponding root toggle.  In the intrinsic **jet gauge**,

```text
(S,beta):(ell,nu) |-> (T_S ell,beta nu).                (22d)
```

Since `2T_S(ell)+1=chi_S(2ell+1) mod N`, the same action in the compiler's
`alpha` coordinate is

```text
(ell,alpha) |-> (T_S ell,chi_S beta alpha).             (22e)
```

Alternatively one may declare an alpha-gauge product action that keeps
`alpha` fixed under root toggles; the two product gauges differ by this
character cocycle.  One must not keep both `alpha` and `nu` fixed under a
partial toggle.

The involution induced by Gaussian conjugation is

```text
(ell,alpha) |--> (-1-ell,-alpha).                       (23)
```

It is the true antipode of the normal atlas.  Since `C_(-1-t)=C_t`, it fixes
`nu_N(t)`.  Hence the exact-grade conjugation quotient is

```text
(F_2^k/<1>) x (Z/NZ)^x                                 (23a)
```

as a torsor after the same gauges, with `2^(k-1)phi(N)` points.  Forgetting
the normal unit collapses `phi(N)` genuinely different compiled rays over
each fixed parent allocation.  In alpha gauge the same quotient is
`(F_2^k x (Z/NZ)^x)/<({all},-1)>`; the jet gauge turns the involution into
`({all},1)` and gives the displayed direct product.

### 4A. Harmonic mass is normalized normal-state count

THM-3359 gives an exact analytic readout of this finite atlas.  Let
`B subset R_(N^2) x (Z/NZ)^x`, and let `a_b in {1,...,N^2}` be the positive
parameter residue represented by `b` under (22b2).  The positive index set

```text
H_B={t>=1:t mod N^2 is represented by some b in B}       (23b)
```

has Dirichlet series and residue

```text
D_B(s)=N^(-2s) sum_(b in B) zeta(s,a_b/N^2),
Res_(s=1)D_B(s)=|B|/N^2.                                (23c)
```

Thus one oriented compiler residue has harmonic coefficient `1/N^2`, one
antipodal parent pair has coefficient `2/N^2`, and the full exact-grade shell
has

```text
delta_exact=2^k phi(N)/N^2.                             (23d)
```

The conjugation quotient has half as many parent labels, but its pullback to
the integer-time carrier is still the union of both residues, so its mass is
not halved a second time.

Equation (23c) is a normalized state count, not a compiler decoder.  Sets of
labels with the same cardinality have the same harmonic residue even when
their prime toggles, normal units, Gaussian allocations, or Berggren words
differ.  Weighting the same ray by its geometric hypotenuse instead gives
`1/C_t=O(t^(-2))` and a convergent series.  The index carrier and the retained
label sidecars are therefore load-bearing.

## 5. Every exact-grade normal ray has a Berggren compiler

Choose the standard representative `1<=alpha<N` of a unit modulo `N`.  The
spinor

```text
(N+alpha,alpha)                                         (24)
```

is primitive, opposite-parity, and positive.  Therefore it has a unique
Berggren word `P_(N,alpha)` satisfying

```text
B(P_(N,alpha))(2,1)=(N+alpha,alpha).                    (25)
```

For `alpha=1`, this prefix is especially simple:

```text
P_(N,1)=D^((N-1)/2).                                    (26)
```

For a representation/gate pair `(g,X)`, put

```text
r_(X,alpha)=ell_X+alpha N,
t_(X,alpha)(s)=N^2s+r_(X,alpha).                        (27)
```

Then for every integer `s>=1`,

```text
B(P_(N,alpha) U^(s-1) S_X)(2,1)
 = (1/N) H_X(t_(X,alpha)(s)+1,t_(X,alpha)(s)).          (28)
```

The proof is a two-line matrix identity.  First,

```text
B_U^(s-1)(N+alpha,alpha)=(Ns+alpha,N(s-1)+alpha).
```

Since the second column of `M_X` is `c_X` and `M_X(1,1)=h`, the left side
of (28) equals

```text
(Ns+alpha)h-Nc_X.                                      (29)
```

On the other hand, (5), (11), and `H_X(1,1)=h` turn the right side into the
same expression.  This proves (28) without primality or squarefreeness.

When U-spine source notation is desired, impose only the explicit threshold

```text
s>=s_(X,alpha,0)=max(1,ceil((1-r_(X,alpha))/N^2)),       (30)
```

so `t>=1`.  The left side of (28) is a Berggren word from the root, hence its
parameter pair is primitive, opposite-parity, and in the positive chamber.

Equations (15), (21), and (28) prove the advertised exhaustion: all
`2^k phi(N)` exact-grade residue rays receive an explicit
fixed-prefix / unary / fixed-suffix transducer.  For every fixed tuple
`(N,alpha,g,X)`, a finite controller recognizes the source exponent modulo
`N^2`, compresses the moving `U` block by `N^2:1`, and emits the fixed words
`P_(N,alpha)` and `S_X`.  No one finite automaton uniform in `N` is asserted.

## 6. Gaussian content, gauges, and CRT locality

The rational orthogonal action `H_X/N` in (28) changes the selected block of
Gaussian prime-power allocations, modulo the same unit/global-conjugation
gauge as THM-3353.  Because (21) makes the complementary norm coprime to `N`,
THM-3336's two primitive Gaussian compositions have exact unordered contents

```text
{N,C_t/N}.                                              (31)
```

This is the folded parent weight, not an orientation of the toggle cube.

For `alpha=1`, write `r_X=ell_X+N`.  The two gates attached to one
representation satisfy

```text
r_A+r_D+1=2N-N^2,                                      (32)
ell_A+ell_D+1=-N^2.                                    (33)
```

Thus the raw compiler section is translated away from conjugation by the
universal amount `2N`, while centering restores the literal antipode modulo
`N^2`.  More generally (23) reverses the normal coordinate as well as the
root.

Centering is also CRT-local.  Let `N=QM` with `gcd(Q,M)=1`, where `Q` is an
admissible unitary divisor, and restrict the Gaussian allocation to `Q`,
choosing its unit/gate gauge compatibly with the one at `N`.
Then

```text
ell_N=ell_Q mod Q^2.                                   (34)
```

For the selected `alpha=1` raw offsets this becomes the cofactor shear

```text
r_N=r_Q+Q(M-1) mod Q^2.                                (35)
```

The centered root, rather than the raw branch intercept, is therefore the
coordinate that factors prime-power by prime-power.

This locality compiles more than the all-`N` diagonal.  Let `t` have exact
`N`-grade and let `Q` range over the admissible unitary divisors `Q>1` of
`N`, together with `Q=1` as the separately declared identity.  Then `t` has
exact `Q`-grade, with

```text
nu_Q(t)=(N/Q) nu_N(t) mod Q.                            (35a)
```

The cofactor is a unit modulo `Q`.  For all sufficiently large `t`, applying
(28) at each `Q` therefore gives a lawful Berggren address for the parent
obtained by conjugating exactly the `Q`-primary Gaussian block.  Varying `Q`
compiles all `2^k` vertices in the oriented `F_2^k` star on the selected
prime-power blocks.  It does not supply the edges between the non-U-spine
outputs.  On the parent quotient these vertices remain distinct when
`C_t/N` has an outside prime coordinate; if `N` is the full grade, complementary
subsets are global-conjugate and only `2^(k-1)` classes remain.  Partial
powers `p^f`, `0<f<e`, are excluded: they leave the same prime in the
complement and can destroy primitive output.  The words and unary scales
depend on `Q`; they are endpoint compilers, not a commuting action or cubical
embedding inside the Berggren tree.

### 6A. The normal bundle is affine-quadratic; the word compiler is not

The Hensel-normal mechanism is more general than the U-spine.  Take a
primitive direction `u`, an integral base point `x_0`, and

```text
A=||u||^2,       h=x_0.u,       C=||x_0||^2,
c=det(x_0,u)!=0,
Q(n)=||x_0+n u||^2=An^2+2hn+C.                          (35b)
```

Then `AC-h^2=c^2`.  Let `M` be admissible with `r=omega(M)` and

```text
gcd(M,2Ac)=1.                                           (35c)
```

Every one of the `2^r` roots of `Q(n)=0 mod M` has a unique lift `ell`
modulo `M^2`.  Every residue above the coarse root is uniquely
`n=ell+alpha M`, and

```text
Q(n)/M = 2alpha(Aell+h) mod M,
gcd(M,Q(n)/M)=gcd(M,alpha).                             (35d)
```

Here `Aell+h` is a unit because its square is `-c^2 mod M`.  Thus the exact
`M`-grade residues admit the intrinsic coordinates

```text
Ret_M(n)=n-Q(n)[Q'(n)]^(-1) mod M^2,
nu_M(n)=Q(n)/M mod M,                                  (35d1)
```

with inverse `(ell,nu) |-> ell+M nu [Q'(ell)]^(-1)`.
They therefore form a torsor

```text
F_2^r x (Z/MZ)^x                                       (35e)
```

of size `2^r phi(M)`.

For a CRT sign `chi_S` modulo `M^2`, the exact root toggle is

```text
T_S(ell)=A^(-1)(chi_S(Aell+h)-h) mod M^2.               (35e1)
```

In the jet coordinate `nu=Q(n)/M mod M`, `(S,beta)` acts as
`(ell,nu) |-> (T_S ell,beta nu)`; in an alpha coordinate the corresponding
factor is `chi_S beta`.  This is the same two-gauge boundary as (22d)--(22e).
The affine reflection

```text
n |-> -2h A^(-1)-n mod M^2                             (35f)
```

negates the local root branch, reverses `alpha`, preserves `Q mod M^2`, and
therefore fixes the intrinsic normal coordinate `Q(n)/M mod M`.  Its quotient
has `2^(r-1)phi(M)` points.  This involution is canonical on the residue atlas;
it is a literal integer symmetry of the enumeration only when `A|2h`.

For an admissible unitary divisor `D|M`, compatible root gauges obey the
strict restriction laws

```text
ell_M=ell_D mod D^2,
alpha_D=(M/D)alpha_M mod D,
nu_D=(M/D)nu_M mod D.                                  (35g)
```

Thus the local products glue with a cofactor shear.  Absolute height, the
integral representative, and the determinant-shell parabolic residue are
not retained by this finite normal atlas.  Since `Q'(ell)=2(x_ell.u)`, the
jet records the directional-derivative trivialization along the affine
parameter line; it is not a geometric normal vector to the determinant shell.

Equations (35b)--(35g) are the maximal arithmetic mechanism: a tame cyclic
lattice line has the same root-cube normal bundle.  The explicit Berggren
identity (28) additionally uses the consecutive U-spine, its odd cusp, and
the `A/D` collision.  A general affine shell has no inherited Berggren word,
and THM-3356's norm-85 carrier hostile rules out a physical LRC transport.

### 6B. Determinant charge is the exact unary-compiler obstruction

There is a sharper reason the word compiler does not automatically extend.
For arbitrary fixed Berggren words `P,S`, put

```text
w=B(P)(2,1),       delta=w_1-w_2,
h_S=B(S)(1,1),
Y_s=B(S)B_U^s w.                                         (35h)
```

Since `B_U^s w=w+s delta(1,1)`, one has

```text
Y_s=Y_0+s delta h_S,
|det(Y_s,h_S)|=|delta|.                                 (35i)
```

Thus every fixed-prefix / unary-`U` / fixed-suffix compiler advances by
exactly one parabolic transvection step: its step coefficient equals its
determinant-shell charge.  The direction `h_S` is primitive odd/odd.

Now apply a **primitive allocation matrix** `G` to an affine shell (35b):
`G` is the determinant-matched integral Gaussian matrix attached to a full
primitive norm-`M` allocation, so

```text
G^T G=M^2 I,       det(G)=+/-M^2,       gcd(entries of G)=1.       (35j0)
```

The last condition is load-bearing: `G=M I` satisfies the norm identity but
is not an allocation matrix and never preserves primitive output.  Let `r` be
an integer representative of a selected coarse root modulo `M` for which
`Gx_r=0 mod M`, chosen in the exact-grade normal class
`Q(r)/M in (Z/MZ)^*`.  Then

```text
Z_s=(1/M)G x_(r+M^2s)                                  (35j)
```

is integral.  Tameness makes `Gu` primitive: a prime outside `M` is excluded by
invertibility of `G`, while a prime dividing `M` would put `u` on an isotropic
kernel and force that prime to divide `A`.  Directly,

```text
Z_(s+1)-Z_s=M Gu,
|det(Z_s,Gu)|=M|c|.                                    (35k)
```

Suppose, up to a fixed index shift and a fixed sign/unit gauge, that this ray
had a one-step word identity

```text
Z_s=B(S)B_U^sB(P)(2,1).                                 (35k0)
```

Both directions are primitive, so (35i) and (35k) force
`h_S=+/-Gu` and `|delta|=M`; their determinant charges then agree only when
`M=M|c|`.  Hence comparison proves the no-go:

```text
a one-step unary-U compiler for (35j) requires |c|=1.   (35l)
```

The U-spine has `c=1`, so (28) lies exactly on the equality boundary.  For
`|c|>1`, splitting `s` into its `|c|` residue classes replaces the step by
`|c|M Gu` and removes the coefficient obstruction; these are precisely
THM-3356's parabolic-orbit sidecars, but the matching alone does not prove a
word exists.

The charge obstruction is independent of parity.  With
`u=(1,1),x_0=(2,-1)`, one has `c=3`, `Q(n)=2n^2+2n+5`, and the tame grade
`M=5`; source parity and the odd cusp are compatible, yet (35k) has step
coefficient `5` and charge `15`.  Conversely, if `u` is opposite-parity then
an odd-grade Gaussian matrix preserves that parity type, so `Gu` cannot equal
the odd/odd direction `h_S` even when `|c|=1`.  For example
`u=(2,1),x_0=(1,0),M=13` is tame and hits this separate cusp-parity no-go.

## 7. Boundaries and decisive hostiles

The hypotheses in (2)--(3) are load-bearing.

- `N=10=3^2+1^2` is even: the cusp is `(14,2)`, not primitive odd, and
  `(N-1)/2` is not an integer.
- `N=45=6^2+3^2` has only a nonprimitive displayed representation and an
  inert prime; `R_N` is empty.
- At the admissible prime power `N=125`, the representation `10^2+5^2` is
  invalid, while the primitive representation `11^2+2^2` is lawful.
- Prime powers themselves are not a defect: `N=25,125,169` satisfy the full
  construction.
- At `N=5`, the selected raw offsets are `{1,8} mod 25`, whose antipodes are
  not one another.  Centering by `5` gives `{21,3}`, and
  `-1-21=3 mod 25`, exactly as (33) requires.
- A nonunit `alpha` does not belong to the exact-grade shell; (20) measures
  its extra valuation rather than silently discarding it.
- In the affine extension, a prime dividing `c` is discriminant-ramified and
  a prime dividing `A` removes the inverse in (35f); neither belongs to the
  Boolean tame atlas.  For example `u=(1,0),x_0=(0,5),M=5` gives
  `Q(n)=n^2+25`, with one coarse root and five square-modulus lifts, while
  `u=(1,2),x_0=(1,0),M=5` makes `A=5` and again destroys the two-root atlas.
  The prime `2` is the separate parity boundary.
- A partial prime-power toggle is genuinely unsafe: at `N=25,t=428`, using
  only `Q=5` produces parameter pair `(600,85)` of content `5`.  This is why
  the Boolean star uses full unitary blocks `p^e`, never `p^f` with `f<e`.

The parent operation in (31) is not a rooted automorphism of the Berggren
tree.  The source word is `U^(t-1)`, while the target word in (28) has fixed
prefix and suffix but a different unary scale.  Their path distance grows as

```text
(N^2+1)s+O_(N,alpha,g,X)(1).                            (36)
```

Indeed the source word is all `U`, whereas the target begins with the fixed
word `P_(N,alpha)`.  This prefix cannot be a pure `U` word:
`B(U^j)(2,1)` has coordinate difference one, while (24) has difference
`N>1`.  The longest common prefix is therefore bounded independently of `s`,
and the two depths have leading terms `N^2s` and `s`.

No ancestry-current, LRC owner, phase, physical mass, tournament orientation,
or globally uniform transducer follows.

## 8. Pell grades give unbounded multi-toggle compilers

THM-3346 constructs square-triangular Pell grades

```text
H=C_n                                                     (37)
```

whose prime factors all split in `Z[i]` and whose `omega(H)` is unbounded.
Applying (15), (22), (28), and the unitary-divisor construction after (35a)
at `N=H` gives

```text
2^omega(H) phi(H)                                       (38)
```

explicit exact-grade residue branches, organized into all oriented Gaussian
allocations and normal units.  On every sufficiently high such source, the
unitary-divisor compilers give the entire oriented `H`-support vertex star;
the `Q=H` branch is its full diagonal.  For large `s`, that diagonal has
folded content distance `log H`, while (36) has ancestry cost on the order of
`H^2s`.  This is a precise Pell/triangular-to-multi-toggle compiler and an
equally precise no-cost-transfer boundary.

It does not identify the prime-toggle cubical topology with the Berggren tree:
the former retains commuting allocation squares and conjugation, while the
latter records source-dependent ancestry and has zero graph `H^1`.

## 9. Exact companion and audit scope

The exact companion compares two routes:

1. the odd-cusp/matrix compiler (4)--(13), followed by an independent
   three-dimensional Berggren evaluation; and
2. direct CRT roots modulo `N`, Newton lifting by (16), and the full normal
   atlas (17)--(21).

The bounded universe covers every admissible odd `N<=1200`, explicit prime
powers, and the rank-four grade `N=32045`.  It checks every representation,
both gates, every normal unit, (28), Gaussian contents (31), conjugation and
CRT shears, the general affine normal bundle, the bounded universal unary-word
identity, the unitary-divisor star, all hostiles in Section 7, and
normal/optimized transcript identity.  Its `1,106,432` compiled word rows and
`87,846` independent unary-word rows are exact integer checks.  The script and
stored output hashes are frozen in the frontmatter.

This establishes the arithmetic/compiler theorem, not a commuting Berggren
cube, an arbitrary-affine word compiler, an LRC owner/current, or LRC(14).
