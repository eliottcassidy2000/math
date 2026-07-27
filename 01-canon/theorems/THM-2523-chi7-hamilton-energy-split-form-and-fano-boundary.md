---
id: THM-2523
title: "Chi7 Hamilton-energy real split form and the Fano orientation boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The multiplicative chi_7
  contrast of the six Hamilton-cycle energies on the real F_13 augmentation
  module is the quadratic form of a symmetric circulant operator A_tau.  Its
  six distinct cyclotomic eigenvalues have minimal polynomial
  t^6-39t^4+299t^2-325, each occurs twice, and the form is nondegenerate of
  real signature (6,6).  Dilation by five is an exact anti-isometry.  Every
  centred point mass is nevertheless isotropic, and there are explicit
  three-dimensional rational fixed and anti-fixed totally isotropic spaces.
  Over Q the Witt index is exactly five, leaving one anisotropic binary
  plane.  The canonical derived partner q=A_tau p gives the coercive scalar
  beta_tau(p,q)=||A_tau p||^2>(13/10)||p||^2 on augmentation.  For the
  THM-2521 collision profile this exceeds (169/10)D_13 whenever D_13>0,
  although q is signed rather than one Boolean ancestry event.  The Fano
  slope product A_1 A_2 A_4 is 5 times the symmetric chi_13 Paley operator;
  because chi_13(-1)=+1 it is an undirected signed graph, not a tournament or
  an orientation of the physical self-cospan.  No Boolean owner charge,
  semantic sheet identification, row exclusion, or LRC(14) proof follows.
source: codex-2026-07-27-chi7-hamilton-split-form
depends_on:
  - THM-2514-cyclic-k14-factor-chart-and-six-phase-ordinary-degree-reconstruction
  - THM-2519-last-digit-collision-drift-and-k13-dirichlet-boundary
  - THM-2521-k13-drift-k14-potential-module-bridge
related:
  - THM-1260-placed-fork-chi7-surjectivity-guardrail
  - THM-2148
  - THM-2168
  - THM-2524-translated-chi7-hamilton-polarization-inversion
script: 04-computation/lrc14_chi7_hamilton_split_thm2523.py
output: 05-knowledge/results/lrc14_chi7_hamilton_split_thm2523.out
script_sha256: 9ef9e16d423a0bdf3bc5eba834c9416025cc5329fc935962dc27ac6accf504b6
output_sha256: e5a61f9bb91979f0eae9ae708ba3604bd38c327915199eec234857d026bf787d
independent_script: 04-computation/lrc14_chi7_hamilton_witt_referee_thm2523.py
independent_output: 05-knowledge/results/lrc14_chi7_hamilton_witt_referee_thm2523.out
independent_script_sha256: 1c79b135582d06d6e09d31ad50c14c3ac4394125aac9f0a1147da3030b42b92f
independent_output_sha256: 32ce81b72f11de09d508d19688cbd2eaff1be3376931b793215c55eb80bb5dc4
hash_basis: working-tree bytes (LF)
---

# THM-2523 -- the `chi_7` Hamilton contrast is real-split, not oriented

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2521 identifies the six finite columns of the cyclic `K_14` chart with
the six antipodal collision-gap energies of the thirteen predecessor roots.
Their sum is positive whenever the last-digit collision drift is positive,
and all six additive septimal characters survive over the rationals.  The
remaining multiplicative probe is

```text
sum_(s=1)^6 chi_7(s) E_s.                                  (1)
```

It is tempting to expect (1) to orient the six columns by the Fano split

```text
{1,2,4} versus {3,5,6}.                                    (2)
```

The exact object says something subtler.  It is a nondegenerate symmetric
quadratic form, but it is maximally indefinite over the reals.  In
particular, nondegeneracy of the polarized form does not make its diagonal
nonzero on a nonzero predecessor profile.  Centred point masses -- the
sharp hostile already present in THM-2521 -- are isotropic.

The same calculation explains the orientation failure conceptually.  The
product of the three residue-slope operators in (2) is the **symmetric**
quadratic-character operator at `13`.  The Fano selector has landed in a
signed undirected Cayley graph, not in a Paley tournament.  An oriented
half-set, owner loop, or chronological sheet is extra data.

## 1. The signed Hamilton operator

Let

```text
U_R={p=(p_x)_(x in F_13) in R^13 : sum_x p_x=0}.             (3)
```

Everything below also holds coefficientwise for vectors in a real Hilbert
space.  Put

```text
chi(s)=+1,  s in {1,2,4},
chi(s)=-1,  s in {3,5,6}.                                   (4)
```

Fix a chart slope `tau in F_13^*`.  For `s=1,...,6`, use the THM-2521
Hamilton energy

```text
E_(tau,s)(p)
 =sum_(h in F_13)||p_(h-tau s)-p_(h+tau s)||^2.              (5)
```

Define

```text
Q_tau(p)=sum_(s=1)^6 chi(s) E_(tau,s)(p).                   (6)
```

Let `T_d` be the cyclic translation `(T_d p)_x=p_(x+d)`.  The
corresponding self-adjoint circulant operator is

```text
A_tau
 =-sum_(s=1)^6 chi(s)(T_(2tau s)+T_(-2tau s)),              (7)

(A_tau p)_x
 =-sum_(s=1)^6 chi(s)
      [p_(x+2tau s)+p_(x-2tau s)].                          (8)
```

Indeed,

```text
E_(tau,s)(p)
 =2||p||^2-<p,T_(2tau s)p>-<p,T_(-2tau s)p>.                (9)
```

Since `sum_s chi(s)=0`, the norm terms cancel, and hence

```text
Q_tau(p)=<p,A_tau p>.                                      (10)
```

The matrix of `A_tau` is symmetric, has zero diagonal, and has zero row
sum.  Thus constants are in its kernel and it descends to (3).  The zero
diagonal is already a warning: every coordinate ray can be isotropic even
when the restriction to the augmentation module is nondegenerate.

## 2. Cyclotomic spectrum and exact characteristic polynomial

Let

```text
zeta=exp(2 pi i/13),
p_hat(k)=sum_x p_x zeta^(-kx).                              (11)
```

The Fourier eigenvalue of (7) at `k` is

```text
lambda_(tau,k)
 =-sum_(s=1)^6 chi(s)
    [zeta^(2k tau s)+zeta^(-2k tau s)].                     (12)
```

It is real and satisfies

```text
lambda_(tau,-k)=lambda_(tau,k),
lambda_(tau,0)=0.                                          (13)
```

Parseval gives the exact diagonalization

```text
Q_tau(p)
 =1/13 sum_(k in F_13^*) lambda_(tau,k)|p_hat(k)|^2.        (14)
```

Multiplication of `k` by `tau` only permutes the eigenvalues, so their
multiset is independent of the chosen affine slope gauge.

The complete algebra is unusually small.  Work in

```text
Z[X,X^(-1)]/(X^13-1)
```

and put

```text
L=-sum_(s=1)^6 chi(s)(X^(2s)+X^(-2s)).                      (15)
```

In the representatives `0,...,12`,

```text
L
 = X+X^3+X^6+X^7+X^10+X^12
   -X^2-X^4-X^5-X^8-X^9-X^11.                              (16)
```

Exact cyclic convolution gives

```text
L^6-39L^4+299L^2-325
 =-25(1+X+...+X^12).                                       (17)
```

At every nontrivial thirteenth root, the right side vanishes.  Therefore
each nontrivial eigenvalue is a root of

```text
f(T)=T^6-39T^4+299T^2-325.                                 (18)
```

The polynomial (18) is Eisenstein at `13`: every nonleading coefficient is
divisible by `13`, while `325=13*25` is not divisible by `13^2`.  Hence it
is irreducible over `Q`.  The six Galois conjugates of one eigenvalue in the
real cyclotomic field `Q(zeta+zeta^(-1))` are exactly the six values indexed
by `F_13^*/{+/-1}`.  They are distinct, and every one occurs twice on the
twelve-dimensional real augmentation module.  Consequently

```text
char_(A_tau|U_R)(T)=f(T)^2.                                (19)
```

In particular `A_tau|U_R` is invertible.  Equation (17) also gives the
explicit inverse

```text
A_tau^(-1)
 =(A_tau^5-39A_tau^3+299A_tau)/325                         (20)
```

on `U_R`.  THM-2524 uses (20) to invert the translated polarization bank;
the present theorem only identifies the diagonal form and its boundary.

For the real inertia, write

```text
g(Y)=Y^3-39Y^2+299Y-325,
f(T)=g(T^2).                                                (21)
```

The exact sign table is

```text
g(1)=-64,    g(2)=125,
g(8)=83,     g(9)=-64,
g(29)=-64,   g(30)=545.                                    (22)
```

Thus `g` has one positive root in each of `(1,2)`, `(8,9)`, and
`(29,30)`.  These are all of its roots.  Therefore `f` has three positive
and three negative real roots.  Each has multiplicity two in (19), proving

```text
signature(Q_tau|U_R)=(6,6).                                (23)
```

For a finite-dimensional real Hilbert coefficient space `H`, the signature
on `U_R tensor H` is `(6 dim H,6 dim H)`.  For the physical
`H=L^2(G(z)dz)`, equations (7)--(20) remain bounded operator identities on
the finite predecessor index.  The word **split** in this theorem refers to
the real signature (23); it is not a claim that the twelve-dimensional
rational form has Witt index six.

## 3. Dilation by five is an anti-isometry

Define the orthogonal coordinate dilation

```text
(D_5 p)_x=p_(5x).                                          (24)
```

On the six antipodal Hamilton directions, multiplication by five induces

```text
pi_5=(1 5)(2 3)(4 6),

(pi_5(1),...,pi_5(6))=(5,3,2,6,1,4).                       (25)
```

Here each product is reduced to its representative in `{1,...,6}` modulo
sign in `F_13`.  Directly from (4),

```text
chi(pi_5(s))=-chi(s)                                       (26)
```

for all six `s`.  Reindexing (5) by `x -> 5x` now gives

```text
E_(tau,s)(D_5p)=E_(tau,pi_5(s))(p),                         (27)

Q_tau(D_5p)=-Q_tau(p).                                     (28)
```

Equivalently,

```text
D_5^* A_tau D_5=-A_tau.                                    (29)
```

This gives a second conceptual proof that the positive and negative inertia
are equal once nondegeneracy is known.  Notice also

```text
D_5^2=D_(-1).                                               (30)
```

Thus the anti-isometry squares to converse/reflection.  It exchanges the two
sign classes; it does not choose one as the forward class.

## 4. Explicit rational isotropic hostiles

For `a in F_13`, put

```text
v_a=e_a-(1/13)1.                                           (31)
```

Every cycle direction sees the same centred point mass:

```text
E_(tau,s)(v_a)=2                    for every s=1,...,6.     (32)
```

It follows that

```text
Q_tau(v_a)=2sum_s chi(s)=0.                                (33)
```

These are thirteen rational, nonzero isotropic profiles, and their
differences span `U_R`.  Hence the isotropic cone is not a negligible
degeneracy.  The form is nondegenerate, yet isotropic profiles span its
entire underlying vector space.

The anti-isometry supplies larger rational families.  Multiplication by
five has three nonzero four-cycles

```text
O_1=(1,5,12,8),
O_2=(2,10,11,3),
O_3=(4,7,9,6).                                              (34)
```

For each `i`, define

```text
f_i=1_(O_i)-4e_0,                                          (35)
```

and, writing `O_i=(a,5a,-a,-5a)`, define

```text
g_i=e_a-e_(5a)+e_(-a)-e_(-5a).                             (36)
```

Then

```text
D_5 f_i=f_i,
D_5 g_i=-g_i.                                               (37)
```

Polarizing (29) shows

```text
<f_i,A_tau f_j>=0,
<g_i,A_tau g_j>=0                                           (38)
```

for all `i,j`.  Thus the fixed and anti-fixed spaces are two explicit
three-dimensional rational totally isotropic spaces.  For `tau=1`, the
cross-pairing matrix in the displayed bases is

```text
[[-20, 32, 16],
 [-16, 12, 24],
 [-16,  8, 12]],                                           (39)
```

whose determinant is

```text
-4160=-2^6*5*13.                                           (40)
```

So these two hostile spaces pair nondegenerately with one another; (38) is
genuine split geometry, not a hidden kernel.

### The rational Witt boundary is one anisotropic binary plane

The preceding six-dimensional hyperbolic piece is exactly the
reversal-even half.  Indeed, `D_5^2=D_(-1)`, so both the `D_5`-fixed and
`D_5`-anti-fixed spaces in (35)--(38) are fixed by reversal.  They have
dimension three each and pair nondegenerately by (39)--(40).  Hence the
reversal-even rational summand is three hyperbolic planes.

For `tau=1`, use the reversal-odd basis

```text
o_r=e_r-e_(-r),                              r=1,...,6.
```

One half of its Gram matrix is

```text
G_odd=
 [ 1  0  0  2 -2 -2]
 [ 0  1  2 -2  0  0]
 [ 0  2 -1  0  0  2]
 [ 2 -2  0  1  2 -2]
 [-2  0  0  2 -1  2]
 [-2  0  2 -2  2 -1],                       det G_odd=-325.   (40a)
```

The two independent vectors

```text
v=(-1,-1,-1,0,-1,0),
w=(-1, 0,-1,0, 0,0)                                         (40b)
```

satisfy

```text
v^T G_odd v=w^T G_odd w=v^T G_odd w=0.                      (40c)
```

Thus the odd Witt index is at least two.  It is not three: a
six-dimensional rational hyperbolic form has determinant squareclass `-1`,
whereas (40a) has squareclass `-13`.  Orthogonality of the reversal halves
therefore gives

```text
Witt_index_Q(A_tau on augmentation)=3+2=5.                  (40d)
```

The remaining rational anisotropic part is a binary plane.  Multiplication
by `tau` permutation-conjugates the form, so the conclusion holds for every
nonzero chart slope.  This sharpens the word **split**: the form is real
split, but not rationally hyperbolic.

## 5. What polarization recovers -- and what the diagonal loses

The symmetric polarization of (10) is

```text
beta_tau(p,q)
 =1/2[Q_tau(p+q)-Q_tau(p)-Q_tau(q)]
 =<p,A_tau q>.                                              (41)
```

By (19),

```text
p!=0  implies  beta_tau(p,q)!=0 for some q.                 (42)
```

But (42) is an existential statement about a **second profile**.  It does
not imply `Q_tau(p)!=0`; equations (31)--(33) are the exact obstruction.
Using polarization physically requires either a lawfully coupled second
ancestry profile or a translated copy whose cross-correlation is actually
available.  It cannot be manufactured from the one diagonal scalar.

Algebraically there is nevertheless one canonical derived partner:

```text
q=A_tau p,

beta_tau(p,A_tau p)
 =<p,A_tau^2p>
 =||A_tau p||^2>0                         for p!=0.          (42a)
```

This is not extra input: it applies the known signed circulant once more.
It is also quantitatively coercive.  With `g(Y)` from (21),

```text
g'(Y)=3Y^2-78Y+299>0                 on [0,13/10],
g(13/10)=-13/1000<0.                                      (42b)
```

Since the three roots of `g` are the squared nonzero eigenvalues of
`A_tau`, none lies in `[0,13/10]`.  Consequently

```text
||A_tau p||^2>(13/10)||p||^2            on augmentation.    (42c)
```

The first row of `A_1^2` is

```text
(12,-5,-1,-1,3,-5,3,3,-5,3,-1,-1,-5).                    (42d)
```

Thus (42a) is a rational, translation-covariant two-step signed stencil.
Its positivity does not Booleanize that stencil or supply a chronological
second ancestry event.

The centred deltas make the lost coordinate explicit.  Since `A_tau`
annihilates constants,

```text
beta_tau(v_x,v_y)=(A_tau)_(x,y).                            (43)
```

At `tau=1`, this is the signed Cayley edge table

```text
positive differences: +/-{1,3,6},
negative differences: +/-{2,4,5},
diagonal: 0.                                                (44)
```

The diagonal `beta(v_x,v_x)` forgets the edge sign completely.  Supplying a
second point exposes it.  The sign in (44) is even under `d -> -d`, so the
polarized object is an undirected signed complete graph, not a tournament.

THM-2524 proves that THM-2519's full translated collision table supplies a
canonical bank of such cross-pairings and that (20) makes that bank lossless
on every nontrivial thirteenth-root mode.  The companion referee below
independently cross-checks its constants and Fourier multiplier.  That
repair remains even under translation converse and therefore still does not
orient the physical self-cospan.

## 6. The Fano slope product lands at the symmetric `13`-Paley boundary

The positive selector

```text
D={1,2,4} subset F_7                                      (45)
```

is a cyclic `(7,3,1)` difference set: the six ordered nonzero differences
between its three elements are exactly `1,...,6`.  Its translates are the
seven lines of the cyclic Fano plane.  Thus (6) is a genuine Fano-labelled
comparison of three Hamilton energies against the complementary three.

Now keep one chosen `F_13` coordinate and consider the three commuting
operators

```text
A_1, A_2, A_4.                                             (46)
```

Let `chi_13` be the quadratic character modulo `13`, extended by
`chi_13(0)=0`, and define the signed Paley operator

```text
(C_13p)_x=sum_(d in F_13^*)chi_13(d)p_(x+d).                (47)
```

Exact cyclic multiplication gives the Fano-slope identity

```text
A_1 A_2 A_4=5C_13.                                         (48)
```

For clarity, the first row of the right side divided by five is

```text
(0,+1,-1,+1,+1,-1,-1,-1,-1,+1,+1,-1,+1).                  (49)
```

The standard quadratic-character correlation gives

```text
C_13^2=13I-J,                                               (50)
```

where `J` is the all-ones operator.  Hence on the augmentation module,

```text
C_13^2=13I.                                                 (51)
```

Equations (48)--(51) provide another exact nondegeneracy witness for the
three-slope product.

The congruence class of the prime is the orientation boundary:

```text
chi_13(-1)=+1                         because 13 == 1 mod 4. (52)
```

Consequently

```text
chi_13(-d)=chi_13(d),
C_13^*=C_13.                                                (53)
```

So (48) is a symmetric signed graph operator.  In contrast,

```text
chi_7(-1)=-1                         because 7 == 3 mod 4,  (54)
```

and a `chi_7` Cayley half-set can orient the Paley tournament on seven
points.  The Fano slope product has crossed from the skew/oriented prime
`7` to the symmetric prime `13`.  It retains the quadratic-character split
but destroys the direction.

An oriented cyclic half-set on `F_13` could restore a direction, but choosing
one is an extra gauge.  Neither the even energies (5), the collision table,
nor the affine predecessor torsor chooses it.  This is why the formal
similarity to a Fano/Paley tournament does not orient the LRC self-cospan.

## 7. Physical consequence and exact stopping boundary

For THM-2521's Hilbert-valued predecessor profile

```text
p_r(z)=q_r(z)-qbar(z) in L^2(G(z)dz),                        (55)
```

equation (6) is exactly the multiplicative `chi_7` contrast of the six
integrated Hamilton collision energies.  The theorem proves:

```text
the polarized signed Hamilton operator is invertible,
its diagonal has real signature (6,6),
and the diagonal can vanish on a nonzero rational predecessor profile.   (56)
```

THM-2521 normalizes the same profile by

```text
integral G(z)||p(z)||^2 dz=13D_13.
```

Applying (42a)--(42c) coefficientwise in its real Hilbert space gives the
canonical polarized detector

```text
R_chi=integral G(z)||A_tau p(z)||^2 dz
      >(169/10)D_13>0                         if D_13>0.      (56a)
```

The centred-delta hostile kills the one-application diagonal `Q_tau(p)`,
not this derived norm.  The boundary remains semantic: `A_tau p` is a
signed predecessor stencil, not a Boolean owner/source/deep event.

Thus positive collision drift does **not** force a nonzero scalar
multiplicative contrast.  The centred-delta control is already a pointwise
rational positive-drift witness to that failure.

The result also does not identify (2) with the Fano incidence that appears
in the scalar `5+3` LRC residual.  THM-1260 proves that one fully placed
sharp handoff fork is locally `chi_7`-surjective.  A physical Fano consumer
must retain incidence among several located handoff occurrences, seam
digits, and owner reuse.  The six static Hamilton labels in (5) do not
supply that chronology.

Specifically, THM-2523 does **not**:

1. turn the signed centred potentials into Boolean owner or arrival events;
2. identify predecessor, source, arrival, old-deep, and future-owner sheets;
3. choose `tau` against `-tau` or orient a collision edge;
4. make the original one-application diagonal value nonzero from positive
   `D_13`, or realize the positive derived norm (56a) as one Boolean event;
5. exclude a scalar row or prove LRC(14).

The exact gain is a faithful algebraic boundary.  The multiplicative probe
is neither absent nor automatically positive: it is a lossless polarized
operator whose one-profile diagonal is an indefinite quotient.  The correct
next object is the translated/coupled cross-pairing, followed by a lawful
Boolean and chronological transplant, not another attempt to force a sign
from (1).

## 8. Exact dependency-free referee

Run

```bash
python3 04-computation/lrc14_chi7_hamilton_split_thm2523.py
python3 -O 04-computation/lrc14_chi7_hamilton_split_thm2523.py
```

Both executions reproduce

```text
05-knowledge/results/lrc14_chi7_hamilton_split_thm2523.out
```

byte-for-byte.  The script uses only exact integer/rational arithmetic and
the finite field `F_547`, which contains primitive seventh and thirteenth
roots.  It verifies:

- all twelve `13 x 13` Hamilton matrices, their exact energy expansion,
  symmetry, row-sum kernel, and rational rank `12`;
- `2,028` coefficientwise dilation anti-isometry identities;
- the group-ring identity (17), the `13`-Eisenstein test, and the exact
  rational root-isolation signs (22);
- the Fano difference-set law, (25)--(26), the slope product (48), and the
  Paley square (50);
- the centred-delta energy and its twelve nonzero translated cross-values;
- the fixed and anti-fixed totally isotropic spaces, including determinant
  `-4160` in (40); and
- as an independent forward cross-check for THM-2524, `468` exact
  position-space translated-polarization identities and `144` finite-field
  Fourier multiplier identities.

The finite-field run checks identities and supplies a hostile reduction.  It
is not the proof of real nonvanishing or signature: those follow from the
integer group-ring identity, Eisenstein irreducibility, and the rational sign
table above.

No computation in the referee assumes Boolean emission, source/deep
alignment, an oriented chart, or an LRC row exclusion.

An independent line audit rederived the energy/operator sign and Fourier
normalization, the Eisenstein minimal polynomial and squared characteristic
polynomial, the real signature, the dilation anti-isometry, both rational
isotropic spaces and their cross pairing, and the Paley product and square.
It also reproduced the normal and optimized transcripts byte-for-byte and
confirmed the recorded hashes.

The rational Witt and coercive additions have a separate dependency-free
referee:

```bash
python3 04-computation/lrc14_chi7_hamilton_witt_referee_thm2523.py
python3 -O 04-computation/lrc14_chi7_hamilton_witt_referee_thm2523.py
```

Both runs byte-match

```text
05-knowledge/results/lrc14_chi7_hamilton_witt_referee_thm2523.out
```

It independently reconstructs `A_tau` over `Q`, verifies the orthogonal
six-plus-six reversal decomposition, the three-dimensional even and
two-dimensional odd hyperbolic pieces, determinant squareclasses `-1` and
`-13`, rational Witt index five, all `2,028` slope-conjugacy entries, and
the exact positive inertia of `10A_tau^2-13I` on augmentation.  It also
checks the `169/10` THM-2521 normalization on a finite rational control.
No Boolean or semantic input is assumed.

**QED.**
