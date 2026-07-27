---
id: THM-2574
title: "Oriented tooth-component holonomy and fixed-frequency descent"
status: >
  PROVED + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.  A pulled-back proper interval under
  x -> kx+s/13 has k physical teeth whose labels acquire the monodromy
  n -> n-1 after one lift s -> s+13.  The rth component character selects
  exactly the physical frequencies X congruent to r modulo k.  Multiplying
  by its compensating connection produces a genuine C_13-periodic complex
  field with target phase (X-r)/k.  For a left frequency X and right
  frequency X+M, the plus-sign target DFT lands at the Euclidean carry
  c=(r+M-t)/k.  Standard residues give only floor(M/k) and floor(M/k)+1.
  If 13 divides M and k is a 13-unit, every carry is zero modulo 13 exactly
  on the resonance k|M; off resonance at least one nonzero carry class has
  a nonzero danger-to-safe fixed-X amplitude.  The same component character
  diagonalizes THM-2573's whole-layer boundary tooth weights.  A bare
  danger/safe gate has uniform weights, so its Abel normal vanishes when
  k does not divide M and, for M=kh, is
  k cos(pi hL/7)/pi^2 times zeta^(-hs); the live 13-divisible bank therefore
  gives only q=0.  A nonzero live normal must see a nonuniform
  (M mod k)-tooth coefficient of the actual total-layer handoff measure.
  The connection is complex and gauge-chosen: changing an integer character
  lift changes target colour while leaving the unconnected component field
  fixed.  No positive Boolean carrier, lawful THM-2334 charge, live-row
  noncancellation, row exclusion, or LRC(14) conclusion follows.
source: lrc-semantic-frontier-2026-07-28-oriented-tooth-holonomy
depends_on:
  - THM-2573-logarithmic-abel-normal-and-common-endpoint-jump-pairing
related:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2507-truncated-radon-toothpick-tomography-and-nonaffine-root-boundary
  - THM-2564-six-tooth-doubly-centered-tomography-and-section-holonomy-boundary
  - THM-2568-full-x-transition-annihilation-and-refined-pair-drift-boundary
  - THM-2569-stationary-diagonal-conditioned-paired-corner-and-frozen-future-role-boundary
script: 04-computation/lrc14_oriented_tooth_holonomy_thm2574.py
output: 05-knowledge/results/lrc14_oriented_tooth_holonomy_thm2574.out
script_sha256: 7f020a32937d35249c7e2e77c645cad3e4112e608528f1049f97a044400fd6a2
output_sha256: 6c3f024d13416f0384fedbdcc6a77adfdeb0f8d73c625eea541261ec9bc35dbc
hash_basis: working-tree bytes (LF)
---

# THM-2574 -- orient one physical tooth before summing frequency

**PROVED + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.**

THM-2568 shows that the completed danger-to-safe endpoint vanishes after
the full physical-frequency sum.  THM-2573 extracts a possible singular
normal from its total-layer boundary handoffs.  Both operations forget one
finite cyclic coordinate: which of the `k` physical teeth supplied a local
piece.

Keeping that coordinate exposes an exact local system:

```text
physical tooth label n in Z/k
  -> component character r in (Z/k)^hat
  -> monodromy around the thirteen target lifts
  -> compensating connection
  -> fixed-X target character = an integer Euclidean quotient.       (1)
```

For two endpoints the quotient difference is a carry.  For the Abel normal,
the same character is the discrete Fourier mode of the actual whole-layer
handoff weights.  This identifies a precise escape coordinate, but also its
sharp type boundary: the connected fields are complex and their target
colour depends on an oriented choice of connection lift.

## 1. The physical tooth atlas must be formed on the cover

Put

```text
p=13,                         zeta=exp(2 pi i/p),
k>=1.                                                              (2)
```

Let `J=(alpha,beta)` be a lifted real interval with

```text
0<beta-alpha<1.                                                   (3)
```

Endpoint conventions are immaterial.  For integers `n,s`, first define on
the universal cover

```text
I~_(n,s)(x)=1_J(kx+s/p-n),                    x in R.              (4)
```

It obeys

```text
I~_(n+k,s)(x+1)=I~_(n,s)(x).                                    (5)
```

Thus one should not truncate a single lift component to a permanently fixed
`[0,1)` domain.  Its physical image is the projected arc

```text
I_([n],s)=pi((J+n-s/p)/k) subset T,             [n] in Z/k,       (6)
```

and (5) makes (6) independent of the representative of `[n]`.  The `k`
arcs are disjoint up to endpoints.  Equivalently, every weighted sum below
can be formed first as a locally finite sum of (4); it is one-periodic and
therefore defines a function on `T`.

For the physical Fourier convention

```text
f_hat(X)=integral_T f(x)exp(-2 pi i Xx)dx,       X in Z,           (7)
```

put

```text
C_J(xi)=integral_J exp(-2 pi i xi y)dy.                         (8)
```

Parametrizing the projected arc in (6) by

```text
x=(y+n-s/p)/k
```

gives the exact component coefficient

```text
(1_(I_([n],s)))_hat(X)
 =1/k C_J(X/k)
    exp(2 pi i Xs/(pk))exp(-2 pi i Xn/k).                       (9)
```

Changing `n` by `k` multiplies the right side by `exp(-2 pi i X)=1`,
as required.  The target lift has nontrivial component monodromy:

```text
I_([n],s+p)=I_([n-1],s).                                      (10)
```

The complete Boolean interval pullback is periodic in `s mod p`; its
individual tooth labels are not fixed by that loop.

## 2. Component characters diagonalize both frequency and holonomy

For an integer character lift `r`, define

```text
F_(r,s)
 =sum_([n] in Z/k) exp(2 pi i r n/k)1_(I_([n],s)).              (11)
```

Only `r mod k` enters (11).  Summing (9) over the complete component cycle
gives

```text
F_(r,s)_hat(X)
 =C_J(X/k)exp(2 pi i Xs/(pk)),       X congruent to r mod k,

 =0,                                otherwise.                 (12)
```

Thus a tooth character is exactly a physical-frequency residue projector.
Equation (10) gives its holonomy

```text
F_(r,s+p)=exp(2 pi i r/k)F_(r,s).                              (13)
```

Choose the compensating connection

```text
G_(r,s)=exp(-2 pi i r s/(pk))F_(r,s).                          (14)
```

Equations (13)--(14) make `G_(r,s)` genuinely `p`-periodic in `s`.  On the
selected physical residue, write

```text
X=r+kQ.                                                       (15)
```

Then

```text
G_(r,s)_hat(X)=C_J(X/k)zeta^(Qs).                             (16)
```

With the plus-sign target transform used in THM-2568,

```text
Hhat(q)=1/p sum_(s in F_p)H(s)zeta^(qs),                      (17)
```

a pure phase `zeta^(Qs)` lands at `q=-Q`.  This sign will reverse once the
right endpoint is conjugated.

## 3. A paired fixed-frequency target colour is a carry

Take lifted intervals `J_L,J_R` for two copies of the same physical speed
`k`.  Fix an integer physical offset `M` and a left frequency `X`.  Choose
the standard representatives

```text
r in {0,...,k-1},             r congruent to X mod k,

t in {0,...,k-1},             t congruent to X+M mod k,       (18)
```

and define the integer carry

```text
c=(r+M-t)/k.                                                   (19)
```

The connected fixed-frequency pair is

```text
K_X(s)
 =G^L_(r,s)_hat(X)
   conjugate(G^R_(t,s)_hat(X+M)).                              (20)
```

Writing

```text
Q_L=(X-r)/k,                  Q_R=(X+M-t)/k=Q_L+c             (21)
```

in (16) gives

```text
K_X(s)
 =C_(J_L)(X/k)conjugate(C_(J_R)((X+M)/k))zeta^(-cs).          (22)
```

Consequently the transform (17) is supported at exactly

```text
q=c mod p,                                                     (23)
```

provided the two displayed interval amplitudes are nonzero.  The sign in
(23) is load-bearing: the left phase minus the conjugated right phase is
`-c`, and the plus-sign DFT lands at `+c`.

Now use Euclidean division

```text
M=ak+d,                         0<=d<k.                       (24)
```

As `r` runs through the standard residues, (19) gives exactly

```text
c=a,                            k-d times,

c=a+1,                          d times,                       (25)
```

where the second line is absent when `d=0`.  The target colour is therefore
the wrap bit of the addition arrow

```text
r + (M mod k) -> t                                                (26)
```

together with the coarse quotient `floor(M/k)`.  This is the same
nonaffine toothpick carry which appears geometrically in THM-2507 and
THM-2564, now represented directly at a physical Fourier pair.

## 4. The exact live-bank resonance split

Assume

```text
p|M,                         gcd(k,p)=1.                       (27)
```

These are the relevant arithmetic types when `M=m c_3` and `k` is one of
the six guard/unit speeds: the deepest blocker is `13`-divisible and the
guard/unit speeds are `13`-units.

If `k|M`, then (24) has `d=0` and

```text
c=M/k congruent to 0 mod p.                                  (28)
```

Conversely, if `k` does not divide `M`, both adjacent values in (25) occur.
Two adjacent integers cannot both vanish modulo `p`, so at least one
physical residue has

```text
c not congruent to 0 mod p.                                  (29)
```

Hence

```text
all standard carry classes are target-trivial
  iff k|M;                                                     (30)

k does not divide M
  -> some standard physical residue has nonzero target colour. (31)
```

Equation (31) is only a character statement.  The interval amplitudes in
(22) still have to avoid their sinc zeros.

## 5. Danger and safe amplitudes cannot kill every lift

For `L in {1,2}`, use the lifted danger and safe intervals

```text
D_L=(-L/14,L/14),                 |D_L|=L/7,

U_L=(L/14,1-L/14),                |U_L|=(7-L)/7.             (32)
```

Their fractional-frequency transforms are

```text
C_(D_L)(xi)
 =sin(pi L xi/7)/(pi xi),

C_(U_L)(xi)
 =exp(-pi i xi)sin(pi(7-L)xi/7)/(pi xi),                    (33)
```

with their positive interval lengths at `xi=0`.  Therefore

```text
C_(D_L)(X/k)=0
  iff X!=0 and 7k divides LX,

C_(U_L)(Y/k)=0
  iff Y!=0 and 7k divides (7-L)Y.                            (34)
```

Fix any residue `r mod k` and put

```text
X_n=r+nk,                    Y_n=X_n+M,       n mod 7.       (35)
```

Because both `L` and `7-L` are units modulo seven,

```text
gcd(Lk,7k)=gcd((7-L)k,7k)=k.                                (36)
```

Each divisibility condition in (34), if soluble in `n`, therefore removes
exactly one residue modulo seven.  The exceptional zero physical frequency
has a nonzero interval transform and can only remove less.  At least

```text
7-1-1=5                                                       (37)
```

choices of `n mod 7` make both amplitudes in (22) nonzero.

Combining (31) and (37) proves the sharp qualitative control:

> If `p|M`, `gcd(k,p)=1`, and `k` does not divide `M`, there are standard
> component residues `r,t` and an integer `X` for which the connected
> danger-to-safe fixed-frequency pair (20) is nonzero and lies in a
> nonzero target colour.

This is an existence theorem for the isolated complex component pair.  It
does not say that the coefficient remains nonzero after convolution with a
complete word, owner, stationary-return packet, or other endpoint factors.

## 6. The same local system diagonalizes the Abel handoff teeth

THM-2573 says that a nonnegative disjoint whole-layer endpoint has
logarithmic Abel normal

```text
N_s(M)=1/(2 pi^2) integral_T exp(2 pi i Mx)dnu_s(x),          (38)
```

where `nu_s` is its positive **total-layer** common-boundary handoff measure.
Suppose a lawful whole-layer orbit has handoff weights on the `2k` boundary
teeth of the `k`-speed danger gate.  Write

```text
x_(j,epsilon,s)
 =(j+epsilon L/14-s/p)/k mod 1,

j in Z/k,                       epsilon in {+1,-1},           (39)

w_(j,epsilon,s)
 =nu_s({x_(j,epsilon,s)})>=0.                                (40)
```

Other common jumps, if present, form a separate part of `nu_s`.  Coincident
factor boundaries must already have been merged into the total-layer weight
(40), as required by THM-2573.

Lawful target covariance around the lifted target loop is exactly

```text
w_(j,epsilon,s+p)=w_(j-1,epsilon,s).                         (41)
```

Define the tooth transform and its connected version

```text
W_(r,epsilon)(s)
 =sum_(j mod k)exp(2 pi i rj/k)w_(j,epsilon,s),

H_(r,epsilon)(s)
 =exp(-2 pi i rs/(pk))W_(r,epsilon)(s).                      (42)
```

Equations (41)--(42) give the same holonomy and connection as
(13)--(14), so `H_(r,epsilon)` descends to `C_p`.

For the standard split

```text
M=r+hk,                         0<=r<k,                       (43)
```

the gate-boundary contribution to (38) is exactly

```text
N^gate_s(M)
 =zeta^(-hs)/(2 pi^2)
   sum_(epsilon=+-1)
     exp(2 pi i M epsilon L/(14k))H_(r,epsilon)(s).          (44)
```

Thus the physical residue `r=M mod k` selects the matching tooth character,
while the integer quotient `h=(M-r)/k` shifts its target spectrum.  Formula
(44) is the common interface between THM-2573's normal and the fixed-`X`
component atlas.  It replaces the vague instruction “retain an oriented
boundary” by one exact obligation:

```text
prove a nonzero r=(M mod k) tooth coefficient
of the actual covariant total-layer handoff weights.         (45)
```

### The uniform bare gate is a sharp hostile

For the bare pair

```text
P_s(x)=d_L(kx+s/p),                  R_s=1-P_s,               (46)
```

every boundary in (39) has handoff weight one.  Hence all nontrivial tooth
characters vanish.  Directly summing (38) gives

```text
N_s(M)=0,                                      k does not divide M,

N_s(kh)
 =k cos(pi hL/7)/pi^2 zeta^(-hs),             k divides M.  (47)
```

The cosine in (47) never vanishes for integral `h` and `L in {1,2}`: a zero
would require the even integer `2hL` to be congruent to `7` modulo `14`.
Under the plus-sign DFT, the surviving colour is

```text
q=h mod p.                                                       (48)
```

Under the live arithmetic hypotheses (27), `k|M` forces `p|h`.  Therefore
the pure danger/safe gate has

```text
Nhat(q;M)=0                         for every q!=0             (49)
```

throughout the allowed deep bank.  Off resonance its teeth cancel; on
resonance its only colour is zero.  A nonzero live Abel normal must therefore
come from nonuniform whole-layer handoff weights in (45), from other common
jumps, or from a different sidecar.  Merely exposing the bare gate boundary
does not defeat THM-2568.

## 7. Gauge dependence and the Boolean boundary

The component character in (11) depends only on `r mod k`, but the
connection does not.  For every integer `a`,

```text
F_(r+ak,s)=F_(r,s),

G_(r+ak,s)=zeta^(-as)G_(r,s).                                (50)
```

For a paired current, changing the left and right lifts by `ak` and `bk`
changes

```text
c -> c+a-b.                                                   (51)
```

A common change `a=b` is harmless; a relative change shifts the target
colour.  The standard section `0<=r,t<k` is a concrete oriented
trivialization, not information contained in the unconnected fields
`F_(r,s),F_(t,s)` themselves.  This is the sharp gauge hostile: the same
physical component characters acquire different target colours under
different relative connection lifts.

Forgetting tooth labels keeps only the trivial character

```text
F_(0,s)=sum_([n] in Z/k)1_(I_([n],s)),                         (52)
```

which is the original Boolean `k`-tooth pullback and has physical Fourier
support only on `kZ`.  Simultaneous left support at `X` and right support at
`X+M` then forces `k|M`; under (27) its carry is zero modulo `p`.  Therefore
the nonzero off-resonant modes proved in Section 5 necessarily use a
nontrivial complex component character.  They are not positive subevents of
the Boolean gate.

Changing the integer lift of the base interval `J` merely relabels teeth and
multiplies a component character by a constant phase.  It does not repair
the relative connection ambiguity in (50)--(51).

## 8. Exact consequence for the live LRC seam

The new structural ledger is

```text
full Boolean endpoint + ordinary full-X sum
  -> zero in every target colour;                      THM-2568

covariant whole-layer endpoint + logarithmic Abel normal
  -> total-layer handoff measure;                      THM-2573

oriented tooth component + fixed X
  -> exact Euclidean-carry target colour;              this theorem

covariant weighted handoff teeth
  -> exact connected tooth transform (44);             this theorem

uniform bare gate in the live deep bank
  -> zero in every nonzero target colour.               this theorem (53)
```

THM-2569 supplies positive common ancestry, an old dangerous `k_a` gate, a
later safe `k_a` repair, a deepest probe, and both moving dipoles on one
packet.  It does **not** yet supply either of the two additional facts needed
to use (22) or (44):

1. a covariant target orbit of the complete old selector and future packet;
2. a nonzero connected tooth coefficient of the resulting total-layer
   handoff weights at the residue selected by the allowed deep offset.

Multiplying (11) by the remaining endpoint factors convolves physical
frequencies, so the isolated nonzero control of Section 5 is not inherited
formally.  Freezing the target-informed selector can also create an auxiliary
colour which is not a THM-2334 relation charge.  The next positive theorem
must build the lawful whole-layer orbit first and then prove (45), or prove a
compatible fixed-`X` coefficient after all word/owner/deep factors are
inserted.

No equality with the relation residue `eta.u`, semantic arrival, aggregate
noncancellation, scalar-row exclusion, or LRC(14) conclusion is asserted.

## 9. Exact companion

Run

```bash
python3 04-computation/lrc14_oriented_tooth_holonomy_thm2574.py
python3 -O 04-computation/lrc14_oriented_tooth_holonomy_thm2574.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_oriented_tooth_holonomy_thm2574.out.
```

The dependency-free referee uses integer phase numerators and exact
divisibility only.  It checks:

- lifted-component monodromy and closure modulo `k`;
- component-character orthogonality and physical residue selection;
- the compensating `C_13` connection, plus-sign DFT convention, and integer
  lift gauge shift;
- the paired carry sign and both exact multiplicities in (25);
- the live `13`-adic resonance equivalence (30);
- the exact danger/safe sinc-zero criterion and the five-of-seven control;
- the connected tooth-weight/Abel-normal phase in (44), the pure-gate
  support, and the nonvanishing cosine; and
- thousands of off-resonant nonzero-colour fixed-`X` controls.

The universal formulas and all quantifiers are proved symbolically above;
the finite census is a hostile sign, gauge, zero-set, and boundary audit, not
a replacement for those proofs.

**QED.**
