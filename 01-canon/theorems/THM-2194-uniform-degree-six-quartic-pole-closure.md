---
id: THM-2194
title: "Uniform degree-six quartic pole closure"
status: >
  PROVED. Let a planar Keller component be
  P=V^2z^4+Vbz^3+gamma z^2+delta z+epsilon after THM-2180, and suppose
  the reduced mate has fibre degree six. Then V divides 4gamma-b^2,
  without any square or nonsquare assumption on V. At a hypothetical bad
  place, the boundary/flux triple first forces the exact square cage of
  THM-2189. An exact five-row Faber phase bank then eliminates the degree-five
  and degree-three odd seeds. The residual E_6+uE_2+vE_1 bank has six
  valuation supports; a cubic-cusp identity and a two-flux comparison make
  each support impossible. Thus the canonical quadratic approximate root is
  polynomial for reduced degrees two and six on either deck. Every remaining
  finite-pole survivor is split, has leading coefficient W^4, and has reduced
  mate degree at least ten. The terminal nonmonic square-prefix step remains
  open.
source: codex-2026-07-24-split-degree-six-pole
depends_on:
  - THM-2129
  - THM-2158
  - THM-2180
  - THM-2189
related:
  - THM-2181
script: 04-computation/jc2_quartic_degree6_phase_bank_codex_20260724.py
output: 05-knowledge/results/jc2_quartic_degree6_phase_bank_codex_20260724.out
script_sha256: af27371539b4cf4cbd41751736c118df03103a4cee6d4039a2b1f003a2c7e83f
output_sha256: a7111755ee9783490e7b02145efcdc2fc00919cc6880ac1e1fe3e19beef5fe80
hash_basis: working-tree bytes (LF)
---

# THM-2194 -- the degree-six phase bank closes the quartic pole

Let

```text
P=V^2z^4+beta z^3+gamma z^2+delta z+epsilon
```

belong to a planar polynomial Keller pair over `C`. Reduce the mate by
polynomial target shears and suppose its remaining fibre degree is

```text
n=6.                                                   (1)
```

THM-2180 gives `beta=Vb`, with `b in C[x]`. Put

```text
D=4gamma-b^2.                                          (2)
```

Then, with no assumption on whether `V` is a square in `C(x)`,

```text
V divides D.                                           (3)
```

By THM-2158, this is exactly the polynomiality of

```text
H_0=Vz^2+(b/2)z+D/(8V).                                (4)
```

## 1. The square cage is uniform

Fix an irreducible `pi|V` and suppose for contradiction that, for the
extended `pi`-valuation `nu`,

```text
e=nu(V)>0,                d=nu(D)<e.                   (5)
```

The square-cage argument in THM-2189, Section 2, does not need deck parity
at its initial negative scale. The degree-six top Faber term is uniquely
deeper than every term of degree at most five; boundary regularity and
constancy of both fluxes therefore force the same three consecutive
boundary coefficients to vanish. THM-2129's classification gives the square
face exactly as before.

Write `f=nu(b)` and `ell_0=nu(delta)`. The resulting inequalities are

```text
2f<d<e,              nu(gamma)=2f,
ell_0+e>3f.                                           (6)
```

Choose `U^2=V`, on either the field or split deck, and put

```text
a=b/(2U),                 A=-nu(a)=e/2-f>0,
C=d-2f>0,                 H=e-d>0.                    (7)
```

Then

```text
A=(C+H)/2,                         0<C<2A.             (8)
```

In the canonically scaled reverse variable `X=au`, THM-2189 gives

```text
T_hat
 =1+2X+(1+c)X^2+lX^3+mX^4
 =K^2+Lambda X^3+Omega X^4,                           (9)

K=1+X+(c/2)X^2,
c=D/b^2,
l=8delta V/b^3,
m=16epsilon V^2/b^4,

Lambda=l-c,
Omega=m-c^2/4.                                       (10)
```

The exact valuations needed below are

```text
nu(c)=C,                  nu(Lambda)>0 or Lambda=0,
nu(m)>=4A,
nu(Omega)=2C.                                             (11)
```

In particular, if

```text
mu=Omega+c^2/4=m,                                      (12)
```

then `nu(mu)>=4A`. This last identity retains information which would be
lost by treating `Omega` as an arbitrary coefficient of valuation `2C`.

## 2. The exact degree-six Faber bank

For every reduced degree `j`, define its normalized boundary and two flux
coefficients by

```text
B_j=[X^j]T_hat^(j/4),
F_j=[X^(j+1)]T_hat^(j/4),
G_j=[X^(j+2)]T_hat^(j/4)
       +(1/2)[X^(j+1)]T_hat^(j/4).                    (13)
```

The translation point is `h=a/2`. Thus the exact contribution of `E_j` to
the polynomial boundary, first flux divided by four, and second flux divided
by four is respectively

```text
a^j B_j,                 a^(j+1)F_j,
a^(j+2)G_j.                                             (14)
```

The coefficient recurrence from

```text
T_hat (T_hat^(j/4))'=(j/4)T_hat' T_hat^(j/4)
```

gives the following exact bank:

```text
B_1=1/2,
F_1=(2c-1)/8,
G_1=Lambda/4;                                           (15)

B_2=c/2,
F_2=Lambda/2,
G_2=-(Lambda-2Omega)/4;                                (16)

B_3=(12Lambda+6c-1)/16,
F_3=3(-16Lambda+32Omega+4c^2-4c+1)/128,
G_3=-3Lambda(2c-1)/32;                                 (17)

B_5=(80Lambda c-40Lambda+160Omega
       +60c^2-20c+3)/256,

F_5=5(32Lambda^2-32Lambda c+16Lambda
       +64Omega c-32Omega
       +8c^3-12c^2+6c-1)/1024,

G_5=-5Lambda(16Lambda-32Omega
       +4c^2-4c+1)/512;                                (18)

B_6=(3Lambda^2+6Omega c+c^3)/8,
F_6=-3Lambda(Lambda-2Omega)/8,

G_6=3((Lambda-2Omega)^2
       +Lambda^2(1-2c))/32.                            (19)
```

Because `c,Lambda,Omega` have positive valuation, `B_3,F_3,B_5,F_5`
are units, while

```text
nu(G_1)=nu(G_3)=nu(G_5)=nu(Lambda)                    (20)
```

when `Lambda!=0`. Formula (19), rather than only its first valuation, is the
sidecar which separates the two odd-seed cancellation phases.

Normalize the nonzero coefficient of `E_6` in the Faber normal form to one:

```text
Q=J(P)+E_6+k_5E_5+k_3E_3+k_2E_2+k_1E_1,
J in C[T],                       k_i in C.             (21)
```

Terms of degree four are target shears and have already been absorbed into
`J(P)`.

## 3. Two valuation lemmas

Put

```text
L=nu(Lambda),
R=nu(Lambda-2Omega),                                  (22)
```

with the convention that the valuation of zero is infinity.

### Cubic-cusp lemma

If

```text
nu(B_6)>3C                                             (23)
```

including `B_6=0`, then

```text
2L=3C,
R=L,
nu(F_6)=nu(G_6)=3C.                                   (24)
```

Indeed, (12) and (19) give the exact identity

```text
8B_6=3Lambda^2-(1/2)c^3+6mu c.                        (25)
```

The last term has valuation at least

```text
4A+C>3C                                               (26)
```

by (8). Therefore (23) forces the first two terms of (25) to have equal
valuation, so `2L=3C`. In particular `L=3C/2<2C=nu(Omega)`, whence
`R=L`. The leading residues of `Lambda-2Omega` and `Lambda` agree. The
formulas for `F_6,G_6` in (19) now give (24); the two leading squares in
`G_6` add with coefficient two.

### Equal boundary/first-flux lemma

Suppose, for some

```text
g in {1,3,5},
```

that

```text
nu(B_6)=nu(F_6)=gA<infinity.                           (27)
```

Let `S=nu(G_6)`. Then

```text
S<L+gA.                                                (28)
```

Moreover the top second-flux term is strictly deeper than the degree-two
second-flux term:

```text
nu(a^8G_6)<nu(a^4G_2).                                (29)
```

From `F_6` in (19),

```text
L+R=gA.                                                (30)
```

If `L<R`, formula (19) gives `S=2L`; if `R<L`, it gives `S=2R`.
If `L=R<2C`, then `Lambda-2Omega` and `Lambda` have the same residue, so
again `S=2L`. The case `L=R>2C` is impossible because their difference
has valuation `2C`. Finally, if `L=R=2C`, equations (19) and (25) give

```text
nu(B_6)=3C,                    nu(F_6)=4C,             (31)
```

contrary to (27). These cases prove (28).
They also show `S<=L+R=gA`, so the top second-flux valuation
`S-8A` is strictly negative for every allowed `g`.

For (29), use `nu(G_2)=R`. If `L<=R`, then `L<=gA/2` and

```text
nu(a^8G_6)-nu(a^4G_2)
 =3L-(g+4)A
 <=(g/2-4)A<0.                                       (32)
```

If `R<L`, the same difference is

```text
R-4A<(g/2-4)A<0.                                     (33)
```

This proves both assertions.

## 4. The degree-five and degree-three seeds vanish

Suppose first that the highest nonzero odd coefficient in (21) is

```text
k_j,                    j in {5,3}.
```

Put

```text
g=6-j,                         g in {1,3}.             (34)
```

By (17)--(18), both `B_j` and `F_j` are units. After multiplication by the
powers of `a` in (14), this term is uniquely deepest among all lower Faber
degrees in the boundary and first flux. Boundary regularity and first-flux
constancy can therefore hold only if

```text
nu(B_6)=nu(F_6)=gA.                                   (35)
```

The equal boundary/first-flux lemma applies. The highest odd second-flux
term has valuation

```text
L-(j+2)A=L-(8-g)A.                                    (36)
```

The difference between the top and this term is

```text
[S-8A]-[L-(8-g)A]=S-L-gA<0                            (37)
```

by (28). Every smaller odd term is shallower by an additional even multiple
of `A`. Equation (29) also places the top strictly below the possible
degree-two term. Thus the top second-flux pole is unique, contradicting
`Psi_Q in C`.

Consequently

```text
k_5=k_3=0.                                            (38)
```

This is the first point where the split deck is genuinely handled: odd
seeds are not discarded by parity but are forced out by the ordered
boundary/first-flux/second-flux bank.

## 5. Exhaust the residual support

Write the remaining reduced part as

```text
Q-J(P)=E_6+uE_2+vE_1,               u,v in C.         (39)
```

Nonzero constants have valuation zero. We now exhaust their four support
patterns, with the final pattern split by the ordering of `C` and `A`.

### Neither lower term is present

If `u=v=0`, boundary regularity forces

```text
nu(B_6)>=6A.                                          (40)
```

Since `6A>3C`, the cubic-cusp lemma gives `nu(F_6)=3C<7A`.
The top first flux is therefore polar and has no lower term to cancel it.

### Only the linear term is present

If `u=0` and `v!=0`, the degree-one boundary and first-flux coefficients
are units. Their only possible cancellations force

```text
nu(B_6)=nu(F_6)=5A.                                   (41)
```

The equal boundary/first-flux lemma with `g=5` makes the top second flux
strictly deeper than the degree-one term, a contradiction.

### Only the quadratic term is present

If `u!=0` and `v=0`, boundary regularity forces

```text
nu(B_6)=C+4A>3C.                                      (42)
```

The cubic-cusp lemma gives

```text
L=3C/2,                         nu(F_6)=3C.            (43)
```

The difference between the top and degree-two first-flux valuations is

```text
[3C-7A]-[3C/2-3A]
 =3C/2-4A<0.                                          (44)
```

The top pole is unique.

### Both lower terms are present

Assume `u!=0` and `v!=0`.

If `C<A`, the degree-two boundary term, of valuation `C-2A`, is strictly
deeper than the degree-one term, of valuation `-A`. Thus (42)--(43) hold
again. The top first flux is deeper than the degree-two term by (44), and
deeper than the degree-one term by

```text
[3C-7A]-[-2A]=3C-5A<0.                               (45)
```

If `C=A`, the two lower boundary poles tie at `-A` and may cancel each
other. To prevent a unique top pole below them one must have

```text
nu(B_6)>=5A>3C.                                       (46)
```

The cubic-cusp lemma now gives

```text
L=3A/2,                         nu(F_6)=3A.            (47)
```

The top first flux has valuation `-4A`; the degree-two and degree-one terms
have valuations `-3A/2` and `-2A`, respectively. Again the top pole is
unique. This is the apparent self-cancellation resonance; the first flux
removes it.

It remains that `C>A`. The degree-one boundary term is now strictly deeper
than the degree-two term, so boundary regularity forces

```text
nu(B_6)=5A.                                           (48)
```

If `L<=A`, then `2C>2A>=2L`, so `R=L` and
`nu(F_6)=2L`. The top first flux is deeper than the degree-two and
degree-one terms by, respectively,

```text
L-4A<0,                         2L-5A<0.              (49)
```

If `L>A`, the degree-one term is strictly deepest among the lower
first-flux terms. Cancellation therefore forces

```text
nu(F_6)=5A.                                           (50)
```

Equations (48) and (50) invoke the equal boundary/first-flux lemma with
`g=5`. Its top second flux is strictly deeper than both the degree-one term,
by (28), and the degree-two term, by (29). This is the final contradiction.

All supports in (39) are impossible. Therefore the assumption (5) was
false at every divisor of `V`, proving (3).

## 6. Exact consequence and boundary

THM-2189 proves (3) for every nonsplit quadratic deck and, uniformly on
either deck, for reduced mate degree two. The present theorem adds degree
six on either deck. Therefore every remaining finite-pole survivor must
satisfy

```text
V=W^2 in C[x],                 V^2=W^4,
reduced mate degree n>=10.                                  (51)
```

This closes the pole, not the planar Jacobian conjecture. Even after (3),
the polynomial approximate root (4) has leading coefficient `V`, and the
terminal nonmonic square-prefix/quadratic-member descent is still required.
For split degrees at least ten, both that terminal issue and the finite-pole
congruence remain open. No statement about general JC(2) or DC(2) is
claimed.

## 7. Exact referee

Run

```bash
python3 04-computation/jc2_quartic_degree6_phase_bank_codex_20260724.py
```

The first path independently derives all fifteen entries in (15)--(19) from
the coefficient recurrence and verifies the cubic-cusp identity (25). The
second path uses an exact `Fraction` Laurent solver on the explicit universe

```text
1<=A<=6,             1<=C<2A,             1<=E<=4A,

T_hat=1+2X+(1+t^C)X^2+lambda t^E X^3,

lambda in {-2,-1,-1/2,0,1/2,1,2}.                    (52)
```

All `4508` hostile rows reject simultaneous boundary and two-flux
regularity. The split degree-six family with `A=2,C=E=3,lambda=1` is a
named negative control. The exact square `T_hat=(1+X)^2` is a positive
algebraic control. The output labels this sweep `FINITE_REFEREE_ONLY`:
neither its bounds nor its single-monomial coefficients supply any
quantifier used in the proof above. QED.
