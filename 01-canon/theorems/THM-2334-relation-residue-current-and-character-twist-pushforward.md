---
id: THM-2334
title: "Relation-residue current and character-twist pushforward"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. The
  27-factor word/deepest-comb/bare Fourier expansion has a canonical
  exact-address pushforward: for every fixed relation address its whole
  gauge-orbit coefficient is absolutely convergent even without Abel
  damping. At positive Abel radius the pushforward is l1 on the relation
  lattice and its total mass is the full marked mixed current. Reduction
  modulo any N is the inverse finite Fourier transform of coordinate-
  translated copies of that current. For the THM-2309 target quotient this
  gives exactly 169 twists, every full-word target-fibre Abel limit exists,
  and at least one fibre survives. A nonzero target vector survives iff
  those 169 twisted currents are nonconstant; the squared nonzero-target
  mass is their exact variance. The transported word is target-charge
  neutral modulo 13 but active modulo 7, so the all-91-unit target mask is
  a genuinely coupled CRT projector. A strictly positive two-factor
  hostile has nonzero full current and six nonzero terms in one exact
  address orbit whose sum is zero. No nonzero target gain, all-unit
  aggregate, visible address, terminal phase, scalar-row exclusion, or
  LRC(14) closure is claimed.
source: codex-2026-07-25-relation-residue-pushforward
depends_on:
  - THM-2302-same-label-expiration-dichotomy-and-pure-terminal-shell-no-go
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2327-two-colour-marked-unit-c3-triangle
related:
  - THM-2059-crt-fiber-product-phase-packet
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2319-crt-unit-bispectrum-needle-and-mixed-polarization-no-go
  - THM-2321-prescribed-root-character-bispectrum-slice-positivity
  - THM-2325-prescribed-target-gain-full-lattice-91-unit-needle-bank
  - THM-2331-two-sided-septimal-address-embedding-in-marked-current
  - THM-2333-abel-target-fibre-sum-landing-and-zero-fibre-boundary
script: 04-computation/lrc14_relation_residue_pushforward_thm2334.py
output: 05-knowledge/results/lrc14_relation_residue_pushforward_thm2334.out
script_sha256: ce220ede12175b6851810782c880f95048fe9f4643cc4f52f47a7f4d8dcb7b0c
output_sha256: d2d9b49db9ef3eabf7e3ae17cea247da554a9f8df2abfc3907243d317b21fec1
hash_basis: working-tree bytes (LF)
---

# THM-2334 -- the current is a measure on relation-address orbits

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2331 proves that every prescribed all-`91`-unit target address occurs
as a nonzero term in the marked mixed current. That is a support theorem.
The simultaneous relation-lattice shift used there also shows why raw term
occurrence is non-discriminating: once one term occurs, a large gauge orbit
of decompositions carries the same address.

The faithful coefficient is the sum over that entire orbit. This theorem
constructs it and then Fourier-transforms the address variable:

```text
27 atomic endpoint harmonics
  -> exact relation address r
  -> absolutely convergent orbit coefficient C(r)
  -> residue current on K_N
  -> coordinate-twist character transform.                         (1)
```

For the two-dimensional target quotient, the transform has only `169`
values. Its nonconstant part is exactly the total squared coefficient mass
on nonzero target vectors. This turns the remaining noncancellation problem
into a finite functional statement, rather than another address-support
problem.

## 1. The 27-factor marked current

Let

```text
w=(w_0,...,w_(d-1)) in Z^d
```

have nonzero coordinates and be primitive. Put

```text
Lambda=ker(w:Z^d->Z),
e_c=(0,...,0,1,0,...,0),                            (2)
```

where `c` is the deepest-speed coordinate. Let `R>=1` be an integer.
For each coordinate let

```text
I_i,J_i:T->R
```

be bounded real functions of bounded variation. In the LRC application,
the nine `I_i` are the interval/complement factors of the exclusive owner
rectangle, the nine `J_i` are the factors of the terminal word, and

```text
R=13^(lambda_j+1).                                  (3)
```

Write

```text
a_i(n)=(I_i)_hat(n),
b_i(n)=(J_i)_hat(n),

a(u)=prod_i a_i(u_i),
b(beta)=prod_i b_i(beta_i).                         (4)
```

The present and transported rectangles are

```text
E(t)=prod_i I_i(w_i t),
W(t)=prod_i J_i(R w_i t),
E_Q(t)=E(t)W(t).                                    (5)
```

Thus `E_Q` has eighteen atomic factors and the conjugate endpoint `E` has
nine: these are the 27 endpoint factors. The deepest danger comb is a
separate fixed Fourier leg. Let its base coefficient at multiplier `m` be

```text
delta_hat(m),
```

and fix an ordinary frequency triangle

```text
Y=X+m w_c.                                         (6)
```

For `0<rho<1`, Abel-smooth every **base** coefficient:

```text
a_(i,rho)(n)=rho^|n| a_i(n),
b_(i,rho)(n)=rho^|n| b_i(n),
delta_hat_rho(m)=rho^|m| delta_hat(m).              (7)
```

Use `a_rho(u)` and `b_rho(beta)` for the products. The left and right
regularized coefficients are

```text
L_rho(X)
 =sum_((u+R beta).w=X) a_rho(u)b_rho(beta),

F_rho(Y)
 =sum_(v.w=Y) a_rho(v).                             (8)
```

Every sum in (8) is absolutely convergent. The regularized marked current
is

```text
J_rho(X,m)
 =delta_hat_rho(m)L_rho(X)conjugate(F_rho(Y)).      (9)
```

## 2. Exact addresses are gauge orbits

One atomic term in (9) has an exact relation address

```text
r=u+R beta+m e_c-v in Lambda.                       (10)
```

For fixed `r`, define its grouped coefficient

```text
C_rho(r;X,m)
 =delta_hat_rho(m)
  sum_((u+R beta).w=X)
    a_rho(u)b_rho(beta)
    conjugate(
      a_rho(u+R beta+m e_c-r)
    ).                                             (11)
```

The right index in (11) automatically has dot product `Y` with `w`.

There is an exact gauge interpretation. Put

```text
G_R(w)
 ={(p,q) in Z^d x Z^d:(p+R q).w=0}.                (12)
```

It acts on decompositions by

```text
(u,beta,v)
 ->(u+p,beta+q,v+p+R q).                            (13)
```

This preserves both frequency equations and the address (10). Conversely,
two decompositions have the same address exactly when their difference is
of the form (13). Hence the decomposition fibre over every `r` is one free
transitive `G_R(w)`-torsor, and (11) is its weighted orbit sum.

THM-2331 uses the subgroup

```text
(p,q)=(p,0),                 p in Lambda,           (14)
```

to move one term away from all septimal Fourier zeros. It proves that many
weights in a prescribed orbit are nonzero. It does not determine the orbit
sum (11).

## 3. A fixed exact-address coefficient is absolutely convergent at rho=1

Bounded variation gives constants `A_i,B_i` such that

```text
|a_i(n)|<=A_i/(1+|n|),
|b_i(n)|<=B_i/(1+|n|).                              (15)
```

The key two-dimensional summability fact is:

> **Three-form lemma.** For fixed nonzero `R` and integer `d`,
>
> ```text
> sum_(p,q in Z)
>  1/[(1+|p|)(1+|q|)(1+|p+R q+d|)]
>  <infinity.                                       (16)
> ```

To prove it, sort the three affine forms

```text
p, q, p+R q+d
```

into dyadic sizes `2^j,2^k,2^ell`. Any two determine `(p,q)` up to a
bounded number of choices, so the number of lattice points in that block
is at most a constant times

```text
min(2^(j+k),2^(j+ell),2^(k+ell)).
```

After multiplying by the reciprocal weight, the block contributes at most
a constant times

```text
2^(-max(j,k,ell)).
```

There are `O(M^2)` triples with maximum `M`, and
`sum_M M^2 2^(-M)<infinity`. This proves (16).

In (11), put

```text
d_i=m(e_c)_i-r_i.
```

Drop the single frequency constraint and apply (16) independently in all
`d` coordinates. The absolute majorant factorizes:

```text
sum_(u,beta)
 |a(u)b(beta)a(u+R beta+m e_c-r)|

 <=prod_i A_i^2 B_i
   sum_(p,q)
    1/[(1+|p|)(1+|q|)(1+|p+R q+d_i|)]
 <infinity.                                        (17)
```

Therefore the undamped exact-address coefficient

```text
C(r;X,m)
 =delta_hat(m)
  sum_((u+R beta).w=X)
    a(u)b(beta)
    conjugate(a(u+R beta+m e_c-r))                 (18)
```

is absolutely convergent. Dominated convergence gives

```text
lim_(rho->1-)C_rho(r;X,m)=C(r;X,m).                (19)
```

Thus an exact address has a canonical ordinary coefficient. Abel
regularization is needed only to sum over infinitely many different
addresses, not to sum one address orbit.

## 4. The Abel relation current is l1 and has the marked current as mass

For fixed `rho<1`, every atomic coefficient sequence is in `l1`. The map

```text
r -> v=u+R beta+m e_c-r
```

is a bijection from `Lambda` onto `{v:v.w=Y}` for fixed `(u,beta)`.
Consequently

```text
sum_(r in Lambda)|C_rho(r;X,m)|<infinity.           (20)
```

Absolute rearrangement is lawful, and every triple `(u,beta,v)` occurs
once. Hence

```text
sum_(r in Lambda)C_rho(r;X,m)
 =J_rho(X,m).                                      (21)
```

If the `I_i,J_i` are interval factors, their Poisson smoothings stay
between zero and one and converge in `L1`. Composition by a nonzero
integer speed preserves Haar `L1`, and a telescoping product estimate
handles all eighteen left factors and nine right factors. Therefore

```text
lim_(rho->1-)J_rho(X,m)
 =delta_hat(m)(E_Q)_hat(X)conjugate(E_hat(Y)).
                                                               (22)
```

For the THM-2327 triangle, the right side of (22) is nonzero.

Equations (18), (21), and (22) define the precise boundary:

```text
one address:
  ordinary absolutely convergent coefficient;

all addresses:
  canonical Abel pushforward whose total mass is the full current.    (23)
```

## 5. Every finite relation-residue quotient is a twist transform

Fix `N>=2` and put

```text
K_N={y in (Z/NZ)^d:y.w=0 mod N}.                    (24)
```

Primitivity supplies a Bezout vector `z.w=1`. If `y in K_N`, choose an
integer lift with `y.w=N h`; then

```text
r=y-N h z
```

lies in `Lambda` and reduces to `y`. The kernel of reduction is
`N Lambda`. Thus

```text
Lambda/N Lambda isomorphic to K_N,
|K_N|=N^(d-1).                                     (25)
```

Define the residue pushforward

```text
C_(rho,N)(y)
 =sum_(r in Lambda;r=y mod N) C_rho(r;X,m).         (26)
```

It is absolutely convergent for `rho<1`.

Write

```text
e_N(s)=exp(2*pi*i*s/N).
```

The characters of `K_N` are represented by

```text
ell in (Z/NZ)^d / <w>,

chi_ell(y)=e_N(ell.y).                              (27)
```

Define the twisted endpoint coefficients

```text
L_(rho,N)^ell(X)
 =sum_((u+R beta).w=X)
   a_rho(u)b_rho(beta)e_N(ell.(u+R beta)),

F_(rho,N)^ell(Y)
 =sum_(v.w=Y)a_rho(v)e_N(ell.v).                   (28)
```

The transform of the relation current factorizes exactly:

```text
H_(rho,N)(ell)
 :=sum_(r in Lambda)e_N(ell.r)C_rho(r;X,m)

 =delta_hat_rho(m)e_N(m ell_c)
   L_(rho,N)^ell(X)
   conjugate(F_(rho,N)^ell(Y)).                    (29)
```

Changing `ell` to `ell+s w` multiplies the three factors by

```text
e_N(sX), e_N(s m w_c), e_N(-sY),
```

whose product is one by (6). Thus (29) is representative-independent.
Finite Fourier inversion gives every residue coefficient:

```text
C_(rho,N)(y)
 =1/N^(d-1)
  sum_(ell mod <w>)
   e_N(-ell.y)H_(rho,N)(ell).                      (30)
```

This is the promised finite reformulation. It replaces an infinite sum
over exact relations by finitely many coordinate-twisted ordinary Fourier
coefficients.

### The coordinate translations

Let `I_(i,rho),J_(i,rho)` be the base Poisson smoothings. The two sums in
(28) are ordinary frequency coefficients of

```text
prod_i
 I_(i,rho)(w_i t+ell_i/N)
 J_(i,rho)(R w_i t+R ell_i/N),                     (31)

prod_i I_(i,rho)(w_i t+ell_i/N),                   (32)
```

respectively. The deep phase in (29) is the same translation on its base
factor. Every twisted current has a boundary limit by the same bounded
`L1` product argument as (22). Hence every residue-class Abel limit

```text
C_N(y;1-)=lim_(rho->1-)C_(rho,N)(y)                (33)
```

exists by the finite inverse transform. Equation (33) is an Abel aggregate
over an infinite residue class; it is not asserted to equal an absolutely
convergent sum of the individual values (18).

## 6. The exact 169-twist target theorem

Now specialize to a positive strict shallow-owner word stratum in
THM-2327. There are nine coordinates, and

```text
R=13^(lambda_j+1),
13|R,
13|X,Y,c_3.                                        (34)
```

Use THM-2309's decomposition

```text
K_13=L_13 direct_sum span(e_a,e_b),
G=K_13/L_13 isomorphic to F_13^2.                  (35)
```

For `q in G`, let

```text
A_rho(q)
 =sum_(r in Lambda;pi(r mod 13)=q) C_rho(r;X,m).   (36)
```

The dual group is

```text
G^=L_13^perp/<w>,
|G^|=169.                                          (37)
```

Equations (29)--(30) descend to

```text
A_rho(q)
 =1/169 sum_(ell in G^)
   e_13(-ell.q)H_(rho,13)(ell).                    (38)
```

Every boundary value

```text
A(q):=lim_(rho->1-)A_rho(q)                        (39)
```

exists. Since (38) is a finite partition of the full current,

```text
sum_(q in G)A(q)
 =delta_hat(m)(E_Q)_hat(X)conjugate(E_hat(Y))
 !=0.                                              (40)
```

At least one full-word target-fibre aggregate survives. Unlike the
zero-mode slice in THM-2333, this statement retains all nonconstant
harmonics of the actual positive terminal word. It still does not prove
that the surviving `q` is nonzero.

## 7. Parseval makes the nonzero-target obstruction exact

Write

```text
H(ell)=lim_(rho->1-)H_(rho,13)(ell),
H_bar=(1/169)sum_ell H(ell)=A(0).                  (41)
```

Finite Parseval applied to (38) gives

```text
sum_(q!=0)|A(q)|^2
 =1/169 sum_(ell in G^)|H(ell)-H_bar|^2.           (42)
```

Therefore the following are equivalent:

```text
some nonzero target vector q has A(q)!=0;

the 169 coordinate-twisted marked currents H(ell) are not constant;

the variance on the right side of (42) is positive.                 (43)
```

This is a necessary-and-sufficient noncancellation theorem for the next
target-vector step. It also gives the cheapest decisive computation:
evaluate the 169 shifted interval products in (31)--(32), not the enormous
bank of relation addresses, and test the exact consequence object (42).

The trivial twist `H(0)` is the nonzero THM-2327 current. If all twists
are equal to it, Fourier inversion puts the entire aggregate in `q=0`.
THM-2333's exact rational group-algebra hostile shows that this zero-only
boundary is possible for arbitrary full-support quotient weights.

## 8. Why the transported word is target-neutral but septimally active

In the target transform,

```text
e_13(ell.R beta)=1                                 (44)
```

for every word harmonic `beta`, because `13|R`. Equivalently, the word
translation in (31) is

```text
R ell_i/13=13^(lambda_j)ell_i in Z,
```

so it is trivial on the circle. Thus the complete terminal word remains
inside the left coefficient, but target characters do not act on it.
They act on the present owner factors, the bare endpoint, and the deep
leg.

Modulo seven the same dilation is a unit:

```text
R mod 7 in {1,-1}.                                  (45)
```

Hence the word harmonics are active in the septimal component of a
mod-`91` character twist. This is the exact CRT asymmetry behind the
post-THM-2331 obstruction.

To state it without hiding the all-unit requirement, let

```text
M_q(y)
 =[every coordinate of y is a unit modulo 91]
  [pi(y mod 13)=q],             y in K_91.          (46)
```

Define the all-unit target aggregate

```text
B_rho(q)
 =sum_(y in K_91)M_q(y)C_(rho,91)(y).              (47)
```

Its exact finite projector is

```text
B_rho(q)
 =1/91^8 sum_(ell mod <w>)
   M_q_hat(ell)H_(rho,91)(ell),

M_q_hat(ell)
 =sum_(y in K_91)M_q(y)e_91(-ell.y).               (48)
```

Every `B_rho(q)` has an Abel boundary limit. Under CRT the arithmetic mask
in (46) splits into a septimal all-unit mask and a mod-thirteen all-unit
target fibre. The transformed **current** in (48) does not split: its
mod-seven word charge is active while its mod-thirteen word charge is
trivial. Thus the `169`-twist variance (42) is the right test for an
unrestricted target aggregate, but not for an all-`91`-unit aggregate.

THM-2325 proves that every nonzero `M_q` has enormous support, and THM-2331
puts nonzero terms on that support. The still-missing statement is

```text
B(q)!=0 for some q!=0,                              (49)
```

or a stronger visible/phase-controlled version. Equation (48) is its exact
finite formulation.

## 9. The target action is not the root-character action

The character in (29) acts through **independent base-coordinate
translations**

```text
I_i(w_i t+ell_i/13).                               (50)
```

THM-2321's root character instead weights the thirteen physical
`T`-preimage branches of a Perron-normalized function. These are different
representations:

```text
relation-address action:
  translate the nine interval coordinates independently, modulo the
  owner-packet annihilator;

root-character action:
  take the cyclic Fourier transform of one physical predecessor-sheet
  index.                                           (51)
```

The common abstract `F_13^2` label and THM-2321's projective map do not
supply an intertwiner between (50) and the predecessor action. Equation
(44) sharpens the mismatch: target twists leave every transported word
harmonic neutral, whereas the root bispectrum is built from the word's
physical root fibres.

Consequently THM-2321's positive root-character rows cannot simply be
inserted into the right side of (42). A successful composition needs
either:

1. an explicit action intertwiner preserving the marked current;
2. a direct positive lower bound for the coordinate-twist variance; or
3. a target-specific phase/uniqueness sidecar that makes one orbit sum
   sign-coherent.

## 10. Exact positive factors can still cancel in one address orbit

Even positivity of every coordinate function does not make (18) positive.
Take

```text
d=2, w=(1,1), X=Y=m=0,
J_1=J_2=1,             delta_hat(0)=1.             (52)
```

Let the first endpoint factor have Fourier coefficients

```text
a_1(0)=1,
a_1(+/-1)=1/90,
a_1(+/-2)=-10/81,
a_1(+/-3)=10/81,                                  (53)
```

and no others. Its trigonometric polynomial is strictly positive, since

```text
a_1(0)-2 sum_(n=1)^3 |a_1(n)|
 =196/405>0.                                       (54)
```

Let the second factor be the positive Poisson kernel with parameter
`9/10`, so

```text
a_2(n)=(9/10)^|n|.                                 (55)
```

On the frequency-zero hyperplane put

```text
A_k=a_1(k)a_2(-k).
```

The nonzero values from `k=-3` through `3` are

```text
9/100,-1/10,1/100,1,1/100,-1/10,9/100.            (56)
```

They sum to one, so the full endpoint current is

```text
|sum_k A_k|^2=1.                                   (57)
```

For the exact relation address

```text
r=(1,-1),
```

equation (18) becomes

```text
C(r)=sum_k A_k A_(k-1)=0.                          (58)
```

All six terms in (58) are nonzero:

```text
-9/1000,-1/1000,1/100,1/100,-1/1000,-9/1000.
```

This is an exact factorized hostile with real strictly positive coordinate
functions and a nonzero full current. It is not claimed to be a canonical
LRC interval rectangle; it proves that factorization, positivity, and
termwise address occurrence alone do not prevent orbit cancellation.

## 11. Scope and exact remaining object

The current hierarchy is now:

```text
THM-2327:
  one nonzero marked word/deep/bare frequency triangle;

THM-2331:
  every prescribed target address occurs termwise;

THM-2334:
  exact addresses are absolutely convergent orbit coefficients;
  every finite residue pushforward is a coordinate-twist DFT;
  one full-word target aggregate survives;
  nonzero-target survival is exactly the variance (42);

still open:
  prove that variance positive, or prove the all-unit projector (49)
  nonzero with visible height and terminal phase.                   (59)
```

No scalar profile is excluded. The exact scalar ledger remains `165`,
the repeated-first and alternative resonance branches remain outside the
THM-2327 input, and LRC(14) remains open.

## 12. Exact companion

The companion checks primitive reduction
`Lambda/N Lambda=K_N` for every `2<=N<=13`, a finite-support 27-factor
pushforward, the exact character factorization, representative
independence, and residue inversion. It verifies the complete
`169`-character orthogonality kernel, the mod-thirteen neutrality and
mod-seven activity of the word dilation, and every rational value in the
positive-factor hostile (53)--(58).

Reproduce with

```bash
python3 04-computation/lrc14_relation_residue_pushforward_thm2334.py
python3 -O 04-computation/lrc14_relation_residue_pushforward_thm2334.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_relation_residue_pushforward_thm2334.out
```

byte-for-byte after LF normalization.
