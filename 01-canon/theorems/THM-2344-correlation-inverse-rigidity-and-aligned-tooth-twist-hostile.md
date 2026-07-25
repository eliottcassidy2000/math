---
id: THM-2344
title: "Correlation-inverse rigidity and the aligned-axis transported-word hostile"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. In the
  centered interval/complement specialization of THM-2334, the phase-free
  target response satisfies K(-ell)=conjugate(K(ell)); equivalently its
  inverse target transform is real. If THM-2343's zero-only boundary
  K=c chi_(-p) occurs, then c is nonzero real and the endpoint residue
  arrays are shifted convolution inverses. On the odd target group, this
  bad response is real exactly on the annihilator of p, so one real
  detecting twist, or global evenness/real-valuedness of K, excludes it.
  More generally, central phase-conjugate covariance about `h` can meet
  the constant boundary only at `h=p`; translated endpoints reduce the
  aligned case to perfect periodic autocorrelation. Reflection and
  physical-factor positivity alone do not. More sharply,
  with only one active translated coordinate and 13 dividing the word
  dilation, an arbitrary nonconstant transported word changes the
  endpoint amplitude but cannot change its target character. A centered
  safe factor times one transported danger factor gives a positive
  12-tooth example with K(s,t)=c zeta^(-t) and all 169 full twists
  H(s,t)=c. Thus at least two independently translated base coordinates
  or genuine cross-coordinate mixing are necessary. The hostile is a
  sharp local factorized control, not a canonical nine-coordinate
  terminal-word row. On a typed nine-factor control, 334,825 exact
  two-sided lifts split as 167,404 positive and 167,421 negative terms,
  so the unweighted sign law is not coherent. No scalar row, grouped-sum
  cancellation, word-matching component, visible all-unit aggregate,
  terminal phase, or LRC(14) closure is proved.
source: codex-2026-07-25-correlation-inverse-rigidity
depends_on:
  - THM-2331-two-sided-septimal-address-embedding-in-marked-current
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2343-deep-comb-affine-target-catalyst
related:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2327-two-colour-marked-unit-c3-triangle
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
  - THM-2340-owner-word-anova-target-landing
script: 04-computation/lrc14_correlation_inverse_hostile_thm2344.py
output: 05-knowledge/results/lrc14_correlation_inverse_hostile_thm2344.out
script_sha256: ce2513ecd6e8290677d69bb5a00b8cb64216035331ebd10887353a34463271fd
output_sha256: c5e17223240b8a383c8a17a6869bb500e363e4010c2a286f62fb4310e11a4bed
secondary_script: 04-computation/lrc14_two_sided_lift_sign_bank_thm2344.py
secondary_output: 05-knowledge/results/lrc14_two_sided_lift_sign_bank_thm2344.out
secondary_script_sha256: 2f37db2a2a9150a83e31bc75cd090f79c8151eb064b35137b5912aab6ac63f23
secondary_output_sha256: ace877a717d5cd2746e9ba88f71bd40c70392d90b80579a77e3504d7b0435d2e
hash_basis: working-tree bytes (LF)
---

# THM-2344 -- reflection narrows the hostile but does not kill it

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The independent audit rederived the quotient-safe reflection law, both
finite Fourier signs in the endpoint cross-correlation, the shifted-unit
equivalence, and the odd-order real-value separator.  It checked the
danger and safe aligned-tooth models directly, including the positive
amplitude and all `169` cancelling full twists.  It separately rederived
the arbitrary same-axis transported-word congruence, checked the
nonconstant twelve-tooth geometry including null boundaries, and verified
its positive amplitude.  Ordinary and optimized replays match the stored
transcript byte-for-byte after LF normalization, and both recorded hashes
match.

THM-2343 identifies one sharp failure line for the full target current:

```text
K(ell)=c chi_ell(-p),               p!=0.           (1)
```

The canonical factors have more structure than an arbitrary complex
function on the `169` target characters: every danger interval, wider
guard, and complement is real and even.  That structure makes `K`
Hermitian.  It forces the scalar `c` in (1) to be real, and it gives a
cheap odd-order escape test.  It does **not** contradict (1).

The exact surviving hostile is an aligned Fourier residue class.  Its
endpoint translation contributes the inverse of the deepest-comb
character, so the two phases cancel and every full twist is equal.  Even
an arbitrary nonconstant transported word on the same axis cannot alter
that residue character.  Thus the next proof needs a genuine
multi-coordinate non-alignment mechanism, not positivity, reflection, or
word nonconstancy alone.

## 1. Reflection of the canonical target twists

Retain THM-2334/2343's notation on a positive strict shallow-owner word
stratum:

```text
G isomorphic to F_13^2,
H(ell)=chi_ell(p)K(ell),
p=(0,m)!=0,                                        (2)

K(ell)
 =d_hat(m)L^ell(X)conjugate(F^ell(Y)).
```

Every base factor `I_i,J_i` of the canonical source rectangle and
terminal word is a centered interval indicator or its complement.  Hence
it is real and even on the circle.  Its Poisson smoothing has the same
two properties.  Since

```text
R=13^(lambda_j+1),
```

the transported-word shift `R ell_i/13` is integral and therefore
invisible on the circle.  At Abel radius `0<rho<1`, write the two physical
products behind the twisted coefficients as

```text
l_(rho,ell)(t)
 =prod_i I_(i,rho)(w_i t+ell_i/13)
          J_(i,rho)(R w_i t),

f_(rho,ell)(t)
 =prod_i I_(i,rho)(w_i t+ell_i/13).                (3)
```

Evenness gives the pointwise identities

```text
l_(rho,-ell)(t)=l_(rho,ell)(-t),
f_(rho,-ell)(t)=f_(rho,ell)(-t).                   (4)
```

The functions in (3) are real.  With the canon's Fourier convention,
(4) implies

```text
L_(rho)^(-ell)(X)=conjugate(L_(rho)^ell(X)),
F_(rho)^(-ell)(Y)=conjugate(F_(rho)^ell(Y)).        (5)
```

The centered danger coefficient `d_hat_rho(m)` is real.  Therefore

```text
K_rho(-ell)=conjugate(K_rho(ell)),
H_rho(-ell)=conjugate(H_rho(ell)).                 (6)
```

The ordinary `L^1` Abel boundary theorem in THM-2334 lets `rho` tend to
one in the finite family (6):

```text
K(-ell)=conjugate(K(ell)),
H(-ell)=conjugate(H(ell)).                         (7)
```

This statement is quotient-safe.  Negation descends to `G^`, while
THM-2334's representative phases cancel in `H`; division by the
well-defined character `chi_ell(p)` makes `K` well-defined as well.

## 2. The hidden object is a real cross-correlation

Define the normalized inverse target transforms of the endpoint twists:

```text
U(x)=1/169 sum_ell conjugate(chi_ell(x))L^ell(X),
V(y)=1/169 sum_ell conjugate(chi_ell(y))F^ell(Y).
                                                               (8)
```

These endpoint functions descend to `G^`: in this specialization
`13|X,Y`, so changing a representative by a multiple of `w` contributes
the trivial phases `e_13(sX)` and `e_13(sY)`.

Fourier inversion gives

```text
L^ell(X)=sum_x chi_ell(x)U(x),
F^ell(Y)=sum_y chi_ell(y)V(y).                     (9)
```

Equation (5), followed by the Abel limit, says exactly that `U` and `V`
are real-valued arrays.  Multiplying (9) and inverting on `G` yields

```text
A_K(r)
 =d_hat(m) sum_(x-y=r) U(x)V(y).                  (10)
```

Thus the phase-free target coefficient is a real endpoint
cross-correlation, not an arbitrary complex amplitude.  If

```text
V^vee(z)=V(-z),
```

then (10) is the group-algebra identity

```text
A_K=d_hat(m)(U*V^vee).                             (11)
```

Combining (1), finite inversion, and (11) gives the exact
correlation-inverse boundary:

```text
K=c chi_(-p)

iff A_K=c delta_(-p)

iff U*V^vee=(c/d_hat(m))delta_(-p).                (12)
```

Because `c!=0` and `d_hat(m)!=0`, (12) makes both endpoint arrays units
of the real group algebra.  Equivalently,

```text
L^ell(X)!=0 and F^ell(Y)!=0 for every ell,          (13)
```

and one endpoint array is a shifted scalar convolution inverse of the
other.  A zero at even one endpoint twist therefore certifies escape from
the bad line.

If `U,V` happen to be nonnegative, singleton support in (12) forces both
supports to be singletons: every difference in
`supp(U)-supp(V)` must equal `-p`.  This rigidity still does not exclude
the boundary; two aligned positive point masses realize it exactly.
Canonical physical factors being nonnegative does not itself imply that
their Fourier-residue arrays are nonnegative.

## 3. What odd-order reflection really buys

Evaluate (7) at the trivial character.  Since THM-2343 gives `K(0)!=0`,

```text
K(0) is nonzero real.                              (14)
```

If the bad line (1) holds, then `c=K(0)`, so

```text
c belongs to R minus {0}.                          (15)
```

This is the full unconditional consequence of reflection.  Indeed, for
real `c`,

```text
c chi_(-p)(-ell)
 =conjugate(c chi_(-p)(ell)),                      (16)
```

so every bad inverse character already has the required Hermitian
symmetry.

There is nevertheless a sharp conditional separator.  On an odd-order
group, a character value is real only when it is `1`.  Hence under (1)
and (15),

```text
K(ell) is real
 iff chi_ell(p)=1
 iff ell belongs to p^perp.                        (17)
```

For `G=F_13^2` and `p=(0,m)`, the annihilator is the single row

```text
p^perp={(s,0):s in F_13}.                          (18)
```

Consequently any one of the following excludes the zero-only boundary:

```text
K(ell) is real for one ell not in p^perp;

K(-ell)=K(ell) for one ell not in p^perp;

K is real-valued on G^;

K is even on G^.                                   (19)
```

The middle equivalences use (7).  Globally, an even bad character would
give `chi_ell(2p)=1` for every `ell`, hence `2p=0`; odd order and `p!=0`
make that impossible.

This identifies a much cheaper possible analytic target than generic
variance: prove one detecting twist is reflection-fixed or real.  No
current theorem supplies that target.  Physical positivity fixes the
sign of selected untwisted factors, not the phase of coordinate-shifted
Fourier coefficients.

## 4. Odd-centre covariance and the perfect-autocorrelation boundary

The preceding Hermitian law is centred at the origin of character space.
A more flexible conditional separator lives directly in correlation space.
For arbitrary complex arrays `U,V` on a finite abelian group `G`, put

```text
R_(U,V)(s)=sum_x U(x)conjugate(V(x+s)),

H(chi)=chi(p)U_hat(chi)conjugate(V_hat(chi)).       (20a)
```

The sign convention of Section 2 gives

```text
A_H(q)=R_(U,V)(p-q).                                (20b)
```

Assume that `G` has odd order and that for some `h in G` and `|eta|=1`,

```text
R_(U,V)(2h-s)=eta conjugate(R_(U,V)(s))
                                                   for every s. (20c)
```

If `H` is a nonzero constant, (20b) says
`R_(U,V)=c delta_p`. Its singleton support must be fixed by
`s -> 2h-s`, so `2h-p=p`. Multiplication by two is injective on an
odd group; therefore

```text
h=p.                                                (20d)
```

Thus covariance about any proved centre `h!=p` certifies nonconstancy.
Equivalently, it is enough to prove that for some `theta`,

```text
exp(-i theta) chi(h)
 U_hat(chi)conjugate(V_hat(chi))
```

is real for every character and that `h!=p`: at a nonzero constant
boundary the displayed value would be a nonzero real scalar times
`chi(h-p)`, while character separation on an odd group supplies a
nonreal value.

The aligned boundary is exact. Define

```text
(T_h U)(x)=U(x-h).
```

If `V=alpha T_h U`, with `alpha!=0` and `U` nonzero, then

```text
H(chi)
 =conjugate(alpha)chi(p-h)|U_hat(chi)|^2.           (20e)
```

For `h!=p` this is nonconstant. For `h=p` it is constant exactly when
`|U_hat(chi)|` is constant, equivalently when the periodic
autocorrelation

```text
R_(U,U)(s)
```

vanishes for every `s!=0`. The residual is therefore a perfect periodic
autocorrelation array, not unspecified cancellation. The point-mass
hostile in the next section attains this boundary.

## 5. Exact aligned danger/safe hostile

Let

```text
d(t)=1_(||t||<1/14),
g(t)=1-d(t),
zeta=exp(2*pi*i/13).                               (20)
```

For nonzero integer `n`,

```text
d_hat(n)=sin(pi*n/7)/(pi*n),
g_hat(n)=-d_hat(n).                                (21)
```

In particular

```text
d_hat(1)>0,
d_hat(2)>0,
g_hat(1)<0,
g_hat(2)<0.                                        (22)
```

Take the primitive three-coordinate carrier

```text
w=(1,13,26),                 c_3=26.                (22a)
```

Modulo thirteen its relation kernel is the two-dimensional target plane
on the last two coordinates. Let `s` act on the inert constant target
factor and let `t` translate the active safe factor. Take

```text
m=1,
X=c_3,
Y=2c_3=X+m c_3.                                   (23)
```

Direct Fourier expansion of `g(c_3 tau+t/13)` gives

```text
L^(s,t)(X)=g_hat(1)zeta^t,
F^(s,t)(Y)=g_hat(2)zeta^(2t).                      (24)
```

The deepest danger leg has base coefficient `d_hat(1)` and character
`zeta^t`.  Therefore the phase-free and full responses are

```text
K(s,t)
 =d_hat(1)g_hat(1)g_hat(2)zeta^(-t)
 =a zeta^(-t),

H(s,t)=zeta^t K(s,t)=a,                            (25)
```

where

```text
a=d_hat(1)^2 d_hat(2)>0.                           (26)
```

All `169` full twists are the same positive real number.  The inverse
phase-free transform is `a delta_(0,-1)` and the full target transform is
`a delta_0`.  At the correlation level, the endpoint arrays are the two
same-axis point masses

```text
U=g_hat(1)delta_(0,1),
V=g_hat(2)delta_(0,2),                             (27)
```

whose difference is `(0,-1)=-p`.

Replacing the endpoint safe factor `g` by the danger factor `d` changes
both endpoint signs and leaves (25)--(26) unchanged.  In that variant

```text
U=d_hat(1)delta_(0,1),
V=d_hat(2)delta_(0,2)
```

are themselves nonnegative.  Thus even nonnegativity of both residue
arrays permits the singleton correlation; it only rigidifies it to this
aligned point-mass form.

This hostile satisfies every unconditional conclusion above:

```text
K(-s,-t)=conjugate(K(s,t)),
c=a>0,
K is real exactly when t=0,
U and V are real one-point convolution inverses up to shift.          (28)
```

It uses nonnegative centered physical functions `d`, `g`, and the inert
constant.  It is a genuine one-tooth instance of the factorized
coordinate-translation formula, with no Fourier truncation.  It is
**not** a canonical nine-coordinate shallow-owner terminal-word row:
the other factors are inert and the transported word is constant.
Accordingly it proves the logical insufficiency of reflection,
centeredness, factor positivity, an untwisted positive scalar, and
one-tooth support.  It does not prove that a canonical row is bad.

## 6. Arbitrary same-axis transported words preserve the bad character

The constant transported word in the first hostile is not the essential
reason for failure.  Normalize the active deepest coordinate by

```text
x=c_3 tau.
```

Let `I,W` be arbitrary bounded-variation functions, let `13|R`, and
consider the one-axis left and bare endpoint factors

```text
l_t(x)=I(x+t/13)W(Rx),
f_t(x)=I(x+t/13).                                  (29a)
```

For integers `a,b`, Abel expansion followed by the ordinary boundary
limit gives

```text
L_t(a)
 =sum_beta I_hat(a-R beta)W_hat(beta)
             zeta^((a-R beta)t)
 =zeta^(a t)L_0(a),                                (29b)

F_t(b)=zeta^(b t)I_hat(b).                         (29c)
```

The second equality in (29b) is the entire mechanism: every word
harmonic changes the base harmonic by a multiple of `R`, hence by zero
modulo `13`.  The word may be nonconstant and may have arbitrarily many
harmonics; it changes `L_0(a)` but not the target character.

If the physical triangle is aligned with this coordinate,

```text
b=a+m,
```

then, whenever the displayed endpoint coefficients are nonzero,

```text
K(s,t)
 =d_hat(m)L_0(a)conjugate(I_hat(b))zeta^(-m t),

H(s,t)=zeta^(m t)K(s,t)
      =d_hat(m)L_0(a)conjugate(I_hat(b)).           (29d)
```

Thus **every** full target twist is constant for an arbitrary same-axis
transported word.  Nonconstancy of the word alone cannot break the
correlation-inverse boundary.

There is an exact centered interval instance with a genuinely
nonconstant word.  In (29a), take

```text
I=g,       W=d,       R=13,       a=1, b=2, m=1.   (29e)
```

Modulo null endpoints,

```text
g(x)d(13x)
```

is exactly the union of the twelve noncentral `13`-adic teeth.  Each has
length `1/91` and centre `k/13`, `1<=k<=12`.  Therefore

```text
C:=L_0(1)
  =sin(pi/91)/pi sum_(k=1)^12 zeta^(-k)
  =-sin(pi/91)/pi <0.                              (29f)
```

Since `g_hat(2)=-d_hat(2)<0` and `d_hat(1)>0`, (29d) becomes

```text
K(s,t)=d_hat(1)C g_hat(2)zeta^(-t),

H(s,t)=d_hat(1)C g_hat(2)>0                        (29g)
```

for all `169` target twists.  This example contains an actual
nonconstant transported danger factor and still lies on the bad line.
It remains a one-active-coordinate local model rather than a canonical
nine-coordinate row.  Its precise consequence is a necessity statement:
the next escape must use at least two independently translated base
coordinates, or a cross-coordinate coupling that prevents all endpoint
harmonics from sharing one target residue.

## 7. A typed nine-factor sign-incoherence control

The one-tooth hostile deliberately makes the other coordinates inert.
The genuine nine-factor sign law still does not supply a coherent
termwise sidecar. Use THM-2331's typed control in the labelled order

```text
w=(13,13^3,2*13^5,1,14,27,40,53,66),

r=(-27,-27,-27,20110798,-41,-27,-27,-27,38),

X=13,                 m=1,                 i_3=2.  (28a)
```

Then `r.w=0` and every coordinate of `r` is a unit modulo `91`. Put

```text
d=m e_(i_3)-r
```

and exhaust

```text
U_7={
 u_bar in F_7^9:
 u_bar.w=X,
 u_bar_i!=0,
 u_bar_i+d_i!=0 for every i
}.                                                  (28b)
```

For each residue vector, take centred representatives and correct the
speed-one guard coordinate by the unique multiple of seven needed to
make the left frequency exactly `X`. This is THM-2331's exact Bezout
lift, and the right endpoint is `v=u+d`.

For the `c_1`-owner rectangle, seven ordinary complements and the guard
complement occur at **each** endpoint. Their sixteen fixed minus signs
cancel. The remaining signs are the exact centred danger-arc signs on
the non-guard right endpoint, the two wider guard-arc signs, and the
positive deepest-comb sign. Independent dynamic programming and literal
enumeration give

```text
|U_7|=334,825,

positive terms=167,404,
negative terms=167,421.                             (28c)
```

Both signs therefore occur in nearly equal numbers in one exact typed
address bank. This proves no cancellation: magnitudes differ, the whole
address orbit is infinite, and an unweighted sign count is not a grouped
coefficient. The control word has the correct strict valuation type but
is not asserted to satisfy the scalar cover. Its consequence is only
that septimal support and interval signs do not furnish a sign-coherent
fibre sum.

## 8. Refined remaining problem

The boundary now has three equivalent descriptions:

```text
twist space:
  K is the real inverse character c chi_(-p);

target space:
  A_K is the real singleton c delta_(-p);

endpoint group algebra:
  U and V^vee are shifted scalar convolution inverses.               (29)
```

Each representation suggests a different decisive attack:

```text
reflection attack:
  find one real or reflection-fixed twist off p^perp;

support attack:
  prove A_K has support away from -p, preferably in THM-2343's
  word-matching affine mask;

factor attack:
  prove the two canonical endpoint residue arrays cannot be shifted
  convolution inverses;

phase attack:
  retain one all-91-unit atomic address or terminal component whose
  phase cannot participate in the aligned inverse character.         (30)
```

The hostile says the missing input must use interaction among multiple
independently translated coordinates, word/base coupling across distinct
axes, owner asymmetry, or a visible-address/terminal-phase sidecar.  A
nonconstant same-axis word is explicitly insufficient.  A tournament is
not intrinsic here: the carrier is a real group-algebra correlation and
its obstruction is an affine shift.  Orienting endpoint residues would
discard amplitudes, convolution, and the distinguished translation `p`.

No scalar profile is excluded.  The ledger remains `165`; repeated-first
rows and alternative resonances remain outside THM-2327; and LRC(14)
remains open.

## 9. Exact companions

The companion verifies target-character orthogonality, the
correlation/product transform with exact rational point masses, the
singleton support shift, all `169` aligned-tooth phase cancellations,
Hermitian reflection, the thirteen-point real-value annihilator, the
odd-order evenness obstruction, and the danger/safe coefficient sign
ledger, including both the positive-array danger variant and the safe
variant.  It also checks the general same-axis residue congruence over
`32,760` hostile parameter choices, the exact twelve-tooth geometry, the
cyclotomic identity `sum_(k=1)^12 zeta^(-k)=-1`, the positive
nonconstant-word amplitude, and all `169` resulting phase
cancellations.  Every load-bearing check raises explicitly under
ordinary and optimized Python.

Reproduce with

```bash
python3 04-computation/lrc14_correlation_inverse_hostile_thm2344.py
python3 -O 04-computation/lrc14_correlation_inverse_hostile_thm2344.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_correlation_inverse_hostile_thm2344.out
```

byte-for-byte after LF normalization.

The separately scoped nine-factor sign bank is reproduced by

```bash
python3 04-computation/lrc14_two_sided_lift_sign_bank_thm2344.py
python3 -O 04-computation/lrc14_two_sided_lift_sign_bank_thm2344.py
```

Both runs must equal

```text
05-knowledge/results/lrc14_two_sided_lift_sign_bank_thm2344.out
```

It checks the exact relation and all-unit address, counts the finite-field
universe by a path independent of literal enumeration, and verifies every
term sign with ordinary and optimized Python.

## 10. Independent audit

The independent audit rederived the Hermitian law including transported-word
neutrality, quotient descent, every DFT/correlation sign, the shifted inverse
criterion, the odd annihilator and covariance tests, and the translated
perfect-autocorrelation boundary. It found and repaired one typing issue in
the hostile: the deepest tooth must lie in the mod-thirteen kernel, now
witnessed by `w=(1,13,26)`. The auditor independently checked the aligned
danger/safe amplitudes and the `334,825=167,404+167,421` sign bank. Both
companions, hashes, dependency slugs, documentation checks, and ordinary and
optimized transcripts passed.
