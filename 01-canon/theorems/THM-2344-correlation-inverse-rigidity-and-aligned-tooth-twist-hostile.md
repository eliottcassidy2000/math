---
id: THM-2344
title: "Correlation-inverse rigidity and the aligned-tooth twist hostile"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. In the
  centered interval/complement specialization of THM-2334, the phase-free
  target response satisfies K(-ell)=conjugate(K(ell)); equivalently its
  inverse target transform is real. If THM-2343's zero-only boundary
  K=c chi_(-p) occurs, then c is nonzero real and the endpoint residue
  arrays are shifted convolution inverses. On the odd target group, this
  bad response is real exactly on the annihilator of p, so one real
  detecting twist, or global evenness/real-valuedness of K, excludes it.
  Reflection and physical-factor positivity alone do not: a centered
  danger/safe one-tooth model has positive scalar amplitude,
  K(s,t)=c zeta^(-t), and all 169 full twists H(s,t)=c. The hostile is a
  sharp local factorized control, not a canonical nine-coordinate
  terminal-word row. No scalar row, word-matching component, visible
  all-unit aggregate, terminal phase, or LRC(14) closure is proved.
source: codex-2026-07-25-correlation-inverse-rigidity
depends_on:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2343-deep-comb-affine-target-catalyst
related:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2327-two-colour-marked-unit-c3-triangle
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
  - THM-2340-owner-word-anova-target-landing
script: 04-computation/lrc14_correlation_inverse_hostile_thm2344.py
output: 05-knowledge/results/lrc14_correlation_inverse_hostile_thm2344.out
script_sha256: 41e8f674c4f8ccc059e3ead0673ff50423115a994f25d71b186792a2276a710b
output_sha256: 69c350342bac1814ac5622d443bc514571ab9371a8d7fced17362da6ed9b0ad7
hash_basis: working-tree bytes (LF)
---

# THM-2344 -- reflection narrows the hostile but does not kill it

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2343 identifies one sharp failure line for the full target current:

```text
K(ell)=c chi_ell(-p),               p!=0.           (1)
```

The canonical factors have more structure than an arbitrary complex
function on the `169` target characters: every danger interval, wider
guard, and complement is real and even.  That structure makes `K`
Hermitian.  It forces the scalar `c` in (1) to be real, and it gives a
cheap odd-order escape test.  It does **not** contradict (1).

The exact surviving hostile is an aligned Fourier tooth.  Its endpoint
translation contributes the inverse of the deepest-comb character, so the
two phases cancel and every full twist is equal.  Thus the next proof
needs a genuine non-alignment mechanism, not positivity or reflection
alone.

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

## 4. Exact aligned danger/safe hostile

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

## 5. Refined remaining problem

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
coordinates, the nonconstant transported word, owner asymmetry, or a
visible-address/terminal-phase sidecar.  A tournament is not intrinsic
here: the carrier is a real group-algebra correlation and its obstruction
is an affine shift.  Orienting endpoint residues would discard amplitudes,
convolution, and the distinguished translation `p`.

No scalar profile is excluded.  The ledger remains `165`; repeated-first
rows and alternative resonances remain outside THM-2327; and LRC(14)
remains open.

## 6. Exact companion

The companion verifies target-character orthogonality, the
correlation/product transform with exact rational point masses, the
singleton support shift, all `169` aligned-tooth phase cancellations,
Hermitian reflection, the thirteen-point real-value annihilator, the
odd-order evenness obstruction, and the danger/safe coefficient sign
ledger, including both the positive-array danger variant and the safe
variant.  Every load-bearing check raises explicitly under ordinary and
optimized Python.

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
