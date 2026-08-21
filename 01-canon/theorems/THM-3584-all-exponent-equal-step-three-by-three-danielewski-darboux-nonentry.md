---
id: THM-3584
title: "All-exponent equal-step three-by-three Danielewski Darboux nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  exponent N>=2 and every squarefree Sigma of degree at least two, a reduced
  polynomial Darboux pair on c^N e=Sigma(b) cannot have both weight supports
  equal to length-three arithmetic progressions with the same positive step.
  All five scalar-row alignments are closed.  For odd N, this includes both
  the gamma=2 middle-family resonance and the additional endpoint family.
  No unequal-step result, higher-support result, Keller-pair construction, or
  consequence for the planar Jacobian conjecture is claimed.
source: kps-s188 / delegated all-exponent Darboux attack, 2026-08-21
audit: >
  An independent reconstruction checked the bracket shift, arm-unit channel,
  both off-central operator factorizations, 17,737,925 admissible central
  profiles, the parity-sensitive gcd and square obstruction, and every odd
  endpoint subfamily.  Normal and optimized stdout agree with the stored
  output.  The audit also found the shorter regularity contradiction recorded
  after (47); the terminal gcd argument remains as an independent closure.
depends_on: []
related:
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
  - THM-3576-higher-exponent-belyi-keller-collision-tower
  - THM-3579-equal-step-three-by-three-danielewski-darboux-nonentry
script: 04-computation/jc_all_exponent_equal_step_three_by_three_danielewski_nonentry_thm3584.py
output: 05-knowledge/results/jc_all_exponent_equal_step_three_by_three_danielewski_nonentry_thm3584.out
script_sha256: 94dcd278a2acfaa945508996a2117a5f3be62cdff312f62956ad4bdd7d1a3949
output_sha256: 4a580c807b8c11b61fd46e7acc7fd1748a56c1de231b65e95c719882beceb639
hash_basis: LF-normalized bytes
---

# THM-3584 -- all-exponent equal-step three-by-three Danielewski Darboux nonentry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This is an
exponent-uniform replacement for the `N=2` argument of THM-3579.  The odd
exponents are not a formal substitution: they create two genuine central
resonances, both closed below.

All rings are over `C`.  Fix an integer `N>=2`, let `Sigma in C[b]` be
squarefree with `deg Sigma>=2`, and put

```text
A_(N,Sigma)=C[b,c,e]/(c^N e-Sigma(b)).                 (1)
```

Equip this smooth surface with the Poisson bracket

```text
{b,c}=c^N,        {c,e}=-Sigma'(b),
{b,e}=-N c^(N-1)e,                                      (2)
```

and the grading

```text
wt(b)=0,            wt(c)=1,            wt(e)=-N.       (3)
```

Suppose, after scalar constants have been subtracted, that `P,Q` each have
exactly three nonzero homogeneous pieces and that their supports are
arithmetic progressions with one common positive step `d`:

```text
P=c^r F(b,c^d),                  Q=c^s G(b,c^d),         (4)

F=f_0+f_1 t+f_2 t^2,            G=g_0+g_1 t+g_2 t^2,
t=c^d.                                                   (5)
```

Thus every retained weight-zero coefficient is nonconstant; if it were a
constant, the preliminary scalar subtraction would remove that support
point.  The theorem is

```text
{P,Q} != 1.                                              (6)
```

## 1. Regularity and the five-row compiler

On the chart `c!=0`, every homogeneous rational function of weight `u` is
uniquely `c^u f(b)`.  The exact extension criterion to `(1)` is

```text
c^u f(b) in A_(N,Sigma)
iff Sigma^ceil(-u/N) divides f             when u<0.     (7)
```

Indeed, a weight-`u` monomial containing `e^j` becomes `c^u Sigma^j` on the
chart, and the least allowable `j` is `ceil(-u/N)`.  This proves both
directions of `(7)`.

For homogeneous pieces define

```text
W_(u,v)(f,g)=v f'g-u f g'.                              (8)
```

A direct calculation from `(2)` gives

```text
{c^u f,c^v g}=c^(u+v+N-1) W_(u,v)(f,g).                (9)
```

Equivalently, with subscripts denoting partial derivatives,

```text
{P,Q}=c^(r+s+N-1) B_(r,s,d)(F,G),                     (10)

B_(r,s,d)(F,G)
 =s F_b G-r F G_b+d t(F_b G_t-F_t G_b).               (11)
```

The five coefficients of `B` have convolution multiplicities

```text
1,2,3,2,1.                                             (12)
```

If `(6)` failed, there would be a unique
`kappa in {0,1,2,3,4}` satisfying

```text
r+s+N-1+kappa d=0,             B=t^kappa.              (13)
```

The `kappa` row is the scalar one and the other four rows vanish.  The rest
of the proof excludes the five values of `kappa`.

## 2. Three local algebra lemmas

### 2.1 The simple-arm scalar gate

Let `beta` be a root of the squarefree polynomial `Sigma`, and put
`z=b-beta`.  Consider one scalar channel, so its weights satisfy

```text
u+v+N-1=0.                                              (14)
```

If the local coefficient orders are `m,n`, the first possible term of its
Wronskian has order and multiplier

```text
m+n-1,                       v m-u n.                  (15)
```

At least one of `u,v` is negative.  If both are negative, `(7)` gives
`m,n>=1`, so no scalar unit is possible.  Otherwise write the negative
weight as `-A`; the other weight is `A-(N-1)`.  A nonzero unit requires the
negative coefficient to have order one and the other to have order zero.
Regularity then gives

```text
N-1<=A<=N.                                              (16)
```

For `A=N-1` the multiplier in `(15)` is zero.  For `A=N` it is nonzero.
Consequently the only scalar channel surviving at a simple arm is

```text
(-N,1),                         (1,-N),                 (17)
```

with local coefficient orders `(1,0)` and `(0,1)`, respectively.

### 2.2 Extreme signs and common powers

Suppose an extreme row vanishes.  At an arm, a negative-weight coefficient
has positive order by `(7)`.  Pairing it with a positive-weight coefficient
gives a nonzero first local Wronskian term, because the two summands in its
multiplier have the same sign.  Pairing it with weight zero can vanish only
when the weight-zero coefficient is constant, which is excluded by the
reduction preceding `(4)`.  Thus a vanishing extreme cannot have mixed
strict signs.  In every interior scalar-row case the low output weight is
negative, so its two endpoint weights are necessarily both negative.  No
blanket nonnegativity claim for the high extreme is needed: when high signs
matter below, the scalar gate lists both possibilities and the mixed one is
excluded directly.

Write the low weights as `-R,-T`, with `R,T>0`, and their coefficients as
`f,g`.  The low equation is

```text
R f g'=T f'g.                                           (18)
```

Let `delta=gcd(R,T)`.  Taking valuations at every irreducible factor in
`C[b]`, or equivalently comparing logarithmic derivatives, gives

```text
f=A h^(R/delta),             g=B h^(T/delta),          (19)
```

where `A,B` are nonzero constants and `h` is a nonconstant polynomial.
Every arm is a root of `h`, because both low coefficients satisfy `(7)`.
The identical common-power statement applies at the high extreme.

### 2.3 Two degree-rigid operators

For polynomials `H,u`, define

```text
E_N(H,u)=H'u+N H u',             J_N(H,u)=H u'+N H'u.  (20)
```

If `H,u` are nonzero and `deg H+deg u>0`, then

```text
deg E_N(H,u)=deg J_N(H,u)=deg H+deg u-1,               (21)

lc E_N=(deg H+N deg u)lc(H)lc(u),
lc J_N=(deg u+N deg H)lc(H)lc(u).                      (22)
```

Both leading multipliers are positive integers.  In particular, neither
operator can equal one when one input is divisible by `Sigma` and hence has
degree at least two.

## 3. The one-channel rows `kappa=0,4`

For `kappa=0`, the scalar row is a single homogeneous bracket.  The gate
`(17)` reduces it, up to the bracket-preserving reflection
`(P,Q)->(Q,-P)`, to

```text
E_N(H,u)=1,                   Sigma divides H.          (23)
```

Equations `(21)--(22)` contradict `deg H>=2`.  Reflecting the supports gives
the same argument for `kappa=4`.  Hence

```text
kappa notin {0,4}.                                      (24)
```

## 4. The lower off-central row `kappa=1`

By the extreme-sign lemma write the low weights as `-R,-T`.  Equation `(13)`
gives

```text
R+T=d+N-1,                                              (25)
```

and the supports are

```text
supp(P)={-R,T-(N-1),R+2T-2(N-1)},
supp(Q)={-T,R-(N-1),2R+T-2(N-1)}.                      (26)
```

The scalar row has two channels.  By `(17)`, one can contribute a unit only
if `R=N` or `T=N`.  Apply `(P,Q)->(Q,-P)` if necessary and take `R=N`.  Formula
`(19)` and the required simple arm then force

```text
delta=N,             T=Nk,             f=A h,
g=B h^k,             Sigma divides h,              k>=1. (27)
```

Thus the normalized supports are

```text
supp(P)={-N,N(k-1)+1,N(2k-1)+2},
supp(Q)={-Nk,1,Nk+2}.                                  (28)
```

Let `a,q` denote the middle coefficients in `P,Q`.  The scalar row factors
exactly as

```text
E_N(h,Aq-kB h^(k-1)a)=1.                               (29)
```

Since `Sigma|h` and `deg Sigma>=2`, `(21)--(22)` rule out `(29)`.  Therefore
`kappa=1` is impossible.

## 5. The upper off-central row `kappa=3`

Let the high endpoint weights be `U,V`.  Equation `(13)` gives

```text
U+V=d-(N-1).                                             (30)
```

Applying gate `(17)` to the two scalar channels and then, if necessary,
applying `(P,Q)->(Q,-P)` gives exactly two high-endpoint possibilities:

```text
(U,V)=(d-N,1)            or            (d+1,-N).       (31)
```

The second pair has mixed strict signs, so its high extreme cannot vanish by
Section 2.2.  In the first pair, `d<N` is mixed-sign for the same reason.  If
`d=N`, then `U=0`; the high extreme is
`W_(0,1)(F,G)=F'G=0`, forcing the retained weight-zero coefficient `F` to be
constant, contrary to the reduction.  Therefore `d>N`.  Put `n=d-N>=1`.
The supports are

```text
supp(P)={-n-2N,-N,n},
supp(Q)={-2n-2N+1,-n-N+1,1}.                           (32)
```

The high common-power equation gives

```text
F=L K^n,                         G=M K.                 (33)
```

For the middle coefficients `a,q`, the scalar row is

```text
J_N(K,Ma-nL K^(n-1)q)=1.                               (34)
```

Both `a` and `q` have negative weight, so `Sigma` divides the second input
of `J_N`.  If that input is zero, the left side is zero; otherwise
`(21)--(22)` give positive degree.  Thus `kappa=3` is impossible.

## 6. Complete classification of the central row `kappa=2`

Again write the low weights as `-R,-T`.  Then

```text
R+T=2d+N-1.                                             (35)
```

Put

```text
p=(T-R-N+1)/2.                                          (36)
```

The three scalar channels have weight pairs

```text
(-R,R-N+1),       (p,-p-(N-1)),       (T-N+1,-T).      (37)
```

The gate says that the middle channel can survive only for `p=1` or
`p=-N`, equivalently

```text
|R-T|=N+1.                                              (38)
```

An endpoint channel can survive only for `R=N` or `T=N`.  Combining its
simple-arm requirement with `(19)` forces, up to `(P,Q)->(Q,-P)`,

```text
R=N,                    T=Nk.                           (39)
```

The integrality of `d=(Nk+1)/2` says exactly that `N` and `k` are odd.  The
middle and endpoint alternatives cannot overlap, since `N` cannot divide
`N+1`.  Therefore the central census is complete:

```text
(M)  |R-T|=N+1                         for every N>=2;

(E)  {R,T}={N,Nk},  N odd, k odd       only for odd N.  (40)
```

There are no omitted negative-weight or zero-weight scalar channels by
Section 2.1.

## 7. Central middle family `(M)`

Reflecting `(P,Q)` by `(P,Q)->(Q,-P)` if needed, take

```text
T=R+N+1,                         d=R+1.                (41)
```

The supports and coefficients are

```text
P=c^(-R)f+c a+c^(R+2)F,
Q=c^(-R-N-1)g+c^(-N)q+c^(R-N+1)G.                     (42)
```

The scalar channel `(a,q)` has weights `(1,-N)`.  Since neither endpoint
alternative in `(40)` overlaps this family, the scalar identity forces `a`
to be an arm-local unit (`a(beta)!=0` at each arm root `beta`) and `q` to be
simple at every arm.  No assertion that `a` is a constant unit of `C[b]` is
used.

### 7.1 The lower bridge fixes the arm order

Set

```text
delta=gcd(R,N+1),       alpha=R/delta,
beta=(R+N+1)/delta.                                      (43)
```

The low extreme gives `f=A h^alpha`, `g=B h^beta`.  The row below the scalar
row is

```text
R f q'-N f'q-(R+N+1)a'g-a g'=0.                        (44)
```

At an arm, let `ord(h)=ell`.  The first two terms in `(44)` have exact order
`alpha ell` and initial multiplier

```text
alpha(delta-N ell).                                     (45)
```

This multiplier cannot vanish: `delta|(N+1)`, hence `gcd(delta,N)=1`.  The
last two terms have exact order `beta ell-1`, because `a` is a unit.  Thus
cancellation requires

```text
((N+1)/delta)ell=1.                                     (46)
```

Consequently

```text
delta=N+1,       ell=1,       R=(N+1)m,
f=A h^m,         g=B h^(m+1),                         m>=1. (47)
```

In particular, the common arm polynomial `h` is simple at every root of
`Sigma`.

There is already a short contradiction here.  The coefficient `f` has
weight `-R=-(N+1)m`, so regularity `(7)` requires

```text
ord_beta(f) >= ceil((N+1)m/N) >= m+1.                 (47a)
```

at every simple root `beta` of `Sigma`.  But `(47)` and the simple-arm
conclusion give `ord_beta(f)=m`.  Thus family `(M)` is empty.  The next two
subsections retain the independent first-integral/terminal closure because
it isolates the odd-exponent resonance that was absent at `N=2`.

### 7.2 The `(N+1)`-power first integral

After substituting `(47)`, equation `(44)` becomes

```text
mA((N+1)h q'-N h'q)
 =(m+1)B h((N+1)h a'+h'a).                             (48)
```

Let `lambda=(m+1)B/(mA)` and set `w=q-lambda h a`.  Then

```text
(N+1)h w'-N h'w=0,

(w^(N+1)/h^N)'
 =w^N h^(-N-1)((N+1)h w'-N h'w)=0.                    (49)
```

If the first integral were a nonzero constant, a simple arm would give
`(N+1)ord(w)=N`, impossible.  Hence

```text
q=lambda h a.                                           (50)
```

### 7.3 The parity-sensitive terminal obstruction

The high weights are

```text
U=R+2,                  V=R-N+1.                        (51)
```

Put

```text
gamma=gcd(U,V)=gcd(2,N+1),       D=(N+1)/gamma,
u=U/gamma,                       v=V/gamma.             (52)
```

Thus `gamma=1` for even `N`, while

```text
gamma=2                         for odd N.               (53)
```

The high extreme gives `F=L K^u`, `G=M K^v`.  Substitution of `(50)` into
the row above the scalar row, with

```text
A_0=Mv-lambda L u h K^D,          Y=a A_0,              (54)
```

reduces it exactly to

```text
gamma K Y'-K'Y=0,

(Y^gamma/K)'=Y^(gamma-1)K^(-2)(gamma K Y'-K'Y)=0.       (55)
```

Therefore

```text
(a A_0)^gamma=C K.                                      (56)
```

The polynomial `A_0` is nonconstant and is coprime to `K`, since
`A_0 mod K=Mv!=0`.  If `K` is constant, the two sides of `(56)` already have
different degrees.  If `K` is nonconstant, `C=0` would force `aA_0=0`, while
`C!=0` would force every factor of `A_0` to divide `K`.  Both contradict the
preceding facts.  This closes family `(M)` for all `N`.

For odd `N`, `(56)` is the load-bearing square equation

```text
(a A_0)^2=C K.                                           (57)
```

It is not covered by the `N=2` cube argument and is one of the two odd-
exponent central families required by the theorem.

## 8. The odd-exponent endpoint family `(E)`

Now let `N,k` both be odd and orient `(39)` as `R=N,T=Nk`.  Define

```text
d=(Nk+1)/2,             S=(Nk-1)/2,
p=(N(k-2)+1)/2,         U=N(k-1)+1.                    (58)
```

The supports are

```text
supp(P)={-N,p,U},                 supp(Q)={-Nk,-S,1}.   (59)
```

The low common powers and the endpoint scalar gate give

```text
f=A h,                    g=B h^k,                      (60)
```

with `h` simple at every arm.  The row below the scalar row factors as

```text
N h X'-S h'X=0,
X=Aq-kB h^(k-1)a.                                      (61)
```

Its first integral is

```text
(X^N/h^S)'=X^(N-1)h^(-S-1)(N hX'-S h'X)=0.             (62)
```

Because `2S=Nk-1`, one has `N` not dividing `S`.  A simple arm therefore
rules out a nonzero constant in `(62)`, and

```text
q=(kB/A)h^(k-1)a.                                      (63)
```

At the high extreme write `F=L K^U`, `G=M K`.  The row above the scalar row
is the exact derivative equation

```text
M(a/K^p)'-U L(q K^S)'=0.                               (64)
```

### 8.1 The nonminimal odd endpoint, `k>=3`

Here `p>0`.  Integrating `(64)` and using `(63)` gives

```text
a D_0=C K^p,

D_0=M-U L(kB/A)h^(k-1)K^(N(k-1)).                      (65)
```

The polynomial `D_0` is nonconstant and `D_0 mod K=M!=0`.  If `K` is
constant, degrees contradict `(65)`.  Otherwise, `C=0` is impossible in the
domain and `C!=0` contradicts `gcd(D_0,K)=1`.  This closes every odd
`k>=3` endpoint resonance.

### 8.2 The minimal odd endpoint, `k=1`

Now both supports equal

```text
{-N,-(N-1)/2,1}.                                       (66)
```

Equation `(63)` says `Aq=Ba`.  The middle scalar channel then vanishes, and
the full scalar row collapses to

```text
(AM-BL)(K h'+N h K')=1.                                (67)
```

If the constant prefactor is zero, `(67)` is immediate nonsense.  If it is
nonzero, the differential polynomial has degree `deg h+deg K-1>=1` and
leading multiplier `deg h+N deg K>0`, because `Sigma|h` and
`deg Sigma>=2`.  Thus `(67)` is also impossible.

This closes family `(E)`, the second odd-exponent central family, and
completes the proof of `(6)`.  **QED.**

## 9. Sharp `N=3` hostiles

The smallest odd exponent already contains every new phenomenon.  The exact
minimal support profiles are

```text
middle:       (-4,1,6)       | (-8,-3,2),
endpoint k=1: (-3,-1,1)      | (-3,-1,1),
endpoint k=3: (-3,2,7)       | (-9,-4,1).              (68)
```

They reduce respectively to the following three terminal gates:

```text
(a[M-(6BL/A)hK^2])^2=C K,                              (69)

(AM-BL)(K h'+3hK')=1,                                 (70)

a[M-(21BL/A)h^2K^6]=C K^2.                            (71)
```

These are hostile profiles, not Darboux pairs.  Equations `(69)--(71)` show
why both odd families and the `gamma=2` parity split must appear explicitly.

Five boundary controls are sharp.

1. A both-negative scalar channel such as `(-1,-1)` has arm order at least
   one and cannot make a unit.
2. The zero boundary `(-2,0)` has multiplier zero even at orders `(1,0)`.
3. The formal middle equation admits `h=z^4,w=z^3`, since
   `4hw'-3h'w=0`; the lower bridge's forced simple arm is indispensable.
4. The endpoint equation admits `h=z^3,X=z`, since `3hX'-h'X=0`; the
   endpoint scalar gate's simple arm is indispensable.
5. The hypotheses `deg Sigma>=2` and polynomial regularity cannot be
   discarded: `h=b,K=1` makes `Kh'+3hK'=1`, while on the Laurent chart

   ```text
   {b,-1/(2c^2)}=1.                                    (72)
   ```

More generally, the last hostile is

```text
{b,-c^(1-N)/(N-1)}=1                                   (73)
```

for every `N>=2`; its pole is exactly what prevents polynomial entry.

## 10. Scope and verification contract

The proved statement is precisely

```text
PROVED:
  for every N>=2 and every squarefree Sigma with deg Sigma>=2,
  every reduced equal-step 3 x 3 weight-support cell is empty;

OPEN / NOT CLAIMED:
  unequal-step 3 x 3 cells;
  2 x 4, 4 x 2, and larger support cells;
  the proposed THM-3576 higher-exponent collision tower;
  existence or nonexistence of arbitrary Darboux pairs on A_(N,Sigma);
  any conclusion for JC(2).                             (74)
```

Reproduce the exact companion with

```bash
python3 04-computation/jc_all_exponent_equal_step_three_by_three_danielewski_nonentry_thm3584.py
python3 -O 04-computation/jc_all_exponent_equal_step_three_by_three_danielewski_nonentry_thm3584.py
```

It uses only the Python standard library and no optimization-sensitive
`assert`.  It checks the Laurent compiler and all five rows; regularity and
the simple-arm gate for `2<=N<=32`; complete central resonance censuses on a
large exact rectangle; both parity branches of the middle family; every
odd-endpoint identity; exact polynomial gcd controls; and all `N=3`
hostiles above.  Normal and `-O` output are byte-identical.

The finite ranges are exact identity and hostile controls.  The theorem for
all `N`, all support parameters, and all degrees is supplied by the
valuation, UFD, first-integral, coprimality, and leading-degree arguments in
Sections 2--8; it is not an extrapolation from the finite scan.
