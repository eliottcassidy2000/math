# Planar JC: the tournament descendant is a Plucker tie-depth response, not an orientation

**Session status (2026-08-23).**  The five-vertex skew-gain identities and the
equality-seam three-case law below are **PROVED algebraically,
VERIFIED-EXACT, and INDEPENDENTLY AUDITED** as necessary statements.  The
session promoted the iterated parity passport
[THM-3894](../01-canon/theorems/THM-3894-cusp-residual-all-degree-gauge-kummer-parity-filtration.md),
closed the first open even terminal at `n=4` in
[THM-3896](../01-canon/theorems/THM-3896-cusp-residual-degree-four-equality-seam-emptiness.md),
and converted the finite `f=0` work into the all-degree closure
[THM-3897](../01-canon/theorems/THM-3897-f-zero-residual-all-degree-global-emptiness.md)
through the high-`y` cutoff
[THM-3895](../01-canon/theorems/THM-3895-f-zero-quartic-covariant-and-high-y-degree-emptiness.md).
All four have status **PROVED + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED**.  The first surviving nonzero-sidecar tariff is now
[THM-3899](../01-canon/theorems/THM-3899-nonzero-sidecar-y-degree-tariff-and-equianharmonic-equality-colors.md),
with the same proved status.  THM-3900--3902 have now also been promoted as
**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED**: THM-3900
classifies the generic `f=0` chart, while THM-3901/3902 give necessary
osculating response laws in the two live `f!=0` degree regions.  None of
these results proves `JC(2)`; the planar Jacobian conjecture remains **OPEN**.

Companion:

- [exact scout](../04-computation/jc2_cusp_residual_plucker_tie_depth_and_two_arm_scout_20260823.py);
- [frozen output](../05-knowledge/results/jc2_cusp_residual_plucker_tie_depth_and_two_arm_scout_20260823.out);
- [degree-six cubic-tau scout](../04-computation/jc2_cusp_residual_cubic_tau_degree6_scout_20260823.py);
- [degree-six frozen output](../05-knowledge/results/jc2_cusp_residual_cubic_tau_degree6_scout_20260823.out).

## 1. Inheritance pass and concept board

The closest proved tournament mechanism is
[THM-2294](../01-canon/theorems/THM-2294-anchored-plucker-tournament-and-kakeya-address-bank.md):
an anchored packet is first a valued Plucker edge field; only after choosing a
tie-free real sign gauge does it cast a tournament shadow.  The canonical
hostile is the nonreduced `(6,4)` cusp, whose reduced support is one point and
whose nilpotent directions are exactly what an orientation erases.

Two corrections are load-bearing:

- [MISTAKE-212](../01-canon/MISTAKES.md) retracts the implication that
  noncancellation is equivalent to tournament transitivity.  Strict dominance
  is sufficient, not necessary.
- [MISTAKE-214](../01-canon/MISTAKES.md) separates Vandermonde node values from
  tournament score exponents.  Repeated nodes call for divided differences or
  derivative rows, not an arbitrary regular tournament.

The least-used relevant sidecars are the first marked response jet of
[THM-3377](../01-canon/theorems/THM-3377-path-colour-deletion-compiler-and-skew-current.md)
and the iterated third-order response suggested at the `ABBA/BAAB` kernel in
[THM-3380](../01-canon/theorems/THM-3380-hamiltonian-deletion-layer-monoid-semiring-and-small-order-boundaries.md).

The live concept board was:

1. the THM-3881 rank-two residual module;
2. an intrinsic skew Plucker gain rather than an ordered sign;
3. depth of a tied edge as a discrete natural-number address;
4. arm deletion and the first response missed by a leading quotient;
5. confluent/Hermite response at the length-six cusp.

The main advance came from comparing every item with item 3: the old degree
filtration was already measuring the first vanishing order of one exact edge.

## 2. Exact connection contract

Retain the notation of
[THM-3881](../01-canon/theorems/THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization.md):

```text
a=x+1,             L=9x+4,
K=y^2-15x^2-15x-4, P=aL^2,
Delta=a^2P-K^2,
r=aT+Kf,           A=KT+aPf,
B=Pf^2-T^2.
```

The tournament-to-JC transfer is lawful only in the following richer form.

| Contract field | Exact content |
|---|---|
| source | THM-3881's rank-two pair `(T,f)` and its two residual sidecars |
| target | a five-vertex complete skew-gain graph over `k[x,y]` |
| vertices | `e_T=(1,0)`, `e_f=(0,1)`, `g=(K,-a)`, `h=(aP,-K)`, `v=(T,f)` |
| binary observable | `omega(p,q)=det(p,q)` |
| preserved predicate | the full residual `S(T,f)` being a polynomial square |
| preserved data | `Delta,r,A,T,f`, hence also `B` and `S` |
| destroyed by orientation | magnitudes, zero minors, field structure, gauge debt, and cancellation |
| needed sidecar | the full gain labels plus the `A` transport edge and the first nonzero tie layer |
| cheapest hostile | algebraically closed `k` has no intrinsic positivity, and the equality seam has `omega(g_top,v_top)=0` |

The upper-triangular edge table is

```text
        e_f       g       h       v
e_T      1       -a      -K       f
e_f              -K     -aP      -T
g                         Delta    r
h                                  A
```

The four-vertex Plucker relations are not decorative.  They are exactly the
matrix-factorization reconstruction laws

```text
aA-Kr       = Delta f,
aP r-KA     = Delta T,
Delta-a^2P+K^2 = 0.                                      (1)
```

Moreover, the elementary vertex slide `v -> v+qg` is precisely the old lift
gauge:

```text
(T,f) -> (T+Kq,f-aq),
r     -> r,
A     -> A-Delta q.                                      (2)
```

Thus `r=omega(g,v)` is an invariant edge, while `A=omega(h,v)` records the
transport debt.  The complete weighted star reconstructs the residual:

```text
S = L^4
  + 2(3A+3P+r^2)L^2 f
  + (8A+6P+3r^2)(Pf^2-T^2).                              (3)
```

This is the exact tournament descendant.  At `(x,y)=(0,0)`, a chosen real
normalization gives the positive cyclic weights `(1,4,1)` on
`e_T -> e_f -> g -> e_T`, but that cycle is only a real sign shadow: vertex
rescaling switches it, and an algebraically closed field supplies no canonical
order.

## 3. Breakthrough: THM-3884 is a tie-depth theorem

Let `n=deg f>=1` and suppose the equality case of
[THM-3884](../01-canon/theorems/THM-3884-cusp-residual-total-degree-leading-gauge-filtration.md)
holds:

```text
deg T=n+1,
f_n=xq,                  T_(n+1)=-K_2q,
K_2=y^2-15x^2,           deg q=n-1.                       (4)
```

For

```text
g_top=(K_2,-x),          v_top=(T_(n+1),f_n),
```

equation `(4)` says exactly

```text
omega(g_top,v_top)=K_2f_n+xT_(n+1)=0.                     (5)
```

The degree filtration has therefore reached a tied Plucker edge, not an
undecided arc.  Since `(5)` cancels the degree-`n+2` part of `r`, put

```text
t=deg r <= n+1,          with t=-infinity when r=0.        (6)
```

The first nonzero depth of this tie is the right discrete address.  It is a
natural-number rank after reversing the finite degree scale; no forced
orientation or sign is needed.

### Exact leading packet

On `(4)`, the unique leading forms of the other two sidecars are

```text
A_(n+4)  =81x^5q,
B_(2n+3) =81x^5q^2.                                      (7)
```

Indeed, `KT` is one degree below `aPf`, and `T^2` is one degree below `Pf^2`.
Expanding `(3)` into seven macro-terms gives the complete degree ledger:

| term | degree ceiling |
|---|---:|
| `8AB` | `3n+7` |
| `3r^2B` | `2t+2n+3` |
| `6AL^2f` | `2n+6` |
| `6PB` | `2n+6` |
| `2r^2L^2f` | `2t+n+2` |
| `6PL^2f` | `n+5` |
| `L^4` | `4` |

Only the first two can lead, and their exact highest forms are

```text
8A_(n+4)B_(2n+3)    =52488x^10q^3,
3r_t^2B_(2n+3)      =243x^5q^2r_t^2.                     (8)
```

Their degree difference is exactly `2t-(n+4)`, and `52488/243=216`.

### The three cases

Suppose `S` is a square.

1. If `2t>n+4`, the second form in `(8)` is the unique leader.  Its valuation
   at the prime `x` is

   ```text
   5+2v_x(q)+2v_x(r_t),
   ```

   which is odd.  This case is impossible.

2. If `2t<n+4`, the first form in `(8)` is the unique leader.  UFD parity of
   `x^10q^3` forces every irreducible valuation of `q` to be even.  Over the
   algebraically closed field,

   ```text
   q=s^2.                                                  (9)
   ```

   Hence `n` is odd.  This condition is necessary, not sufficient.

3. If `2t=n+4`, then `n` is even and the common degree `3n+7` is odd.  The
   entire top form must cancel:

   ```text
   243x^5q^2(r_t^2+216x^5q)=0.                            (10)
   ```

   Domain and UFD parity give

   ```text
   q=x s^2,
   r_t=rho x^3s,              rho^2=-216,
   deg s=(n-2)/2.                                         (11)
   ```

   In particular `f_n=x^2s^2` is itself a square.

Over a non-algebraically-closed characteristic-zero field, `(9)` and `(11)`
retain scalar squareclasses.  Characteristic zero is also needed so that
`216` and `243` do not vanish.

### Incoming THM-3886 is the first-response specialization

The concurrently promoted THM-3886 writes the next edge symbol as

```text
R=r_(n+1)=K_2f_(n-1)+xT_n-y^2q.                           (12)
```

If `R!=0`, then `t=n+1`, so the comparison governing the three cases becomes

```text
2t-(n+4)=n-2.                                             (13)
```

This identifies all three of its regimes without a new degree table:

- `n>=3` is supercritical, hence impossible unless `R=0`.  Reducing `R=0`
  modulo `x` and using `gcd(x,K_2)=1` gives
  `f_(n-1)=q+xp` and `T_n=-K_2p+15xq`, exactly THM-3886's second gauge jet.
- `n=2` is the critical case `(11)`, with `q=xu^2` and
  `R=rho x^3u`.  A separate full-universe replay in this session verifies
  `v_x(S_12)=8` with `[x^8y^4]S_12=-648u^6`, followed by
  `v_x(S_11)=3` with `[x^3y^8]S_11=8u^6`.  The equations
  `S_12=g_6^2`, `S_11=2g_6g_5` are incompatible, for both Kummer signs.
- `n=1` is subcritical and only says that the constant `q` is a square up to
  scalar, which is automatic over `k`.  The independent square-root replay
  then forces successively `c=q`, `v=0`, and `u=15q`, leaving the constant
  pure gauge `(T,f)=(-Kq,aq)`.  THM-3881's divisibility gate `a|q` excludes it.

Thus the scout audits the complete leading trichotomy and stable-range
two-jet conclusion, and the separate homogeneous replays discharge both
exceptional obligations.  This is the independent hostile audit recorded in
the promoted THM-3886 status.

### Half-depth gauge staircase: the incoming signal iterates

The same argument does not stop at two jets.  It is the independently proved
and hostile-audited
[THM-3894](../01-canon/theorems/THM-3894-cusp-residual-all-degree-gauge-kummer-parity-filtration.md).

Let `q_i` be homogeneous of degree `n-1-i` and put

```text
W_j=q_0+...+q_(j-1),
Ttilde_j=T+K W_j,              ftilde_j=f-a W_j.
```

The first equality-seam jet is `W_1=q_0=q`.  Suppose after `j>=1` jets

```text
deg Ttilde_j<=n+1-j,           deg ftilde_j<=n-j.
```

Gauge invariance gives

```text
r=a Ttilde_j+K ftilde_j,
deg r<=n+2-j.
```

No square equation is asserted for `(Ttilde_j,ftilde_j)`.  Every use of the
tie-depth law is on the original square residual, with its original `q_0`,
`A_top`, and `B_top`; only the invariant edge `r` is computed through the
gauge-subtracted pair.

If its top component is nonzero, then `t=n+2-j`, and the tie-depth comparison
is exactly

```text
2t-(n+4)=n-2j.
```

Whenever `n>2j`, this is the impossible supercritical regime.  Hence the top
component vanishes:

```text
x(Ttilde_j)_(n+1-j)+K_2(ftilde_j)_(n-j)=0.
```

Since `gcd(x,K_2)=1`, there is a homogeneous `q_j` with

```text
(ftilde_j)_(n-j)=xq_j,
(Ttilde_j)_(n+1-j)=-K_2q_j.
```

Here later `q_j` are allowed to be zero; “homogeneous” includes that nominal
graded component.

Adding `q_j` to `W_j` produces the next full gauge jet.  Induction forces

```text
J=ceil(n/2)
```

leading jets for every equality-seam survivor with `n>=3`.  Explicitly, with
negative-index `q_i` interpreted as zero,

```text
f_(n-i)=xq_i+q_(i-1),
T_(n+1-i)=-K_2q_i+15xq_(i-1)+4q_(i-2),   0<=i<J.
```

The terminal parity passport is forced rather than guessed:

- for odd `n`, the terminal bound is strictly subcritical, so `q_0` is a
  square up to scalar;
- for even `n`, a subcritical terminal would require `n` odd and is therefore
  impossible.  The terminal edge must have degree `n/2+2` and be critical:
  `q_0=xs^2`, `r_top=rho x^3s`, `rho^2=-216`.

This is the natural-number encoding promised by the tournament analogy: the
discrete address is the number of consecutive zero Plucker responses before
the terminal parity carrier, namely `ceil(n/2)`.  The exact companion freezes
the recurrence and all odd/even boundary inequalities through 1,188 active
gates.  It remains an associated-graded classification; it does not make
gauge peeling preserve the square equation.

This is a sharp new necessary trichotomy, but not a closure theorem.  In the
critical case the next degree already depends on

```text
8(A_(n+4)B_(2n+2)+A_(n+3)B_(2n+3))
+3(r_t^2B_(2n+2)+2r_t r_(t-1)B_(2n+3)),                  (14)
```

which is the first marked response after the tied leading edge.  Formula
`(14)` is where the THM-3377 lesson becomes literal: the scalar leading
response has vanished, so its first jet is load-bearing.

## 4. Niche result: both exact arms close the `f=0` lane through degree six

The inherited theorem
[THM-3885](../01-canon/theorems/THM-3885-cusp-residual-f-zero-arm-dichotomy-and-quadratic-closure.md)
closes every `f=0` candidate through total degree three.  A forgotten
deletion/intersection idea gives a genuinely new three-degree extension on the
zero branch of both exact arms.

**Theorem (THM-3893).**  Let

```text
S(T,0)=L^4-6aL^2T^2-8KT^3-3a^2T^4,
T(0,0)=0.
```

If `S(T,0)` is a square,

```text
T(-1,y)=0,              T(-4/9,y)=0,                     (15)
deg T<=6,
```

then `T=0`.

**Proof.**  The two arm conditions give `aL | T`; write `T=aLH`.  Exact
substitution yields

```text
S=L^3 [ L(1-6a^3H^2-3a^6H^4)-8a^3KH^3 ].                (16)
```

Modulo `L`, the bracket is

```text
-8(5/9)^3(y^2-8/27)(H mod L)^3.                          (17)
```

If `L` did not divide `H`, `(16)` would have odd `L`-valuation three.  Hence
`H=LH_1` and

```text
T=aL^2H_1,              H_1(0,0)=0.                      (18)
```

Since `deg(aL^2)=3`, the degree-six universe has `deg H_1<=3`.  At `x=0`, put

```text
T(0,y)=tau,               tau=py+qy^2+ry^3.              (19)
```

The specialized residual is

```text
F=256-96tau^2-8(y^2-4)tau^3-3tau^4.                      (20)
```

After choosing the sign of a square root, write

```text
G=16+g_2y^2+...+g_6y^6.
```

Coefficients in degrees two through six determine every `g_i` recursively;
the recursion divides only by `32`.  The six remaining equations in degrees
seven through twelve have an exact 15-element grevlex Groebner basis containing

```text
r^4,
q^4-6p^2r^2,
p^5+16p^2r+16pq^2+(448/5)r^3.                             (21)
```

There is no saturation or division by `p,q,r`.  At every field-valued zero,
`(21)` gives successively `r=0`, `q=0`, and `p=0`.  Hence `H_1(0,y)=0` and

```text
H_1=xJ,                  deg J<=2.                        (22)
```

Now split by `d=deg_y J`.

- If `d=0`, then `T` depends only on `x`.  The missing `y` term and nonzero
  value `S(0,0)=256` contradict `[y^2]S=-8T^3` unless `T=0`.
- If `d=1`, let `t_1(x)=[y]T`.  The unique highest `y` term is

  ```text
  [y^5]S=-8t_1^3,
  ```

  so `S` has odd `y`-degree.
- If `d=2`, the coefficient `b=[y^2]J` is a nonzero constant and

  ```text
  t_2=[y^2]T=baL^2x,
  [y^8]S=-t_2^3(8+3a^2t_2).                              (23)
  ```

  At the prime `a=x+1`, `v_a(t_2)=1`, while the second factor in `(23)` is
  `8 mod a`.  Thus the leading `y` coefficient has odd `a`-valuation three
  and cannot be a square in `k[x]`.

All cases force `T=0`.  In degrees four and five the same proof has an
elementary shorter form: the `x=0` specialization is only quadratic, with
hostile Groebner basis `<z,q^2>` for `z=p^2`; after it vanishes, a mixed term
has unique odd `y`-degree five, and the remaining `T` is univariate.

Consequently the first possible nonzero zero/zero-arm candidate has
`deg T>=7`, equivalently `deg H_1>=4`.  The degree-seven cell is open.

The concurrent
[THM-3888](../01-canon/theorems/THM-3888-f-zero-equianharmonic-jacobian-and-two-section-integrality.md)
reframes the generic `f=0` quartic as a rank-six rational elliptic surface
with a two-boundary-section integrality problem.  The lemma here is a
complementary finite-degree **descent** checkpoint: generic Mordell--Weil
geometry does not see simultaneous evaluation at `a=0` and `L=0`, whereas
those two labelled fibres plus the first `L`-adic response close every global
polynomial point through degree six.  An effective elliptic height argument
should carry these evaluation maps as sidecars rather than expecting generic
two-section integrality alone to enforce `k[x,y]` descent.

## 5. Breakthrough: the first open even Kummer terminal is empty

THM-3894 turns the tied Plucker edge into an integer depth and forces the
even `n=4` terminal data

```text
q_0=xs^2,                    r_4=rho*x^3s,
deg s=1,                     rho^2=-216.                 (24)
```

The forgotten THM-3377 lesson is to keep the first response after this tie.
After the two gauge jets, put `W=q_0+q_1+q_2` and

```text
bar f=f-aW,                  bar T=T+KW,
D=Delta W^2,                Q=8Delta W+3r^2.             (25)
```

The Kummer equation is exactly `Q_8=0`.  Every residual term of degree at
least seventeen then collapses to the single marked product `DQ`:

```text
S_18=D_11Q_7,
S_17=D_11Q_6+D_10Q_7.                                  (26)
```

Write `s=alpha*x+beta*y`.  If `beta!=0`, `S_18=g_9^2` forces
`v_x(g_9)=4`, but

```text
[x^3y^14]S_17=8beta^6,
```

so `v_x(S_17)=3`, impossible for `2g_9g_8`.  If `beta=0`, write
`s=sigma*x`.  The degree-eighteen square forces

```text
Q_7=x^3h_2^2,            g_9=9sigma^2x^7h_2.
```

The next equation forces `h_2|Q_6` and determines `g_8`; the following shell
then has the universal coefficient

```text
[x^4y^12](S_16-g_8^2)=2sigma^6/81.                       (27)
```

Its left side has `x`-valuation four, while `2g_9g_7` is divisible by
`x^7`.  Both branches are empty.  The positive control

```text
Q_7=-8x^3K_2^2,
Q_6=6(40+rho*z)x^4K_2
```

genuinely lifts the degree-eighteen and degree-seventeen equations before
`(27)` fails.  Thus [THM-3896](../01-canon/theorems/THM-3896-cusp-residual-degree-four-equality-seam-emptiness.md)
is a response obstruction, not a disguised leading-term contradiction.

## 6. Breakthrough: the whole polynomial `f=0` lane closes

The degree-six arm result was only the finite shadow.  THM-3895 first proves

```text
G^2=S(T,0)  ==>  deg_y T<=2.                              (28)
```

Its quartic covariant gives a complete generic proof.  An independent
forgotten-arm audit gives a shorter global mechanism: with `s^2=-3`,

```text
C=8K^2-9a^3L^2,
G_0=s(aT^2+4KT/(3a)-C/(9a^3)),
(G-G_0)(G+G_0)
 =(C^2+27a^6L^4-24a^2KCT)/(27a^6).                      (29)
```

For `d=deg_yT>2`, this forces `deg_y(G-G_0)=6-d`.  The `a=0` arm says every
positive-`y` coefficient `t_i(x)` of `T` is divisible by `a`.  But then

```text
[y^4]G_0/s
 =a sum_(i+j=4)t_i t_j+4(t_2-Ft_4)/(3a)-8/(9a^3),       (30)
```

whose last term is an uncancellable `a^(-3)` pole.  This independently
proves `(28)` and identifies global descent as the load-bearing coordinate.

THM-3897 closes the three remaining channels.  If
`T=u(x)y^2+v(x)y+w(x)`, its leading coefficient in the quadratic case is

```text
-u^3(8+3a^2u).
```

The two factors are coprime.  Square parity gives `u=r^2` and
`s^2=8+3a^2r^2`; factoring the latter over `sqrt(3)` makes both factors
units of `k[x]`, hence forces `r=0`.  The linear channel has unique odd
degree five, and the constant channel's missing linear coefficient conflicts
with the address value `S(0,0)=256`.  Therefore

```text
f=0, S(T,f)=G^2, T(0,0)=0  ==>  T=0, G=+/-L^2.          (31)
```

Every remaining residual survivor must pay a genuinely nonzero sidecar.

[THM-3900](../01-canon/theorems/THM-3900-f-zero-generic-y-polynomial-root-color-response-classification.md)
also classifies the generic `k(x)[y]` chart without importing the provisional
elliptic-height claims.  At each root of the normalized quadratic `u`, the
square root has intrinsic color `+1` or `-1`.  Same-color roots force exactly
the rational hostile `T_*=-2K/(3a^2)`.  Opposite colors require their first
derivatives to vanish; the second-derivative response then forces

```text
A=-8/3, B=0, C^2=1/8, f^2=169/512.
```

This split-color family exists as a sharp constant-parameter hostile, but the
actual `f^2=F^2/(a^3L^2)` is nonconstant.  Thus the only generic coordinates
are `T=0,T_*`, and global `a`-integrality deletes the latter.  Root color plus
marked second response is therefore sufficient here; root color alone is not.

## 7. The live nonzero-sidecar response

THM-3899 gives the first all-degree tariff.  If
`m=deg_yT`, `n=deg_yf`, and `f!=0`, then

```text
m>=n.                                                     (32)
```

For `m=n>=1`, with leading coefficients `u,v,g`, write `g=vh`.  Then

```text
h^2=3(aL^2v^2-u^2),
(h-du)(h+du)=3aL^2v^2,        d^2=-3.                    (33)
```

The two color valuations at `a=x+1` sum to `1+2ord_a(v)`.  This is a lawful
binary observable, but the canonical payment `(h-du,h+du)=(a,3L^2)` proves
that color parity alone is not a no-go.

The exact [THM-3902](../01-canon/theorems/THM-3902-nonzero-sidecar-equality-color-two-jet-response.md)
restores the deleted response.  With `epsilon=y^(-1)`, `T=y^nU`, `f=y^nV`,
`G=y^(2n+2)VH`, and `C_+/-=H+/-dU`, it gives

```text
C_-C_+=3aL^2V^2+2L^2epsilon^nV
 +6epsilon^2(kappa+aU/V)(aL^2V^2-U^2)  mod epsilon^3.    (34)
```

Thus the forgotten `K^2f^3` term enters at response depth `n`: first jet at
`n=1`, second at `n=2`, and later for `n>=3`.  At `n=1` the first response
has `a`-valuation exactly one less than the leading color product; the marked
color loses its simple `a`-factor.  Explicit address-compatible data can
still lift two jets, so `(34)` is a sharp necessary response law rather than
a closure.

For `m>n`, write `delta=m-n`.  The leading fan has only three cells:

```text
delta=1:   -3u^2v^2,
delta=2:   -3u^2(au+v)^2,
delta>=3:  -3a^2u^4.                                    (35)
```

All are already squares over `k`; the only cancellation wall is
`delta=2, au+v=0`.
[THM-3901](../01-canon/theorems/THM-3901-nonzero-sidecar-osculating-response-fan.md)
uses the sign-matched root `sTr` and the exact defect

```text
S-(sTr)^2=S+3T^2r^2                                     (36)
```

to expose the first response fan and its divisibility passports.  In the
formerly automatic high-gap cell it subtracts the second osculating root
`s(Tr-a^2L^2f^2/2)` and obtains another exact three-cell fan.  These are
proved necessary laws after independent audits; their value is the exact
location of the next information, not a claim that `JC(2)` has closed.

## 8. Wildcard: confluence, not a tournament, at the sextic cusp

For the principal `(6,4)` finite-exact face, the local algebra is

```text
R_0=k[u,v]/(6v^2-u^3, u^2(2u+3v), u^4)
```

with basis `1,u,u^2,u^3,v,uv`.  Its reduced support is one point, so there is
no collection of geometric branches to orient.  In `m/m^2`,

```text
(alpha u+beta v)^2
  = alpha^2u^2+2alpha beta uv                 mod m^3.    (37)
```

Thus `kv` is the unique square-null tangent line.  But the named complementary
coordinate `ku` is not intrinsic: for `t!=1/2`,

```text
u -> u+tv,
v -> (1-2t)v-(t/4)u^2                                    (38)
```

preserves all three defining relations and has invertible tangent determinant
`1-2t`; at `t=1` it moves `ku` to `k(u+v)`.  This refutes a tournament built
from the named `u/v` coordinate pair.  It does not claim that every possible
intrinsic flag of `R_0` has been classified.

The lawful Vandermonde descendant is therefore confluent.  Compute the six
Macaulay/Hermite dual functionals of `R_0`, retain the full bucket-selected
response and its first Rees jet, and contract it with the intrinsic skew
Jacobian fibre form.  Only if that first skew response vanishes should one
compute the quadratic Kuranishi obstruction and the iterated third-order
response suggested by THM-3380.

A forgotten reconstruction key supplies a useful negative control.  The pair
`(OCF polynomial, det(I+S))` is only finite-exact at order six, not universal:
THM-3380's nonisomorphic ordered words `ABBA` and `BAAB` have the same
commutative OCF response, and
`det(I+S)=det(C+2I)` is just an evaluation of the same spectral response.
Another scalar determinant therefore cannot replace the ordered/transverse
current that the confluent cusp calculation needs.

The old open Jacobian Plucker face complex makes the same warning at global
coefficient level:

```text
[x^w]Jac(P,Q)
  = sum_(u+v=w+(1,1)) det(u,v) p_u q_v.                   (39)
```

An edge list alone misses the four-exponent, three-perfect-matching Plucker
circuit.  Any decoupled channel quotient must restore actual minor values,
zeros, output fibres, and coefficient Segre constraints before it can feed a
proof.

## 9. Generated next tasks

1. **Next even Kummer terminal.**  THM-3896 closes `n=4`; repeat its `DQ`
   response construction at `n=6`, retaining enough shells to distinguish
   the line factors of the quadratic terminal `s`.
2. **Odd terminal square lift.**  Substitute `q_0=s^2` after all forced jets,
   retain `r` rather than gauge-peeling it away, and locate the first degree at
   which the `A` transport debt becomes a transverse obstruction.
3. **Transplant the generic root response.**  THM-3900 proves that the
   `k(x)[y]` `f=0` chart has only `T=0,T_*`.  Test whether its same/split
   color and second-derivative dichotomy survives the first nonzero sidecar.
4. **Nonzero-sidecar response completion.**  Extend THM-3901/3902 one jet
   beyond their payable controls.  The primary
   targets are the wall `delta=2, au+v=0` and the equal-degree `n>=3` source
   which enters only after response depth two.
5. **Confluent sextic response.**  Build the Macaulay dual on the displayed
   six-element basis, then compute tangent space, cokernel, first nonzero
   obstruction order, and the fixed descent character.  A scalar resultant is
   not the target.
6. **Four-exponent hostile.**  In a small coefficient support with one tied
   output fibre, compare the full Plucker/Segre packet with its orientation and
   gain-only quotients.  Freeze the minimal false positive before using the
   complex at higher degree.

The conceptual shift is small but decisive: map a discrete object to the
natural number measuring **how long an intrinsic alternating edge remains
tied**, while retaining its first response sidecar.  In this JC residual that
number is `deg`-depth of `r`; it produces exact algebra immediately, whereas
choosing an orientation at the tie destroys the very cancellation that makes
the theorem work.

## 10. Reproduction

```bash
python3 04-computation/jc2_cusp_residual_plucker_tie_depth_and_two_arm_scout_20260823.py
python3 -O 04-computation/jc2_cusp_residual_plucker_tie_depth_and_two_arm_scout_20260823.py
python3 04-computation/jc2_cusp_residual_cubic_tau_degree6_scout_20260823.py
python3 -O 04-computation/jc2_cusp_residual_cubic_tau_degree6_scout_20260823.py
```

Each pair of streams must byte-match its correspondingly named frozen output.
The degree-six Groebner replay is intentionally exact and may take several
minutes per run.
