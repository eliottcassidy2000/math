# Planar JC: the tournament descendant is a Plucker tie-depth response, not an orientation

**Session status (2026-08-23).**  The five-vertex skew-gain identities and the
equality-seam three-case law below are **PROVED algebraically,
VERIFIED-EXACT, and INDEPENDENTLY AUDITED** as necessary statements.  The
two-arm `f=0` closure through total degree six is now canonized as
[THM-3893](../01-canon/theorems/THM-3893-cusp-residual-f-zero-two-arm-degree-six-closure.md),
with status **PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED**.  The
companion scripts have normal/optimized byte matching as their reproducibility
gate.  During the session,
[THM-3886](../01-canon/theorems/THM-3886-cusp-residual-equality-seam-second-layer-trichotomy.md)
was independently promoted on `origin/main` as a proof candidate.  This
session then completed the requested hostile audit: the tie-depth theorem
independently verifies its leading filtration envelope and stable range, while
separate full homogeneous replays verify both exceptional `n=1,2` closures.
THM-3886 is now `PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED`.
None of these results proves `JC(2)`; the planar Jacobian conjecture remains
**OPEN**.

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

The same argument does not stop at two jets.  This is an independently proved
and hostile-audited algebraic packet for the currently
`RESERVED / UNPROVED EMPTY STUB`
[THM-3894](../01-canon/theorems/THM-3894-cusp-residual-all-degree-gauge-kummer-parity-filtration.md);
the reserved file is not promoted here.

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

## 5. Wildcard: confluence, not a tournament, at the sextic cusp

For the principal `(6,4)` finite-exact face, the local algebra is

```text
R_0=k[u,v]/(6v^2-u^3, u^2(2u+3v), u^4)
```

with basis `1,u,u^2,u^3,v,uv`.  Its reduced support is one point, so there is
no collection of geometric branches to orient.  In `m/m^2`,

```text
(alpha u+beta v)^2
  = alpha^2u^2+2alpha beta uv                 mod m^3.    (24)
```

Thus `kv` is the unique square-null tangent line.  But the named complementary
coordinate `ku` is not intrinsic: for `t!=1/2`,

```text
u -> u+tv,
v -> (1-2t)v-(t/4)u^2                                    (25)
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
  = sum_(u+v=w+(1,1)) det(u,v) p_u q_v.                   (26)
```

An edge list alone misses the four-exponent, three-perfect-matching Plucker
circuit.  Any decoupled channel quotient must restore actual minor values,
zeros, output fibres, and coefficient Segre constraints before it can feed a
proof.

## 6. Generated next tasks

1. **Even terminal Kummer lift.**  After the forced half-depth staircase,
   parameterize the terminal `q_0=xs^2`, `r_top=rho x^3s` for `n=4,6`.
   Test whether THM-3886's consecutive-valuation contradiction at `n=2`
   scales or whether a genuine higher Kummer carrier first appears.
2. **Odd terminal square lift.**  Substitute `q_0=s^2` after all forced jets,
   retain `r` rather than gauge-peeling it away, and locate the first degree at
   which the `A` transport debt becomes a transverse obstruction.
3. **Degree-seven elliptic/descent cell.**  Put `H_1` quartic in `(18)`.  The `x=0`
   specialization now has degree four, while mixed `y^3`/`y^4` channels can
   tie rather than die by a unique odd degree.  Express the same shell in
   THM-3888's elliptic coordinates and retain the full bivariate `(a,L)`-adic
   response instead of growing a fixed arm-jet bank.
4. **Off-equality degree region.**  Rebuild the seven-macro-term degree ledger
   for `deg T>=deg f+2`, using the weighted gain star to identify the first
   genuine tie rather than applying another coarse total-degree maximum.
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

## 7. Reproduction

```bash
python3 04-computation/jc2_cusp_residual_plucker_tie_depth_and_two_arm_scout_20260823.py
python3 -O 04-computation/jc2_cusp_residual_plucker_tie_depth_and_two_arm_scout_20260823.py
python3 04-computation/jc2_cusp_residual_cubic_tau_degree6_scout_20260823.py
python3 -O 04-computation/jc2_cusp_residual_cubic_tau_degree6_scout_20260823.py
```

Each pair of streams must byte-match its correspondingly named frozen output.
The degree-six Groebner replay is intentionally exact and may take several
minutes per run.
