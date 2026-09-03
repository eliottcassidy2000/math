---
id: THM-4345
title: "LRC(14) half-periodic anchor-strip Euclidean remainder and the sharp current-only cubic envelope"
status: >
  PROVED ANALYTIC + VERIFIED-EXACT; LRC(14) OPEN. For every integrable
  half-periodic observable pulled back by an odd common dilation, restriction
  to the danger strip of an even anchor is exactly a quotient term plus one
  residual nested-wall integral indexed by the Euclidean division D=7q+r.
  This gives exact restricted signed-current covariances, sharp observable
  discrepancy bounds, and an exact current-energy sidecar. Separately, the
  nonlinear q-cubic rebate has a sharp lower envelope depending only on the
  absolute sheet current. Exact primitive denominator-complete controls show
  that the full Brownian pair kernel does not determine anchor-strip energy
  and that its second-moment rebate remains negative even when the strip is
  known exactly. No new common-dilation LRC family is claimed: that consequence
  is already subsumed, much more strongly, by THM-616.
source: anchor_current / LRC14 continuation session, 2026-09-02
depends_on:
  - THM-616-one-tightener-useless-all-m
related:
  - THM-4346-lrc14-halfturn-current-brownian-kernel-and-cubic-moment-boundary
  - THM-616-one-tightener-useless-all-m
  - THM-859-hamming-six-dilation-conjugacy-and-order-one-gate
  - THM-4228-common-gcd-two-outsider-periodic-observable-haar-ray
  - THM-4233-pair-specific-primitive-observable-oscillation-haar-charts
artifacts:
  - 04-computation/lrc14_anchor_current_strip_probe_anchor_current_20260902.py
  - 05-knowledge/results/lrc14_anchor_current_strip_probe_anchor_current_20260902.out
script_sha256: e6c0260bbce54e31c1f7282c90eca7bd8cc7617d563397bc1fadff344d5842f0
output_sha256: 09c8368547f251893796eba24b53f4e0b3fc6b1fe9f7c80ff86615c08c249dec
hash_basis: raw LF bytes
audit: >
  PASS. A literal Fraction wall sweep checks the identity on the constant,
  current-energy, signed cross-current, q-cubic, free-sheet-count, and
  symmetric one-sheet Bonferroni observables. Ten odd multipliers cover every
  residue r=0,...,6, both quotient parities, nontrivial gcd reductions, exact
  independence at r=0, and the r=1/r=6 narrow-wall/complement controls. A
  separate 16,800-case rational audit checks the general radius/phase
  transducer. The same program checks every absolute-current envelope value
  through depth twelve and the exact physical hostile margins. The analytic
  proof is not finite-range.
---

# THM-4345 -- the half-periodic anchor-strip remainder operator

**PROVED ANALYTIC + VERIFIED-EXACT. LRC(14) REMAINS OPEN.**

## 1. Inheritance and the exact operator

The closest proved mechanism is THM-4228's periodic-observable endpoint
cocycle. The canonical hostile is common odd dilation: it leaves the full
current kernel unchanged while moving energy into and out of the fixed
anchor strip. The corrected near miss is therefore to subtract an
"independent `1/7`" of the full energy without an address. The least-used
sidecar is the Euclidean-remainder wall below, including the parity of its
quotient.

Put, up to null endpoints,

```text
A(x)=1_{||x||<1/14}.                                  (1)
```

Let `h,m` be positive integers with `m` odd, and let `f` be any integrable
function on `R/Z` satisfying

```text
f(x+1/2)=f(x) a.e.                                    (2)
```

Define

```text
c=gcd(m,2h)=gcd(m,h),   D=m/c,   H=2h/c,
D=7q+r,                 0<=r<=6.                      (3)
```

Here `c` and `D` are odd, `H` is even, and `gcd(D,H)=1`. On the ordinary
half-base `[0,1/2)`, put

```text
F(f)=integral_0^(1/2) f(x)dx,
R_(H,q,r)={x: ||Hx+q/2||<r/14},
R(f)=integral_(R_(H,q,r)) f(x)dx.                     (4)
```

For `r=0`, the residual set is empty. In every case

```text
|R_(H,q,r)|=r/14.                                     (5)
```

> **Exact half-periodic anchor-strip identity.**
>
> ```text
> I_(h,m)(f)
> :=integral_0^(1/2) f(mt) A(2ht)dt
>  =[q F(f)+R(f)]/D.                                  (6)
> ```

Consequently the complementary anchor-safe core satisfies

```text
J_(h,m)(f)
:=integral_0^(1/2) f(mt)[1-A(2ht)]dt
 =[(D-q)F(f)-R(f)]/D,                                 (7)

I_(h,m)(f)-F(f)/7=[R(f)-(r/7)F(f)]/D,
J_(h,m)(f)-6F(f)/7=-[R(f)-(r/7)F(f)]/D.               (8)
```

Thus an apparently three-frequency restriction has only one new coordinate:
the integral of the primitive observable on a single nested residual wall.
The quotient `q`, the remainder `r`, and the parity of `q` specify that wall.

### Proof

Both factors in the integrand of `(6)` are half-periodic in `t`: `(2)` and
oddness of `m` handle the first, while `H` is even for the anchor factor.
It is therefore enough to compute twice `(6)`, on the full circle. The
degree-`m` change of variables `x=mt` gives

```text
2I_(h,m)(f)=integral_0^1 f(x) Q(x)dx,

Q(x)=(1/m) sum_(j=0)^(m-1) A(2h(x+j)/m).              (9)
```

The multiset of shifts in `(9)` consists of each point of a `D`-grid exactly
`c` times. Since `gcd(D,H)=1`, multiplication by `H` permutes that grid, so

```text
Q(x)=(1/D) sum_(j=0)^(D-1) A(Hx/D+j/D).               (10)
```

The following grid-count identity is the load-bearing step. For every real
`y`, away from the null set of endpoints,

```text
sum_(j=0)^(D-1) A(y+j/D)
 =q+1_{||Dy+q/2||<r/14}.                              (11)
```

One elementary proof counts integers in a sliding open interval of length
`D/7=q+r/7`. A useful sign audit comes from Fourier series. Averaging the
Fourier series of `A` over the grid kills every frequency not divisible by
`D`, leaving coefficients

```text
sin(pi nD/7)/(pi nD)
 =(-1)^(nq) sin(pi nr/7)/(pi nD).                     (12)
```

These are exactly the coefficients of
`[q+1_{||Dy+q/2||<r/14}]/D`, including the half-shift when `q` is odd.
Substitute `y=Hx/D` in `(11)`, then use the half-periodicity of `f` and of
the residual indicator (`H` is even). Equations `(9)--(11)` become `(6)`.
The residual indicator has duty cycle `r/7` on the full circle and hence
mass `r/14` on the half-base, proving `(5)`. Subtracting `(6)` from `F(f)`
gives `(7)--(8)`.

### The nested-wall transducer

The same proof resolves the iteration suggested by `(11)`. For
`s in {0,...,6}` and a phase bit `epsilon in {0,1}`, set

```text
A_(s,epsilon)(x)=1_{||x+epsilon/2||<s/14}.             (11a)
```

For any positive odd `D`, divide

```text
sD=7Q+R,                 0<=R<=6,
epsilon' =D epsilon+Q                         (mod 2). (11b)
```

Then, almost everywhere,

```text
sum_(j=0)^(D-1) A_(s,epsilon)(y+j/D)
 =Q+A_(R,epsilon')(Dy).                                (11c)
```

Thus repeated strip pushforwards form a finite transducer on the fourteen
states `(s,epsilon)`: multiply the radius numerator by `D`, retain its
remainder modulo seven, and add the quotient parity to the phase bit. The
rule is associative. Indeed, if a second degree `E` gives
`RE=7Q_2+R_2`, then the sequential total quotient is `QE+Q_2`, exactly the
quotient in `sDE=7(QE+Q_2)+R_2`, and the sequential phase is

```text
E(D epsilon+Q)+Q_2
 =DE epsilon+(EQ+Q_2)                         (mod 2). (11d)
```

This proves that the nested-wall sidecar does not proliferate under further
odd common dilations; it stays one radius numerator and one phase bit.

## 2. Sharp discrepancy bounds

Equation `(8)` immediately gives, for an arbitrary signed integrable `f`,

```text
|I_(h,m)(f)-F(f)/7|
 <=eta_r/D integral_0^(1/2)|f(x)|dx,                  (13)

eta_0=0,
eta_r=max(r,7-r)/7                    (1<=r<=6).       (14)
```

Indeed, the multiplier of `f` in `(8)` is
`1_R-r/7`, whose two absolute values are `1-r/7` and `r/7`.

If `0<=f<=M` and `F=F(f)`, then the exact residual formula gives the
one-sided estimates

```text
qF/D <= I_(h,m)(f) <= (q+1)F/D,                       (15)

-rF/(7D) <= I_(h,m)(f)-F/7
 <=r(M/14-F/7)/D,                                     (16)
```

and the range-refined sharp interval

```text
[qF+max(0,F-M(7-r)/14)]/D
 <=I_(h,m)(f)
 <=[qF+min(F,Mr/14)]/D.                               (17)
```

The word `sharp` in `(17)` is for the unrestricted measurable-observable
class with prescribed `F` and range `[0,M]`: it is attained by putting as
much mass as possible inside, or outside, the residual set. In `(15)`, the
lower equality means `f=0` almost everywhere on the residual wall; the upper
equality means `f=0` almost everywhere off it. These are not asserted to be
simultaneously attainable by a physical runner observable.

When `r=0`, equations `(6)` and `(8)` give exact independence:

```text
I_(h,m)(f)=F(f)/7,            J_(h,m)(f)=6F(f)/7.      (18)
```

## 3. Restricted signed-current covariance and energy

For an odd positive speed `w`, define the sheet current

```text
sigma_w(x)=A(wx)-A(w(x+1/2)).                          (19)
```

For odd `u,v`, the product `f=sigma_u sigma_v` is half-periodic. Let

```text
K(u,v)=integral_0^1 sigma_u(x)sigma_v(x)dx.            (20)
```

The exact full-circle Brownian kernel computes `K(u,v)`. Applying `(6)` gives
the exact restricted covariance formula

```text
integral_0^(1/2) A(2ht)sigma_(mu)(t)sigma_(mv)(t)dt
 =q K(u,v)/(2D)
  +(1/D) integral_(R_(H,q,r)) sigma_u(x)sigma_v(x)dx. (21)
```

This displays precisely what the full pair kernel forgets: one signed
residual-wall covariance per pair.

Equation `(21)` is not restricted to a globally common-dilate family. For
any two actual odd speeds `a,b`, take

```text
m=gcd(a,b),             u=a/m, v=b/m,                 (21a)
```

and form `c,D,H,q,r` from this `m` as in `(3)`. Then `(21)` computes their
anchor-strip covariance. Therefore every arbitrary odd family `W` has the
exact finite decomposition

```text
V_strip(h;W)=sum_((a,b) in W^2) Lambda_h(a,b),         (21b)
```

where `Lambda_h(a,b)` is the right side of `(21)` after the reduction
`(21a)`. The full Brownian term is scalar, while the remaining debt is a
labelled bank of primitive-pair residual walls. Different pairs generally
have different `D,H,q,r`, so they cannot be collapsed to one unlabelled
energy statistic.

More compactly, let `B` be any finite odd-speed family and put

```text
C_B(x)=sum_(w in B) sigma_w(x),
V_B=integral_0^(1/2) C_B(x)^2 dx.                     (22)
```

Since `C_(mB)(t)^2=C_B(mt)^2`, the anchor-strip energy of the common-dilate
family is

```text
V_strip(h;mB)
 =[q V_B+integral_(R_(H,q,r)) C_B(x)^2dx]/D.          (23)
```

In particular,

```text
q V_B/D <=V_strip(h;mB)<=(q+1)V_B/D,                  (24)

(D-q-1)V_B/D <=V_core(h;mB)<=(D-q)V_B/D.              (25)
```

The inequalities are pointwise-optimal absent more information about where
the energy lies. They are not an assertion that current energy can realize
the arbitrary equality cases.

The residue extremes make the nested-wall structure concrete:

- If `D=7q`, the strip takes exactly `1/7` of every half-periodic observable.
- If `D=7q+1` with `q` even, the residual is the original narrow
  radius-`1/14` anchor wall, at the reduced anchor `H`.
- If `D=7q+6` with `q` odd, the residual radius is `3/7` around the half
  phase, hence it is the complement of that same narrow wall.

So the Euclidean remainder does not merely bound the strip; it says which
coarser or complementary wall must be remembered.

## 4. The sharp current-only envelope for the q-cubic rebate

Let

```text
p(n)=1-n+binom(n,2)-binom(n,3),
q_*=min(a,b),       C=a-b.                            (26)
```

The values of `p` are nonincreasing on the nonnegative integers, with the
plateau `p(1)=p(2)=p(3)=0`. Hence

```text
p(min(a,b))=max(p(a),p(b))
 =[p(a)+p(b)+|p(a)-p(b)|]/2.                          (26a)
```

Thus the exact q-cubic identity writes its nonlinear rebate as
`|p(a)-p(b)|/2`. If only the absolute current `d=|C|` is retained, the sharp
pointwise lower envelope over all nonnegative integer depths is

```text
g(0)=g(1)=g(2)=0,       g(3)=1/2,
g(d)=d(d^2-6d+11)/12                    for d>=4.      (27)
```

Equivalently,

```text
p(min(a,b))
 >=[p(a)+p(b)]/2+g(|a-b|).                            (28)
```

For a sheet-symmetric region `E`, if

```text
B_3(E)=integral_E p(a)=integral_E p(b),                (29)
```

then

```text
integral_E p(min(a,b))
 >=B_3(E)+integral_E g(|C|).                          (30)
```

This is the strongest possible pointwise bound depending on the two
one-sheet cubic terms and `|C|` alone.

To prove `(27)`, put

```text
delta_n=p(n)-p(n+1).
```

Then

```text
delta_0=1, delta_1=delta_2=0,
delta_n=binom(n-1,2)                         (n>=3).   (31)
```

For fixed `d`, half the difference `p(x)-p(x+d)` is half the sum of `d`
successive `delta` values. The minimum is zero for `d=0,1,2`, and is `1/2`
for `d=3`. For `d>=4`, shifting the block right from `x=0` replaces
`delta_x` by `delta_(x+d)` and cannot decrease the sum; it increases already
at the first shift. Thus the minimum is at `(a,b)=(0,d)` up to swap, and

```text
[p(0)-p(d)]/2=d(d^2-6d+11)/12.                        (32)
```

This proves sharpness and `(28)--(30)`. The earlier quadratic-current bound

```text
g(d)>=max(d^2-4,0)/12                                 (33)
```

is valid but can lose the decisive signal by collapsing the current
histogram to one second moment.

## 5. Exact physical controls and the second-moment hostile

Take

```text
h=420,
B=(1,3,5,7,9,11,13,15,17,19,21,45).                  (34)
```

The rows `{840} union mB` for `m=1,13,127` are primitive and satisfy the
inherited denominator-completeness conditions. Their full half-current
energy is the same:

```text
V_B=35666/15015.                                      (35)
```

The anchor-strip energies are nevertheless different. If

```text
S=15361/45220                                         (36)
```

is the `m=1` strip energy, then `(23)` gives the exact relations

```text
V_strip(h;B)=S,
V_strip(h;13B)=(2V_B-S)/13,       # D=13=7+6
V_strip(h;127B)=(18V_B+S)/127.    # D=127=7*18+1      (37)
```

The `m=49` nonprimitive control has `c=7`, `D=7`, and hence exact
independence `V_strip=V_B/7` as required by `(18)`. Thus common dilation
preserves the full Brownian pair energy but not the anchor-strip energy.
Formula `(23)`, rather than the full kernel alone, carries the missing
coordinate.

There is also a decisive hostile to the second-moment completion. On the
primitive denominator-complete `m=1` row, the exact anchor-core values are

```text
V_core=39490603/19399380,
B_3=-14821441/31951920,

B_3+[V_core-12/7]/12
 =-71225183/162954792<0,                              (38)

B_3+(1/12) integral_core max(C^2-4,0)
 =-3392354543/9777287520<0.                           (39)
```

So even exact knowledge of the strip and even the exact positive-part
quadratic current do not close this row. But the sharp current-only envelope
does:

```text
B_3+integral_core g(|C|)
 =2572403/33948915>0,                                 (40)

integral_core p(min(a,b))
 =526429/5819814>0.                                   (41)
```

The positive signal therefore lives above the second moment but need not
retain the full total-depth coordinate `a+b`: the absolute-current histogram
already supplies a positive, though nonexact, certificate here.

## 6. Novelty boundary and scope

The closest general mechanism is THM-4228's periodic-observable endpoint
cocycle. It gives a component-count/BV discrepancy bound for a general
periodic target. Equation `(6)` is a different specialization: the target is
the exact radius-`1/14` even-anchor strip, half-periodicity is retained, and
the discrepancy is evaluated exactly by one Euclidean-remainder wall rather
than bounded by an oscillation norm.

THM-4233 also contains a `7q+r` law, but for the density and primitive
oscillation of a varying two-comb pair. It is not the grid average `(11)` or
the observable transport `(6)`. THM-859 proves exact common-dilation
conjugacy and identifies ramification degree for component actions; it does
not compute intersection with an undilated anchor strip. The new result is
the synthesis `(6)--(8)` and its covariance/current specializations.

Applying `(7)` to the nonnegative free-sheet observable does give a cofinal
positivity estimate. If `M=max B`, settled LRC through thirteen and the
Lipschitz interval give `F>=1/(182M)`, while `(5)` gives the conservative
bound

```text
J>=6F/7-r/(14D)>0                   when D>91M.        (42)
```

This is **not new family-level content**. After dividing by `c`, the physical
row is `DB union {H}` with `gcd(D,H)=1`; whenever `D>1`, THM-616 already gives
the stronger margin at least `1/13`, with no size threshold. Equation `(42)`
is retained only as a normalization check and an illustration of the
operator. Neither `(6)` nor `(30)` proves the remaining arbitrary minority-
anchor branch, and LRC(14) remains open.

## 7. Connection contract and reproduction

```text
source:       half-periodic primitive observable f and row {2h} union mB
target:       its anchor-danger strip and anchor-safe core integrals
map:          degree-m pushforward, gcd quotient, D=7q+r grid count
preserved:    exact dt normalization, sheet swap, common-dilate observable
destroyed:    pointwise address inside each D-fibre after averaging
sidecar:      one nested wall R_(H,q,r), including q parity
hostile:      the full current energy is dilation-invariant while strip
              energy varies; quadratic rebate remains negative on (34)
decisive test: residual-wall integral in (6), or g(|C|) beyond variance.   (43)
```

Reproduce from the repository root:

```bash
python3 -B 04-computation/lrc14_anchor_current_strip_probe_anchor_current_20260902.py
```

The frozen output is
`05-knowledge/results/lrc14_anchor_current_strip_probe_anchor_current_20260902.out`.
