---
id: THM-4056
title: "Divisor phase compiler and Duffin--Schaeffer LCM clock"
status: >
  PROVED elementary identities + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED, with CITED Khinchin and Duffin--Schaeffer inputs. Every
  cyclic N-clock is the disjoint union of primitive exact-denominator packets,
  with an explicit inverse and Ramanujan transform. A second natural-edge
  lift clock and a reduced-residue phase clock have the same finite
  Duffin--Schaeffer mean but different information loss: the lift is divisor-
  invertible, while the phase clock retains only radical-class aggregates.
  Pointwise named-number firewalls, a finite digit-geometric-mean exponent
  theorem, reciprocal-rigidity, and exact e+pi proof gates are included.
  Nothing here proves e+pi irrational, identifies Khinchin's constant at a
  named number, reproves Duffin--Schaeffer, or proves LRC(14).
source: codex-khinchin-ds-rational-tournament-20260824
audit: >
  The primary exact-denominator companion checks the cyclic compiler,
  inverse, gcd strata, unit orbits, rational weights, Ramanujan transforms,
  prefix filters, prime-power jumps, reflections, Pythagorean content, and
  overlap hostiles. The secondary and independent companions separately
  audit reduction fibres, radical versus lift clocks, Möbius inversion,
  metric normalization, continued-fraction consequences, reciprocal loss,
  and named-point counterexamples, in normal and optimized modes.
related:
  - HYP-2211-rational-shadow-carrier-obstruction
  - HYP-2212-quadratic-discriminant-carrier-pi-e
  - HYP-2247-recursive-fourth-face
  - HYP-2628-lrc14-euler-copy-squarefree-profile
  - THM-2000-support-harmonic-abel-dini-figurate-surface
  - THM-2041-frobenius-stability-of-exact-period-projectors
  - THM-873-ramanujan-fourier-expansion-of-interval-core-good-sets
  - THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone
  - THM-3509-reduced-fraction-harmonic-k4-face-and-fibonacci-unit-cassini-ray
  - THM-3744-pell-prefix-loneliness-constant-carry-exact-formula
  - THM-3756-odd-square-ordinal-berggren-affine-descent
  - THM-4042-prime-sector-ap-cover-exact-clock-and-holonomic-law
  - THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge
script: 04-computation/divisor_phase_duffin_schaeffer_clock_thm4056.py
output: 05-knowledge/results/divisor_phase_duffin_schaeffer_clock_thm4056.out
script_sha256: a15fb134839f84b5e1a4f07131f893af790adcc2c6b8189780c851a9e90d77e9
output_sha256: bce0b8ccc22db995a0d96560aa66c7707b581b45e419a02589920a4249b28947
secondary_script: 04-computation/divisor_phase_duffin_schaeffer_thm4056.py
secondary_output: 05-knowledge/results/divisor_phase_duffin_schaeffer_thm4056.out
secondary_script_sha256: 9b9f234a2f9b49b4c254063d04a86d37324906c00c3c7789b6603ba0c7fe66a9
secondary_output_sha256: fa1c2c019c5c8e31d0496a2f6a0229885305c113dd9fbddac772a2fe1a09f601
independent_audit_script: 04-computation/divisor_phase_duffin_schaeffer_thm4056_independent_audit.py
independent_audit_output: 05-knowledge/results/divisor_phase_duffin_schaeffer_thm4056_independent_audit.out
independent_audit_script_sha256: efb303a62ed272111b23f56c61befc6e5ba940eae9bc6fc5fe0090ed09a85cd4
independent_audit_output_sha256: 444a040cb2f9679043c2c30d058d5a200d55cc3ae59d5a6fc1bd6e9871d33ed1
hash_basis: raw bytes
---


# THM-4056 -- primitive rational addresses compile exactly into one LCM clock

**PROVED + VERIFIED-EXACT.**

This theorem sharpens the exact-period packet observation in
[HYP-2628](../../05-knowledge/hypotheses/HYP-2628-lrc14-euler-copy-squarefree-profile.md).
The identity `sum_(d|N) phi(d)=N` was already present there.  The new content
is the explicit labelled bijection and inverse, its direct-system and
prime-power laws, the arbitrary weighted/Ramanujan transform, the exact
Duffin--Schaeffer block and prefix dictionaries, and the primitive
Pythagorean subpacket.  No novelty claim is made for Euler's identity or for
the cited Duffin--Schaeffer theorem.

## 1. The divisor-phase compiler

Write

```text
C_N = Z/NZ,
U_q = (Z/qZ)^*                         (q>1),
U_1 = {0}.                              (1)
```

The convention at `q=1` is load-bearing: the one residue modulo one is the
identity of the trivial unit group and represents the zero phase.

### Theorem 1.1 (explicit compiler and inverse)

For every `N>=1`, the map

```text
iota_N : disjoint_union_(q|N) U_q -> C_N,
iota_N(q,a) = (N/q)a mod N                              (2)
```

is a bijection.  If `x` is represented by `0<=x<N`, put

```text
g_N(x)=gcd(x,N),             d_N(x)=N/g_N(x).           (3)
```

Its inverse is

```text
x |-> (d_N(x), x/g_N(x)),                               (4)
```

where the second coordinate is read as `0 in U_1` when `x=0`.
Consequently

```text
#{x in C_N : d_N(x)=d}=phi(d)             for d|N,      (5)
sum_(d|N) phi(d)=N.                                      (6)
```

Moreover, multiplication by `U_N` preserves every `d_N` stratum and acts
transitively on it.  Thus the exact-denominator packets are precisely the
`U_N`-orbits in the additive clock.

### Proof

For `(q,a)` in the source of `(2)`,

```text
gcd((N/q)a,N)=N/q,                                      (7)
```

because `gcd(a,q)=1`.  Therefore `(3)--(4)` recover `(q,a)`, proving
injectivity.  Conversely `(4)` is reduced because
`gcd(x/g_N(x),N/g_N(x))=1`, and substituting it into `(2)` returns `x`.
This proves surjectivity and `(5)--(6)`.

If `d|N`, reduction `U_N -> U_d` is onto.  Prime by prime, lift a unit modulo
`p^e|d` to any congruent unit modulo the larger `p^E|N`, and use `1` at primes
not dividing `d`; the Chinese remainder theorem joins the lifts.  Hence
`U_N` is transitive on the image of `U_d` in `(2)`.

The rational meaning is literal:

```text
x/N reduced = a/d_N(x),       x=(N/d_N(x))a.            (8)
```

So `d_N(x)` is both the additive order of `x in C_N` and the exact
denominator of its circle point.

## 2. Every denominator weight has an exact Ramanujan transform

Let `W` take values in any commutative ring (or, additively, any abelian
group).  Grouping `C_N` by `(2)` gives

```text
sum_(x in C_N) W(d_N(x))
  =sum_(d|N) phi(d)W(d).                                (9)
```

For the Fourier refinement, let `W` be complex-valued (or take values in a
`Q(zeta_N)`-algebra with the displayed exponential interpreted as `zeta_N`).
For an integer frequency `k`, let

```text
c_d(k)=sum_(a in U_d) exp(2*pi*i*k*a/d)                 (10)
```

be the Ramanujan sum, with `c_1(k)=1`.  The full discrete Fourier transform
is

```text
sum_(x in C_N) W(d_N(x)) exp(2*pi*i*k*x/N)
  =sum_(d|N) W(d)c_d(k).                                (11)
```

Indeed, substituting `x=(N/d)a` in the left side produces `(10)` packet by
packet.  Equation `(11)` is the character-side refinement of `(9)` and the
finite-clock analogue of the exact-period projectors in
[THM-2041](THM-2041-frobenius-stability-of-exact-period-projectors.md).
It also explains why the primitive-centre transforms in
[THM-873](THM-873-ramanujan-fourier-expansion-of-interval-core-good-sets.md)
are Ramanujan sums: they are the Fourier images of the same denominator
orbits.

## 3. Duffin--Schaeffer blocks and prefixes

For a nonnegative function `psi`, the Duffin--Schaeffer layer at denominator
`q` uses the reduced centres

```text
a/q,                  a in U_q,                         (12)
```

with circle radius `psi(q)/q`.  Its formal two-radius mass is

```text
2phi(q)psi(q)/q.                                        (13)
```

Koukoulopoulos--Maynard prove that divergence of the infinite sum of the
half-masses `phi(q)psi(q)/q` forces almost-everywhere infinite approximation;
the exact cited statement and boundary are recorded in the
[continued-fraction/Duffin--Schaeffer primary pin](../../05-knowledge/reference/continued-fractions-khinchin-duffin-schaeffer-pins.md).

### Corollary 3.1 (finite divisor-complete block)

For every `N`, the complete block `q|N` becomes one labelled cyclic sum:

```text
sum_(q|N) 2phi(q)psi(q)/q
 =sum_(x in C_N) 2psi(d_N(x))/d_N(x).                  (14)
```

More strongly, `(2)` preserves each centre, exact denominator, radius, cyclic
order, and pairwise circle distance.  It therefore preserves the full finite
labelled overlap graph of the approximation arcs, not merely `(14)`.

If `0<=psi(q)<=1/2`, the arcs inside one fixed denominator layer have
disjoint interiors, so `(13)` is that layer's union measure.  Arcs from
different denominator layers can overlap.  Therefore `(14)` is a
multiplicity/first-moment identity and is **not** the measure of the whole
block union.

### Corollary 3.2 (exact denominator prefix)

For an integer `Q>=1`, let

```text
L_Q=lcm(1,2,...,Q).                                    (15)
```

Every `q<=Q` divides `L_Q`, and

```text
{Duffin--Schaeffer centres with q<=Q}
 <-> {x in C_(L_Q) : d_(L_Q)(x)<=Q}.                   (16)
```

The number of compiled prefix centres is

```text
sum_(q<=Q) phi(q).                                     (17)
```

The filter in `(16)` cannot be omitted.  The full clock contains all exact
denominators `d|L_Q`, often including many `d>Q`; for example `L_5=60`
contains `6,10,12,15,20,30,60` in addition to the prefix denominators.

## 4. Prime powers are exactly the clock-growth events

If `N|N'`, the phase-preserving injection is

```text
C_N -> C_(N'),        x |-> (N'/N)x.                   (18)
```

It preserves the rational point `x/N` and its exact denominator.  Ordinary
residue reduction is the wrong variance: the half phase is `3 in C_6` and
`6 in C_12`, while reducing `6 mod 6` gives zero.

Put `N_Q=L_Q`.  For `Q>=2`, `N_Q=N_(Q-1)` unless `Q=p^e` is a prime power.  At a
prime-power jump,

```text
N_Q=pN_(Q-1),
#(C_(N_Q) \ image(C_(N_(Q-1))))=(p-1)N_(Q-1),          (19)
```

and the new phases are exactly

```text
{x : v_p(d_(N_Q)(x))=e}.                               (20)
```

Only the `p`-adic denominator depth changed, so a divisor of `N_Q` fails to
divide `N_(Q-1)` exactly when its `p`-valuation is the new exponent `e`.
This proves `(19)--(20)`.

The four clocks from the preceding session are therefore

```text
L_3=6,        L_5=60,        L_7=420,        L_11=27720. (21)
```

They are the first four odd-prime samples of one valuation-depth direct
system, not four unrelated coincidences.  Apéry-style denominator clearing
uses the same `L_Q` (and powers `L_Q^m`) in the opposite direction: it clears
all small denominators instead of enumerating their primitive centres.  This
shared clock does not make a convergent sequence an irrationality proof;
the error must still beat the cleared height, as required by MISTAKE-209 and
the [Apéry framework](../../05-knowledge/reference/apery-style-irrationality-framework.md).

## 5. The primitive Pythagorean/Berggren subpacket

Let `N=2^s m` be even with `m` odd.  Under `(4)`, retain the nonzero phases
whose reduced pair `(a,d)` has opposite parity.  Then

```text
(a,d) |-> (d^2-a^2, 2ad, d^2+a^2)                     (22)
```

is the usual injective Euclid parametrization of primitive ordered
Pythagorean triples with `0<a<d` and `d|N`.  Their number is

```text
B(N)=N-(m+1)/2.                                        (23)
```

For even `d`, every unit `a` is odd, contributing `phi(d)`.  For odd `d>1`,
the involution `a |-> d-a` pairs one even and one odd unit, so exactly
`phi(d)/2` units have parity opposite to `d`.  Hence

```text
B(N)
 =sum_(d|N, 2|d)phi(d)+(1/2)sum_(d|m,d>1)phi(d)
 =(N-m)+(m-1)/2,
```

which is `(23)`.  On the master clocks this gives

| `N` | odd part `m` | primitive Pythagorean phases |
|---:|---:|---:|
| `6` | `3` | `4` |
| `60` | `15` | `52` |
| `420` | `105` | `367` |
| `27720` | `3465` | `25987` |

This is an LCM-height slice of the Berggren tree, not a ternary-depth level.
It retains denominator divisibility but forgets Berggren ancestry; THM-4057
supplies the orthogonal Stern--Brocot-depth grading.

## 6. Symmetry and loss ledger

The native circle reflection is

```text
x |-> -x in C_N,
a/d |-> (d-a)/d       (d>1),                           (24)
```

The unique `d=1` point is fixed modulo one.  Negation preserves `d_N`, every
denominator weight, the approximation radius, and the circle metric.
Multiplication by any unit of `C_N` preserves `d_N`, weights, and radii, but
need not preserve the circle metric or overlap geometry.  In contrast,
reciprocal
edge reversal

```text
a/d |-> d/a                                           (25)
```

usually changes the denominator and leaves the `[0,1]` chart.  It is the
global Stern--Brocot mirror used in THM-4057, not the Duffin--Schaeffer circle
reflection.

The transfer contract is:

| source | target/map | preserved | destroyed / sidecar |
|---|---|---|---|
| primitive fractions `a/d`, `d|N` | `C_N` by `(2)` | centre, exact denominator, cyclic distance, radius | none inside the finite divisor block |
| prefix `d<=Q` | `C_(L_Q)` | every prefix centre | extra divisor denominators unless filter `(16)` is retained |
| denominator weight | phase function `W(d_N(x))` | total and all Fourier modes `(9)--(11)` | numerator-specific weights |
| DS layers | labelled arcs on `C_N` | finite centres and overlap geometry | infinite limsup and second-moment independence |
| clock phase | reciprocal rational edge | projective rational | denominator/radius; keep homogeneous height |
| opposite-parity phase | primitive Euclid triple | ordered primitive triple | LCM height is not Berggren ancestry |

Koukoulopoulos--Maynard's hard theorem controls correlations between infinitely
many denominator layers via second moments and GCD graphs.  Equations
`(14)--(16)` repackage the finite labelled input; they do not supply that
correlation estimate and do not reprove Duffin--Schaeffer.

## 7. Hostiles, reproduction, and scope

The exact companion includes the following mandatory controls:

- `U_1={0}`; omitting it loses the zero phase and breaks bijectivity;
- `q|L_Q` is not the same condition as `q<=Q`;
- the injection of clocks is `(18)`, not residue reduction;
- cyclic negation preserves exact denominator, reciprocal reversal need not;
- at denominators `2` and `4` with `psi=1/2`, formal layer mass is `1` but
  the cross-layer union has measure `3/4`;
- coprimality without opposite parity gives content two in `(22)`.

Run

```bash
python3 -B 04-computation/divisor_phase_duffin_schaeffer_clock_thm4056.py
python3 -B -O 04-computation/divisor_phase_duffin_schaeffer_clock_thm4056.py
```

Normal and optimized transcripts are byte-identical.  The audit includes
exact cyclotomic polynomial remainders for `(11)` rather than floating-point
Fourier comparisons.

This theorem proves a finite compiler, its symmetries, transforms, and exact
subpackets.  It proves no irrationality statement about Khinchin's constant,
`e+pi`, or any other constant; no Duffin--Schaeffer correlation estimate; no
union-coverage theorem; no classification of the Berggren tree; and no
LRC(14) result.

## 8. Complementary natural-edge, metric, and named-point layer

The cyclic exact-denominator compiler above and the two time sequences below
are compatible but not interchangeable. The first keeps labelled
exact-denominator packets inside one clock. The reduced-residue phase sequence
forgets prime-power exponents, while the natural-edge lift retains a
divisor-invertible scale coordinate. Their equality of complete-period means
therefore does not imply pointwise equality or equal overlap geometry.

**PROVED elementary identities and implications; CITED metric theorems;
VERIFIED-EXACT and INDEPENDENTLY HOSTILE-AUDITED.** No literature-priority or
global-novelty claim is made. The theorem gives an exact discrete model for
the first-moment quantity in Duffin--Schaeffer theory and proves several
pointwise consequences. It does not reprove the Duffin--Schaeffer theorem,
control interval overlaps, or settle the arithmetic nature of `e+pi` or
Khinchin's constant.

### 1. Inheritance, conventions, and closest hostiles

[THM-2000](THM-2000-support-harmonic-abel-dini-figurate-surface.md) separates
support from indexed multiplicity and forbids promoting fast convergence to
irrationality. [THM-3744](THM-3744-pell-prefix-loneliness-constant-carry-exact-formula.md)
shows that ordered continued-fraction recurrence and carry, rather than a
digit scalar, are load-bearing. The valid `e,pi` carrier in the mixed-status
[HYP-2212](../../05-knowledge/hypotheses/HYP-2212-quadratic-discriminant-carrier-pi-e.md)
is

```text
S=e+pi,        P=e*pi,        D=e-pi,        D^2=S^2-4P.       (1)
```

No two of `S,P,D` are algebraic, but their owner is not identified.

Use the Koukoulopoulos--Maynard normalization

```text
|x-a/q| <= psi(q)/q,       a in Z, gcd(a,q)=1,       q>=1,      (2)
```

and write `W'(psi)` for the subset of `R` satisfying (2) infinitely often.
This is the period-one extension of the standard `[0,1]` limsup set. For the
finite clocks below fix `Q>=2`, put

```text
L_Q=lcm(1,2,...,Q),
Psi_Q=sum_(q=2)^Q 2 phi(q) psi(q)/q.                         (3)
```

`Psi_Q` is raw interval length per unit interval, counted with multiplicity.
It is not the measure of the union when intervals overlap.

### 2. Natural edges and the exact reduction fibre

Let

```text
E_N={(a,b):1<=a<b<=N}
```

be the arcs of the natural-order transitive tournament, and reduce an arc by

```text
rho(a,b)=(a/g,b/g),       g=gcd(a,b).                         (4)
```

For every primitive `1<=p<q<=N`,

```text
rho^(-1)(p,q)={(kp,kq):1<=k<=floor(N/q)}.                    (5)
```

Thus the fibre has size `floor(N/q)`. There are `phi(q)` primitive types in
the denominator shell `q`, so exact double counting gives

```text
C(N,2)=sum_(q=2)^N phi(q) floor(N/q).                        (6)
```

Taking the first difference in `N`, a primitive type `p/q` acquires one new
lift exactly when `q|N`. Equivalently, the `N-1` new arcs ending at `N` split
as one copy of every reduced `p/q` with `q|N`, and

```text
N-1=sum_(q|N,q>=2) phi(q).                                  (7)
```

Equations (5)--(7) concern the complete edge set before reduction. Deleting
noncoprime pairs first leaves an incomplete coprimality graph, not a
tournament. Forgetting `k` in (5) loses scale multiplicity.

### 3. Two exact LCM clocks with one mean

#### 3.1 Reduced-residue phase clock

Let `U_q` be the reduced residue classes modulo `q`. Define the periodic
phase clock

```text
P_Q(t)
 =2 sum_(q=2)^Q psi(q) sum_(a in U_q) 1_[t=a mod q]
 =2 sum_(q=2)^Q psi(q) 1_[gcd(t,q)=1].                      (8)
```

The identity

```text
1_[gcd(t,q)=1]=sum_(d|gcd(t,q)) mu(d)
               =sum_(d|q) mu(d)1_[d|t]                    (9)
```

compiles (8) into the signed divisor clock

```text
P_Q(t)=sum_(d=1)^Q c_d 1_[d|t],
c_d=2 mu(d) sum_(2<=q<=Q,d|q) psi(q).                      (10)
```

Formula (8) proves `P_Q(t)>=0`, although the coefficients in (10) are signed.
The period divides `L_Q`; cancellation can make it smaller. Every unit class
modulo `q` occurs `L_Q/q` times in one complete period, hence

```text
(1/L_Q)sum_(t=1)^L_Q P_Q(t)=Psi_Q.                         (11)
```

The phase compiler is not injective on denominators. Both
`1_[gcd(t,q)=1]` and `phi(q)/q` depend only on `rad(q)`. Therefore it retains
only the radical-class aggregates

```text
A_r=sum_(rad(q)=r) psi(q).                                 (12)
```

Moving unit weight from `q=2` to `q=4` leaves (8), (10), and (11) unchanged.
Prime-power denominator depth is the missing sidecar.

#### 3.2 Scale-lift divisor clock

Give each lifted edge of primitive type `p/q`, `q<=Q`, the pulse amplitude
`2psi(q)`, and let `W_N` be the total amplitude in `E_N`. By (5),

```text
W_N=2 sum_(q=2)^Q psi(q)phi(q)floor(N/q).                   (13)
```

Adding vertex `t` gives the positive lift clock

```text
R_Q(t)=W_t-W_(t-1)
      =2 sum_(q=2)^Q psi(q)phi(q)1_[q|t].                  (14)
```

This clock is invertible as a divisor sum. If

```text
b(q)=2phi(q)psi(q),       b(1)=0,
```

then

```text
R_Q(n)=sum_(q|n)b(q),
b(n)=sum_(d|n)mu(n/d)R_Q(d).                              (15)
```

Its complete-period mean is again

```text
(1/L_Q)sum_(t=1)^L_Q R_Q(t)=Psi_Q.                        (16)
```

The two clocks are not pointwise equal. At the smallest hostile,
`Q=2,psi(2)=1`, their values on one period are

```text
P=(2,0),        R=(0,2).                                  (17)
```

The `q=2` versus `q=4` radical hostile is sharper: the phase clocks and raw
masses coincide, but the lift clocks on four ticks are respectively

```text
(0,2,0,2),      (0,0,0,4).                                (18)
```

For a fixed edge `a/q`, its phase ticks `a+jq` match its lifted copies
`((j+1)a,(j+1)q)` over a complete period. This bijection depends on the edge,
so it proves equality of first moments, not synchronized overlap data. For an
arbitrary prefix `T`, (13) gives only

```text
|W_T/T-Psi_Q| <= (1/T)sum_(q=2)^Q 2psi(q)phi(q).           (19)
```

### 4. Imported Duffin--Schaeffer theorem and the pointwise firewall

Koukoulopoulos--Maynard prove on `[0,1]` (and hence period-one on `R`) that,
for arbitrary nonnegative `psi`,

```text
sum_q phi(q)psi(q)/q=infinity  ==>  W'(psi) has full measure. (20)
```

This is **CITED**, not reproved by the clocks. The clocks recompile exactly
the partial sums in (20); they do not supply the overlap theorem.

The following pointwise statements are elementary and **PROVED**.

#### Pointwise theorem

1. If `psi(q)->0`, every member of `W'(psi)` is irrational.
2. For `psi_0(q)=1/q`,

   ```text
   W'(psi_0)=R\Q.                                          (21)
   ```

3. For

   ```text
   alpha=(sqrt(5)-1)/2,       psi_*(q)=1/(4q),             (22)
   ```

   the series in (20) diverges, so `W'(psi_*)` has full measure, but
   `alpha` is not in `W'(psi_*)`.

For a reduced rational `A/B` and any different reduced `a/q`,

```text
|A/B-a/q|=|Aq-aB|/(Bq)>=1/(Bq).                           (23)
```

When `psi(q)->0`, (23) excludes all but finitely many approximants; the exact
fraction itself occurs at only one reduced denominator. This proves part 1.
Every irrational has infinitely many reduced convergents with error
`<1/q^2`, while part 1 excludes rationals, proving (21).

For (22), let `beta=-(sqrt(5)+1)/2`. If `r=p/q` lies in `[0,1]`, then

```text
|(r-alpha)(r-beta)|=|p^2+pq-q^2|/q^2>=1/q^2,
|r-beta|<4,                                                (24)
```

so `|r-alpha|>1/(4q^2)`. Rationals outside `[0,1]` are farther than `1/4`.
On the other hand,

```text
sum_q phi(q)/(4q^2)=infinity,                              (25)
```

already by the prime terms. This explicit named quadratic irrational is the
hostile against reading a full-measure theorem pointwise.

### 5. The lawful Khinchin--Duffin--Schaeffer bridge

Let an irrational

```text
x=[a_0;a_1,a_2,...]
```

have a finite digit-geometric-mean limit

```text
(a_1...a_n)^(1/n) -> G < infinity.                         (26)
```

Then its irrationality exponent is exactly

```text
mu(x)=2.                                                    (27)
```

Indeed, with `S_n=sum_(j<=n)log(a_j)`, (26) gives

```text
log(a_(n+1))/n
 =((n+1)/n)(S_(n+1)/(n+1))-S_n/n ->0.                     (28)
```

The convergent denominators satisfy `q_n>=F_(n+1)`, hence
`log(a_(n+1))/log(q_n)->0`. The standard continued-fraction bounds, together
with Legendre's criterion, give

```text
mu(x)=2+limsup_n log(a_(n+1))/log(q_n),                    (29)
```

and prove (27). Thus, in the normalization (2),

```text
x in W'(q |-> 1/q),
x notin W'(q |-> q^(-1-epsilon))       for every epsilon>0. (30)
```

This is the valid infinite connection between a finite Khinchin limit and
metric approximation. A single finite digit mean has no such consequence:
any finite prefix admits both exponent-two and Liouville continuations.

The classical Khinchin theorem says that (26) holds with
`G=2.685452...` for almost every `x`; it does not describe that named
constant's own continued fraction. Euler's continued fraction supplies a
sharp atypical control:

```text
e=[2;1,2,1,1,4,1,1,6,1,...],
product_(first 3n positive digits)=2^n n!,                 (31)
```

so the geometric mean for `e` diverges. Nevertheless `a_n=O(n)` and
`q_n>=F_(n+1)` give `mu(e)=2` by (29); the converse to (26) is false.

### 6. Reciprocal reflection and its exact losses

Reflection sends a primitive edge `p->q` to `q->p`. It preserves
coprimality, the unordered endpoints, maximum height, and the scale-fibre
size

```text
floor(N/max(p,q)).                                         (32)
```

It does not preserve the distinguished denominator shell or a general
Duffin--Schaeffer decoration. For example, with `psi(n)=1/(n+1)`, the raw
lengths of `3->5` and `5->3` are `1/15` and `1/6`.

More generally, the raw edge length

```text
ell(p->q)=2psi(q)/q                                        (33)
```

is reversal-invariant for every primitive pair if and only if

```text
psi(n)=c n                                                 (34)
```

for one constant `c`: use the primitive edge `1<->n` to force
`psi(n)/n=psi(1)`, and the converse is immediate. If `psi(n)->0`, or if
`psi(n)<=1/2` for every `n`, (34) forces `c=0`. Thus no nontrivial standard
decaying or bounded approximation function is globally reciprocal-even.

There is also a continued-fraction normalization trap:

```text
3/5=[0;1,1,2],       5/3=[1;1,2].                         (35)
```

Deleting only a leading zero gives the same **projective Euclidean
coefficient word** `(1,1,2)`. The standard finite Khinchin word always
deletes the integer part `a_0`, giving `(1,1,2)` versus `(1,2)` in (35).
Therefore the standard finite Khinchin mean is not exactly reciprocal-even.
For an infinite limit, deleting finitely many initial digits does preserve
the limit and irrationality exponent.

### 7. Exact `e+pi` proof gates

The metric theorem cannot name `e+pi`. Two elementary reductions identify
what a genuine proof would still need.

#### Trace-only recurrence gate

Every symmetric polynomial in variables `x,y` is uniquely `F(S,P)` with
`S=x+y,P=xy`. If `S=e+pi` were algebraic, then `P=e*pi` would be
transcendental: otherwise `e,pi` would both be roots of an algebraic
quadratic. Consequently every specialized coefficient `F(S,P)` which remains
nonconstant in `P` is transcendental. It cannot be denominator-cleared from
the hypothesis `S in Qbar` alone.

Thus a symmetric recurrent-integral attack must satisfy at least one exact
gate:

```text
(i) all coefficients and initial denominators descend to Q(S);
(ii) every P-dependence cancels after specialization;
(iii) a separate P- or D-sidecar supplies denominator control.              (36)
```

For example `e^2+pi^2=S^2-2P`, so merely symmetrizing independent recurrences
usually fails (36).

#### Synchronized-linear-form criterion

Suppose integers `U_n,V_n,R_n,W_n` give

```text
L_n=U_n+V_n e,        M_n=R_n+W_n pi.
```

Then

```text
N_n=W_nL_n+V_nM_n
   =(W_nU_n+V_nR_n)+V_nW_n(e+pi).                         (37)
```

If `N_n` is nonzero infinitely often and tends to zero, then `e+pi` is
irrational: under `e+pi=A/B`, every nonzero `N_n` lies in `(1/B)Z`. A useful
sufficient rate packet is

```text
|L_n|<=exp(-an), |V_n|<=exp(bn),
|M_n|<=exp(-cn), |W_n|<=exp(dn),       a>d, c>b,           (38)
```

plus nonvanishing. Existing repo work does not provide this synchronized
cross-decay/nonzero packet.

### 8. Validity boundary and generated frontiers

- The Koukoulopoulos--Maynard theorem and classical Khinchin theorem are
  **CITED**. The clock, pointwise firewall, exponent implication, reciprocal
  rigidity, and proof gates are **PROVED** here.
- `e+pi` irrationality remains **OPEN**. The criterion (37) is a target
  format, not a constructed sequence.
- The arithmetic nature of Khinchin's constant is not settled by any cited
  result here. Typical digit means concern other real numbers.
- Equality of clock means loses cross-denominator overlap synchronization;
  the radical-class collision (18) is the minimal witness.
- Promising next work is a trace-purity audit of every repo irrationality
  recurrence, a paired growth/decay/nonvanishing ledger for (37), and a
  congruence-restricted version of the phase clock retaining prime-power
  depth.

### 9. Literature pins

- Koukoulopoulos--Maynard,
  [*On the Duffin--Schaeffer conjecture*](https://arxiv.org/abs/1907.04593),
  proves (20).
- Koukoulopoulos--Maynard--Yang,
  [*An almost sharp quantitative version...*](https://arxiv.org/abs/2404.14628),
  uses the finite mass `Psi_Q` in its quantitative counting theorem.
- Cellarosi--Hensley--Miller--Wellens,
  [arXiv:1402.0208](https://arxiv.org/abs/1402.0208), records the classical
  Khinchin digit-mean theorem and the atypical Euler expansion.
- Sondow,
  [*Irrationality Measures, Irrationality Bases, and a Theorem of Jarnik*](https://arxiv.org/abs/math/0406300),
  records formula (29); the short proof above also identifies its mechanism.

### 10. Replay

From the repository root:

```text
python -B 04-computation/divisor_phase_duffin_schaeffer_thm4056.py
python -B -O 04-computation/divisor_phase_duffin_schaeffer_thm4056.py
python -B 04-computation/divisor_phase_duffin_schaeffer_thm4056_independent_audit.py
python -B -O 04-computation/divisor_phase_duffin_schaeffer_thm4056_independent_audit.py
```

Both normal/optimized pairs are byte-identical. The independent path uses a
totient sieve, direct clock arrays, recursive divisor inversion, and the
finite-Khinchin normalization hostile. **QED.**

