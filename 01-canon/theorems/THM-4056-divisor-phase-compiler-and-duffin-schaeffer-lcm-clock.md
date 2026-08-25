---
id: THM-4056
title: "Divisor phase compiler and Duffin--Schaeffer LCM clock"
status: >
  PROVED + VERIFIED-EXACT.  Every cyclic N-clock is the disjoint union of
  primitive rational-centre packets indexed by d|N; the explicit bijection,
  inverse, unit-orbit structure, arbitrary weighted identity, and Ramanujan
  transform are exact.  For N=L_Q this losslessly compiles every
  Duffin--Schaeffer denominator q<=Q after retaining the exact-denominator
  filter, while the full clock compiles all q|L_Q.  Clock growth occurs only
  at prime powers.  The opposite-parity packet of every even clock is an
  exact primitive-Pythagorean/Berggren slice.  This is a finite carrier and
  first-moment/Fourier theorem; it does not prove Duffin--Schaeffer,
  irrationality, interval coverage, or LRC(14).
source: codex-khinchin-ds-rational-tournament-20260824
audit: >
  The companion checks the bijection and inverse, gcd strata, unit orbits,
  rational weights, cyclotomic Ramanujan transforms, prefix filters,
  prime-power jumps, reflections, Pythagorean content, and hostile overlap
  controls in normal and optimized modes.  An independent hostile audit
  identified and verified the q=1 convention, prefix/divisor distinction,
  multiplicity/coverage boundary, and direct-system map.
related:
  - HYP-2628-lrc14-euler-copy-squarefree-profile
  - THM-2041-frobenius-stability-of-exact-period-projectors
  - THM-873-ramanujan-fourier-expansion-of-interval-core-good-sets
  - THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone
  - THM-4042-prime-sector-ap-cover-exact-clock-and-holonomic-law
  - THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge
script: 04-computation/divisor_phase_duffin_schaeffer_clock_thm4056.py
output: 05-knowledge/results/divisor_phase_duffin_schaeffer_clock_thm4056.out
script_sha256: a15fb134839f84b5e1a4f07131f893af790adcc2c6b8189780c851a9e90d77e9
output_sha256: bce0b8ccc22db995a0d96560aa66c7707b581b45e419a02589920a4249b28947
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
