---
id: THM-2043
title: Period-14 parity-Hasse jets are complete local coordinates but are globally magnitude-blind
status: >
  PROVED / exact algebraic theorem and sharp LRC carrier no-go; not an LRC(14)
  proof.  The characteristic-seven parity-Hasse map is an isomorphism on
  reduced period-14 functions.  Nevertheless an infinite family of strict
  17/41 witnesses has exactly the tight AP's complete period-14 owner-residue
  data and q-threshold.  No fixed finite 7-adic lift-height truncation repairs
  this.  An owner-labelled mod-13 x mod-14 chart restores every positive
  height at most 181; universally one needs exact height or a resolved
  cross-denominator certificate.
source: codex-2026-07-21-NC2-transfer
related:
  - THM-2041
  - THM-2000
  - THM-671
  - HYP-8800
  - HYP-2963
  - HYP-2979
  - HYP-3036
  - THM-404
script: 04-computation/lrc14_parity_hasse_jet_packet_audit_codex_20260721.py
output: 05-knowledge/results/lrc14_parity_hasse_jet_packet_audit_codex_20260721.out
---

# THM-2043 -- Period-14 parity-Hasse completeness and its exact limit

This theorem resolves HYP-8800's first characteristic-seven question.  Hasse
jets are the correct complete coordinates for one period-14 function after
reduction modulo seven.  They do **not** repair information lost before that
reduction, and in particular they cannot see how high a speed lies above its
residue class.  The smallest valid LRC carrier therefore needs a height or
multi-modulus sidecar.

No claim here proves LRC(14).

## 1. Complete parity-Hasse coordinates

Let

```text
A = F_7[X]/(X^14-1).
```

In characteristic seven,

```text
X^14-1 = (X^2-1)^7 = (X-1)^7 (X+1)^7.                 (1)
```

The two factors on the right are coprime.  Hence the Chinese remainder
theorem gives

```text
A ~= F_7[X]/(X-1)^7 x F_7[X]/(X+1)^7.                 (2)
```

For `F(X)=sum_(a=0)^13 f(a)X^a`, write `D^(j)` for the `j`th Hasse
derivative.  The map

```text
J_14(f) = (D^(j)F(1), D^(j)F(-1))_(0<=j<=6)           (3)
```

is an `F_7`-linear isomorphism from the fourteen-dimensional space of
functions on `Z/14Z` to `F_7^14`.

### Proof

Taylor expansion with Hasse derivatives identifies the first factor in (2)
with the seven values `D^(j)F(1)` and the second with the seven values
`D^(j)F(-1)`.  Thus (2) is exactly (3).  Equivalently, the explicit `14 x 14`
matrix

```text
J_(epsilon,j),a = binomial(a,j) epsilon^(a-j),
epsilon in {1,-1}, 0<=j<=6, 0<=a<=13,                 (4)
```

has rank fourteen over `F_7`.  The companion audit performs both the CRT
factorization and exact row reduction.

The word "complete" is intentionally local: (3) recovers the coefficient
vector `f mod 7`.  It does not recover an integral lift of `f`, nor data that
was already forgotten when a speed set was pushed to `Z/14Z`.

## 2. The infinite magnitude-blind collision

For a thirteen-speed set `S` and `a in Z/14Z`, define three exact integer
phase functions:

```text
N_S(a) = #{v in S : av = 0 mod 14},
W_S(a) = 1[min(av mod 14,14-av mod 14) >= 1 for every v in S],
B_S(a) = #{v in S : min(av mod 14,14-av mod 14) = 1}.  (5)
```

Let

```text
S_0 = {1,2,...,13},
S_h = {1,2,...,11,13,12+14h},       h>=1.             (6)
```

Then, for every `h`, the complete owner-indexed residue data of `S_h` modulo
fourteen equals that of `S_0`: the changed runner remains congruent to `12`.
Consequently

```text
(N_(S_h),W_(S_h),B_(S_h)) = (N_(S_0),W_(S_0),B_(S_0))                 (7)
```

over the integers, and all forty-two parity-Hasse coordinates of these three
functions agree after reduction modulo seven.

Nevertheless the LRC behavior differs sharply.

```text
M(S_0) = 1/14.                                                        (8)
```

Indeed `t=1/14` gives the lower bound.  For arbitrary `t`, two of the fourteen
points `0,t,...,13t` on the circle have distance at most `1/14`; their
difference is `kt` for some `1<=k<=13`, giving the upper bound.

If `6` does not divide `h`, then `t=1/12` is a strict witness for `S_h`:

```text
(12+14h) mod 12 = 2h mod 12 != 0,                                    (9)
```

and none of `1,...,11,13` is divisible by twelve.  Every runner is therefore
at distance at least `1/12 > 1/14` from an integer.

This proves an infinite collision between a boundary-tight row and strict
witness rows inside one complete period-14 Hasse packet.  There is a stronger
collision that also defeats the usual `q_threshold` sidecar.  Put

```text
T_n = {1,2,...,11,13,96+3444n},       n>=0.                          (10)
```

Here `3444=lcm(12,14,41)`.  The replacement is congruent to `12` modulo
fourteen, so the full period-14 owner-residue and Hasse packets still equal
the AP packet.  For each `q=2,...,11,13`, the retained speed `q` blocks the
phase `1/q`; the replacement is divisible by twelve; and no speed is
divisible by fourteen.  Thus `T_n` and `S_0` have the same blockedness ledger
through thirteen and the same `q_threshold=14`.

On the other hand, at `t=17/41` the residue distances of the retained speeds
`1,...,11,13` are

```text
(17,7,10,14,3,20,4,13,11,6,18,16),
```

while `96+3444n = 14 mod 41` and its distance is eight.  Therefore

```text
min_(v in T_n) ||17v/41|| = 3/41 > 1/14,
14*3-41 = 1.                                                        (11)
```

The exact `+1` is a genuine strict LRC exit, not another stable shadow.

## 3. Single-period no-go

Any proposed certificate that factors only through the owner-indexed residues
modulo fourteen has the same value on `S_0` and every `S_h`.  This includes,
without owner heights:

- every unlabelled period-14 phase statistic;
- the complete integer functions in (5);
- their Fourier, Ramanujan, or characteristic-seven Hasse transforms;
- any nonlinear function of the complete Hasse packet.

Therefore no such certificate can by itself distinguish boundary tightness
from strict loneliness, classify the global proof route, or supply the missing
LRC seed-to-exit implication.  Deeper Hasse order cannot recover information
absent from the input function.

The finite audit makes the obstruction visible beyond one example.  Replace
each AP speed `d` in turn by `d+14h`, with `1<=d<=13` and `1<=h<=12`.  All 156
rows have the same complete period-14 owner-residue packet.  Their first
denominator `q<=14` with no zero runner has exact census

```text
q:       8   9  10  11  12  13  14
count:   9  11  10  11  10  12  93.                                 (12)
```

Thus 63 of the 156 rows already have a direct strict witness with `q<=13`.
An exact primitive-phase scan tests the integer slack

```text
min_(v in S) (14 min(av mod q,q-av mod q)-q) > 0                       (13)
```

for every coprime numerator `a` and `q<=42`.  It finds a strict certificate
for all 156 aliases, with first strict-denominator histogram

```text
{8:9, 9:11, 10:10, 11:11, 12:10, 13:12, 15:10, 16:1,
 17:21, 19:26, 21:2, 22:10, 23:15, 25:4, 27:2, 41:2}.                (14)
```

This is an exact bounded audit, not a claim about the full HYP-2963 bank.  It
also shows why `q_threshold` alone is not a complete repair: 93 strict aliases
remain on its `q_threshold=14` fiber.  The resolved `(q,a,positive slack)`
ledger separates this finite alias family.

The eleven named HYP-2979 controls give an independent route check.  Their
three exact functions in (5) form only three fibers.  One seven-row fiber is

```text
AP, GW 12->24, residue liar 12->26, near 12->36,
petal 10->20, P10+GW, P10+K33.                                     (15)
```

It mixes two AP/GW labels, one direct `q` witness, two petal labels, and two
K33/state-lift labels.  Every Hasse truncation from depth zero through depth
six has the same three-fiber census and the same mixed seven-row fiber.  The
route labels here are the previously audited HYP-2979 labels; this computation
does not assert the still-open global HYP-2963 classification theorem.

## 4. Fixed finite 7-adic height also fails

The obstruction is not removed by attaching any fixed number of base-seven
digits of lift height.  Fix `k>=0`, choose a positive integer `u_k` with

```text
7^k u_k = 1 mod 41,
```

and set, for `m>=0`,

```text
x_(k,m) = 12 + 84*7^k*(u_k+41m).                                   (16)
```

Then simultaneously

```text
x_(k,m) = 12 mod 14,
12 divides x_(k,m),
x_(k,m) = 14 mod 41,
(x_(k,m)-12)/14 = 0 mod 7^k.                                      (17)
```

Replacing `12` in the AP by `x_(k,m)` therefore preserves its period-14
packet, its blockedness ledger through thirteen, its `q_threshold`, and its
first `k` lift-height digits, but the strict certificate (11) still has margin
one.  For `k=0,1` this gives replacements `96,3540`; at `k=1` the lift height
is `252=0 mod 7`, exactly the AP height digit.

Hence neither raw owner labels, `q_threshold`, nor any fixed finite 7-adic
height truncation is a faithful global carrier for the unbounded speed
problem.

## 5. The constructive repair: a multi-modulus owner chart

The incoming THM-2000 support-measure viewpoint makes the defect transparent.
Pushing the owner-labelled atomic support along

```text
Z -> Z/14Z                                                        (18)
```

forgets the kernel direction `v -> v+14`.  The Hasse isomorphism (3) changes
coordinates only after this pushforward and therefore cannot undo it.

Add the owner-labelled residue modulo thirteen.  Since `gcd(13,14)=1`,

```text
v |-> (v mod 13, v mod 14)                                          (19)
```

is injective modulo `182`.  In particular it recovers every positive speed
`1<=v<=181` exactly.  The companion audit checks all 181 heights with zero
collisions.  For the sharp pair `S_0,S_1` from (6), even the period-12 phase
section already distinguishes the rows and exhibits the strict witness.

For an unbounded family such as (10), the clean positive sidecar is instead
the resolved phase certificate

```text
C_(q,a)(S) = min_(v in S) (14 min(av mod q,q-av mod q)-q).           (20)
```

It is an exact exit whenever it is positive; here `C_(41,17)(T_n)=1` for
every `n`.  Exact owner height is another faithful, deliberately maximal,
repair.

This is the positive conclusion of the theorem:

> On bounded banks, parity-Hasse jets can be used as a complete local chart,
> but they must be glued to an owner-labelled coprime-modulus or exact-height
> chart.  On unbounded families the modulus must grow adaptively, a resolved
> positive phase certificate must be produced, or exact height must remain
> attached.

The word **owner-labelled** is load-bearing.  Two unlabelled residue multisets
at different moduli do not specify how residues are paired runner by runner.

## 6. A characteristic-three scalarization guardrail

Let `w` be the weak-safe mask of the AP at period fourteen.  Direct exact
projection in `F_3[C_14]` gives

```text
e_14 w = (0,2,1,2,1,2,1,0,1,2,1,2,1,2) != 0,                       (21)
E_14(w) = 0 mod 3.                                                   (22)
```

Thus THM-2041 can preserve the nonzero projected vector even when the scalar
Ramanujan energy vanishes: over a finite field the quadratic pairing may be
isotropic.  The characteristic-three program must propagate a vector or
owner-current, not merely a scalar energy.

## 7. General bad-prime cyclic template

The local algebra is not special to `14`.  Let `N=M p^a` with `p` prime and
`p` not dividing `M`.  In characteristic `p`,

```text
F_p[C_N]
 = F_p[X]/(X^N-1)
 = F_p[X]/((X^M-1)^(p^a)).                                         (23)
```

The polynomial `X^M-1` is square-free because its derivative is
`M X^(M-1)`.  Therefore, over a splitting field,

```text
k[C_N] ~= product_(zeta^M=1) k[X]/(X-zeta)^(p^a),                   (24)
```

and the Hasse values `D^(j)F(zeta)`, `0<=j<p^a`, are complete local
coordinates.  Over `F_p` itself, group the roots into irreducible Frobenius
orbits.  Equations (1)--(3) are the case `(p,a,M)=(7,1,2)`.

This gives a reusable bad-prime packet for repeated-root cyclic codes,
periodic signals, and group-algebra transport.  THM-2043's LRC collision is
also the universal guardrail: these jets recover the reduced cyclic function,
not the pre-pushforward object that produced it.

## 8. Consequence for the LRC program

The raw characteristic-seven target in HYP-8800 is now decided negatively.
The surviving target is narrower:

```text
endpoint-owner Hasse module
  + exact height or adaptive resolved-denominator chart
  + a sign/duality implication
  -> safe phase or a named AP/GW, petal, K33, or covering exit.        (25)
```

Equation (25) still lacks its general final implication.  THM-2043 proves which
coordinates are faithful and which are impossible; it does not create the
THM-671 `B5` seed or the Fejer/Toeplitz exit.

## 9. Tournament Analysis

The audit uses proof carriers rather than runners as vertices.  The pairwise
observable is the lexicographic retention vector

```text
(strict-witness status, endpoint labels, magnitude, full phase, Hasse depth).
```

The carrier tournament is transitive, with score histogram
`{0:1,1:1,...,7:1}`, zero directed three-cycles, eight singleton SCCs, and one
Hamiltonian path.  This preserves the hierarchy of proof obligations and
destroys raw circle geometry.  The challenged assumption is exactly the one
refuted by (6)--(17): increasing jet depth or retaining finitely many lift
digits cannot recover the exact height that the period-14 pushforward erased.
