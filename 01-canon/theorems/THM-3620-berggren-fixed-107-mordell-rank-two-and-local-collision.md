---
id: THM-3620
title: "Berggren fixed-107 rank-two 3-saturated subgroup and local denominator collisions"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE / PENDING INDEPENDENT AUDIT.
  On E:Y^2=X^3+107^3-2, the two Berggren/two-cube points
  P=(232,3703) and Q=(4960,349321) generate a rank-two, 3-saturated
  subgroup of E(Q).  Their difference and sum exhibit exact denominator
  collisions at 3 and 197.  This proves neither rank equality nor a complete
  integral-point classification.
source: kps-s189 / Berggren fixed-summand Mordell wildcard, 2026-08-21
audit: PENDING -- exact companion and proof require independent hostile audit
depends_on:
  - THM-3370-berggren-two-cube-biquadratic-norm-collision
related:
  - THM-3375-berggren-positive-two-cube-pell-ray
  - THM-3580-berggren-positive-cube-slope-atlas-completion-through-101
  - THM-3594-berggren-positive-cube-slope-atlas-through-201
script: 04-computation/berggren_fixed_107_mordell_rank_two_thm3620.py
output: 05-knowledge/results/berggren_fixed_107_mordell_rank_two_thm3620.out
script_sha256: d04425f7b9c0590f6e456476eb76cdd846286750b82602fdf2ba806fe617fc80
output_sha256: 74e21293e69d3d3058810011dd1d2f6a7b556bb2b0321a4d96800fb6e6458396
hash_basis: LF-normalized bytes
---

# THM-3620 -- fixed `107` supports a rank-two 3-saturated subgroup

**RESERVED / PROVISIONAL PROOF CANDIDATE / PENDING INDEPENDENT AUDIT.**

THM-3370 found two positive Berggren/two-cube intersections sharing the
summand `107`.  This theorem studies that fixed-summand fibre as one elliptic
curve.  Its new content is a rank lower bound and a saturated subgroup, not a
rediscovery of the previously audited points.

## 1. Statement and translation

Put

```text
B=107^3-2=1,225,041=3*408,347,
E: Y^2=X^3+B,
P=(232,3,703),                  Q=(4,960,349,321).       (1)
```

Then `408,347` is prime, `P,Q in E(Q)`, and

```text
H=<P,Q> is isomorphic to Z^2 and is 3-saturated in E(Q). (2)
```

Here 3-saturated means

```text
3T in H, T in E(Q)                 ==> T in H.           (3)
```

The Berggren translation is exact.  A fixed-summand representation

```text
107^3+X^3=(2r+1)^2+2                                  (4)
```

is an integral point `(X,Y)=(X,2r+1)` on `E`.  Thus `P,Q` encode precisely
the already known depths

```text
r=1,851,                         r=174,660.              (5)
```

No assertion is made that `(5)` exhausts this Mordell fibre.

## 2. Torsion is trivial

The displayed integral model has discriminant

```text
Delta=-432B^2=-2^4*3^5*408,347^2.                       (6)
```

Consequently `7,13,19` are primes of good reduction.  Exact enumeration gives

```text
#E(F_7)=4,                #E(F_13)=12,              #E(F_19)=27. (7)
```

At a good prime greater than `2`, reduction on rational torsion is injective.
In particular its order divides both `4` and `27`; hence

```text
E(Q)_tors=0.                                             (8)
```

The coprime controls in `(7)` are load-bearing.  A single reduction count
would not prove `(8)`.

## 3. Two independent classes modulo `3E(Q)`

For a finite group `G`, write `3G={3R:R in G}`.  Direct group-law enumeration
on the two useful reductions gives

```text
3E(F_13)={O,(1,0),(3,0),(9,0)},
aP+bQ in 3E(F_13)       iff a+b=0 mod 3;                (9)

3E(F_19)={O,(0,4),(0,15)},
aP+bQ in 3E(F_19)       iff a-b=0 mod 3.               (10)
```

The two kernels in `(Z/3Z)^2` meet only at `(0,0)`.  If
`aP+bQ in 3E(Q)`, reduction at both good primes therefore forces

```text
a=b=0 mod 3.                                             (11)
```

Since `(8)` removes the torsion contribution,

```text
E(Q)/3E(Q) is isomorphic to (Z/3Z)^rank(E(Q)).           (12)
```

Equations `(9)--(11)` exhibit two independent classes, proving
`rank E(Q)>=2` and `H isomorphic to Z^2`.

They also prove saturation.  If `3T=aP+bQ`, `(11)` writes
`a=3a_0,b=3b_0`, so

```text
3(T-a_0P-b_0Q)=O.
```

Equation `(8)` then gives `T=a_0P+b_0Q`, which is `(3)`.

Each prime separately leaves a nonzero coefficient line: at `13`,
`(1,2),(2,1)` survive, while at `19`, `(1,1),(2,2)` survive.  Thus the
two-prime intersection is essential rather than cosmetic redundancy.

## 4. Local collisions become rational denominators

The exact rational group law gives

```text
Q-P=(3448/9,-204659/27),                                (13)
P+Q=(94164361/788^2,1062189113125/788^3),
788=4*197.                                               (14)
```

The corresponding chord factorizations are

```text
(224X-63077)^2-9(X^3+B)
  =-(X-232)(X-4960)(9X-3448),                           (15)

(57603X-10445932)^2-788^2(X^3+B)
  =-(X-232)(X-4960)(620944X-94164361).                  (16)
```

They record two distinct local coordinate collisions:

```text
P=(1,1)=Q              on the singular cubic modulo 3,  (17)
P=(35,157)=-Q                in E(F_197).                (18)
```

Thus subtracting the points exposes powers of `3` in `(13)`, while adding
them exposes powers of `197` in `(14)`.  These are exact examples of a useful
denominator heuristic; they are not asserted to classify all denominators in
`H`.  Since `3` divides the discriminant, `(17)` is deliberately not written
as an equality in an elliptic-curve group and plays no role in the rank or
saturation argument.

## 5. Fixed-fibre congruence support

For any integral point in the Berggren form `Y=2r+1`, reduction modulo `7`
forces

```text
X^3=1, Y=0 mod 7,                    hence r=3 mod 7.     (19)
```

Modulo `9` it forces `X^3=1,Y=+-4`, hence

```text
r=2 or 6 mod 9.                                          (20)
```

The CRT gives the sharper fixed-`107` support

```text
r=24 or 38 mod 63.                                       (21)
```

Both points `(5)` lie in class `24`.  This is a necessary local sieve only.
It neither proves that class `38` occurs nor that class `24` contains only
the two known points.

## 6. Exact-computation boundary

The companion independently verifies `(1),(6)--(21)` using only integer and
rational arithmetic.  It also scans the explicit finite interval

```text
-107 <= X <= 10,000,000                                 (22)
```

and finds only `P,Q` with positive odd ordinate.  This is finite evidence,
not a height theorem and not a complete computation of `E(Z)`.  In
particular this theorem proves none of:

- `rank E(Q)=2`;
- `H=E(Q)`;
- a complete integral-point list;
- finiteness or infinitude of fixed-`107` Berggren depths; or
- a new global census replacing THM-3370, THM-3580, or THM-3594.

Reproduce the exact certificate with

```bash
python3 04-computation/berggren_fixed_107_mordell_rank_two_thm3620.py
python3 -O 04-computation/berggren_fixed_107_mordell_rank_two_thm3620.py
```

Both streams must agree with one another and with the stored output after LF
normalization.  The companion uses 271 explicit optimization-safe gates and
pins THM-3370 by LF hash
`9abf73a45d789fe5a804977235a3bce3e415c3ad6da865a628d8f338159fa53a`.
