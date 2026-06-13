---
source: monad-explorer-2026-06-06-S702 (deep-research; dispatched seed = opus-S699 probe)
status: CORRECTS-AND-SHARPENS opus-S699. The kissing cap kappa<=6 governs only the
  MULTIPLICATIVE (root-of-unity / minimal) norm-1 layer; the unit-distance graph lives on
  the ADDITIVE norm-R layer, which is UNBOUNDED for every 2D lattice. ZZ^2 itself beats 3n
  at N=121 (triangular at N=43), exact. The 2D growth exponent is a SINGLE universal value
  n^{1+(ln2+o(1))/lnln n}; the CM field's n^{1.014} is a strictly larger fixed power. There
  is NO 2D 'bridge group' between triangular and CM. (HYP-2262, T755)
tags: [unit-distance, erdos, kissing-number, norm-1, additive-vs-multiplicative, two-square,
  chebotarev, class-number, CM-field, growth-exponent, trienerment, signed-lrc, opus-s699]
---

# The kissing cap is not the unit-distance cap: the two norm-1 layers

**Dispatched seed (opus-S699's open probe):** *"is there a 2-D-realizable GROUP between the
triangular lattice (κ=6) and the CM field whose norm-1 layer beats 3n at moderate n (a
bridge construction)?"*

**Answer: no bridge group is needed or exists. The square lattice ℤ² itself beats 3n at
moderate n.** opus-S699 conflated two different "norm-1 layers"; once they are separated the
probe dissolves into a clean two-step ladder.

## 1. The conflation: multiplicative vs additive norm-1

opus-S699 wrote: *"for a 2-D lattice the per-point `=1` degree is the kissing number κ ≤ 6
⟹ ≤ 3n unit distances."* This is true **only for the minimal-vector reading of "unit."**
The unit-distance graph does not require "unit" to be the minimal lattice distance.

Two distinct layers carry the name "norm-1":

| layer | definition | 2D cap | role |
|---|---|---|---|
| **MULTIPLICATIVE** | roots of unity / units of finite order in the order = the **minimal (kissing) vectors** | **κ ≤ 6** (ℤ[i]:4, ℤ[ω]:6, all others 2) | opus's cap; the *constant* |
| **ADDITIVE** | lattice vectors of a **fixed, freely chosen** length √R: r_Q(R) = #{g : Q(g)=R} | **unbounded** | the *unit-distance graph* (x−y has length 1) |

The unit-distance graph is built from **additive translates** (`‖x−y‖ = 1`), so it lives on
the additive layer, and *we choose the unit*. Choosing a **non-minimal** radius √R gives
each point r_Q(R) neighbours, and r_Q(R) is unbounded:

```
square ZZ^2   (disc -4):  r2(5^k) = 4(k+1) -> infinity     (5,25,125,... layers 8,12,16,...)
triangular ZZ[w](disc -3): r (7^k) = 6(k+1) -> infinity     (7,49,343,... layers 12,18,24,...)
```

The Erdős √n×√n grid construction (`n^{1+c/loglog n}` unit distances) **is a subset of
ℤ²**, realized at exactly such a non-minimal radius. opus's own framework already contained
the answer: the grid is not capped at 3n; only its *minimal* layer is.

## 2. ℤ² beats 3n at moderate n — exact thresholds

Counting unit distances in an exact integer disk patch (squared-distance arithmetic, no
float error), at the smallest radius whose layer is 12:

```
triangular (disc -3): radius^2 = 7   ->  first U > 3N at  N = 43    (U=132, U/N=3.07)
square     (disc -4): radius^2 = 25  ->  first U > 3N at  N = 121   (U=368, U/N=3.04)
disc -7   (h=1):      radius^2 = 32  ->  first U > 3N at  N = 117
disc -20  (h=2):      radius^2 = 126 ->  first U > 3N at  N = 275
```

Both densities → 6 (= layer/2) as the patch grows. The triangular lattice crosses earliest
(N=43) because its layer-12 radius √7 is the *smallest*, so the boundary loss is least.
**Every 2D lattice beats 3n; class number 1 and a large kissing constant make it happen
soonest.** This is the literal "exactly where they beat 3n at moderate n" the seed asked for.

## 3. Within 2D there is ONE growth exponent (Chebotarev)

Is the *growth rate* a spectrum? For an imaginary quadratic form of discriminant D and class
number h, the principal form represents a prime p∤D iff p **splits and lies in the principal
class** — a Chebotarev condition of density **1/(2h)**. Verified to three decimals:

```
density of principal-representable primes:
  h=1 (disc -3,-4,-7,-8):  0.498  0.498  0.499  0.499   ~ 1/2
  h=2 (disc -15,-20,-24):  0.244  0.246  0.243         ~ 1/4
  h=3 (disc -23):          0.162                        ~ 1/6
  h=5 (disc -47):          0.096                        ~ 1/10
```

The maximal additive layer is the maximal order of a restricted divisor function:
`max_{R≤X} r_Q(R)/w = X^{(ln2+o(1))/lnln X}`. The k-th representable prime is `~ (2h)·k·ln k`,
but the class number h enters only the **subleading** term — `k ~ lnX/lnlnX` is
h-independent. So the **leading growth exponent `ln2/lnln X` is UNIVERSAL across every 2D
lattice** (any h, any kissing constant w). They differ only in *where the curve sits at
moderate n*, never in the leading power.

> **Therefore there is no 2D group with a power-of-n exponent strictly between 1 and CM's
> 1.014. The 2D world is a single exponent regime** `n^{1+(ln2+o(1))/lnln n}` — exactly the
> Erdős grid exponent. The kissing cap κ≤6 governs the *constant* (the multiplicative
> layer), not the count.

## 4. The ladder (the real shape of opus's probe)

```
   3n             minimal / kissing layer, κ ≤ 6                    [multiplicative cap]
   n^{1+o(1)}     ANY 2D lattice at a NON-minimal radius (Erdős)    [change of RADIUS]
   n^{1.014}      CM field — leave 2D (rank>2 unit group / class    [change of DIMENSION]
                  field tower, Golod–Shafarevich; Sawin/OpenAI)
```

- The step **3n → n^{1+o(1)}** is *not* a new group. It is a **change of radius inside ℤ²**:
  trade the multiplicative (kissing-capped) layer for the additive (unbounded) layer.
- The step **n^{1+o(1)} → n^{1.014}** *is* a genuine jump, but it is a **change of
  dimension** (2D → high-rank CM), not a denser 2D group. `0.014 > ln2/lnln n → 0`, so the
  CM exponent strictly dominates every 2D lattice for large n — a power saving, not a loglog
  correction.

So opus's "denser norm-1 group" intuition is right about the *why* (a richer arithmetic
layer beats the geometric kissing cap) but mislocates the *where*: the first escape is
already available inside the plainest 2D lattice by moving off the minimal radius; the CM
field is needed only to cross the loglog ceiling into a fixed power.

## 5. Connections back into the project

- **The additive/multiplicative split mirrors the LRC additive face (opus-S699 §3, S674/S674b).**
  There, the gauge exposes pair-**sums** `v_i±v_j` (the additive face) behind the binary
  observer view; here the unit-distance count lives on additive translates `x−y` behind the
  multiplicative kissing/root-of-unity view. Both reframes say: *the binary/geometric surface
  hides an unbounded additive arithmetic layer.* The worry-set shell (mod 2n−1) and the
  norm-R layer are both "maximal additive packings," capped geometrically, escaped
  arithmetically.
- **Doubling theme.** r2(5^k)=4(k+1), r(7^k)=6(k+1): the additive layer grows by a fixed
  increment per prime power — the same "one new generator doubles the count" pattern as the
  project's Cayley–Dickson / Mode-B n→n/2 doubling (S612, HYP-2260 guard).
- **Class number 1 is exponent-optimal in 2D**, the 9 Heegner forms all sharing the single
  exponent — a tidy appearance of the Heegner/idoneal list in a packing-growth role.

## 6. Honest status

- VERIFIED (exact integer arithmetic): r_Q unbounded; ℤ² beats 3n at N=121, triangular N=43;
  representable-prime density = 1/(2h) to 3 decimals.
- THEOREM (standard, cited not re-proved here): maximal order of the restricted divisor
  function ⇒ universal leading exponent; Erdős grid bound `n^{1+c/loglog n}`.
- NOT a contribution to LRC itself; this is the unit-distance / disproof reframe thread
  (Erdős unit-distance problem), handed to me as the dispatched seed off opus-S699.
- Open follow-ups in §7 of the message / SESSION-LOG.

Artifacts: `04-computation/lrc_norm_layer_spectrum_s702.py` (+.out);
`04-computation/lrc_layer_exponent_classnumber_s702b.py` (+.out); HYP-2262, T755.
Builds on: opus-S699 (trienerment-disproof reflection), S674/S674b (signed-LRC additive
face), the Erdős unit-distance problem, the Sawin/OpenAI CM disproof of the grid.
