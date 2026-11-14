# Rule 30's center column as a continued fraction and as a 2-adic number

mac-mini-2026-11-13-cuts · **All three Rule 30 prize problems remain OPEN.** Status of what follows:
- **ELEMENTARY** (§1, the reformulation);
- **FINITE-EXACT** (§2–§3, with positive and negative controls);
- **SURVEY** (§4, the repository's own Rule 30 program).

Scripts: [`center.c`](scripts/center.c) (packed single-seed evolution, checked bit by bit against direct simulation
for 200 steps), [`bm.c`](scripts/bm.c) (packed Berlekamp–Massey, matching the reference
[`bmref2.py`](scripts/bmref2.py) on 65,536 bits), [`raa.py`](scripts/raa.py) and [`raa_run.py`](scripts/raa_run.py)
(Klapper–Goresky rational approximation; recovers −5/13 exactly).

## 1. Problem 1 as an irrationality problem

Write c_0 c_1 c_2 … = 1101110011000101… for the center column from a single black cell. Attach two numbers to it:
- the power series α(X) = Σ c_t X^t in F_2[[X]];
- the 2-adic integer β = Σ c_t 2^t in Z_2.

A sequence is eventually periodic iff α is a rational function, and iff β is a rational number. So **Prize 1 asks
whether α ∉ F_2(X), equivalently whether β ∉ Q.** Each has an exact rational-approximation profile:
- for α, the linear complexity L(n), computed by Berlekamp–Massey. Its jumps are the degrees of the partial quotients
  of the continued fraction of α in F_2((1/X)).
- for β, the 2-adic complexity log2 Φ(n), computed by Klapper–Goresky. It is the bit size of the smallest rational
  p/q (q odd) that agrees with the first n digits.

For a random sequence, both profiles follow n/2 with small fluctuations. The partial-quotient degrees are then
geometric, P(deg = k) = 2^(−k), and there are about n/4 quotients.

## 2. The linear complexity profile (FINITE-EXACT)

| n | L(n) | max over m ≤ n of abs(2L(m) − m) | ones among c_0..c_(n−1) |
|---|---|---|---|
| 2^20 = 1,048,576 | 524,288 = n/2 | 19 (at m = 937,375) | 524,976 (+688, 0.67σ) |
| 2^22 = 4,194,304 | 2,097,152 = n/2 | 25 (at m = 4,033,763) | 2,097,861 (+709, 0.35σ) |

There are 261,724 partial quotients up to 2^20 (n/4 = 262,144 expected). Their degree histogram is

1: 130,572 · 2: 65,720 · 3: 32,588 · 4: 16,342 · 5: 8,186 · 6: 4,138 · 7: 2,095 · 8: 1,029 · 9: 536 ·
10: 271 · 11: 124 · 12: 65 · 13: 32 · 14: 11 · 15: 8 · 16–18: 2 each · 19: 1,

halving at every step, as for a random series. The large quotients (degree ≥ 13) sit at scattered positions with no
dyadic pattern: 19,520, 35,022, 46,593, …, 420,933 (degree 18), 937,375 (degree 19).

**Consequences.**
- *Prize 1 at finite scale.* An eventually periodic sequence with onset t_0 and period p has L(n) ≤ t_0 + p for every
  n. Hence **t_0 + p ≥ 2,097,152** for any eventual periodicity of the center column.
- *Prize 3 in the LFSR model.* For every n ≤ 2^22, no linear feedback shift register shorter than n/2 − 12.5 outputs
  c_0 … c_(n−1). In this model the center column is as incompressible as a random sequence at every scale tested.

**Controls.**
- A pseudo-random file from Python's Mersenne Twister has linear complexity **19,938** at n = 65,536. That is
  19,937 + 1: Berlekamp–Massey finds MT19937's GF(2)-linear state, a 19,937-bit register (19,937 is a Mersenne
  exponent).
- The same file has full 2-adic complexity (§3). The FCSR model is blind to GF(2)-linearity, and Rule 30 passes in
  both models.
- A periodic control (onset 5, period 7) stops at L = 10 ≤ 12.

### 2.1 Extension to 2^22

There are 1,047,637 partial quotients (n/4 = 1,048,576 expected). Their degree histogram is

1: 523,280 · 2: 262,189 · 3: 130,759 · 4: 65,587 · 5: 32,916 · 6: 16,428 · 7: 8,271 · 8: 4,077 · 9: 2,079 ·
10: 1,040 · 11: 514 · 12: 261 · 13: 120 · 14: 53 · 15: 36 · 16: 8 · 17: 12 · 18: 3 · 19: 2 · 20: 1 · 25: 1.

The single degree-25 quotient (at m = 4,033,763) is the only mild outlier. Among about 10^6 geometric samples, one at
least that large has probability about 6%.

### 2.2 The 2-kernel (automaticity)

By Christol's theorem, α is algebraic over F_2(X) iff c is 2-automatic, iff its 2-kernel {(c_(2^k n + r))_n} is
finite. On the 2^22-bit prefix, **all 32,767 kernel elements with k ≤ 14 are pairwise distinct already on their
first 256 terms** ([`kernel.py`](scripts/kernel.py)). Any automaton that reads t least-significant digit first and
outputs c_t on this range therefore needs at least 32,767 states. This extends the 8,191-element (k ≤ 12) check on
131,072 cells reported in arXiv:2609.25077. It is evidence, not proof: a finite check cannot exclude a larger finite
kernel.

## 3. The 2-adic complexity profile (FINITE-EXACT)

| n | log2 Φ(n) | max over m ≤ n of abs(2 log2 Φ(m) − m) |
|---|---|---|
| 4,096 | 2,045 | 11 |
| 16,384 | 8,194 | 12 |
| 65,536 | 32,768 = n/2 | 21 (at m = 42,397) |

The random control (Mersenne Twister bits) gives 32,768 and 16. So **β has no rational approximation p/q of height
below 2^(n/2 − 10.5) agreeing with its first n digits, for every n ≤ 65,536.** Equivalently, no feedback-with-carry
register with a connection integer of fewer than n/2 − 10.5 bits generates the prefix. Such registers are the
2-adic division machines x ↦ (x − b·q)/2, the same carry dynamics as the repository's Collatz parity numbers
(THM-4634 dictionary).

## 4. Where the repository's Rule 30 program stands (survey of 36 theorems, THM-3456 … THM-4210)

- **Reductions.**
  - Left-permutive trace bijection (THM-3456): the center trace together with the right half-line determines the left
    boundary.
  - Periodicity from a given onset with a given period is exactly the vanishing of one inverse left-boundary word
    (THM-4048).
  - Cesàro balance is equivalent to o(2^m) discrepancy on dyadic shells (THM-4048), and Prizes 1 and 2 reduce through
    the dyadic block-current Cartier tree (THM-4210).
- **Fixed columns are tame.**
  - Every fixed left-front column is eventually periodic with dyadic period, through an affine monodromy clock
    (THM-4047).
  - There are infinitely many isolated eventual-zero columns, the next after depth 87,866 being 1,420,878,968
    (THM-4065).
  - The right edge is a 2-adic odometer (THM-3458).
- **The barrier** (THM-4047's own sentence). The center is the moving diagonal c_t = ℓ_t(t), and the onset bounds of
  the fixed-column theorems do not reach it. Haar-random-row independence (THM-4085) is likewise not a single-seed
  statement.
- **Literature (2026 preprints, read only as abstracts).** [arXiv:2609.25077](https://arxiv.org/abs/2609.25077)
  covers diagonal periods and Newton supports of Rules 30, 86 and 135. [arXiv:2609.09431](https://arxiv.org/abs/2609.09431)
  shows that finite seeds cannot give a constant trace. Both state that the prizes remain open.

**What this note adds.** A third and fourth exact observer, the F_2 continued fraction and the 2-adic approximation
profile, neither of which appeared in the repository before. Both are blind to the moving-diagonal barrier, because
they read the center column alone. Their certificates are finite, but they bound *every* linear (or carry-linear)
mechanism at once, not one onset and period at a time.

## 5. Proof paths and predictions

- **R1 (transcendence route).** Algebraic over F_2(X) means 2-automatic (Christol), which means a finite 2-kernel.
  §2.2 finds no coincidence among 32,767 kernel elements. A proof that the center's 2-kernel is infinite would show
  α transcendental, hence irrational, and so **would settle Prize 1**, since eventually periodic sequences are
  automatic. It is a stronger statement than Prize 1, but it is a statement about a self-similar family, the kernel,
  where THM-4210's dyadic block-current tree already supplies the recursion. The concrete target is a structural
  reason why no two kernel elements (c_(2^k n + r))_n and (c_(2^k′ n + r′))_n can coincide, propagated down the
  kernel tree.
- **R2 (irrationality measure).** The maximal partial-quotient degree up to n is about log2 n: 19 at 2^20 and 25 at
  2^22. A proof that it is o(n) would *not* settle Prize 1. A proof that infinitely many quotients exist would.
- **R3 (Prize 3 in restricted models).** The LFSR and FCSR models are now quantified. The next model is quadratic
  feedback shift registers (Chan–Games quadratic span). Rule 30 is itself quadratic over F_2, so the question is
  whether its center has small quadratic span. **P-R3:** it does not; the quadratic span grows linearly.
