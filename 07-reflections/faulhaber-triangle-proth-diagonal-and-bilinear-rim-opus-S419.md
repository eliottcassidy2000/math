# The Faulhaber triangle, the Proth diagonal, and the bilinear +1 rim

**Instance:** opus-2026-07-20-S419 (owner: the triangle {1},{2,1},{3,3,1},...,{7,21,55,101,99,33,1};
2n+1 and 2^x+1; the n*2^x+1 table; Fibonacci/powers-of-2/Moser; two rational series).
**HYP-8155.** Scripts + frozen outs: `faulhaber_triangle_series_opus_S419.py`, `..._S419b.py`.

## 1. The triangle IS the Faulhaber (power-sum) triangle — the "third perspective" identified

**T(n,k) = Sum_{j=1}^{n-k+1} j^(k-1)** matches 25/28 given cells EXACTLY, including every
boundary law: col 1 = n (Sum j^0), col 2 = triangular (Sum j), col 3 = square pyramidal
(Sum j^2), **penultimate col = Sum_{j<=2} j^(n-2) = 2^(n-2)+1** (the owner's 2^x+1 is
itself a Faulhaber column evaluated at just two terms!), last col = 1. The three
deviations are +1 at the deep-interior cells (6,4), (7,4), (7,5); the flat law
"+1 iff k>=4 and n-k>=2" fits all 28 cells (S419b B1). No integer linear shift-recurrence
exists up to 5 terms over 9 shifts (B2) — consistent, since column k is a degree-k
polynomial in n, so no fixed small shift-stencil can generate the array.

**The third perspective:** triangular numbers sit in three towers — polygonal (fix dim 2,
vary shape), polyhedral/simplicial (fix shape, vary dimension), and now **Faulhaber
(vary the EXPONENT: Sum j^0, Sum j^1, Sum j^2, Sum j^3, ...)**. The owner's triangle
packages the third tower as a bona fide triangle, whose rim is Fermat-shaped (2^x+1).

**Row-8 prediction under the flat-correction law:** 8, 28, 91, 226, 355, 277, 65, 1.
(Falsifiable: the correction law is underdetermined by 3 data points — this is the test.)

**Fibonacci hint concretized:** diagonal sums 1,2,4,7,12,21,37 satisfy
**a(n) = a(n-1) + a(n-2) + a(n-4)** (verified on all three available instances).

**Moser cross-check:** T(7,5) = 99 = R(8) = C(8,4)+C(8,2)+1, i.e. 1+2^4+3^4+1 = 70+28+1.
Charm-grade coincidence — but ON THEME: the session's subject is exactly why low-order
jets of different growth families collide (Moser's 1,2,4,8,16,31 vs 2^n).

## 2. The two rational series = harmonic transforms of lines in the n*2^x+1 table

Define the harmonic transform of a line L(k) in the table: **H[L](n) = Sum_{k<=n} L(k)/k.**

- **Series 1 = H[all-ones row] = H_n** (1, 3/2, 11/6, 25/12, 137/60): EXACT, 5/5.
- **Series 2 = H[diagonal] = Sum (k*2^(k-2)+1)/k = H_n + 2^(n-1) - 1**: increments
  1/k + 2^(k-2) = **Proth numbers over k** (numerators 3, 7, 17, 41 = the table diagonal).
  Matches printed terms 1, 2, 4 exactly (109/12 arbitrates AGAINST the rival
  Erdos-Borwein partial sums Sum(2^k-1)/k = Sum C(n,k)/k [identity verified], which
  gives 103/12); term 3 matches modulo a denominator slip (29/6 vs printed 29/3);
  term 5 mismatch 1037/60 vs printed 1079/60 (Delta = 7/10) — flagged honestly, both
  candidates logged. The hybrid's exact 4th-term hit + the Proth-increment structure
  make it the identification.
- Micro-illusion found en route: the Proth diagonal 3, 7, 17, 41 coincides with the
  companion-Pell sequence (3, 7, 17, 41, 99) for four terms, then breaks (97 vs 99) —
  the session generated its own Moser-style trap, with 99 appearing AGAIN.

Owner's table corners verified: (x,0)=1, (1,n)=2n+1, (x,1)=2^x+1; note (0,n) = n+1,
not n (off-by-one in the prompt as stated — the n*2^0+1 column is n+1).

## 3. The bilinear +1 rim — the repo-wide unification the owner pointed at

"2n+1 and 2^x+1 are key to mathematical clarification": every quantized denominator our
LRC(14) proofs produced has the SAME bilinear shape m*b + e, e in {+1,-1} — the Proth
shape n*2^x+1 is one column of a family the repo has been mining all along:

| object | form |
|---|---|
| Proth numbers / owner's table | n*2^x + 1 |
| odd axis / GW escape | 2N+1 (M(GW) = 2/(2N+1)) |
| S-T ladder rungs | s/(ns+1) |
| (D,s) rung ladder (THM-1269) | D/((N+1)D - 1) |
| ghost packing (THM-1292) | K/(K(m+1) + 1) |
| **the deep well** | **14/183, 183 = Phi_6(14) = 13*14 + 1** |

The deep well's 183 is literally Proth-shaped (13*14+1). The "clarification" reading:
these +-1 offsets are the SAME quantization — a resonance denominator is a lattice count
m*b plus a duty of +-1, and the floor problems we solve are about which duty sign the
geometry can afford. The harmonic crown (THM-1153, currency H_n = series 1) and the
dyadic sheet duties (2^(n-1)-1 = Mersenne ladder inside series 2) are the two AXES of
the same table; series 2's hybrid is the first object we've met that spends both
currencies in one ladder.

Moser's circle B(n-1,4) = 1,2,4,8,16,31,57,99 is the canonical partial-binomial-sum
trap; the repo's own instance is the width-of-G_n formula C(n-2,floor((n-2)/2)) —
exact at n=3..6, FAILS at n=7 (predicted 10, actual 15). Same genus: low-order
binomial jets impersonating exponentials. MISTAKES category E, now with a named
external archetype.

## 4. Handoffs

(a) Row-8 test: obtain/derive the owner's row 8 and test the flat-correction prediction
(8, 28, 91, 226, 355, 277, 65, 1) — if corrections grow instead, the correction triangle
is a new object (candidate: its own Faulhaber cascade). (b) Series-2 terms 3 and 5:
confirm with owner whether 29/3, 1079/60 are exact as intended; if YES, the rule is
genuinely different and non-monotone-increment (increments would include -7/12) —
a new hunt. (c) The bilinear rim as a THEOREM: state and prove the "duty sign" law —
which m*b+-1 denominators are achievable as LRC(14) floors (connects THM-1269's
D = M*s to Proth-primality-style constraints). (d) Fibonacci: locate the
a(n)=a(n-1)+a(n-2)+a(n-4) diagonal-sum law inside the Faulhaber identity (should
follow from the polynomial column structure + rim; 3 verified instances only).

Cross-links: THM-1153 (harmonic crown) . THM-1269 (D = M*s) . THM-1292 (ghost packing)
. S399 synthesis (MISTAKES category E / small-n extrapolation) . S410 rung theory
(duty calculus) . scripts + frozen outs S419, S419b.
