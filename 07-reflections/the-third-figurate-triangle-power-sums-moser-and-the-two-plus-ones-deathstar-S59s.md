# The third figurate triangle — power sums, the Moser break, and the two +1's as one table

**death-star-2026-07-20-S59s** (HYP-8165, THM-1355; owner: the triangle
{1},{2,1},{3,3,1},{4,6,5,1},{5,10,14,9,1},{6,15,30,37,17,1},{7,21,55,101,99,33,1},
its "third perspective on triangular numbers," Fibonacci / powers of 2 / Moser,
the family n·2^x+1, two rational series, and "how odd degree relates to
tournaments"). Method: three parallel repo-mining sweeps + exact computation.
Everything marked verified was computed this session; credits: opus-S317 /
mac-mini-S109 / klein-S313 (the figurate-triangle program this extends),
boxeph-S146 (odd-degree conjecture), HYP-6865 (harmonic = staircase resistance),
THM-466/871/448 (the 2^x+1 stacks), the observer principle (opus-2026-06-30).

## 1. The triangle is the THIRD figurate triangle — and the "third perspective" is literal

The repo already runs a **figurate-triangle program** (opus-S317's Vandermonde
truncation law): two triangles whose columns are figurate families.

- **Polyhedral triangle** (Pascal): column j = simplex numbers C(m+j−1, j);
  **row sums = 2^n**; shallow diagonals = **Fibonacci**.
- **Polygonal triangle**: column j = (j+2)-gonal numbers; **row sums = Moser's
  circle numbers A000127** (1,2,4,8,16,31,57…); shallow diagonals = a Fibonacci
  analog "D".

The owner's triangle is the **third** member: the **power-sum triangle**, column
j = Σ_{k=1}^{m} k^j. Column 0 = the naturals (Σk⁰), column 1 = the **triangular
numbers** (Σk), column 2 = the **square-pyramidal numbers** (Σk²) — all exact
(verified). So triangular numbers appear as column 1 of *all three* triangles:
polygonal (spatial figurate), polyhedral (Pascal / C(n,2)), and power-sum (Σk).
That is precisely the owner's "third perspective on the triangular numbers, in
addition to the polygonal and polyhedral." Its own row sums are 1, 3, 7, 16, 39,
106, 317 — neither 2^n nor Moser: a genuinely distinct third row-sum law.

## 2. The Moser break is the point — and the repo has been teaching it all along

The power-sum columns are exact for j = 0, 1, 2 and then **break at the cubes**:
column 3 is 1, 9, **37, 101** where Σk³ would give 1, 9, **36, 100** (verified).
The triangle looks like the power sums and then deviates — exactly Moser's circle
numbers looking like 2^n (1,2,4,8,16) and breaking to 31. This is not a defect;
it is the lesson. The repo has three independent instances of the same grammar:

- **Moser**: 1,2,4,8,16 → 31 (the polygonal row sums; 2^n truncated at binomial
  column 4; the "32nd region = the x⁵ coefficient", THM-868).
- **The hypotenuse**: H = 2^{n−2}+1 is exact for a single apex flip at *every* n
  (THM-250), but the *full tiling class* value 2^{n−2}+1 is only a small-n
  coincidence and **breaks at n=7** (31≠33; MISTAKE-115/119) into a Tribonacci
  law.
- **The power-sum triangle**: Σk^j exact through j=2, breaks at j=3.

Three "looks-simple-then-breaks" phenomena, one moral — the repo's own
MISTAKES-ledger refrain ("patterns that hold at n=3,4,5 break at 6,7"), now
visible as a *structural* feature of figurate truncation (opus-S317: the
deviation is always Pascal shifted two layers deep).

## 3. Fibonacci becomes cubic

Pascal's shallow diagonals sum to Fibonacci (x² = x+1, φ = 1.618). The owner
triangle's shallow diagonals sum to **1, 2, 4, 7, 12, 21, 37, 65, …**, satisfying
**a(n) = 2a(n−1) − a(n−2) + a(n−3)** (verified), whose growth constant is the
real root of **x³ − 2x² + x − 1 = 0, ≈ 1.7549** — a cubic Pisot analog of the
golden ratio, distinct from both φ and the tribonacci constant 1.839. So the
"Fibonacci of the power-sum triangle" is a degree-3 algebraic number: the
triangle deforms Pascal, and Pascal's quadratic golden ratio deforms with it into
a cubic. (This dovetails with mac-mini-S137's golden JC₂ corner and opus-S410's
θ-flow Fibonacci endpoint — the repo's Fibonacci appears whenever a shallow
diagonal or a worst-approximable slope does; here it appears cubically.)

## 4. n·2^x+1: the two +1's are one table (THM-1355)

The centerpiece. The family f(x,n) = n·2^x+1 (the **Proth numbers**) has the
repo's two master "+1" constants as its axes: **x=1 gives 2n+1** (the observer /
Rédei / LRC-modulus axis) and **n=1 gives 2^x+1** (the hypotenuse / Fermat /
Cayley-Dickson axis). The Fermat numbers are the n=1, x=2^m entries; the Cullen
numbers n·2^n+1 are the diagonal; the whole table (x≥1) is odd. Before today
these were two separate constants with two separate theorem stacks (THM-871 and
THM-466 on one axis; the observer principle and THM-401 on the other); the Proth
table is their common home, and the shared +1 is the OCF vacuum digit α₀ = 1 —
the same +1 the S59p–S59r Jacobian arc tracked as the conserved unit u = 1+xy.
The subdiagonal of the owner's *triangle* is 2, 3, 5, 9, 17, 33 — the 2^x+1
column of this table — so the triangle and the Proth table are stitched together
at the hypotenuse.

## 5. The two rational series are the two axes again

Series 1 is the **harmonic numbers** H_n = 1, 3/2, 11/6, 25/12, 137/60 (exact).
And the repo gives them a load-bearing meaning: HYP-6865 proves H_{n−2} is the
**electrical pole-resistance of the staircase Δ_{n−2}** — the very right-isoceles
triangle whose tilings are the tournaments — with γ (Euler–Mascheroni) as its
growth constant. So the harmonic series is not decoration; it is the staircase's
resistance, the "logarithmic" shadow (H_n ∼ ln n). Series 2 (the owner's second
list, transcription-noisy but a monotone 2-power-weighted harmonic, best fit
Σ(2^k−1)/k = Σ2^k/k − H_n) is the **exponential** shadow (∼ 2^n/n). Two series,
two growth regimes — logarithmic and exponential — mirroring the table's two
axes (the linear/observer 2n+1 and the exponential/hypotenuse 2^x+1). The smooth
(rational-series) world and the discrete (figurate-triangle) world tell the same
two-axis story.

## 6. Odd degree and tournaments (the closing question, answered)

n·2^x+1 is odd for all x ≥ 1, and oddness is the tournament parity invariant:
Rédei's H = 1 + 2·(…) is always odd because its 2-adic bottom digit is α₀ = 1
(THM-466). The Jacobian mirror is boxeph-S146's **odd-degree conjecture**
(Rédei-odd ↔ Keller cover-degrees odd; the verified counterexample's fiber is
1 + 2 = 3), and the observer principle makes A₁ ↔ the 3-cycle, the atomic odd
cycle. So "odd degree relates to tournaments" three ways at once: the family's
oddness *is* the Rédei signature; the Jacobian counterexample's cover-degree is
odd for the same parity reason; and the observer's +1 that makes everything odd
is the marked vertex / the vacuum digit / the conserved unit. The owner's whole
prompt is one object seen from many sides: the +1, and the odd numbers it
generates.

## 7. What is verified vs woven

Verified exactly: the three figurate row-sum laws; the power-sum columns and
their break at the cubes; the cubic-Fibonacci shallow diagonal and its growth
constant; the Proth table's four distinguished slices; harmonic = series 1. Woven
(cited, graded): the "third perspective" framing (real — triangular numbers in
all three triangles), the Moser-lesson unity across three break-phenomena, the
harmonic-as-resistance meaning (HYP-6865), the two-series/two-axis mirror, and
the odd-degree↔Rédei bridge (boxeph-S146). Open thread for a future session
(backlog): the exact closed form / OEIS identity of the owner's power-sum
triangle (columns mimic Σk^j with a small Pascal-deep correction; candidate
A143542, unverified — OEIS blocked this session) and whether its row-sum law
1,3,7,16,39,106,317 has a Vandermonde-truncation derivation like Moser's.

## Cross-links

THM-1355 (the Proth unification) · HYP-8165 · opus-S317 / mac-mini-S109 /
klein-S313 (the figurate-triangle program) · THM-466 (the +1 = α₀ vacuum) ·
THM-871 / THM-448 (Fermat / CD rungs) · THM-250/284 (single-flip H=1+2^d) ·
the observer principle (opus-2026-06-30; S59m §4; S59r) · boxeph-S146
(odd-degree) · HYP-6865 (harmonic = staircase resistance) · mac-mini-S137 /
opus-S410 (Fibonacci corners).
