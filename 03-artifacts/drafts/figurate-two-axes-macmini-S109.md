# The two lifts of the triangular numbers — polygonal vs polyhedral, Moser vs Fibonacci, and the support filtration between them

*mac-mini-2026-07-15-S109; owner directive (locker parity law, no-holes completeness,
polygonal/polyhedral triangles, A000127, the 1,1,2,3,5,8,13,21,33,… diagonal sums, the
difference patterns). Companion canon: THM-865 (locker parity law refuted), THM-866 (axis
level completeness proved). Script: `04-computation/figurate_two_axes_macmini_S109.py` →
`05-knowledge/results/figurate_two_axes_macmini_S109.out` (every claim below machine-exact
unless marked classical).*

*Concurrent-work note (same day, independent): kind-pasteur cont.21 (II) proved the
"rank-2 polygonal law" (diagonals affine in the column index, differences in the
C(j+d,d)-basis with tetrahedral leaders) and opus-S317 the "Vandermonde truncation law"
(polygonal = the k ≤ 1 Vandermonde layers, differences = Pascal-row tails) — both are the
§3/§4 structure in different bases; the A000127 row sums and the G sequence were found by
all three sessions independently (triple confirmation). New here beyond those: the
run/support/skeleton trinity (§3), the T_j tower with the proved first-miss law (§3),
G's GF/recurrence/closed form and OEIS status (§2), the A060488 identity (§4), Brown
completeness and the Fermat-polygonal split (§5), the master array N(s,d,m) (§0), and the
locker-tournament refutation (§6, THM-865).*

---

## 0. The two triangles, aligned

Left-justified triangles, row r ≥ 0, columns k = 0..r; write **m := r − k + 1** for the
index within column k (column k starts at row k).

- **Polyhedral (= Pascal).** Entry C(r,k). Column k, read down, is the k-simplex
  (hypertetrahedral) sequence C(m+k−1, k): ones, naturals, triangulars, tetrahedrals,
  pentatopes, … Rows sum to 2^r; shallow diagonals sum to Fibonacci. This is the triangle
  the repo has long used (tile count m = T_{n−2} is its (s,d) = (3,2) cell).
- **Polygonal.** Column 0 = ones; column k ≥ 1 = the (k+1)-gonal numbers
  Q(r,k) = P(k+1, m) = C(m,1) + (k−1)·C(m,2): ones, naturals, triangulars, squares,
  pentagonals, hexagonals, …

The two triangles **share the first three columns** (ones, naturals, triangulars) **and the
first two entries of every column** (1, then the shape/dimension label) — because, as §3
makes exact, disagreement requires both m ≥ 3 and k ≥ 3.

Both are slices of the **master figurate array**

> **N(s, d, m) = (s−2)·C(m+d−2, d) + C(m+d−2, d−1)**   (m-th d-dimensional s-gonal number)

with N(3,d,m) = the simplex column d and N(s,2,m) = the polygonal column s−1. The two
axes through the triangular numbers T = N(3,2,·) behave differently in kind:

- **shape axis (s): affine.** ∂N/∂s = C(m+d−2, d) is s-free; in the plane,
  P(s+1, m) = P(s, m) + T_{m−1} — adding a side glues one more triangular wedge
  (n² = T_{n−1} + T_n, the two-staircase gluing, is the s = 3→4 instance).
- **dimension axis (d): integral.** N(s, d, m) = Σ_{i≤m} N(s, d−1, i) — raising dimension
  is coning/stacking (partial summation), which is why simplex columns have growing degree
  while every polygonal column stays quadratic.

*The polygonal triangle is the quadratic world; Pascal is what iterated integration builds
over it.*

## 1. Row sums: 2^r vs Moser (A000127)

**Theorem 1.** Row r of the polygonal triangle sums to
**1 + C(r+1,2) + C(r+1,4) = A000127(r+1)** — Moser's circle numbers 1, 2, 4, 8, 16, 31,
57, 99, 163, 256, …

*Proof.* Σ_k C(m,1) = Σ_{j≤r} j = C(r+1,2), and by two hockey sticks
Σ_k (k−1)C(m,2) = Σ_{j≤r} (r−j)C(j,2) = Σ_i Σ_{j<i} C(j,2) = Σ_i C(i,3) = C(r+1,4). ∎

So the polygonal triangle does to A000127 what Pascal does to 2^r — exactly the owner's
observation. The famous "the circle count stops doubling at n = 6" is row r = 5: 31 = 32−1
(§4 locates the missing 1 precisely).

## 2. Diagonal sums: Fibonacci vs G(n)

Shallow (skipped-row) diagonal sums G(n) = Σ_k Q(n−k, k):

> **G = 1, 1, 2, 3, 5, 8, 13, 21, 33, 51, 76, 111, 157, 218, 295, 393, 513, 661, 838,
> 1051, 1301, 1596, 1937, 2333, 2785, …**

(the owner's 14 terms confirmed and extended; **not in OEIS** as of 2026-07-15 — submission
candidate). It is to Fibonacci precisely what Moser is to 2^r:

- **GF:**  G(x) = (1 − 2x + 3x³ − 2x⁴ + x⁶) / ((1−x)⁵(1+x)²).
- **Recurrence (order 7):** G(n) = 3G(n−1) − G(n−2) − 5G(n−3) + 5G(n−4) + G(n−5) − 3G(n−6) + G(n−7).
- **Closed form (exact, quasi-polynomial):**
  **G(n) = n⁴/96 − n³/16 + n²/3 + n/32 + 53/64 + (−1)ⁿ·(11/64 − n/32).**
- **Asymptotics:** G(n) ~ n⁴/96, vs Moser ~ n⁴/24: same degree, quarter constant — the
  diagonal/row ratio echoing F(n) ~ φⁿ vs 2ⁿ.
- **Fibonacci defect** F(n+1) − G(n) = 0, 1, 1, 2, 3, 5, 8, 13, **22**, 38, 68, 122, … —
  itself Fibonacci-like until its own second-order holes open.

## 3. The fundamental object: the run/support filtration

Everything above is one theorem seen twice. By Vandermonde,

> **C(m+k−1, k) = Σ_{s≥1} C(m, s)·C(k−1, s−1)**,

and the s-th summand has three equivalent readings (all machine-verified):

| s-th summand C(m,s)C(k−1,s−1) counts… | language |
|---|---|
| k-multisets from m types with **support exactly s** | multisets |
| k-subsets of [r] (r = m+k−1) with **exactly s maximal runs** | subsets |
| lattice points on the **open (s−1)-faces** of the dilated simplex k·Δ_{m−1} | Ehrhart |

**Theorem 2 (support-2 truncation).** The polygonal entry is the s ≤ 2 truncation:
Q(r,k) = C(m,1) + C(m,2)(k−1) = multisets with ≤ 2 distinct types = subsets with ≤ 2 runs
= lattice points on the **1-skeleton** of k·Δ_{m−1} (m vertices + (k−1) interior points on
each of the C(m,2) edges). Pascal is the untruncated sum (the full simplex body). ∎

This explains the shared frame instantly: support ≥ 3 needs m ≥ 3 **and** k ≥ 3, so the
triangles agree on the first three columns and the first two entries of every column.
It also gives the cleanest row-sum proof: a subset of [r] with s runs ⟺ 2s boundary
marks among the r+1 gaps, so #subsets with ≤ 2 runs = C(r+1,0) + C(r+1,2) + C(r+1,4) —
Theorem 1 again, now bijectively, and visibly the **even-binomial truncation** of
2^r = Σ_i C(r+1, 2i).

**The tower.** Let T_j be the triangle of ≤ j-run subset counts (T_1 ⊂ T_2 = polygonal ⊂ …
→ Pascal). Then (verified j ≤ 4, r ≤ 25; proofs one-line from the above):

- row sums of T_j = Σ_{i≤j} C(r+1, 2i): j = 1 gives 1 + C(r+1,2) (central polygonal /
  lazy-caterer family), j = 2 gives Moser, j → ∞ gives 2^r;
- diagonal sums G_j(n): G_1(n) = 1 + ⌊n²/4⌋ = **A033638** (quarter-squares + 1),
  G_2 = G above, G_3 = 1,1,2,3,5,8,13,21,34,55,89,**143**,228,358,554,841,… (not in
  OEIS), → Fibonacci;
- **first-miss law (proved):** the unique minimal support-(j+1) cell is
  (m,k) = (j+1, j+1), so T_j's rows first miss 2^r at r = 2j+1, deficit exactly 1, and
  G_j first misses Fibonacci at **n = 3j+2**, deficit exactly 1. For j = 2, r = 5:
  31 = 32 − 1 — *Moser's missing circle region is literally the unique 3-run subset
  pattern on 5 elements*; and G's first miss 33 = 34 − 1 at n = 8 is its diagonal twin.

## 4. The difference triangle — the owner's question answered

Δ(m,k) := C(m+k−1,k) − Q(r,k) is the Vandermonde tail:

> **Δ(m,k) = Σ_{j≥2} C(m, j+1)·C(k−1, j)**
> = (# of j-faces of Δ_{m−1}) · (interior lattice points of the dilated j-face), summed —
> the face-stratified Ehrhart tail beyond the 1-skeleton.

Per diagonal (fixed entry-index m), the coefficient row on the basis C(k−1, j), j ≥ 2 is
(C(m,3), C(m,4), …, C(m,m)) — **Pascal's row m with its first three entries deleted**:

| m | coefficients | difference sequence (k = 3, 4, 5, …) | identification |
|---|---|---|---|
| 3 | (1) | 1, 3, 6, 10, 15, 21, … | triangular numbers C(k−1,2) — the owner's first diagonal (with the typo fixed: the entry over octagonal-21 is C(9,2) = **36**, and 36 − 21 = 15) |
| 4 | (4, 1) | 4, 13, 28, 50, 80, 119, 168, … | **= A060488** ("4-block ordered tricoverings of an unlabeled n-set"), closed form (k−1)(k−2)(k+9)/6 — the owner's "subtle" second diagonal, now exact: 4·C(k−1,2) + C(k−1,3) = 4·triangular + tetrahedral |
| 5 | (10, 5, 1) | 10, 35, 81, 155, 265, 420, … | 10C(k−1,2) + 5C(k−1,3) + C(k−1,4); not in OEIS |
| 6 | (20, 15, 6, 1) | 20, 85, 231, 511, 1002, … | Pascal row 6 minus first three entries |

Down a fixed column instead: Δ(m,3) = C(m,3) — the squares column sits a full
**tetrahedral sequence** below the tetrahedral column. So: across the first differing
diagonal the deficit is *triangular*, down the first differing column it is *tetrahedral*,
and in general the deficit scheme is **Pascal re-entering with its first three columns cut
off** — the difference between the two lifts of the triangle is governed by the part of
Pascal beyond what they share. Self-similarity, exactly where the owner pointed.

## 5. No-holes completeness — the two axes of representability

Call a sequence *complete* if every positive integer is a sum of **distinct** terms.
**Brown's criterion:** a nondecreasing sequence with a₁ = 1 is complete iff it never
overshoots: a_{k+1} ≤ 1 + Σ_{i≤k} a_i. (⟸ by the greedy exchange walk downward — take
a_k iff the remainder allows, the no-overshoot inequality keeps the remainder in range;
⟹ because an overshoot makes A_k + 1 permanently unrepresentable — a hole.)

Machine-checked (60 terms each) and structurally clear:

- **2^r (Pascal rows): complete and TIGHT** — equality at every step — the unique
  no-slack case = binary numeration, each integer represented exactly once.
- **Fibonacci (Pascal diagonals): complete with slack** — the F3-style exchange
  F_{k+1} = F_k + F_{k−1} powers the walk; adding the no-consecutive rule restores
  uniqueness (Zeckendorf numeration).
- **Moser rows and G diagonals: still complete** (the truncations lose the numeration
  tightness but never overshoot; polynomial growth keeps them complete forever).
- **The columns FAIL**: squares (4 > 1+1) and triangulars (3 > 1+1) are incomplete as
  distinct sums — their completeness is the *bounded-repetition* kind: Gauss's
  three-triangular theorem and Lagrange's four squares, i.e., **Fermat's polygonal number
  theorem: s terms for the s-gonal column** (spot-verified to 500).

So "no holes" splits along the same two axes: **sum-sequences (rows/diagonals) are
distinct-part complete via a greedy exchange walk; columns are repetition-complete via
Fermat-polygonal**, with the tight/binary case 2^r as the corner. The same
greedy-plus-exact-step shape reappears in today's THM-866 (the F3 tie-splitting walk
proving the axis level lattice has no holes) — the repo's completeness proofs are all
"exchange walks with a conserved step".

## 6. The locker thread (cross-references)

The integer locker (τ(n) odd ⟺ n square, by the involution d ↔ n/d; open-locker gaps =
the odd numbers) selects the **s = 4 column** of the polygonal triangle. Today's THM-865
shows the τ-parity mechanism does **not** transfer to the locker tournament's Hamiltonian
count: H(D_11) = 4027 ≡ 3 (mod 4) refutes THM-853(III)'s conjectured law — parity taming
survives dimension-lifting in Pascal's world, but not the passage from divisor lattices to
tournament path counts. The two results are the same lesson with opposite outcomes: the
axis level lattice (THM-866) has a conserved exact step (+8) and hence no holes; α₁(D_n)
mod 2 has no conserved pairing and hence no law.

## 7. Engineering note (mandate)

Σ_{i≤j} C(n, 2i)-type truncations are Sauer–Shelah/VC growth functions: the T_j row sums
are the growth functions of "≤ j interval unions" concept classes (VC dimension 2j) —
used directly in learning theory for unions-of-intervals hypotheses, and the exact-run
counts C(m,s)C(k−1,s−1) are the standard reliability-theory census for
consecutive-k-out-of-n systems. The T_j/G_j family gives closed-form growth and
"scan-statistics" curves with exact quasi-polynomial forms (§2) — a small clean library
candidate (`figurate_growth.py`: N(s,d,m), T_j entries, row/diagonal sums, Brown checks).
Backlogged.

## 8. Open leads (added to INVESTIGATION-BACKLOG)

1. OEIS submissions: G (§2), G_3, the m = 5 difference row (10, 35, 81, …). Δ-row m = 4
   is A060488 — *why do 4-block ordered tricoverings equal the second-diagonal deficit?*
   A bijection tricoverings ↔ support-3/4 multisets would be new content for that entry.
2. A geometric incidence model for G(n): Moser's rows count circle regions; what planar
   process does the diagonal truncation count? (The tiling reading exists already: G(n) =
   square/domino tilings of an n-strip whose dominoes form ≤ 2 runs — via the multiset
   dictionary; a regions picture would be better.)
3. The G_j "first-miss at 3j+2, deficit 1" law generalizes: second miss, full deficit GF =
   Fibonacci-diagonal sum of the Δ-tail — compute the deficit triangle's own diagonal GF.
