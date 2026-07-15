# The Two Axes of the Triangular Numbers and the Vandermonde Tail

**Instance:** klein-2026-07-15-S313
**Status:** PROVED (all identities verified computationally in
`04-computation/polygonal_polyhedral_vandermonde_tail.py`, output in
`05-knowledge/results/polygonal_polyhedral_vandermonde_tail.out`; proofs below are complete)
**Thread:** additive-basis / Pascal-slope thread (T1083–T1087, HYP-2998–HYP-3003)

---

## 0. Setup: one seed, two axes

The triangular numbers T_j = C(j+1, 2)-ish live at the corner of a 2-parameter family.
Fix position index j ≥ 1 within a sequence and generalize:

- **Polyhedral (dimension) axis** — the k-simplex numbers
  `S_k(j) = C(j+k-1, k)`: column k of Pascal's triangle
  (k=0: all ones; k=1: naturals; k=2: triangular; k=3: tetrahedral; k=4: pentatope; ...).

- **Polygonal (gonality) axis** — the s-gonal numbers
  `P(s, j) = (s-2)·C(j,2) + j`
  (s=3: triangular; s=4: squares; s=5: pentagonal; ...), placed in column k = s−1
  (k=0: ones; k=1: naturals = "2-gonal"; k=2: triangular; k=3: squares; ...).

Two triangles result. Row n of the **polyhedral triangle** is Pascal's row C(n, k).
Row n of the **polygonal triangle** has entries `G(n,k) = P(k+1, n-k+1)` (with G(n,0)=1):

```
1
1  1
1  2  1
1  3  3  1
1  4  6  4  1
1  5 10  9  5  1          <- first deviation: 9 = P(4,3) (square) vs 10 = C(5,3) (tetrahedral)
1  6 15 16 12  6  1
1  7 21 25 22 18  7  1
```

**Agreement zone (proved in §2):** the two tables agree exactly on columns k ≤ 2 and on the
first two entries (j ≤ 2) of every column — 1, then k+1 — and NOWHERE else.

---

## 1. The master identity: polygonal = rank-2 truncation of a Vandermonde expansion

**Lemma 1 (Vandermonde mode expansion of the simplex numbers).**
For k ≥ 1, j ≥ 1:

    S_k(j) = C(j+k-1, k) = Σ_{i≥0} C(k-1, i) · C(j, i+1).

*Proof.* Vandermonde: C(m+n, r) = Σ_i C(m, i) C(n, r−i) with m = j, n = k−1, r = j−1
applied to C(j+k−1, k) = C(j+k−1, j−1):
C(j+k−1, j−1) = Σ_t C(j, j−1−t) C(k−1, t) = Σ_t C(j, t+1) C(k−1, t). ∎

**Lemma 2 (polygonal numbers are the i ≤ 1 modes).**

    P(k+1, j) = C(j, 1) + (k-1) · C(j, 2).

*Proof.* P(s,j) = ((s−2)j² − (s−4)j)/2 = (s−2)·j(j−1)/2 + j with s = k+1. ∎

Geometrically Lemma 2 is the **fan decomposition** of the (k+1)-gon from one vertex into
(k−1) triangles: a polygonal number is one spine of j points plus (k−1) triangular wedges
T_{j−1}. The polygonal axis is *affine in k* — a rank-2 motion in the space spanned by
C(j,1), C(j,2). The polyhedral axis genuinely raises dimension and needs every mode.

**Theorem 3 (the Vandermonde tail).** The difference table Δ_k(j) := S_k(j) − P(k+1, j) is
exactly the tail of the mode expansion:

    Δ_k(j) = Σ_{i≥2} C(k-1, i) · C(j, i+1).

In particular Δ ≥ 0 everywhere, and Δ_k(j) = 0 iff k ≤ 2 or j ≤ 2 (the agreement zone:
C(k−1,i) = 0 for i ≥ 2 needs k ≤ 2; C(j,i+1) = 0 for i ≥ 2 needs j ≤ 2). ∎ (Lemmas 1−2.)

### 1a. The diagonals of the difference table (the requested patterns)

Fixing j (the "j-th diagonal", i.e. the j-th entry of each column) and letting k run:

| j | law in k | values (k = 3,4,5,...) |
|---|----------|------------------------|
| 3 | `C(k-1,2)` | 1, 3, 6, 10, 15, ... (**triangular numbers** — as observed; note 36−21, not 35−21) |
| 4 | `4·C(k-1,2) + C(k-1,3) = (k-1)(k-2)(k+9)/6` | 4, 13, 28, 50, 80, 119, 168, 228, 300, 385, ... (**= OEIS A060488**, 4-block ordered tricoverings) |
| 5 | `10·C(k-1,2) + 5·C(k-1,3) + C(k-1,4)` | 10, 35, 81, 155, 265, 420, ... |
| j | `Σ_{i=2}^{j-1} C(j, i+1) · C(k-1, i)` | coefficient vector = **Pascal row-j tail** (C(j,3), C(j,4), ..., C(j,j)) |

So the "subtle pattern" of the second diagonal is `4·(triangular) + (tetrahedral)`:
4·(1,3,6,10,15) + (0,1,4,10,20) = (4,13,28,50,80). Every deeper diagonal is the same
construction one column further into Pascal: the coefficient rows (1), (4,1), (10,5,1),
(20,15,6,1), ... are the rows of Pascal's triangle truncated to columns ≥ 3.

Dually, fixing the column k and letting j run: Δ_3(j) = C(j,3) (tetrahedral − square =
shifted tetrahedral), Δ_4(j) = 3C(j,3) + C(j,4), Δ_5(j) = 6C(j,3) + 4C(j,4) + C(j,5) —
coefficients C(k−1, 2), C(k−1, 3), ..., the transposed law.

**Theorem 4 (the difference table is smoothed Pascal).** The bivariate generating function is

    Σ_{j,k} Δ_k(j) x^j y^k = x³y³ / [ (1-x)³ (1-y)² (1-x-y) ].

Equivalently: Δ is Pascal's table C(a+b, a) integrated (cumulatively summed) **three times
along the j-axis and twice along the k-axis**, then shifted to start at (3,3). The
polyhedral−polygonal defect is not new data — it is Pascal's triangle again, five
integrations deep. (Proof: sum the geometric series in Theorem 3's expansion; verified
entry-by-entry for j,k ≤ 21.)

---

## 2. Row sums: powers of 2 versus Moser (A000127)

**Polyhedral rows:** Σ_k C(n,k) = 2^n. In mode form, convolving Lemma 1 along a row with the
diagonal Vandermonde identity Σ_a C(a,i)C(n−a, i+1) = C(n+1, 2i+2) gives the **even-binomial
mode decomposition**

    2^n = 1 + Σ_{i≥0} C(n+1, 2i+2)   (mode i contributes the (2i+2)-th even binomial).

**Polygonal rows:** only modes i ≤ 1 survive (Lemma 2), so

    R(n) = 1 + C(n+1, 2) + C(n+1, 4) = Σ_{j=0}^{4} C(n, j) = A000127(n+1)  — Moser's circle numbers,

i.e. 1, 2, 4, 8, 16, **31**, 57, 99, 163, 256, ... The celebrated "looks like 2^n then
breaks" is, from this vantage, literal: Moser IS 2^n with the Vandermonde tail (modes i ≥ 2)
deleted, and the tail identity Σ_{i≥2} C(n+1, 2i+2) = Σ_{j≥5} C(n, j) converts the missing
modes into the missing binomial columns. Both classical forms of A000127 drop out of one
Vandermonde convolution.

---

## 3. Skip-row (shallow-diagonal) sums: Fibonacci versus the new sequence D

**Polyhedral diagonals:** Σ_k C(d−k, k) = F(d+1) (classical). Applying the mode expansion
along the *skew* direction (j = d − 2k + 1) gives the **Fibonacci layer decomposition**

    F(d+1) = 1 + Σ_{i≥0} L_i(d),    L_i(d) = Σ_k C(k-1, i)·C(d-2k+1, i+1),
    GF of L_i:  x^{3i+2} / [ (1-x²)^{i+1} (1-x)^{i+2} ],

and summing the geometric series in i reproduces 1/(1−x−x²) exactly — Fibonacci as a tower
of Vandermonde layers (the skew analogue of the even-binomial decomposition of 2^n).

**Polygonal diagonals (the new object):** keeping modes i ≤ 1 only,

    D(d) = 1 + L_0(d) + L_1(d),   L_0(d) = ⌊d²/4⌋  (quarter-squares A002620),

with values **1, 1, 2, 3, 5, 8, 13, 21, 33, 51, 76, 111, 157, 218, 295, 393, 513, 661, 838,
1051, ...** — Fibonacci through d = 7, then permanently below it. This sequence is **not in
OEIS** (checked 2026-07-15, including interior windows) — submission candidate.

**Theorem 5 (closed form).**

    D(d) = [ 2d⁴ − 12d³ + 64d² + 6d + 159 + (−1)^d (33 − 6d) ] / 192,
    D(2m) = (m⁴ − 3m³ + 8m² + 6)/6,     D(2m+1) = (m⁴ − m³ + 5m² + 7m + 6)/6.

*Proof.* GF: D(x) = 1/(1−x) + x²/[(1−x²)(1−x)²] + x⁵/[(1−x²)²(1−x)³] (the three retained
layers). Partial fractions of the last term (A₅..A₁ = 1/4, 1/4, 3/16, 1/8, 5/64; B₁ = 5/64,
B₂ = 1/32 at (1+x)) and ⌊d²/4⌋ = (2d²−1+(−1)^d)/8 assemble to the stated quasi-polynomial;
verified for d ≤ 40. ∎

**Theorem 6 (D is to Fibonacci what Moser is to 2^n).**

    D(d) = F(d+1) − t(d),   t(x) = x⁸ / [ (1-x²)² (1-x)³ (1-x-x²) ],

and the Fibonacci pole **cancels**: (1−x)³(1−x²)² − x⁸ is divisible by 1−x−x². The
cancellation certificate is pleasingly Fibonacci-flavored: modulo 1−x−x² one has
x^n ≡ (−1)^n(F(n−1) − F(n)x), so x⁸ ≡ 13 − 21x, while (1−x)³(1−x²)² ≡ (5−8x)(1−x)·x ≡
13 − 21x as well. Hence D grows as a quartic quasi-polynomial (degree 4 pole at x=1,
parity ripple from (1+x)²) — the exact analogue of Moser's quartic 1 + C(n+1,2) + C(n+1,4),
with the parity ripple as the fingerprint of row-skipping. ∎

After the cancellation the reduced form is

    D(x) = (1 − 2x + 3x³ − 2x⁴ + x⁶) / [ (1−x)⁵ (1+x)² ],

with **palindromic** numerator Q(x) = x⁶·Q(1/x) — a reciprocal polynomial, so D shares the
self-dual flavor of the h-vector world (Q(1) = 1 fixes the leading quartic coefficient
1/(4!·2²)·... = 2/192, i.e. D(d) ~ d⁴/96).

**Proportion:**  2^n : A000127 :: Fibonacci : D — in both cases "keep Vandermonde modes
i ≤ 1", i.e. row sums / diagonal sums of the polygonal triangle versus Pascal.

---

## 4. Locker parity law (divisor-pairing involution)

**Theorem 7.** For n ≥ 1, τ(n) (number of divisors) is odd iff n is a perfect square.
Consequently, in the locker process (locker n toggled once per divisor of n, starting
closed), locker n ends **open iff n is a square**.

*Proof.* The map ι(d) = n/d is an involution on Div(n). Its orbits partition Div(n) into
2-element orbits {d, n/d} and fixed points; d is fixed iff d² = n, which has at most one
positive solution, existing iff n is a square. So τ(n) = 2·#(pairs) + [n is a square], and
τ(n) ≡ [n is a square] (mod 2). Locker n's final state is the parity of its toggle count,
which is τ(n). ∎

This is the same engine as Rédei's theorem: *an involution with controlled fixed points
computes a parity*. Rédei pairs Hamiltonian paths (path-reversal / arc-flip involutions,
fixed-point analysis mod 2); the locker law pairs divisors. It is the ur-example of the
project's parity method.

## 5. No-holes completeness via the F3 exchange walk

Call a nondecreasing sequence a₁ ≤ a₂ ≤ ... **no-holes complete** if every positive integer
is a sum of distinct terms. ("Hole" = an unrepresentable integer.)

**Lemma 8 (Brown's criterion, with proof).** If a₁ = 1 and a_{m+1} ≤ 1 + Σ_{i≤m} a_i for all
m, the sequence is complete: by induction every N ≤ Σ_{i≤m} a_i is a sum of distinct terms
with indices ≤ m. (If N ≤ Σ_{i≤m−1}, use the induction hypothesis; otherwise
N > Σ_{i≤m−1} ≥ a_m − 1 forces N ≥ a_m, and N − a_m ≤ Σ_{i≤m−1}, recurse.) A "hole" appears
exactly when some term jumps past the reachable interval — hence the name. ∎

**Theorem 9 (F3 exchange walk ⇒ Zeckendorf; constructive completeness of Fibonacci).**
Work with multisets of Fibonacci indices ≥ 2 (F₂ = 1, F₃ = 2, ...), value Σ F_c. The
**F3 exchange walk** applies the 3-term Fibonacci relation as rewrite rules:

    (merge)      {k, k+1} → {k+2}                    [F_k + F_{k+1} = F_{k+2}]
    (duplicate)  {k, k}   → {k+1, k−2}   (k ≥ 4)     [2F_k = F_{k+1} + F_{k−2}]
                 {3, 3}   → {4, 2},  {2, 2} → {3}    [base cases]

Every rule preserves the value. With the integer potential Φ = Σ w(c), w(2) = 8, w(c) = 5c
(c ≥ 3), **every rule strictly decreases Φ** (by 3; ≥10; 1; 2; ≥5 respectively — the weight
8 = "1.6 × 5" is tuned so that the two base cases both drop). Hence the walk terminates.
Terminal states admit no rule: no duplicates and no consecutive indices — i.e. **Zeckendorf
form**. Starting from N·{2} (N ones), the walk therefore *constructs* a Zeckendorf
representation of every N ≥ 1: Fibonacci is no-holes complete. Uniqueness (hence global
confluence of the walk, no Newman's lemma needed): with top index m, the maximal Zeckendorf
value is F_m + F_{m−2} + ... = F_{m+1} − 1 < F_{m+1}, so F_m ≤ N < F_{m+1} forces the top
index; induct downward. ∎

(Verified: walk run from N·{2} and from deliberately messy representations for N ≤ 2000;
potential asserted to drop at every step; terminal state always the greedy Zeckendorf form.
This is the "carry-confluence" endpoint of the additive-basis ladder of HYP-2998/HYP-3000 —
the walk IS Fibonacci carry-propagation, and Φ is an explicit termination certificate.)

**Binary side.** The same walk with the 2-term rule {k,k} → {k+1} (2·2^k = 2^{k+1}) proves
every N has a unique binary representation: powers of 2 are the *extremal* no-holes complete
sequence (Brown holds with equality; uniqueness ⟺ extremality).

**Theorem 10 (completeness survives the Vandermonde-tail truncation).**
D and Moser are also no-holes complete.
*Proof.* For D: t(d) − t(d−1) − t(d−2) = [x^d] x⁸/((1−x²)²(1−x)³) ≥ 0 (all three factors
have nonnegative series), so D(d) ≤ D(d−1) + D(d−2) ≤ 1 + Σ_{i<d} D(i), with the last step
since both terms appear among the (positive) earlier terms; D starts 1, 1 — Brown applies.
(Equality D(d) = D(d−1) + D(d−2) holds iff d ≤ 7.) For Moser: C(n,j) = C(n−1,j) + C(n−1,j−1)
gives R(n) ≤ 2R(n−1) after truncation at j ≤ 4, and Brown again applies. ∎

So the two triangle sum-laws produce the two canonical complete sequences (binary and
Zeckendorf), and their polygonal shadows (Moser, D) remain complete — **no holes open up
when the tail is cut**, even though uniqueness (the extremal/binary property) is lost.

---

## 6. Synthesis — what the fundamental objects are

1. **The kernel.** Everything is the two-variable Pascal kernel 1/(1−x−y). The polyhedral
   table is the kernel itself; the polygonal table is its **rank-2 truncation** (modes
   i ≤ 1 of the Vandermonde expansion); the difference is the **Vandermonde tail**, which is
   again the Pascal kernel, five integrations deep (Theorem 4). The triangle measures its
   own shadow with itself.
2. **Two functionals, one truncation.** Row-sum functional: kernel → 2^n, truncation →
   Moser A000127. Skew functional (skip rows): kernel → Fibonacci, truncation → D (new,
   quartic quasi-polynomial, Theorem 5). The parity ripple (−1)^d(33−6d)/192 in D is the
   trace of the period-2 skew; Moser, an un-skewed sum, has none.
3. **Mode ladders.** 2^n = 1 + Σ_i C(n+1, 2i+2) and F(d+1) = 1 + Σ_i L_i(d) are the same
   decomposition under the two functionals; Moser and D are their common i ≤ 1 prefix.
4. **Parity and completeness.** The locker law (divisor involution) and the F3 exchange
   walk (rewriting to Zeckendorf normal form with an explicit potential) are the two
   verification engines — an involution computing a parity, a terminating confluent walk
   computing a normal form — and completeness is robust to the tail truncation (Theorem 10).
5. **Project resonance.** Affine axis vs dimension axis mirrors Mode A (hypotenuse, fast,
   1-step) vs Mode B (leg removal, dimensional descent) on the staircase; and "truncate a
   Vandermonde/binomial tail, certify what survives" is precisely the LRC(14) density-route
   pattern (explicit-tail transfers). Same algebraic move, opposite goals: there we bound
   the tail, here the tail IS the object.

## 7. Errata to the prompt (for the record)

- First-diagonal pair "35 − 21 = 15" should read **36 − 21 = 15** (C(9,7) = 36, octagonal
  P(8,3) = 21); the pattern (triangular differences) is correct as stated.
- The polygonal skip-row sums as listed (1,1,2,3,5,8,13,21,33,51,76,111,157,218) are
  confirmed exactly.
