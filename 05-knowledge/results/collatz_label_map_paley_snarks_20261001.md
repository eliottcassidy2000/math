# Collatz XX: the owner's label map, a Paley-coded restatement of Collatz, and snarks

**opus S15, twentieth note, 2026-10-01.** Worktree `collatz-functional-uniqueness-20261001`.

- Script: `04-computation/experiments/collatz_label_map_snarks_20261001.py` (+ `.out`, ALL CHECKS PASSED).
- Companion note (the tournament part of the same prompt):
  [`shaved_tournaments_unavoidable_cores_20261001.md`](shaved_tournaments_unavoidable_cores_20261001.md).
- Canon: THM-4527 (Theorems 1–2 and Proposition 3).

**Status.**
- PROVED (elementary): Theorems 1–2, Propositions 3–5.
- FINITE-EXACT: all tables.
- **NO PROOF of Collatz.** The Paley bridge reformulates Collatz exactly but cannot prove it (§2.3).
- Snarks: a dictionary plus Proposition 5 (a re-derivation of Máčajová–Škoviera's "all seven points and at
  least four lines").
- Independent audit: OWED.

The owner's fourth prompt asked for four things:
1. consider the family of snarks, in view of the Paley symmetry split (the Frobenius survives, the
   translations do not);
2. try to prove Collatz with the Paley parity bridge;
3. study the corrected map `F(2N) = 3N`, `F(2N−1) = 2N−1−K_N` and the structure of `K` (the guess "two copies
   of each difference" was "a shot in the dark");
4. the shaved 4-tournament (companion note).

| Question | Answer | Label |
|---|---|---|
| Is the corrected `F` Collatz? | Yes: `F` is the Syracuse map `S(A) = oddpart(3A+1)` in labels `M = (A+1)/2` | PROVED (Thm 1) |
| Structure of `F` | `F(2N) = 3N`, `F(4j+1) = 3j+1`, `F(4n−1) = F(n)`; the last rule is an exact microcosm | PROVED (Thm 1) |
| Two copies of each difference? | False in `K` (mean multiplicity 1.3437). True for ascents + shallow descents: the units `2^1−3 = −1` and `2^2−3 = 1` | PROVED (Thm 2) |
| Multiplicity of `m` in `K` | `#{v ≥ 2 : (2^v − 3) ∣ 6m+1}`, from `6K + 1 = (2^v − 3)·S(A)` | PROVED (Thm 2) |
| Collatz via the Paley bridge | Exact restatement: Collatz ⟺ `7·Φ_T(n) ∈ Z` for all `n ≥ 1`, and then `7Φ_T(n) ≡ −2^σ(n) ∈ NQR_7` | PROVED (Prop 3); no proof of Collatz |
| Why the bridge stalls | It is the `q = 1` shadow; the obstacle is the resonance `2^L ≈ 3^k` and the sign | PROVED (Prop 4) + typed |
| Snarks | 3-edge-colouring = colouring by the single Paley line `{1,2,4}`; snarks need at least 4 translates of it | PROVED (Prop 5) + DICTIONARY |

## 1. The owner's map F

### 1.1 Decoding

Write `S(A) = oddpart(3A+1)` for odd `A`, and label odd numbers by `M = (A+1)/2 ∈ {1, 2, 3, …}`.

**Theorem 1.**
- (a) `F(M) := (S(2M−1) + 1)/2 = (oddpart(3M−1) + 1)/2`. It satisfies `F(2N) = 3N` and
  `F(2N−1) = 2N−1−K_N`, with `K_N = (A − S(A))/2` at `A = 4N−3`. This reproduces the owner's
  `K = 0,2,1,4,2,10,3,9,4,15,5` and the examples `S(1) = 1`, `S(5) = 1`, `S(9) = 7`, `S(13) = 5`. So `F` is
  conjugate to the Syracuse map, and **Collatz ⟺ every `M ≥ 1` reaches 1 under `F`**.
- (b) Three rules define `F` completely:
  - `F(2N) = 3N`;
  - `F(4j+1) = 3j+1`;
  - `F(4n−1) = F(n)`, **the microcosm**.

  In Syracuse language the last rule is `S(4A+1) = S(A)`. The quarter `M ≡ 3 (mod 4)` of the labels is an
  exact copy of the whole map, rescaled by 4. Its fixed point is `M = 1/3`, and the chain `1, 3, 11, 43, …` is
  the trunk `(4^j − 1)/3` (all mapped to 1).
- (c) `K_{2j+1} = j`, `K_{4m} = 5m − 1`, `K_{4m−2} = K_m + 6m − 4`. So `K` is a 2-regular sequence, and its
  positions `≡ 2 (mod 4)` replay the whole sequence plus a linear term.

*Proof.*
- (a) `3(2M−1)+1 = 2(3M−1)`, so `S(2M−1) = oddpart(3M−1)`.
- (b) For `M = 4j+1` (`A = 8j+1`): `S(A) = oddpart(24j+4) = 6j+1`, whose label is `3j+1`. For `M = 2N`:
  `3M−1` is odd, so `F(2N) = (6N−1+1)/2 = 3N`. For `M = 4n−1`: `3M−1 = 4(3n−1)`, so
  `oddpart(3M−1) = oddpart(3n−1)` and `F(4n−1) = F(n)`.
- (c) Follows from (b), using `K_N = M − F(M)` at `M = 2N−1`. ∎

Checked for `M, N ≤ 10^6` (Check A).

### 1.2 The difference spectrum: where the "two copies" really live

Extend `K` to every odd `A` by `K(A) = (A − S(A))/2`. It is negative on the ascending branch `A ≡ 3 (mod 4)`.
Let `v = v_2(3A+1)`.

**Theorem 2.**
- (a) `6K(A) + 1 = (2^v − 3)·S(A)` for every odd `A`.
- (b) Every `m ∈ Z` occurs as `K(A)` exactly `#{v ≥ 1 : (2^v − 3) ∣ 6m+1 and (6m+1)/(2^v−3) > 0}` times:
  - every negative `m` exactly once (`v = 1`, the ascent `F(2N) − 2N = N`);
  - every `m ≥ 0` exactly `#{v ≥ 2 : (2^v − 3) ∣ 6m+1}` times.

  In particular **the owner's `K` (the descents, `A ≡ 1 mod 4`) contains `m` exactly
  `#{v ≥ 2 : (2^v − 3) ∣ 6m+1}` times.** The mean multiplicity is
  `Σ_{v≥2} 1/(2^v − 3) = 1 + 1/5 + 1/13 + 1/29 + 1/61 + … = 1.34367…`, not 2.
- (c) The divisors `2^v − 3 = −1, 1` (`v = 1, 2`) are the only units, and they give **one ascent and one
  shallow descent of every size, plus the single 0 at the fixed point `M = 1`**. That is exactly "two copies
  of each element of N plus one 0", once the ascents are counted. The extra copies are the divisors
  `5, 13, 29, 61, 125, 253, 509, …` of `6m+1`. They are the deep descents, which live in the microcosm
  quarter: `Δ(4n−1) = Δ(n) − 3n + 1` with `Δ(M) = F(M) − M`.

*Proof.*
- (a) `S(A) = (3A+1)/2^v`, so `6K + 1 = 3A − 3S + 1 = 2^v S − 3S`.
- (b) Fix `m` and `v`, and let `s = (6m+1)/(2^v−3)`. Then `s` is odd and must be positive. The value
  `A = (2^v s − 1)/3` is an integer because `2^v s ≡ (2^v − 3)s = 6m+1 ≡ 1 (mod 3)`. It is odd and positive,
  `v_2(3A+1) = v_2(2^v s) = v`, and `K(A) = m`. Conversely `(A, v)` determines `(m, v)`. ∎

**Data.** On `[0, 10^5]` the multiplicity histogram is `{1: 69619, 2: 26617, 3: 3549, 4: 210, 5: 6}` (Check
A). For example:
- 1 occurs once;
- 2 occurs twice (`13 = 2^4 − 3`);
- 24 occurs three times (`145 = 5·29`).

Over the labels `M ≢ 3 (mod 4)`, `|F(M) − M|` takes 0 once and every `m ≥ 1` exactly twice.

**Gersonides again.** The two copies come from `|2^v − 3| = 1`, the same equation `|2^a − 3^b| = 1` that made
the trivial cycle integral (nineteenth note).

**The steps of the progressions.** The steps `2^v − 3` are the denominators of the one-odd-step rational
cycles `x_v = 1/(2^v − 3)`. On the branch with exactly `v` halvings,
`A − S(A) = (2^v − 3)(A − x_v)/2^v`, so the descent measures the distance to that branch's fixed point
(EXPLAINED).

**`2^{K_1}, 2^{K_2}, 2^{K_3} = 1, 4, 2`.** The next values are `16, 4, 1024`. These are forced small values
(the fixed point, its microcosm image, the first shallow descent): **NUMEROLOGY**.

## 2. Trying to prove Collatz with the Paley parity bridge

### 2.1 An exact Paley restatement

Let `T` be `n ↦ n/2, 3n+1`, and let `Φ_T(n) = Σ_j (T^j(n) mod 2)·2^j ∈ Z_2` be its parity code. The trivial
cycle `1 → 4 → 2` has word `100`, so `Φ_T(1) = −1/7` (nineteenth note: its 2-adic code is `NQR_7`).

**Proposition 3.** Collatz holds iff `7·Φ_T(n) ∈ Z` for every `n ≥ 1`. In that case
`7Φ_T(n) = 7·(Σ_{j<σ} b_j 2^j) − 2^σ`, where `σ = σ(n)` is the number of steps to reach 1. So
`7Φ_T(n) mod 7 = −2^σ(n) mod 7 ∈ {6, 5, 3} = NQR_7`, according as `σ ≡ 0, 1, 2 (mod 3)`.

*Proof.*
- A 2-adic integer lies in `(1/7)Z` iff its bits are eventually periodic with period dividing 3.
- `T`-words never contain `11`. So the only possible tails are `0^∞` (impossible for `n ≥ 1`) and rotations
  of `(100)^∞`.
- The parity-vector map is injective (it refines the bijective `T_1` map of Lagarias), so tail `(100)^∞`
  means the orbit enters `{1, 4, 2}`.
- The residue is read off the tail: a tail `−c/7` with `c ∈ {1, 2, 4}` gives residue `−2^σ`. ∎

This is the `T`-version of Lagarias' `Φ_{T1}(N) ⊆ (1/3)Z`.

**Checks (B).**
- For `n ≤ 20000`, the formula matches the actual parity bits to 40 bits past the entry into 1. The residues
  split `{3: 6665, 5: 6715, 6: 6620}`.
- Negative integers have tails of periods 2, 5, 18 (cycles through −1, −5, −17), so their codes have
  denominators `3, 31, 2^18 − 1`.

### 2.2 The bridge is the q = 1 shadow

In the family `x ↦ x/2, (qx+1)/2`, a cycle with parity word `w` (length `L`, `k` ones) sits at
`x_w(q) = c_w(q)/(2^L − q^k)`, where `c_w(q) = Σ_t b_t q^(#ones after t) 2^t`.

**Proposition 4.**
- (a) At `q = 1` the cycle points are `c_w(1)/(2^L − 1)`: the binary readings of the words, i.e. the Paley
  codes. **`QR_7/7` and `NQR_7/7` are the two rational 3-cycles of the `q = 1` map** (words `100` and `110`).
- (b) At `q = 1`, Collatz is trivially true. `c_w(1) < 2^L − 1` for every non-constant word, so no
  non-constant word gives an integer cycle, and the only integer cycles are `{0}` and `{1}`.
- (c) At `q = 3` the same size bound reads `x_min ≥ (3^k − 2^k)/(2^L − 3^k)`. It stops excluding cycles exactly
  near the resonances `2^L ≈ 3^k`. The positive-sheet pairs where the crude bound exceeds 1 are
  `(L, k) = (5,3), (7,4), (8,5), (10,6), (12,7), (13,8), (15,9), (16,10), (18,11), (20,12), (21,13), (23,14), …`
  (Check B).

The Paley bridge carries no information about these resonances, because the codes do not depend on `q`.

**What closes the resonant cases is Diophantine and archimedean.** Steiner (1977) used Baker's theory to
exclude 1-cycles. Simons and de Weger (2005) excluded `m`-cycles for `m ≤ 68`. Hercher (J. Integer Seq. 26,
2023, Article 23.3.5) proved there are no `m`-cycles with `m ≤ 91`. Every integer `T1`-cycle with period
`L ≤ 16` is one of the five known (`0`; `1`; `−1`; `−5`; `−17`), checked over all `2^L` words (Check B).

### 2.3 Why the bridge cannot finish

1. **It is a tail condition read 2-adically.** Proposition 3 is about the tail of `Φ_T(n)`. Any argument
   that looks at finitely many bits of codes (a cylinder condition) treats a positive `n` and a negative
   `n' ≡ n (mod 2^K)` alike. But the codes of negative integers are not in `(1/7)Z`. This is the
   eighteenth note's drift- and sign-blindness, in parity-code form.
2. **It does not see `q`.** The bridge is identical for `3x+1`, for `5x+1` (cycles 1, 13, 17), and for the
   `q = 1` map. The cycle condition `(2^L − 3^k) ∣ c_w(3)` does see `q`, and that is where Collatz lives.
3. **It does not extend past period 3.** A single ×2-orbit is a Paley set only when `L = 3` (nineteenth note,
   audited).

**Verdict: NO PROOF.** A proof through parity codes would need an archimedean input: Baker-type bounds for
the cycles, and a drift or measure statement for divergence. The Paley structure supplies neither.

## 3. Snarks

### 3.1 Dictionary

Use Singer coordinates: the points of the Fano plane are `α^i` (`i ∈ Z_7`) in `F_8 = F_2[α]/(α³+α+1)`.
- **The 7 Fano lines are the translates `D + t` of the Paley set `D = {1,2,4}`**, since
  `α^(1+t) + α^(2+t) + α^(4+t) = 0`.
- The Paley tournament's out-neighbourhoods are exactly these lines: `N^+(x) = x + D`.
- The Frobenius `x ↦ x²` fixes `D` and rotates its points `1 → 2 → 4`; this is the trivial cycle's own
  dynamics.
- A **Fano colouring** of a cubic graph labels the edges by points so that every vertex star is a line. It is
  the same thing as a nowhere-zero `Z_2^3`-flow, so every bridgeless cubic graph has one (Jaeger's 8-flow
  theorem).

### 3.2 What snarks need

**Proposition 5.** A Fano colouring is a 3-edge-colouring in disguise if it uses
- (i) one line,
- (ii) two lines,
- (iii) any number of concurrent lines, or
- (iv) colours avoiding some point `p`.

Hence every Fano colouring of a snark uses all 7 points and at least 4 lines. With exactly 4 lines, they form
**a pencil (three lines through one point) plus one line off it**. The other 4-line type, the four lines
missing a point, avoids that point and so falls under (iv).

*Proof.*
- (i)–(iii): if every star contains the point `a`, the `a`-edges form a perfect matching. On the complementary
  2-factor, neighbouring vertices share an edge whose colour lies in both lines minus `a`. Two lines through
  `a` meet only in `a`, so each cycle keeps one line and alternates its two other points. Hence every cycle is
  even, which gives a 3-edge-colouring.
- (iv): project `F_2^3 → F_2^3/⟨p⟩ ≅ F_2^2`. This gives a nowhere-zero `Z_2^2`-flow, i.e. a 3-edge-colouring.
- Three non-concurrent lines cover only 6 points. ∎

This re-derives Máčajová–Škoviera (Theoret. Comput. Sci. 349 (2005) 112–120). They also prove that 6 lines
always suffice and conjecture that 4 do.

### 3.3 Computations (Check C)

- **Petersen:** exactly **28560** nowhere-zero `Z_2^3`-flows. This equals the flow polynomial at 8 and is
  `170 · |GL(3,2)|`. Each one uses all 7 points. By number of lines:

  | Lines used | Flows |
  |---|---|
  | 4 (pencil + line) | 3360 |
  | 5 | 10080 |
  | 6 | 10080 |
  | 7 | 5040 |

- **Flower snark `J_5` and the two 18-vertex dot products `P·P`** (non-isomorphic: 10 versus 8 five-cycles,
  so these are the two Blanuša snarks): each has no Fano colouring of types (i)–(iv) or the 4-line
  quadrilateral, and each is colourable with a pencil + line.

### 3.4 In the owner's language (DICTIONARY, no implication)

- A 3-edge-colouring is a colouring by the single Paley line `{1,2,4}`. Its colours are rotated by the
  Frobenius, the trivial cycle's dynamics.
- **A snark is exactly a bridgeless cubic graph that cannot do without the translations.** It needs at least
  4 translates of `{1,2,4}`, three of them through a common point.
- So the 7 translations that "have no Collatz counterpart" (nineteenth note) are precisely what snarks force.
- Nothing transfers between the problems.

### 3.5 Collatz graphs are the opposite of snarks

The Collatz graph is a functional graph, so each component is a tree or unicyclic.
- Every non-cycle edge is a bridge.
- With maximum degree 3, it is trivially 3-edge-colourable. (On the 4358 vertices reached from `[1, 2000]`
  the maximum degree is 3.)

Snarks are the opposite: bridgeless, colour-resistant, and connectivity-rich.

The genuine shared pattern is the **minimal-counterexample reduction**:
- A minimal counterexample to the cycle-double-cover or 5-flow conjecture must be a snark.
- A minimal Collatz counterexample is the minimum of its component (eighteenth note).
- But the four-colour-theorem route (a finite unavoidable set of reducible configurations) is closed for
  Collatz. The residue classes mod `2^k` with no coefficient descent in `k` steps have positive density for
  every `k`: `1.25e−1, 6.25e−2, 2.6e−2, 5.8e−3, 6.6e−4` at `k = 5, 10, 20, 40, 80`. The class of `−1` always
  ascends.

## 4. Directions

- **D73.** Is there a Fano/Paley-coded invariant of a cubic graph that sees the resonance `2^L ≈ 3^k`? (Most
  likely not; it is listed only to close the loop.)
- **D74.** The 2-regular sequence `K`: closed form of its summatory function, and the distribution of
  `#{v : (2^v − 3) ∣ 6m+1}` (a divisor problem over `2^v − 3`).
- **D75.** An exact `q`-interpolation. For which odd `q` does the size bound of Proposition 4 fail for
  infinitely many `(L, k)`, and how does that track the known cycles of `qx+1`?

## 5. Reproduction

```bash
python 04-computation/experiments/collatz_label_map_snarks_20261001.py
```

This takes about 3 minutes.
