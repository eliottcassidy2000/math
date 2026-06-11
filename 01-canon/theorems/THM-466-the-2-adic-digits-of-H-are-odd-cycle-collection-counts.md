# THM-466 — The 2-adic digits of H are odd-cycle collection counts (the higher-Rédei tower)

**Type:** Theorem
**Status:** PROVED (parts (i)–(iii), all n; one-line consequences of THM-002/OCF) + VERIFIED
(part (ii) machine-checked with 0 failures on ALL 2,130,016 labeled tournaments n = 3..7 and
200,000 random tournaments at n = 8, by an independent Held-Karp path DP vs.
directed-cycle censuses). The score-determination statements (iv)–(v) are
PROVED-by-counterexample / classical respectively, with full-census verification ranges stated inline.
**Source:** kind-pasteur-2026-06-10-S2, Thread B lab (HYP-2379; upgrades tangent T007).
**Builds on:** THM-002 (OCF, Grinberg-Stanley arXiv:2307.05569; Irving-Omar arXiv:2412.10572), THM-001 (Rédei).
**Script:** `04-computation/ocf_digit_tower_kpo2.py` → `05-knowledge/results/ocf_digit_tower_kpo2.out`
**Tags:** #ocf #2-adic #higher-redei #digit-tower #score-sequences #t007

---

## Statement

Let T be a tournament on n vertices, H = H(T) its number of directed Hamiltonian paths, and
α_k = α_k(Ω(T)) the number of collections of k pairwise vertex-disjoint **directed odd cycles**
of T (α_0 = 1), i.e. the independence numbers of the conflict graph Ω(T) of THM-002.

**(i) Digit lemma (the higher-Rédei tower).** For every m ≥ 1:

```
H  ≡  Σ_{k<m} α_k 2^k   (mod 2^m).
```

H mod 2^m is determined by the collections of **fewer than m** disjoint odd cycles. The
residues of H at all powers of 2 — equivalently the 2-adic digits of H — read off the
disjoint-odd-cycle collection counts.

**Corollaries of (i).**
- (m = 1) H ≡ α_0 = 1 (mod 2): **Rédei's theorem** is the ground floor of the tower.
- (m = 2) H ≡ 1 + 2·α_1 (mod 4), where α_1 = c_3 + c_5 + c_7 + … is the **total number of
  directed odd cycles**.
- (valuation — the answer to T007) v_2(H) = 0 always; the right object is
  ```
  v_2(H − 1) = 1 + v_2(α_1 + 2α_2 + 4α_3 + …),
  ```
  with v_2(H − 1) = ∞ iff T is transitive (α_k = 0 ∀k ≥ 1 ⟺ no 3-cycle ⟺ transitive).
  In particular **v_2(H − 1) = 1 ⟺ α_1 is odd** (odd number of directed odd cycles).

**(ii) Finite identity.** A collection of k pairwise disjoint odd cycles uses ≥ 3k vertices, so
α_k = 0 for k > ⌊n/3⌋. Hence the tower terminates and gives **exact** identities:

```
n ≤ 5:  H = 1 + 2α_1                            (single cycles only)
n ≤ 8:  H = 1 + 2α_1 + 4α_2                     (singles and disjoint pairs)
```

Pair types in α_2: at n = 6,7 only (3,3); at n = 8 both (3,3) and (3,5). (7-cycles enter α_1
at n ≥ 7; the first triple-collection channel (3,3,3) opens at n = 9.)

**(iii) Reversal invariance.** Arc reversal T ↦ T^op induces a graph isomorphism
Ω(T) ≅ Ω(T^op). Hence α_k(T) = α_k(T^op) for every k, I(Ω(T), x) = I(Ω(T^op), x) as
polynomials, and **every 2-adic digit of H is an invariant of the merged class {T, T^op}** —
a node invariant of the merged metagraph G_n/Z_2.

**(iv) The digits are NOT score functions (n ≥ 5).** α_1 mod 2 — equivalently H mod 4 —
is **not** determined by the score sequence. Smallest counterexample: n = 5 (explicit pair
below). Consequently v_2(H − 1) is not a score function either. Non-constant score classes:
1/9 (n = 5), 5/22 (n = 6), 28/59 (n = 7) — full labeled censuses.

**(v) H is score-determined exactly up to n = 4.** For n ≤ 4 the only odd cycles are
3-cycles, so H = 1 + 2c_3, and c_3 is score-determined (Kendall–Babington Smith:
c_3 = C(n,3) − Σ_v C(s_v,2)). For n = 5 the pair in (iv) has equal scores and H = 11 vs 13.

---

## Proofs

**(i).** THM-002 (PROVED, Grinberg-Stanley): H = I(Ω(T), 2) = Σ_k α_k 2^k. Every term with
k ≥ m is divisible by 2^m. ∎

(Valuation corollary: H − 1 = 2(α_1 + 2α_2 + 4α_3 + …) since α_0 = 1, so
v_2(H−1) = 1 + v_2(α_1 + 2α_2 + …), infinite iff all α_{k≥1} = 0, which happens iff T has no
directed 3-cycle, i.e. T transitive — a tournament with any directed cycle has a directed
3-cycle.)

**(ii).** Directed odd cycles have ≥ 3 vertices and the cycles in a collection are pairwise
vertex-disjoint, so a k-collection occupies ≥ 3k vertices: α_k = 0 for 3k > n. For n ≤ 8:
⌊n/3⌋ ≤ 2; for n ≤ 5: ⌊n/3⌋ ≤ 1. The pair-type list is the list of (a,b), 3 ≤ a ≤ b odd,
a + b ≤ n. ∎

**(iii).** C = (v_1 → v_2 → … → v_{2j+1} → v_1) is a directed cycle of T iff the reversed
cyclic sequence (v_1 → v_{2j+1} → … → v_2 → v_1) is a directed cycle of T^op, on the **same
vertex set** and of the same (odd) length. The map C ↦ C^op is a bijection on directed odd
cycles which preserves vertex sets, hence preserves the shares-a-vertex relation: it is an
isomorphism Ω(T) → Ω(T^op). Independent sets correspond, so α_k(T) = α_k(T^op). ∎
(This refines H(T) = H(T^op): not just the count but the entire conflict graph, hence the
entire polynomial I(Ω, x) and every digit, is reversal-invariant.)

**(iv).** Proof by explicit counterexample (machine-verified by brute-force permutation count
AND independent Held-Karp DP; codes in the lexicographic pair-bit convention of the script).
Both tournaments below have score sequence (1,2,2,2,3); they differ **only** by reversing the
single directed 3-cycle on {0,1,3} (a score-preserving move):

```
B (code 12):  0→{3,4}, 1→{0}, 2→{0,1}, 3→{1,2}, 4→{1,2,3}
              c3=4, c5=1, α1=5, α2=0, H=11 ≡ 3 (mod 4), v2(H−1)=1
A (code 41):  0→{1,4}, 1→{3}, 2→{0,1}, 3→{0,2}, 4→{1,2,3}
              c3=4, c5=2, α1=6, α2=0, H=13 ≡ 1 (mod 4), v2(H−1)=2
```

Reversing a cyclic triangle preserves c_3 (score-determined) but **can flip c_5 parity**. ∎

**(v).** At n ≤ 4 odd cycles are exactly 3-cycles (5 > n), so (ii) gives H = 1 + 2c_3, and c_3
is a function of the scores by the classical cyclic-triple identity. ∎

---

## Verification record (`ocf_digit_tower_kpo2.py`, runtime 14 s)

Identity H = 1 + 2α_1 + 4α_2 checked with H from an independent Held-Karp path DP and
α_1, α_2 from directed-cycle censuses (rooted directed-Hamiltonian-cycle enumeration per odd
subset — MISTAKE-023-compliant; e.g. one 5-set can carry 3 directed 5-cycles, max observed
LUT_5 value = 3; max directed 7-cycles on 7 vertices = 24):

| n | tournaments | coverage | failures | extra exact checks (all passed) |
|---|-------------|----------|----------|-------------------------------|
| 3 | 8 | all labeled | 0 | ΣH, Σc3 exact; #(H=1) = 3! |
| 4 | 64 | all labeled | 0 | ΣH, Σc3 exact; α2 ≡ 0; #(H=1) = 4! |
| 5 | 1,024 | all labeled | 0 | ΣH, Σc3, Σc5 = C(5,5)·24·2^5 exact; α2 ≡ 0 |
| 6 | 32,768 | all labeled | 0 | + Σα2 = 10·2^11 exact |
| 7 | 2,097,152 | all labeled (covers all 456 iso classes) | 0 | + Σc7 = 720·2^14, Σα2 = 70·2^17 exact |
| 8 | 200,000 | random sample | 0 | (3,5)-pair channel active in 174,495 samples; means c3/c5/c7/H match theory 14/42/45/315 |

12+12+12+8 random codes at n = 4..7 and 4 at n = 8 additionally cross-validated against a
pure-Python brute force (explicit directed-odd-cycle lists + n!-permutation path scans).
Rédei (H odd) machine-asserted on every tournament computed.

**Distribution of v_2(H − 1), full labeled censuses (T007 data):**

| n | v2=1 | v2=2 | v2=3 | v2=4 | v2=5 | v2=6 | v2=7 | ∞ (transitive) |
|---|------|------|------|------|------|------|------|----------------|
| 4 | .2500 | .3750 | — | — | — | — | — | .3750 |
| 5 | .2969 | .3516 | .2344 | — | — | — | — | .1172 |
| 6 | .3750 | .3223 | .1562 | .0439 | .0806 | — | — | .0220 |
| 7 | .4479 | .2573 | .1442 | .0673 | .0393 | .0264 | .0152 | .0024 |

P(v_2(H−1) = 1) = P(α_1 odd) climbs 0.250 → 0.297 → 0.375 → 0.448 (→ 1/2?).

**Score-class determinism (full labeled enumeration, sorted-score classes = A000571):**

| invariant | n=4 | n=5 | n=6 | n=7 |
|---|---|---|---|---|
| c3 mod 2 (classical) | const | const | const | const |
| α1 mod 2 (= H mod 4) | const | **1/9 non-const** | 5/22 | 28/59 |
| c5 mod 2 | const (0) | 1/9 | 5/22 | 24/59 |
| c7 mod 2 | — | — | — (0) | 20/59 |
| (c5+c7) mod 2 | const | 1/9 | 5/22 | 28/59 |
| α2 mod 2 | const (0) | const (0) | 7/22 | 35/59 |
| v2(H−1) full | const | 1/9 | 7/22 | 34/59 |

The natural sharper conjectures are also REFUTED: c5 (and c5+c7) is NOT always even
(n = 5: 400/1024 tournaments have c5 odd), so α_1 ≢ c_3 (mod 2) in general.

---

## Honest scope

- Parts (i)–(iii) are proved for **all n** (they are formal consequences of THM-002, which is
  proved for all n). No caveat.
- The **computational** verification of (ii) covers n ≤ 8. The α_3 channel (digit
  contributions from triples of disjoint odd cycles) first opens at n = 9 and is untested
  territory: for 9 ≤ n ≤ 11, H = 1 + 2α_1 + 4α_2 + 8α_3 exactly. Next rung: sampled n = 9.
- Score-determination results (iv) are exhaustive for n ≤ 7; at n = 8 sampled evidence only
  (114/166 sampled score classes show both α_1 parities — existence, not exhaustive).
- MISTAKE-028/036/055 note: the "P(v_2=1) → 1/2" reading of the distribution table is a
  CONJECTURE from four data points, not a law.

## T007 upgrade note

Tangent T007 asked (pre-proof) for the 2-adic tower of higher Rédei theorems and a
combinatorial characterization of the 2-adic valuation of H. Both are now closed: the m-th
floor of the tower is H ≡ Σ_{k<m} α_k 2^k (mod 2^m), and v_2(H−1) = 1 + v_2(α_1 + 2α_2 + …).
The valuation is computable from ≤ ⌊n/3⌋ collection counts, is a merged-class (T ↔ T^op)
invariant, and is **not** a score function for n ≥ 5.

## Key references

- Grinberg & Stanley, arXiv:2307.05569 (Rédei-Berge symmetric function; proves OCF).
- Irving & Omar, arXiv:2412.10572 (Corollary 20, matrix-algebra derivation).
- Kendall & Babington Smith 1940 (c_3 is score-determined — the ground floor's digit IS a
  score function; THM-466(iv) shows the tower's higher digits are not).
