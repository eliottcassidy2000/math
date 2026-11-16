# The seed monoid: order-9 banks, banks that keep a line of exteriors alive, arc banks, and the Fibonacci-box descent

Session `mac-mini-2026-11-15-monoid`, 2026-11-15. **The Chan–Pak conjecture is OPEN; nothing here proves it.**
Not independently audited (owner directive of 2026-11-13). Builds on THM-5282 (one head over tape–block–tape).

`H` = number of Hamiltonian paths. A *seed* is a block `C` with a pattern `X` of contacts of one extra vertex (the
head). Attached paper read for this session: H. Lee Price, *The Pythagorean Tree: A New Species*, arXiv:0809.4324.

## 0. Results at a glance

| # | statement | status |
|---|---|---|
| 1 | **Order 9.** Among all 259,280 seeds on the 786 blocks of order 3–9 with count `h <= 45`, no synchronised bank is longer than its block count. Longest by `h`: 2 (`h = 3`), 3 (5), 4 (9), 3 (11), 5 (15), 4 (17), 5 (33), 6 (45). | FINITE-EXACT |
| 2 | **Seed monoid.** Over `R = Z[e]/(e² - 1)` a seed is the matrix `M = [[h, b, c], [0, G, a], [0, 0, h]]`; ordinal sum of blocks is matrix product and `H = h + 2·(even part of c)` for every chain. A tape vertex is a one-vertex seed. A synchronised bank is a set of seeds that differ only in the corner; composition multiplies corner steps by the other count. | PROVED; checked on 40 products and 30 chains |
| 3 | **Fibonacci-box descent.** On a tape the character `e = 1` acts on `(b, c)` by `(2b, c + b)` (non-contact) and `(2b + 1, c)` (contact): Price's rule "one entry doubled, the other bumped". The inverse halves the address at every level and subtracts when it is even; `H = 1 + c⁺ + c⁻`, `c⁻` = number of factors `10`. Price's two other generators double the *count* coordinate, which tournaments can do only at ports. | PROVED (identities); the comparison is exact arithmetic |
| 4 | **Line banks (compensated exteriors).** Seeds with the same `(h, G)` whose signatures `(a, b)` lie on a line of direction `(p, -q)` differ by constants for *every* exterior on the line `pα - qβ = ρ`. Along such lines there are banks of **5** consecutive counts (3 without compensation), e.g. `h = 41`, `α = 2β`: counts 179..187, 14583..14591, 51771..51779, 111579..111587. | PROVED (criterion) / FINITE-EXACT (lengths; literal check) |
| 5 | **Arc banks.** Two arcs can be independent digits of weights 2 and 4 (counts 9, 11, 13, 15, affine). No three arcs of any tournament of order at most 7 give eight consecutive counts (606,480 cubes); the nearest is 9, 11, 13, 15, 17, 19, 23, 27. There is no congruence obstruction: 1,012 cubes permute the odd residues mod 16. | FINITE-EXACT |
| 6 | Readings of Chair44 and of `{4, 2, 5, 3, 7, 13}` (section 6). | ANALOGY / NUMEROLOGY, flagged |

## 1. Order 9

`scripts/cross_bank.c` with `HMAX=45` (`out/cross_bank_upto9_h45.txt`). A bank longer than its count needs a small
count, so blocks with `h > 45` (198,160 of order 9) are skipped. With one tape the best ratios are 2/3 and 3/5;
with two arbitrary tapes the longest bank is still 3 (now also at `h = 9`, key `(9, 672, 672, 127, 118, 218, 209)`).
The self-reproducing regime `L >= 2h` of THM-5282 (D) is not reached. HYP-9832 (a) stays open.

## 2. The seed monoid

For a seed `(C, X)` put `h = H(C)` and, with `e` marking an odd number of path vertices,

    G = G_e + G_o e,   a = Σ_{s ∉ X} (λ_e[s] + λ_o[s] e),   b = Σ_{s ∈ X} (κ_e[s] + κ_o[s] e),   c = Σ_{u ∈ X, v ∉ X} (W_e(u,v) + W_o(u,v) e)

(notation of THM-5282; `W_e` is the `W` of the cut form). Then

    M(C₁ → C₂, X₁ ∪ X₂) = M(C₁, X₁) · M(C₂, X₂),        H(chain + head) = h + 2·[c]_even.

*Proof.* A path of `C₁ → C₂` is a path of `C₁` followed by a path of `C₂`; weights and parities multiply. The four
entries of the product are the four ways a path can start and end relative to the two blocks. ∎

One vertex: `M₀ = [[1,0,0],[0,1+e,e],[0,0,1]]` (non-contact), `M₁ = [[1,e,0],[0,1+e,0],[0,0,1]]` (contact).
`(1+e)² = 2(1+e)` is where the powers of two come from. The master formula of THM-5282 is the product
(upper word)·`M`·(lower word). In this language:

* a **synchronised bank** is a set of seeds equal except for `[c]_even`: the corner is the low digit;
* for a product, the corner is `h₂c₁ + b₁a₂ + h₁c₂`: steps are multiplied by the other count, and the cross term
  `b₁a₂` is the carry between the two factors;
* two unit-count factors commute up to the corner exactly when `a₁(G₂ - 1) = a₂(G₁ - 1)` and
  `b₁(G₂ - 1) = b₂(G₁ - 1)`, and then the corner changes by `b₁a₂ - b₂a₁`; otherwise the order changes the middle
  entries as well.

## 3. Fibonacci boxes

Price's boxes `[q q'; p p']` (`q' = p - q`, `p' = p + q`) carry a ternary tree on coprime pairs of opposite
parity with children `(2q, p+q)`, `(p-q, 2p)`, `(p+q, 2p)` and the unique descent *halve the even entry, replace
the other by the positive difference*; its binary subtree is the Pythagoras family `B(n) → B(2n), B(2n+1)`.

Under `e = 1` the tape letters act on `(b⁺, c⁺)` (first row of the product) by

    non-contact:  (b, c) ↦ (2b, c + b),          contact:  (b, c) ↦ (2b + 1, c),

and under `e = -1` the corner counts the factors `10`. So for a one-apex word with binary value `b`,

    H = 1 + Σ_{k : bit k of b is 0} ⌊b / 2^(k+1)⌋ + #{factors 10},

e.g. `11101001000`: `(b⁺, c⁺) = (1864, 1725)`, `c⁻ = 3`, `H = 1729`. The first letter rule is Price's first
generator; the descent is "halve the address, and when it was even subtract the half from the count". Price's
other two generators double the *second* coordinate. For tournaments that is the dyadic step `n ↦ 2n` of THM-5272,
available only at ports. No family of tournaments can realise all three generators on a pair whose second entry
is the half-count: Price's theorem would then make every half-count a count, including 3 and 10. In this picture
the obstruction to Chan–Pak sits in the two count-doubling generators. (The fixed out-degrees `0, 1, 1, 2` of the four-vertex atom below
are the Fibonacci 4-tuple before `[1, 1, 2, 3]`; no more than that is claimed.)

## 4. Banks that keep a line of exteriors alive

Take the lowest upper vertex a non-contact and the highest lower vertex a contact. Then (THM-5282)

    H = h(1 + 2N_up + 2N_low) + 2[c + aα + bβ + 2Gαβ],        a = a_o + a_e,  b = b_o + b_e,  G = G_e + G_o.

**Criterion (PROVED).** Let seeds `i = 0..L-1` share `(h, G)` and satisfy `a_i = a_0 + p t_i`, `b_i = b_0 - q t_i`.
For every exterior with `pα - qβ = ρ` their counts differ by `2(c_i - c_0 + t_i ρ)`. If these numbers are
consecutive, the `L` tournaments have `L` consecutive odd counts on the whole line of exteriors. The synchronised
case is `t_i = 0`; freezing one exterior is the case of a single point.

`scripts/line_bank.py` (`out/line_bank_h45.txt`), seeds of order 3–9 with `h <= 45`, directions with `p, q <= 4`:

| direction `(p, q)` | (0,0) | (1,0), (0,1) | (1,1) | (1,2), (2,1) | (1,3), (3,1) | (2,3), (3,2) | (1,4), (4,1) | (3,4), (4,3) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| longest bank | 3 | 5 | 5 | 5 | 5 | 4 | 5 | 5 |

Example (literal check, `scripts/line_bank_check.py`): block `0 5 17 7 107 79 143 63` with patterns 187, 219, 69,
37 and block `2 4 25 19 35 79 31 127` with pattern 159; `h = 41`, `G = 1424`, signatures
`(169, 1045) + t(1, -2)`, `t = 0, 6, 333, 327, 9`, cuts 69..73. On the line `α = 2β` the counts are 179..187 (no
tapes), 14583..14591, 51771..51779, 111579..111587; off the line (`α = β = 1`) they are scattered.

The trade is explicit: a synchronised two-tape bank has two free words and length 3; a line bank has one free
parameter and length 5. What compensation buys is length, at the price of one dimension of exteriors.

## 5. Arc banks and the four-vertex atom

The atom of `trees5_tournaments4_erdos_sos_20261005` (the path `1→2→3→4` with `1→4`, the two span-2 arcs free)
realises the four tournaments on four vertices with counts 1, 3, 3, 5: the two free arcs are independent
(second difference 0) but of equal weight. By THM-5272 two contacts of one head are never independent; two arcs
at different places can be. `scripts/arc_bank.c` searches all sets of `k` arcs of all tournaments of an order:

* `k = 2`: complete banks exist and are always the counts 9, 11, 13, 15 at orders 5 and 6 (20 and 40 cubes); 8 and
  16 of them are affine, i.e. two independent binary digits of weights 2 and 4;
* `k = 3`: none through order 7. Among the 178,784 cubes with eight distinct counts the least spread is 9
  (`9, 11, 13, 15, 17, 19, 23, 27`), against 7 for a bank;
* `k = 4`: none at orders 5 and 6.

Residues do not explain it: at order 7, 185,412 cubes have `H mod 4` balanced, 23,016 have each odd residue
mod 8 twice and 1,012 permute the odd residues mod 16. A third independent digit of weight 8 is what is missing,
for arcs as for contacts (HYP-9842).

## 6. Two readings

* **Chair44** (`recursive_interfaces_20261016b.md`): a period is divided by 2 at every level of the unique
  coarsening, so only zero survives. The tape descent of section 3 is the same mechanism with a finite end: the
  address `b` is halved at every level and the count is the sum of what is left at the even levels. What tournaments
  lack is the coarsening of the *count* (ANALOGY).
* **`{4, 2, 5, 3, 7, 13}`** is `{1 + d : d | 12}`. Its odd members are the counts with half-count dividing 6: 3, 5, 7,
  13, of which 7 is the missing one and `13 = 2·7 - 1`. These are also the levels `N` with `(N - 1) | 12`, for which
  `(η(τ)/η(Nτ))^(24/(N-1))` is a Hauptmodul with even exponent (Monster classes 2B, 3B, 4C, 5B, 7B, 13B; recalled,
  not verified here). For the Baby Monster I only recall that it arises from the centraliser `2·B` of a 2A
  element; I found no statement that singles out this set. No structural link to tournament counts is claimed
  (NUMEROLOGY).

## 7. Open

1. A third independent digit: three arcs or three contacts with eight consecutive counts (HYP-9842 (a)).
2. Line banks of unbounded length; planes of exteriors with banks longer than 3 (HYP-9842 (b)).
3. A realisation of a count-doubling generator on a family closed under the tape letters (HYP-9842 (c)).
4. Chan–Pak.

## Files

`scripts/`: `cross_bank.c` (env `HMAX`, `DUMP`), `line_bank.py`, `line_bank_check.py`, `seed_monoid.py`,
`arc_bank.c`, `master_formula_check.py` (shared helpers). `out/`: the outputs quoted above.
