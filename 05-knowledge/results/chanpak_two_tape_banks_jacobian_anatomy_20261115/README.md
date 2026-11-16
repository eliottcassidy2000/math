# One head over tape–block–tape: no three-contact bank, synchronised seed banks, and a two-tape bank whose leading coefficient is a free product

Session `mac-mini-2026-11-15-tapes`, 2026-11-15. **The Chan–Pak conjecture is OPEN; nothing here proves it.**
Not independently audited (owner directive of 2026-11-13). The Jacobian part of this session is in
[JACOBIAN.md](JACOBIAN.md).

`H` = number of Hamiltonian paths. A *base* is `TT_i → C → TT_j`: an upper transitive tape, a block `C` (any small
tournament), a lower tape, in ordinal sum. One new vertex `x` (the head) has contacts `U` on the upper tape,
`X` on `C`, `L` on the lower tape.

## 0. Results at a glance

| # | statement | status |
|---|---|---|
| 1 | **Master formula.** `H = h(1 + 2N_up + 2N_low) + 2[cut_C(X) + a_o(α+u₁) + a_e α + b_o(β+z₁) + b_e β + G_e((α+u₁)(β+z₁) + αβ) + G_o((α+u₁)β + α(β+z₁))]`, with every symbol explicit (section 1). | PROVED; brute-force check (60 random instances) |
| 2 | **No three-contact bank on these bases.** For every block of order at most 8 (all 7,410 classes), every three contacts in the block, every block context and **all tapes of all lengths**, the eight on/off patterns never give eight consecutive odd counts. Two-contact banks exist only at isolated tape data, never along a line. | FINITE-EXACT over blocks, exact in the tapes (linear algebra) |
| 3 | **Synchronised banks.** Seeds (block, pattern) with equal tape invariants `(h, G_e, G_o, a_o, a_e, b_o, b_e)` have counts differing by a constant for *every* pair of tape words. Longest banks through block order 8: 6 consecutive counts with an upper tape whose lowest vertex is free, 4 with any upper tape, **3 with two arbitrary tapes**. | PROVED (invariance) / FINITE-EXACT (lengths) |
| 4 | **An explicit two-tape bank.** For all integers `A, B >= 0` the three numbers `K + 90, K + 92, K + 94` are counts, `K = 432AB + 174A + 303B + 45 - 30(A mod 2) - 15(B mod 2) - 90(J(A) + J(B))`, `J` the pair energy of the binary digits. The leading coefficient `432AB` is a free product. | PROVED; literal check on 40 tape pairs (orders up to 15) |
| 5 | This one bank reaches 11.6 %, 14.3 %, 16.9 % of the odd numbers below `10^5, 10^6, 10^7`; with the upper tape alone it stays at 3.9 %. | FINITE-EXACT |
| 6 | **Composition law.** For blocks in ordinal sum, `h = h₁h₂`, `G = G₁ ⊗ G₂` (by parity) and `cut = h₂ cut₁ + h₁ cut₂ + b_e a_e + b_o a_o`: banks compose like a numerical semigroup with generators `h₁, h₂`, so a bank at least as long as its own block count would reproduce itself. No such bank is known (longest/`h`: 4/9, 6/45). | PROVED (40 random checks) / open question |

Not claimed: a selector for arbitrary targets, any density theorem, completeness of any family.

## 1. The master formula

Heights on the upper tape are counted from the block (height 1 is adjacent to `C`), depths on the lower tape
likewise. Put `Z` = lower non-contacts and

    α = Σ_{a ∈ U, height ≥ 2} 2^(height-2),   u₁ = [height 1 ∈ U],   β = Σ_{b ∈ Z, depth ≥ 2} 2^(depth-2),   z₁ = [depth 1 ∈ Z].

For the block: `h = H(C)`; for a directed path `Q` in `C` let `w(Q) = H(C - V(Q))`; then
`λ_par[s]` = sum of `w(Q)` over paths ending at `s` with `|Q|` of that parity, `κ_par[s]` the same for paths
starting at `s`, `G_par` the sum over all paths (the empty path is even), `W(a,b)` the sum over even paths from `a`
to `b`, and

    cut_C(X) = Σ_{a ∈ X, b ∉ X} W(a,b),   a_par(X) = Σ_{s ∉ X} λ_par[s],   b_par(X) = Σ_{s ∈ X} κ_par[s].

`N_up`, `N_low` are the one-apex half-counts of the two tape words (`N_up = α - J(U)`, `N_low = β - J(Z)`).

*Proof.* The cut form of THM-5272 (A) writes `H(B + x) - H(B)` as twice a sum over pairs (contact `a`, non-contact
`b`) of the even paths `a ⇝ b` weighted by the count of the complement. A path is an initial run in the upper
tape, a path in `C`, a run in the lower tape; the complement of a path keeps both tapes transitive, so its count
is `w(Q)`; a tape vertex at height `h` reaches the block through `2^(h-2)` runs of each parity (`h ≥ 2`). ∎

`G_e + G_o` is the number of covers of `C` by a path and a Hamiltonian path of the rest, the "two-slot path-cover
count" of `variable_reader_coverage_20261116/BLOCKS.md`; `scripts/master_formula_check.py` checks the formula.

## 2. No three-contact bank

In cut units the marginal of a contact `s` of the block (other bank contacts off) is

    Δ_s = m_s + κ_o[s](β + z₁) + κ_e[s] β - λ_o[s](α + u₁) - λ_e[s] α,

`m_s` block-internal, and the pair interaction `ω_st = W(s,t) + W(t,s)` does not involve the tapes. A complete
`k`-bank needs `Σ Δ_s y_s - Σ ω_st y_s y_t` to be a bijection of the cube onto `2^k` consecutive integers; for a
given `ω` the admissible `Δ` are a finite list, and each is a linear system in `(α, β)` with `u₁, z₁ ∈ {0, 1}`.
`scripts/cube_bank_tapes.py` solves all of them exactly.

| block order | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|
| classes | 2 | 4 | 12 | 56 | 456 | 6,880 |
| three-contact banks: isolated tape data / lines | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| two-contact banks: isolated tape data / lines | 0 / 0 | 20 / 0 | 42 / 0 | 66 / 0 | | |

So on these bases a complete three-contact bank does not exist at any count, and two-contact banks occur only at
isolated tapes, e.g. the strong 4-block with contacts 0, 3 and tape data `(α, β) = (3, 4)`. By THM-5272 (B) a
three-contact bank needs its contacts within a bounded window, so "tape–block–tape" is the general shape of a
base with long transitive parts; what is not covered is a base with several separated blocks. HYP-9822 (a) stands.

**Why a cube cannot follow the tapes.** `Δ_s` has slope `(κ[s], -λ[s])` in `(β, α)`, with `κ[s] = U_C(s)` and
`λ[s] = U*_C(s)` the port quantities of THM-5262, both positive. A single contact always drifts as the tapes
grow; only *swaps* that keep `Σ_{s ∈ X} κ[s]` and `Σ_{s ∉ X} λ[s]` fixed are tape-independent. Banks that survive
a varying leading coefficient are therefore made of patterns with equal signatures, not of independent contacts.

## 3. Synchronised banks

By the master formula two seeds `(C, X)`, `(C', X')` with the same `(h, G_e, G_o, a_o, a_e, b_o, b_e)` give counts
differing by `2(cut - cut')` for all tapes. `scripts/sync_bank.c` (patterns on one block) and `scripts/cross_bank.c`
(all 1,823,696 seeds on blocks of order 3–8) find the longest runs of consecutive cut values:

| setting | key | one block | all seeds |
|---|---|---:|---:|
| upper tape, height 1 free | `(h, a_o + a_e)` | 3 | **6** (`h = 45`) |
| upper tape, any word | `(h, a_o, a_e)` | 3 | 4 (`h = 9`) |
| two tapes, height 1 and depth 1 free | `(h, G, a, b)` | 2 | 3 (`h = 9`; 40 such banks) |
| two tapes, any words | all seven | 2 | **3** (`h = 45`; 34 such banks, 404 of length 2) |

The one-tape rows are the regime of `variable_reader_coverage_20261116` (three- and five-count banks); the
two-tape rows are new.

## 4. The explicit two-tape bank

Seeds: block masks `22 44 24 49 34 5` with head pattern 23, and `38 28 24 49 33 6` with patterns 39 and 23
(bit `j` of mask `i` = arc `i → j`; bit `s` of the pattern = the head beats vertex `s`). All three have
`h = 45`, `G_e = G_o = 216`, `(a_o, a_e) = (72, 57)`, `(b_o, b_e) = (144, 114)` and cut values 45, 46, 47
(`H(C + x) = 135, 137, 139`). With `A = 2α + u₁` (the upper contacts read as a binary number, height 1 the unit
bit) and `B = 2β + z₁` (the lower non-contacts), the master formula collapses to row 4:

    K = 432AB + 174A + 303B + 45 - 30(A mod 2) - 15(B mod 2) - 90(J(A) + J(B)),        counts K + 90, K + 92, K + 94.

`K ≡ 3 (mod 6)`, so the three counts are the three odd residues mod 6: one full turn of a low digit in radix 6,
for every pair of binary words. `scripts/two_tape_bank_check.py` builds the tournaments literally.

**What the second tape buys.** With one tape the number of presentations below `X` is at most a constant times
`X` (the capacity bound of `LIMITS.md`); here `K ≈ (432A + 303)(432B + 174)/432`, so there are about
`(X ln X)/432` windows below `X`. `scripts/bank_coverage.py`: 5,787 / 71,340 / 843,577 distinct odd numbers below
`10^5 / 10^6 / 10^7` (6,117 / 77,160 / 931,405 presentations, so collisions cost under 10 %), against a constant
3.9 % with the upper tape alone. The exact low-digit turn is kept while the leading coefficient runs over all
products; what is still missing is the selector: solving `K(A, B) = N'` for a prescribed `N'`, where the pair
energies `J` are the carry.

## 5. Composition

For `C = C₁ → C₂` and patterns `X₁, X₂`: `h = h₁h₂`, `G_e = G₁e G₂e + G₁o G₂o`, `G_o = G₁e G₂o + G₁o G₂e`, and

    cut_C(X₁ ∪ X₂) = h₂ cut_{C₁}(X₁) + h₁ cut_{C₂}(X₂) + b_e(X₁) a_e(X₂) + b_o(X₁) a_o(X₂).

So two synchronised banks of lengths `L₁, L₂` compose to the value set `{h₂ i + h₁ j}`: steps are multiplied by the
other block's count (the reason given in `BLOCKS.md` why a common block destroys the unit step), and
`h₁` steps of size `h₂` equal `h₂` steps of size `h₁`: an exact two-digit carry in mixed radix. For coprime
`h₁, h₂` the composite contains a run of about `h₂L₁ + h₁L₂ - 2h₁h₂` consecutive values, so banks with
`L ≥ 2h` would reproduce themselves at every scale. The longest known have `L/h = 4/9` and `6/45`. Finding a bank
longer than its block count is the sharp form of the bottleneck.

## 6. Open

1. A synchronised bank with `L ≥ h` (HYP-9832 (a)); equivalently long runs among seeds with a common small `h`.
2. A selector for the two-tape bank: for which `N'` is `K(A, B) = N'` solvable (HYP-9832 (b))?
3. Three-contact banks on bases with several separated blocks.
4. Chan–Pak.

## Files

`scripts/`: `master_formula_check.py`, `cube_bank_tapes.py`, `sync_bank.c`, `cross_bank.c` (inputs: lines of
out-masks per order, from `gentourng`), `two_tape_bank_check.py`, `bank_coverage.py`, `jacobian_levels.py`.
`out/`: the outputs quoted above.
