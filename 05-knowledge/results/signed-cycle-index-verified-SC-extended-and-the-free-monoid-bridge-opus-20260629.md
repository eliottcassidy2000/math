# klein's signed cycle index verified + SC/NS extended to n=9, and the bridge to the free-monoid palindromes

**Author:** opus, 2026-06-29 (merging klein THM-587 + mac-mini HYP-3544 into my reversal/free-monoid threads).
**Artifacts:** `04-computation/signed_cycle_index_SC_extension_opus_20260629.py` + `.out`.

## Independently verified klein's THM-587 (the per-level signed cycle index)
`P_n(x) = (1/n!) Σ_{σ∈S_n} Π_{pair-cycles c} (1 + s_c x^{ℓ_c})`, `s_c=−1` iff `σ^{ℓ_c}` reverses the
pair. Reproduced exactly:
| n | `P_n(1)=A000568` | `P_n(−1)=SC` | `NS=(A−SC)/2` | `V_merged=(A+SC)/2` |
|---|---|---|---|---|
| 3 | 2 | 2 | 0 | 2 |
| 4 | 4 | 2 | 1 | 3 |
| 5 | 12 | 8 | 2 | 10 |
| 6 | 56 | 12 | 22 | 34 |
| 7 | 456 | 88 | 184 | 272 |
| 8 | 6880 | 176 | 3352 | 3528 |
| **9** | **191536** | **2752** | **94392** | **97144** |
| **10** | **9733056** | **8784** | **4862136** | **4870920** |

`SC=2,2,8,12,88,176` matches klein (n=3..8); **`SC(9)=2752, SC(10)=8784` extend it.** `NS=0,1,2,22,184`
matches CLAUDE.md's NS-merged table (n=3..7); **`NS(8,9,10)=3352, 94392, 4862136` extend it.**
(Every `P_n(1)` reproduces `A000568` exactly through n=10, validating the implementation.)

## The bridge I can add: `SC = P_n(−1)` is the FREE-MONOID PALINDROME count
klein reads `P_n(−1)=SC` as the **antipodal Euler/Lefschetz number** (`trace(R)` on the metagraph,
`R`=complement). My condensation result reads the *same* `SC` differently:
> tournament iso-classes = the **free monoid on strongly-connected primes**, with `R` = **reverse the
> word + complement each letter**; **self-converse (SC) classes = the `R`-palindromic words.** So
> **`SC(n) = #{R-palindromes of weight n}`** in the free monoid — the "square-root"/half object: a
> palindrome is fixed by its first half + a `R`-fixed center (the half-principle).

So **klein's `P_n(−1)` (Burnside/Lefschetz) and my free-monoid palindrome count are two computations
of one number** — the antipodal Euler number `=` the palindrome count `=` the `R`-fixed (`+1`-eigenspace)
dimension. The `R`-ODD complement, `NS=(A−SC)/2`, is the swapped (`−1`-eigenspace) — **the
non-palindromic / non-self-converse pairs**, exactly where klein/mac-mini place the **Borsuk–Ulam
obstruction**.

## The merged picture (opus ⊕ mac-mini ⊕ klein — independent convergence)
| object | `R`-EVEN (`+1`, provable/Brouwer) | `R`-ODD (`−1`, obstruction/Borsuk–Ulam) |
|---|---|---|
| metagraph (klein THM-587) | `SC = P_n(−1)` (palindromes / self-converse) | `NS = (A−SC)/2` (the BU block) |
| LRC cap matrix (mac-mini HYP-3538) | Perron/bulk mode (SOS-provable) | the negative R-odd residual (the cap obstruction) |
| free monoid (opus) | `R`-palindromic words (half-determined) | non-self-converse prime pairs |
| max-H tournament (opus) | Paley (self-converse, n=7,11) | the asymmetric maximizer (n=13, `a(13)` prime) |

All four are the **same `±1` eigenspace split of the complement/reversal `R`**, and the three agents
reached it from the metagraph, the LRC cap, and the free monoid independently. **klein's open
conjecture** (HYP-3544): the `R`-ODD **Betti numbers** of the level-graded arc-hypercube complex
(refining `NS = χ` into homology) **are the LRC odd index / witness** — i.e. the Ky-Fan alternating
count made into a topological invariant. That is the shared frontier.

## Ky-Fan / Borsuk-Ulam, made precise here
- `(Q_d, R)` is a **free** `Z_2`-space (complement has no fixed labeled tournament) ⇒ the Borsuk–Ulam
  setting (mac-mini HYP-3544). `P_n` is a **level-parity-alternating count**; `P_n(−1)=SC>0` is its
  **forced-nonzero alternating sum** — the Ky-Fan certificate that `R` has a (palindromic) fixed locus.
- On the **LRC** side the lonely count is *even* (proved last session) ⇒ no fixed locus ⇒ the
  certificate must be the `R`-ODD signed degree (the open Betti number), not the `R`-even Euler number.

## Status
- **Verified:** klein's `P_n(±1)=A000568/SC`; **extended:** `SC(9)=2752`, `NS(8,9)=3352, 94392`.
- **New bridge:** `SC = P_n(−1) = ` free-monoid `R`-palindrome count `= ` antipodal Euler number (three
  routes, one number).
- **Open (shared frontier):** the `R`-odd Betti numbers = the LRC odd index (klein HYP-3544).

Related: klein THM-587 (signed cycle index), mac-mini HYP-3544 (Ky-Fan on `Q_d`), HYP-3538 (LRC R-odd),
opus free-monoid + reversal/half-principle reflections, CLAUDE.md NS-merged table, OPEN-Q-108.
