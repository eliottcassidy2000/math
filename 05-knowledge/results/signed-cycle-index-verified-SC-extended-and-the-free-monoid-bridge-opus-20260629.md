# klein's signed cycle index verified + SC/NS extended to n=9, and the bridge to the free-monoid palindromes

**Historical synthesis, corrected 2026-10-29 by propagation of the existing
2026-10-24 canon repair.** The numerical table and reversal-word fixed-class
interpretation survive. SC is a trace and a fixed-basis-object count, not
the full even eigenspace dimension. The old Euler/Betti/LRC identifications
had no specified complex or predicate-preserving map and are withdrawn.
Current source: [THM-587, per-level signed cycle index](../../01-canon/theorems/THM-587-per-level-signed-cycle-index-antipodal-euler-number.md).

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

## The surviving bridge: `SC = P_n(−1)` is the free-monoid fixed-word count

The repaired statement is the graded antipodal trace `trace(R)`, with
`R` the converse involution on the class basis. Condensation gives a second
computation of the same fixed-class count:
> tournament iso-classes = the **free monoid on strongly-connected primes**, with `R` = **reverse the
> word + complement each letter**; **self-converse (SC) classes = the `R`-palindromic words.** So
> **`SC(n) = #{R-palindromes of weight n}`** in the free monoid — the "square-root"/half object: a
> palindrome is fixed by its first half + a `R`-fixed center (the half-principle).

Writing `A=A000568(n)`, the correct linear-space dimensions are
`dim V_+=(A+SC)/2` and `dim V_-=(A-SC)/2`. Each swapped class pair adds
one vector to each sector; each fixed class adds one even vector. At order6,
`A=56`, `SC=12`, `dim V_+=34`, and `dim V_-=22`. Thus the fixed-word count
is not the even-sector dimension.

## Withdrawn transfers and retained provenance

The original June synthesis identified this decomposition with an LRC cap
obstruction, a Borsuk--Ulam certificate, and unspecified Betti numbers.
Those identifications are not established by the counts above. In
particular, freeness on labelled tournaments and fixed isomorphism classes
refer to different actions. A signed trace alone supplies neither a
differential nor a topological fixed-point theorem. Historical HYP-3538
and HYP-3544 were inspiration, not proved dependencies. The old maximum-H
table and its cross-domain equivalence are not revalidated here.

## Status
- **Verified:** klein's `P_n(±1)=A000568/SC`; **extended:** `SC(9)=2752`, `NS(8,9)=3352, 94392`.
- **Surviving bridge:** `SC = P_n(−1) =` the free-monoid `R`-fixed-word count.
- **Not supplied:** an Euler/Betti realization or an LRC witness-preserving map.

The repair and the distinct roles of12 and22 are traced in
[the October29 recovery atlas](twelve_twenty_two_trace_20261029/ARCHAEOLOGY.md).

Related: klein THM-587 (signed cycle index), mac-mini HYP-3544 (Ky-Fan on `Q_d`), HYP-3538 (LRC R-odd),
opus free-monoid + reversal/half-principle reflections, CLAUDE.md NS-merged table, OPEN-Q-108.
