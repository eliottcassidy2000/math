# Whip turnbacks, a fixed-count obstruction, and a wider fold

The physical inspiration is Goriely and McMillen's [traveling-loop model of a cracking whip](https://doi.org/10.1103/PhysRevLett.88.244301). They distinguish motion of the loop from motion of material at the loop and study tension, taper and boundary effects. In a tournament the analogous distinction is between changing the visible tail word and sending an actual edge back into the earlier body. The equations below are tournament identities, not a model of the sound barrier.

## Close a seam rather than repaint a still-ordinal tail

- [Backward seam](BACK_CROSS.md): one reversed ordinal edge raises the ordinary spine count by the exact rank-two product S_A S_B−d_A d_B. A genuine single crossing cannot retain the old h. Two physical arms can pay the same price: at h=41, opposite first tail marks and different back-edges give 891 and 893 with the **same** (h,S,d)=(41,1262,0). The +2 survives every common transitive continuation.
- [Retrofit two actual fold arms](FOLD_FIRST_SEAM.md): allow the first tail vertex every contact subset of the four inherited prepared spines at heights 3,4,5. No first-vertex-only retrofit gives a common eight-arm state. A compensated return involving two of the original arms does work for every j≥3: 173+181·2^(j−1) and 175+181·2^(j−1), at shared h=77 and S=181·2^j−154, d=0. The pair has changed its absolute position, not the original eight-target window. [Its exact boundary](H77_BOUNDARY.md) rules out adding a third neutral arm at h=77, λ=181, even with arbitrary spine order, and bounds that congruence-complete pair bank by (2/647)X+O(log X) distinct counts through X.
- [Exact return spectrum](TURNBACK_FOCUS.md) and [turnback lock](TURNBACK_LOCK.md): read a literal returning tip with all ordered slots. Matching ordinary h can leave an exponentially growing two-slot difference. In another physical pair both growing channels cancel while a constant 48 survives at every depth. An actual cyclic reader converts that into exactly 48 Hamilton paths. A two-slot lock is not full continuation equivalence.

## Why the original eight cannot close at their old h

[Whole fixed-state inversion](INVERSE_STATE_FIBRE.md) works backward through the strong-component bank, recovering the old start by exact division at each step. It allows arbitrary final spines of the prescribed ordinary count, not just edits at the old seam. The original eight odds 150401–150415 cannot share a continuation at **h=15, any order and any start**, even allowing independent count-neutral transitive preparation. The complete four-target subintersection is already empty; the code checks five unit targets. Exact stopping bounds make this stronger than a bounded-edge search. This is not an obstruction at a different h.

## A wider physical window elsewhere

[Six and twelve](SIX_TWELVE.md) gives six strong matched seeds at h=15, starting with 20227–20237 by twos at preparation height 4. Every common transitive tail preserves all six. Neighboring fold addresses then give **twelve consecutive odd counts** at every odd preparation height ≥5, starting with 789799435–789799457 by twos. The two folded terminal states still differ. This is increased absolute coverage, not twelve matched arms.

The useful next step is to choose a new ordinary spine count and *both* seam responses from an absolute target, so that many arms pay the same h cost and their remaining coarser/finer channels match simultaneously. The exact h=41 and h=77 pairs do this for two arms; the old h=15 obstruction makes clear why simply enlarging an edit budget at fixed h cannot close the original window. [Chan–Pak](https://arxiv.org/abs/2308.10214) remains open.

## Random routes and replays

The separate random draws are recorded in the notes: the source-linear/quadratic boundary pair and adjacent-factor precision audit suggested completing a whole fiber before claiming a missing response; the rhombus and adaptive Collatz notes suggested distinguishing a genuine carrier change from a nominal edit and allowing a new partner length; the Farey resonance note suggested matching exact powers and odd parts. These are research prompts, not theorem transfers.

Run the companion Python scripts and compile the C++ programs with C++17. The companion receipts separate literal counts from finite search scope; the all-order statements follow from the stated structural stop rules and transfer laws.
