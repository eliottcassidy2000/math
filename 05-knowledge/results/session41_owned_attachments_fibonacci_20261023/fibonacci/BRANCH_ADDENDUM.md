# Both discrepancy branches: zero-sensitive rescue and an all-horizon larger universe

Status: PROVED for explicitly different grammars; FINITE-EXACT branch receipts. This addendum does not change the frozen canonical model in RESULTS.md.

## Incoming-work credit

Before publication, inspected incoming commit `0429c7d8fdf`, `05-knowledge/results/odd_cell_gluing_20261023/ZERO-SLOT.md`. Its §1 already proves the four-sheet retained-zero fibre theorem; §3 proves arbitrary-subset C→B recolorings and the global historical35 repair, explicitly with occurrence-dependent choices and a shared absolute-position subset. Those results were obtained convergently here and are credited to that incoming source, not claimed as a new Session41 discovery. The sharpening supplied below is the COMPLETE κ-language classification (including independent B→A choices and necessity), its exact prefix counts/exponential rates, and the fixed-dictionary one-step rescue versus full-language distinction. The older D-language theorem remains credited to `reset_20260926_colours.md`.

## 1. Which data and which conventions?

The inherited source `reset_20260926_colours.md`, §4, supplies exactly35 colors, not a36th authoritative color. Its36 analysis is a predicted obstruction, not a supplied observation. The same is true of `crossroads_crossing_20260926_colour.md`, §1. Translated toA,B,C, the supplied35th symbol isB, whereas the canonical wordV hasC there. We preserve both.

Three different mathematical universes must not be identified:

1. Fixed canonical block words: each occurrence ofFn uses the same canonical blockW_n from RESULTS.md; only the outer support can vary.
2. A globally chosen dictionary: choose one possibly noncanonical expansion for each Fibonacci value and reuse it in every occurrence. This lies between the fixed dictionary and the next universe; it is not classified here.
3. Occurrence-dependent recursive choices: every occurrence ofFn can independently use any legal proper smaller-slot decomposition, recursively, and choices may differ between rows. Only the resulting colored rows must satisfy prefix coherence. This is the larger universe classified below. No requirement of one shared expansion per Fibonacci value is silently retained.

Use R for the empty-zero-slot grammar andD for the compressed three-unit grammar. W is the ascending word; V is its reversal. Coherence is adding one unit at the left ofW, equivalently extendingV at the right.

## 2. Exact fixed-block branch: rescue at36, obstruction at37

Let j≥6 and writeQ=W_j. Complete fixed-block row languages at F_j+t, for t=1,2,3, are

- t=1, bothR,D: AQ, BQ, CQ.
- t=2, D: ACQ. InR: ACQ orABQ.
- t=3, bothR,D: AACQ orBACQ.

These follow immediately from the complete fibres in RESULTS.md; j≥6 keeps all low summands separate fromF_j and below the next Fibonacci value. They are languages of words, so multiple supports giving the same word are counted once.

Therefore the branchBQ at F_j+1 dies immediately inD. InR it extends uniquely toABQ at F_j+2, then dies at F_j+3 because neitherAACQ norBACQ has suffixABQ. Forj=9, this is the actual35→36→37 branch:

35: B W9;
36: A B W9, possible only inR;
37: no continuation if all inner blocks remain canonical.

This proves an all-size fixed-block boundary obstruction, not an obstruction for all imaginable recursive interpretations. Complete supports and full words for34–40 are in branch-receipt.json.

## 3. The larger recursive languages factor exactly

Let L_n^R be the set of ascending leaf words obtained by recursively expandingFn through any proper R support. All three marked units remain terminal. Choices at distinct occurrences are independent. DefineL_n^D similarly.

For a letterc put
κ_R(A)={A}, κ_R(B)={A,B}, κ_R(C)={B,C};
κ_D(A)={A}, κ_D(B)={A,B}, κ_D(C)={C}.
Extendκ to a word by independent Cartesian choices at its positions.

**Theorem4 (exact block languages).** For every n≥3,
L_n^R=κ_R(W_n), and L_n^D=κ_D(W_n).

The D statement is already proved in `reset_20260926_colours.md`, §4. Here is a simultaneous boundary-sensitive proof.

For n=3, D has onlyAC; R hasAC orAB. For n=4, D hasAAC orBAC; R hasAAC,BAC,AAB,BAB. These are precisely the asserted Cartesian languages.

For odd n, the two proper R supports agree on all higher indices and end in A+C orA+B. The sole D support ends inA+C. For even n, the two proper supports in either grammar differ only in their final low marked unitA orB, with the next blockF3. After the higher recursive choices, these boundary choices are independent of every inner choice. More succinctly, for both grammars and n≥5,
L_n=L_(n−2)L_(n−1),
as a concatenation of languages with independent choices: the low portion of a properFn expansion is exactly a properF(n−2) expansion, while its highest block isF(n−1). Base cases and induction prove Theorem4. This recurrence is for the larger language, not a requirement that repeated values share one expansion.

## 4. All infinite coherent stripes, in both universes

**Theorem5 (complete global classification under occurrence-dependent expansion).** LetV be the canonical infinite word beginningCABCACAB. An infinite wordx has a legal recursively expanded row of every length, with those rows prefix-coherent, if and only if

x_n∈κ_R(V_n) for every n inR,
or x_n∈κ_D(V_n) for every n inD,
respectively.

Necessity: inspect any Fibonacci-length rowFn. An outer singletonFn followed by recursive expansion, or any proper outer support followed by expansion, gives exactly the block language of Theorem4. Taking arbitrarily largeFn imposes the stated coordinatewise restrictions at every fixed position. Standalone low-row boundary liberties cannot escape this larger-row test.

Sufficiency: use ordinary Zeckendorf support for a row of lengthN and concatenate its recursively expanded blocks in descending order. The canonical blocks giveV[1..N]; Theorem4 permits all prescribed coordinatewise changes inside every block of index≥3. If a terminal ordinary unitC is present, D requires it to stayC. InR the desired terminal letter can instead beB. Replace that outerC slot by the markedB slot. This remains legal because the next selected ordinary index, if present, is at least4. The other blocks are unchanged. Thus every permitted prefix has a valid root support and recursive expansion. Choosing these for allN gives the required coherent rows. No consistency of internal decomposition trees between different rows is required beyond equality of their visible prefixes.

**Consequences.** The old35thB is impossible in any infiniteD stripe under this terminal-unit grammar, even allowing all recursive choices. InR it is completely valid: choose x35=B instead of canonicalC and keep every other position canonical. This one-site alteration has legal expansions at every finite horizon. Its36th letter is necessarilyA; its37th may beB orA. Thus the zero-retained universe repairs the old branch for all horizons, not merely one extra row, once occurrence-dependent inner choices are admitted.

This does not determine the intended rule from35 supplied symbols. It proves that both the canonical branch and the old35B branch inhabit a precisely defined larger universe. It also identifies exactly what extra freedom does the work.

## 5. How much larger?

Let a(N),b(N),c(N) count canonical letters inV[1..N]. The number of lengthN prefixes that extend to infinite coherent stripes is exactly

2^b(N) inD;
2^(b(N)+c(N))=2^(N−a(N)) inR.

These are global-extension counts, not counts of arbitrary standalone row representations. Every free coordinate is independently binary. Consequently both spaces have continuum cardinality, and their exponential prefix growth rates, in natural logarithms per visible unit, are

φ^−3 log2 inD, andφ^−1 log2 inR.

The bounded discrepancy already proved implies the exponents differ from their linear growth terms byO(1). This is a precise form of hidden boundary information producing macroscopic choice. It is not a statement that the original three labels form a cyclic symmetry, nor a Collatz or tournament transport theorem.

## 6. Verification and limits

`check_branches.py` independently builds the recursive languages throughF7=13 and compares them exactly with the Cartesian descriptions. It checks fixed-block branching forj=6,...,20, enumerates full fixed supports/words34–40, and constructs actual recursive parse witnesses for every prefix through300 of the canonical stripe with only35 changed toB. It rejects the36-prefix inD. Normal and python−O outputs match.

The all-size conclusions are proved above, not extrapolated from these checks. The globally fixed noncanonical dictionary in universe2 remains a separate problem; the theorem must not be cited as solving that stricter synchronization requirement.
