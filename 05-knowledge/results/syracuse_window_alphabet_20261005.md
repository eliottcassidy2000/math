# The Syracuse 4-window alphabet: the four 4-tournaments are the truth table of two overlapping 2-step descents, the vortices are the mixed windows, and the densities are 11 : 5 : 5 : 11

2026-10-05, session opus-2026-10-05-S12 (worktree `codex/session-trees-tournaments`), the first "next target" of [the trees-versus-4-tournaments note](trees5_tournaments4_converse_20261005.md) (S11). Owner's directive: "do the Syracuse odd-to-odd 4-window alphabet next".

**PROVED (elementary, all non-degenerate windows):** Theorem 1 (the B2 law under the descent gauge), Theorem 2 (the ascent gauge sees only the two self-converse classes), Proposition 3 (time reversal is the converse of the opposite-gauge reading), Proposition 4 (the four classes have exact natural densities 11/32, 5/32, 5/32, 11/32 among odd integers), Theorem 5 (the 5- and 6-window alphabets are exact for all odd x: 8 of 12 and 25 of 56 classes). **FINITE-EXACT:** the census of all odd x below 2^22 for m = 4, 5, 6. **CITED (repository):** the order-laws note's short-lag genericity and the Redei four-vertex note's Theorem C and involution dictionary, of which the Lemma and Proposition 3 are special cases; the credit note's receipt rule. **NUMEROLOGY (so labelled):** 6 merged 5-window classes = 6 trees on 6 vertices. **OPEN:** nothing universal; Collatz untouched. Audit status: three adversarial subagent auditors (computation with fresh code, proof/typing, repository duplication) attacked eight claims; all eight CONFIRMED, one DOWNGRADED in wording (the credit-episode reading, corrected in section 5), two upgrades adopted (exact densities; exactness of the longer alphabets), one stale docstring fixed (MISTAKE-562). An independent-session audit is OWED.

Script: [syracuse_window_alphabet_20261005.py](../../04-computation/experiments/syracuse_window_alphabet_20261005.py) (imports the class builder of `trees5_tournaments4_converse_20261005.py`); saved output [syracuse_window_alphabet_20261005.out](syracuse_window_alphabet_20261005.out). All assertions pass.

## 1. Inheritance and portfolio

- **Closest proved mechanisms.** The [Redei four-vertex note](procgen_tourn_20260924_four_vertex_sheets_redei.md): a tournament on a Collatz window is read with the map arcs plus numerical order on the remaining pairs ("smaller -> larger"); its Theorem C proves that reversing the word gives the converse tournament, with the remark that the proof works for Syracuse windows, and its involution dictionary gives converse = negation o time reversal. The [order-laws note](collatz_procgen_20260922_order_laws.md), section 3.3: short-lag genericity, every comparison at lag < 17 on the plus sheet equals its word prediction, so every window of at most 12 odd iterates has the generic pattern of its valuation word (PROVED there, computer-assisted). The S11 note, section 9: the compressed-map 4-window alphabet is parity-driven (ooo -> TT4, eoo -> C3 over sink, ooe -> source over C3, the other five words -> STRONG). THM-4512 / the information-dimension bridge: a k-step word descends iff 2^(A_k) > 3^k (coefficient descent).
- **Hostile.** The fixed point 1: windows that reach 1 repeat a value and have no tournament; they are counted (52 below 2^22 at m = 4) and excluded. Windows that merely end at 1 (103 below 2^22, e.g. (17, 13, 5, 1)) are covered by the law.
- **Corrected near miss.** The repository's established gauge (smaller -> larger) is degenerate on odd-to-odd steps (Theorem 2): the drift of the Syracuse map makes "forward chord = descent" the informative convention, and the vortex labels of this note cannot be compared label by label with the compressed-map table.
- **Least-used sidecar.** The valuation word (a_1, a_2, a_3) as the coordinate instead of the parity word.

Portfolio. **Anchor:** the 4-window theorem. **Niche:** what the 5- and 6-windows reach. **Wildcard:** the credit-grammar words of the owner's pasted notes read through the alphabet.

Board: **descent bit / overlapping windows / gauge / time reversal versus converse / exact density / Sturmian thresholds**.

## 2. Definitions and the comparison lemma

For odd x let S(x) = (3x+1)/2^a with a = v_2(3x+1), and write a_1, a_2, a_3 for the valuations of the first three steps, A_k = a_1 + ... + a_k. The 4-window is W(x) = (x, Sx, S^2x, S^3x). Its tournament T(x) has the three time arcs x -> Sx -> S^2x -> S^3x and, on the three non-consecutive pairs {x, S^2x}, {x, S^3x}, {Sx, S^3x}, an order arc:

    gauge DESC: from the larger value to the smaller (a forward chord is a descent);
    gauge ASC : from the smaller value to the larger (the Redei note's convention).

W(x) is degenerate if it repeats a value, i.e. if the orbit reaches 1 within two steps (x = 1, 3, 5, 13, 21, 53, 85, 113, ...). Exact formulas:

    Sx    = (3x + 1) / 2^(a_1),
    S^2x  = (9x + 3 + 2^(a_1)) / 2^(A_2),
    S^3x  = (27x + 9 + 3*2^(a_1) + 2^(A_2)) / 2^(A_3).

**Lemma (the three comparisons).** For every non-degenerate window:

    x  < S^2x  iff  A_2 <= 3;      Sx < S^3x  iff  a_2 + a_3 <= 3;      x < S^3x  iff  A_3 <= 4.

*Proof.* x < S^2x iff (2^(A_2) - 9) x < 3 + 2^(a_1). For A_2 <= 3 the left side is negative. For A_2 = 4, a_1 <= 3, the right side is at most 11 and the left side 7x >= 21 (x = 1 is degenerate). For A_2 >= 5, 2^(A_2) - 9 - 3 - 2^(A_2 - 1) = 2^(A_2 - 1) - 12 > 0. The second comparison is the first one applied to Sx (Sx = 1 is degenerate). x < S^3x iff (2^(A_3) - 27) x < R := 9 + 3*2^(a_1) + 2^(A_2). For A_3 <= 4 the left side is negative. For A_3 = 5 the six words (3,1,1), (2,2,1), (1,3,1), (2,1,2), (1,2,2), (1,1,3) have R = 49, 37, 31, 29, 23, 19, so an exception needs 5x < R, i.e. odd x <= 9, 7, 5, 5, 3, 3; but a_1 fixes x modulo 2^(a_1+1), and the only odd x within these bounds in the right residue class are 5 (word (3,1,1)), 1 (words with a_1 = 2) and 3 (words with a_1 = 1), all degenerate; the actual words of 7 and 9 are (1,1,2) and (2,1,1) with A_3 = 4. For A_3 >= 6 and x >= 3, (2^(A_3) - 27) x >= 3*2^(A_3) - 81 while R <= 9 + (5/4) 2^(A_3), and the difference (7/4) 2^(A_3) - 90 is positive (x = 1 is the equality case 37 = 37 at A_3 = 6). QED. (This is the lag-2 and lag-3 case of the order-laws note's short-lag genericity; the proof is repeated here because the two constants are all the theorem needs.)

The chord pattern determines the class: all three chords forward gives TT4; only {x, S^2x} backward gives the 3-cycle x -> Sx -> S^2x -> x over the sink S^3x; only {Sx, S^3x} backward gives the source x over the 3-cycle Sx -> S^2x -> S^3x -> Sx; every other pattern (the long chord backward, or both short chords backward) is STRONG (a direct check of the eight patterns; the S11 note's created-cycle law gives the same, since a backward long chord closes the Hamiltonian path).

## 3. The theorem: B2 of the two overlapping 2-step descents

Put a = [x > S^2x] = [a_1 + a_2 >= 4] and b = [Sx > S^3x] = [a_2 + a_3 >= 4] (the Lemma).

**Theorem 1 (gauge DESC).** For every non-degenerate window, T(x) is the following function of (a, b):

| (a, b) | class of T(x) | past x | present S^3x | 3-cycles | density |
|---|---|---|---|---|---|
| (1, 1): both 2-step windows descend | TT4 | source (score 3) | sink (score 0) | 0 | 11/32 |
| (0, 1): only the last two steps descend | C3 over a sink (scores 0222) | on the 3-cycle (score 2) | the sink | 1 | 5/32 |
| (1, 0): only the first two steps descend | source over C3 (scores 1113) | the source | on the 3-cycle (score 1) | 1 | 5/32 |
| (0, 0): neither | STRONG (1122) | score 1 if A_3 <= 4, else 2 | score 2 if A_3 <= 4, else 1 | 2 | 11/32 |

*Proof.* In the three cells with a = 1 or b = 1 we have A_3 >= 5, so by the Lemma the long chord {x, S^3x} is forward; the pattern is then all-forward, only-{x,S^2x}-backward, or only-{Sx,S^3x}-backward respectively, and the classification of patterns gives the three classes with the stated marked vertices. In the cell (0,0) both short chords are backward and the pattern is STRONG whichever way the long chord points; the long chord is backward iff A_3 <= 4 (words (1,1,1), (1,1,2), (1,2,1), (2,1,1)) and forward iff the word is (2,1,2), the only word with a_1 + a_2 <= 3, a_2 + a_3 <= 3 and A_3 >= 5. QED. **PROVED**; the census below 2^22 finds no exception among the 2,097,100 non-degenerate windows. The number of 3-cycles of T(x) is the number of failed 2-step descents, 2 - a - b.

So the four 4-tournaments of the S11 chain are the B2 = {{}, {a}, {b}, {a,b}} of the [free-gas reflection](../../07-reflections/the-free-gas-and-the-B2-atom-triangular-numbers-and-the-four-4-tournaments.md), realised with a = "the first two steps descend" and b = "the last two steps descend" instead of a = +1, b = /2. The two bits are the k = 2 instances of the coefficient-descent criterion 2^(A_k) > 3^k on the two overlapping 2-step sub-windows. The converse pair is the pair of mixed windows and the self-converse classes are the symmetric ones: TT4 = both descend, STRONG = neither. In the S11 grading by cycle content, L(T(x)) = 1, 3, 4 according as the window has 2, 1, 0 two-step descents: the cycle content of the window tournament is its non-descent. In the repository's ASC convention the vortices never occur on Syracuse windows (Theorem 2), so the vortex labels here are in the opposite gauge to the compressed-map table of the Redei note (section 4.5) and of S11 section 9, and the two tables must not be compared label by label.

**Theorem 2 (gauge ASC).** Under the Redei note's convention only the two self-converse classes occur: T(x) is TT4 iff A_3 <= 4 (the window ascends over three steps, which forces both 2-step ascents) and STRONG otherwise; the vortices never occur. *Proof.* All chords forward means three ascents; x < S^3x gives A_3 <= 4, hence A_2 <= 3 and a_2 + a_3 <= 3. If the long chord is backward the pattern is STRONG. A vortex would need x < S^3x together with one 2-step descent, impossible since that forces A_3 >= 5. QED. The mechanism: the ascent comparisons are nested (a 3-step ascent forces both 2-step ascents) while the descent comparisons are not (a 3-step descent leaves the two 2-step bits free), so the gauge aligned with the drift of the map is the one that reads the two bits.

**Proposition 3 (time reversal).** Reversing the time arcs of W(x) while keeping the value-determined chords is the same as reversing every arc of the opposite-gauge tournament (a labelled-graph identity, asserted on every window of the census); hence the time-reversed window read with gauge DESC is the converse of T_ASC(x), i.e. TT4 iff A_3 <= 4 and STRONG otherwise. Pure time reversal therefore sends TT4 and both vortices to STRONG and is not a map of the DESC alphabet to itself. The converse of T_DESC(x), which the Redei note identifies as negation o time reversal and, at the word level, as word reversal (its Theorem C), is the swap (a, b) -> (b, a): it exchanges the two vortices and fixes TT4 and STRONG. **PROVED** (an instance of the Redei note's involution dictionary).

**Proposition 4 (exact densities).** The class of T(x) depends only on the capped word (min(a_1, 4), min(a_2, 4), min(a_3, 4)), and the odd x with a given capped word form a union of residue classes modulo 2^(A+1) with A <= 12. Hence the four classes have exact natural densities among the odd integers, equal to the Haar masses of the capped word (P(a = k) = 2^(-k)):

    d(TT4) = P(a = b = 1) = sum_(a_2) 2^(-a_2) P(a_1 >= 4 - a_2)^2 = 1/32 + 1/16 + 1/8 + 1/8 = 11/32,
    d(each vortex) = 1/2 - 11/32 = 5/32,
    d(STRONG) = P(A_3 <= 4) + P(word (2,1,2)) = 10/32 + 1/32 = 11/32,

and under gauge ASC d(TT4) = 5/16, d(STRONG) = 11/16. The two bits each have density 1/2 and are positively correlated through the shared a_2. **PROVED**; the exact rational enumeration of words in the script reproduces both laws, and the census below 2^22 (2,097,100 windows) gives

| class | count | frequency | density |
|---|---|---|---|
| TT4 | 720845 | 0.343734 | 11/32 = 0.343750 |
| C3 over sink | 327680 | 0.156254 | 5/32 = 0.156250 |
| source over C3 | 327679 | 0.156253 | 5/32 = 0.156250 |
| STRONG | 720896 | 0.343759 | 11/32 = 0.343750 |

(the C3-over-sink count is exactly 5 * 2^16; the deviations, at most 1.6e-5, are the degenerate windows and the residue boundary). The marked scores are exactly as in the table of Theorem 1 (STRONG: (1,2) on the 655,360 windows with A_3 <= 4 and (2,1) on the 65,536 windows with word (2,1,2)).

## 4. Longer windows: what the 5- and 6-window alphabets reach

For the m-window the chord {S^i x, S^j x} is forward (a descent) iff A(i, j) = a_(i+1) + ... + a_j >= t_(j-i), where t_k = floor(k log_2 3) + 1 = 4, 5, 7, 8, 10, 12, 13, ... (k = 2, ..., 8) is the least total valuation that makes a k-step word contract. The differences t_(k+1) - t_k in {1, 2} follow the Sturmian word of log_2 3 (the order-laws note's ladder-epoch lemma), so a descent over k steps forces a descent over k+1 steps exactly when t_(k+1) = t_k + 1 (k = 2, 4, 7, ...) and not when the threshold jumps (k = 3, 5, 6, ...: (3,1,1,1) descends over three steps and ascends over four, since 64 < 81).

**Theorem 5 (the 5- and 6-window alphabets are exact).** For m <= 6 the chord patterns realised by non-degenerate windows are exactly the generic patterns of the valuation word, for every odd x. *Proof.* A k-step chord starting at the odd value y with word w, total s and 2^s > 3^k is a descent iff (2^s - 3^k) y > c_w, where c_w is the constant of S^k y = (3^k y + c_w)/2^s; with 2^s < 3^k it is always an ascent. So a non-generic chord needs y < c_w/(2^s - 3^k), whose maximum over words is 11/7, 49/5, 331/47, 1121/13 (k = 2, ..., 5), below 90; every odd y is itself a window start in the census, which checked all y < 2^22 with zero mismatches. (The order-laws note's short-lag genericity gives the same for every window of at most 12 odd iterates.) QED. The alphabets (gauge DESC; census support = exact support):

| m | chord patterns realised | classes occurring | self-converse | converse pairs | modulo converse | missing classes |
|---|---|---|---|---|---|---|
| 4 | 4 of 8 | 4 of 4 | 2 | 1 | 3 | none |
| 5 | 16 of 64 | 8 of 12 | 4 | 2 | 6 | 4, all strong, among them the regular tournament (2,2,2,2,2) |
| 6 | 49 of 1024 | 25 of 56 | 7 | 9 | 16 | 31, all with largest strong component >= 5 (8 with L = 5, 23 with L = 6) |

m = 5 exact densities: the strong class with scores (1,1,2,3,3) 19/64; TT5 13/64; the L = 3 converse pair (1,1,1,3,4) / (0,1,3,3,3) 1/8 each; the L = 4 pair (1,1,2,2,4) / (0,2,2,3,3) 5/64 each; the self-converse (0,2,2,2,4) 3/64 and the strong (1,2,2,2,3) 3/64 (census frequencies within 3e-4 at 2^22; the densities are exact by the residue-class argument of Proposition 4 with cap 8, valid for m <= 6 because 2^8 > 3^5). By L the merged occurring classes are distributed (1, 2, 1, 2) over L = 1, 3, 4, 5. No Syracuse 5-window is regular and every missing class is strong: the realisable patterns are threshold patterns of a monotone word and cannot produce the balanced arc set of the regular tournament; at m = 6 the classes (2,2,2,2,2,5) and (0,3,3,3,3,3), a vertex over or under a regular 5-tournament, are among the missing for the same reason.

Numerology trap, logged so it is not promoted: the 6 merged 5-window classes equal the 6 trees on 6 vertices, but the L-fibres (1,2,1,2) differ from the diameter fibres (1,2,2,1) of the S11 note, and the 6-window count 16 is not the 11 trees on 7 vertices.

## 5. The credit-grammar episodes through the alphabet

The owner's pasted credit notes use G(x) = (9x+5)/8, the valuation word (1,2), and H(x) = (729x+669)/1024, a common-future dependency with the word v = (1,2,1,1,1,2) (U(F_v(x)) = U(H(x)); H is not itself a forward Syracuse word). No two consecutive letters of G, v, GG, or the letter concatenations v.G and G.v sum to 4 or more, so every 4-window inside these letter words has (a, b) = (0, 0): STRONG. But a letter concatenation is not the orbit word of an episode: the credit note's section 5 receipt rule replaces the suffix (b, tail) by (v, b+2, tail), so an H-then-G episode has the orbit word (1,2,1,1,1,2, b+2, 2, ...) with b + 2 >= 3 at the junction. Witness x = 12443 = 155 mod 2048 (the native guard of H), H(x) = 8859 = 11 mod 16 (G legal after H): the Syracuse word of x is (1,2,1,1,1,2,3,2,1) and its 4-windows read STRONG, STRONG, STRONG, STRONG, then C3 over sink (1,2,3), TT4 (2,3,2), source over C3 (3,2,1). The descent the episode pays for shows up exactly at the junction valuation, not inside the letter words. **PROVED for the witness, FINITE-EXACT as a reading**; it is a reading, not a payment statement: the alphabet sees cycle content (non-descent) and nothing about credit. Cheapest decisive test for the credit thread: tabulate, for the bank of checked ROOT words, the B2 profile of its episode orbit words at the junctions (prediction: the junction letters carry the transitive and vortex windows, the letter words only STRONG). BLOCKED until the bank file is in the repository.

## 6. Verdict, traps, next

- **Settled.** The Syracuse 4-window alphabet is exactly the B2 truth table of the two overlapping 2-step coefficient descents (Theorem 1), with the converse pair as the mixed windows and exact densities 11 : 5 : 5 : 11 (Proposition 4). The Redei gauge is the wrong gauge for odd-to-odd steps (Theorem 2). Time reversal is the converse of the opposite-gauge reading, not an involution (Proposition 3). The 5- and 6-window alphabets are exact (Theorem 5) and never contain the regular tournament.
- **What is old and what is new.** Old: the window class is a function of the valuation word, with word reversal = converse (Redei note, Theorem C); short-lag genericity (order-laws note); the compressed-map table. New: the odd-to-odd m = 4 table and its B2 reading, the DESC gauge and the vortex-free ASC statement, the exact densities, the m = 5, 6 supports and densities, the Sturmian threshold structure of the implications, the credit-episode junction reading.
- **Compared with the compressed map.** There the alphabet is parity-driven and STRONG dominates; here the valuation word drives it and the two self-converse classes carry 22/32 of the density with the vortices sharing 10/32. The S11 bridge prediction (a pairwise field, not a degree-one one, forces the global atom) is untouched: the window tournament is a pairwise object by construction, and its class is a two-bit function of the word.
- **Traps.** The 6 = 6 of section 4. Reading "TT4 = fully paid window" as credit: the bits are descents of value, not of the credit account (v is a six-step value descent and yet every one of its 4-windows is STRONG). Treating a letter concatenation as an orbit word (MISTAKE-562).
- **Next.** (1) The m-window alphabet in general: which n-tournament classes are threshold patterns of a monotone word with the Sturmian thresholds t_k; conjecture from m = 5, 6: every missing class is strong with L >= m - 1 and the regular tournament never occurs for m >= 5 (the script's cap must be raised to >= 10 for m = 7). (2) The 3x - 1 sheet: the order-laws note proves short-lag genericity there for lag < 12 with different gate crossings (165, 309, 549 at lag 12), so the m <= 6 alphabets transfer with the same thresholds; record whether the sheet changes any density (prediction: no, the capped-word densities are sheet-independent). (3) The owner's bank, when supplied.

Reproduction:

    python -X utf8 04-computation/experiments/syracuse_window_alphabet_20261005.py --N 22 --save

META-PATTERNS cards used: *Type every analogy and every implication*; *Respect symmetries by searching orbit representatives* (the gauge is a symmetry choice and the two gauges give different alphabets); *Expose the obstruction first, choose the scale second* (the Sturmian threshold jumps are the obstruction to monotone implications between window scales); *Search the statement before the method* (the lemma was already in the order-laws note as short-lag genericity).
