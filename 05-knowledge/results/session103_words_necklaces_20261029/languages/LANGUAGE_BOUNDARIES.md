> **Historical pre-recovery note.** Any missing-definition or open-criterion status below records an earlier stage only. The exact grammar is recovered and the all-order criterion is now proved; see [current synthesis](../README.md) and [actual grammar refinements](ACTUAL_GRAMMAR_REFINEMENTS.md). The free-language examples remain valid controls, not the actual insertion language.

> **RECOVERY UPDATE:** The actual five-letter insertion grammar is now recovered at commit `d3132eafbfa9c43c6d9f1110566b1e9c3a3286e1`. Read [ACTUAL_GRAMMAR_REFINEMENTS.md](ACTUAL_GRAMMAR_REFINEMENTS.md) for the current definitions and results. The missing-definition statements below describe the earlier pre-recovery stage only. These test-alphabet models remain controls, not the actual grammar.

# Insertion words: exact language toolkit and the missing grammar boundary

**PROVED:** automaton, weighted generating-series, first-hit and cyclic-language formulas below for the explicitly defined formal languages. **FINITE-EXACT:** independent brute-force/transfer checks through word length 9. **OPEN / NOT RECOVERED:** the actual insertion alphabet, seed, admissibility rules, order grading and construction-to-isomorphism quotient behind the reported 182,409,919; equality of forbidden-pattern avoidance with nonsingularity of the tournament observer. No tournament has been constructed or tested here.

## 0. Evidence and inheritance

Source checkout `75bd8b27b6dd` contains HYP-9445, `05-knowledge/hypotheses/HYP-9445-run-observer-torsion-and-hamiltonian-cycles-star-law-fermat-routes-and-gaussian-gf-questions.md`. Its item 2 gives nonsingular counts 1,1,3,7,16,36,81 at tournament orders 3..9; 182 at order 10 is still called a prediction there. The parent has newer transcript evidence of 182,409,919; we retain this without claiming the older file verifies it. The manuscript source does not define b/e insertion operations. That missing definition is a real mathematical dependency, not a clerical detail.

Closest inherited mechanism: the scalar recurrence candidate in HYP-9445. Hostile: a matching characteristic polynomial does not fix boundary vectors. Least-used sidecar: initial/terminal automaton state, letter grading and marked construction endpoint. Concept board: admissible insertion grammar; first forbidden adjacency; observer rank versus count; cyclic boundary; isomorphism quotient. No all-order tournament classification follows here.

## 1. A sharp cheap control against guessing the grammar

Take the **test alphabet**, not a recovered alphabet, Σ={a,b,e}, and let L consist of all words containing neither be nor ee. Let f_m count length-m words including the empty word. The two live states are N (empty or last letter a) and H (last letter b or e). With rows the old state and columns the new state, the transition-count matrix is

    T = [[1,2],[1,1]],   f_m = [1,0] T^m [1,1]^t.

Thus

    f_0=1, f_1=3, f_m=2 f_(m−1)+f_(m−2),
    Σ f_m t^m = (1+t)/(1−2t−t²),
    f_m = 1,3,7,17,41,99,239,577,1393,3363,... .

The target 1,3,7,16,36,81,182,409,919 differs already at m=3. Even **if** the actual alphabet is ternary and m=n−4, bare be/ee avoidance is insufficient: an additional admissibility/seed/quotient constraint is essential.

Another tempting guess, excluding eb as a base-grammar rule, gives adjacency

    [[1,1,1],[1,1,0],[1,0,0]]

with characteristic polynomial λ³−2λ²−λ+1. But its unrestricted word counts are 1,3,6,14,..., not 1,3,7,16,...; its GF is (1+t−t²)/(1−2t−t²+t³). Equal characteristic polynomials do not identify languages or initial data. This explicitly corrects the initial scratch suggestion in the lane messages.

The target sequence itself obeys the candidate recurrence and would predict 2065,4640,10426,23427 next. These are only **conditional predictions**, not computed observer counts.

### 1a. Exact nonnegative boundary realization of the target sequence

**PROVED as a marked-walk identity; CONDITIONAL ONLY as a possible model for the reported counts.** Let

    T = [[1,1,1],[1,1,0],[1,0,0]],
    α = [1,2,0],   β = [1,0,0]^t,   h_m = α T^m β.

Interpret the three states as 0,1,2. There is one marked initial copy of state 0, **two distinguishable marked initial copies** of state 1, and none of state 2; the terminal state must be 0. Then h_m counts length-m walks with these boundary conditions. The mark on either initial copy of state 1 persists as a counting multiplicity even when their underlying state walks coincide. At m=0 only the initial copy of state 0 contributes, so h_0=1. Direct multiplication gives h_1=3 and h_2=7. Since

    det(λI−T)=λ³−2λ²−λ+1,

Cayley–Hamilton proves h_m=2h_(m−1)+h_(m−2)−h_(m−3) for m≥3, and therefore

    Σ_(m≥0) h_m t^m = (1+t)/(1−2t−t²+t³),
    h_m = 1,3,7,16,36,81,182,409,919,2065,4640,... .

This is an exact nonnegative weighted-walk realization of the candidate sequence, **not a recovered insertion grammar**. A walk of m edges has m+1 states, so it also cannot be silently equated with an unrestricted length-m word. The same bulk matrix T, with unrestricted word boundaries, instead gives the 1,3,6,14,... sequence above: a transfer polynomial is only bulk data; initial multiplicities, terminal selection and grading are indispensable sidecars.

No map from these marked walks to unique-Hamiltonian-cycle tournaments is provided, and no equality with observer nonsingularity is asserted. A separately verified chord-row encoding or other tournament constructor must be compared through an explicit map that preserves counts, fibres and the observer predicate; matching scalar series is insufficient.

## 2. Weighted free-language formula

Give letters a,b,e weights x,y,z, respectively. The two-state live matrix is

    T(x,y,z) = [[x,y+z],[x,y]].

The sum of all safe-word weights, with t marking length, is

    F(t;x,y,z) = (1+t z)/(1−t(x+y)−t² xz).

**Proof.** The first state is N; a returns to N from either state, b enters H from either state, e enters H only from N. Summing the geometric matrix series gives `[1,0](I−tT)^−1[1,1]^t`, the displayed rational function. Equivalently multiplying by its denominator gives a polynomial recurrence refined by the counts of each letter. All empty/start/end cases are included. ∎

This is a letter-content refinement, not a claim that letters encode tournament run sizes.

## 3. Product construction for the actual grammar, once supplied

Suppose the actual admissible construction words form a deterministic finite automaton with state set S, initial row α, terminal column β, and transition matrices A_a,A_b,A_e. Unambiguous finite automata work equally well; ambiguous automata instead count accepting paths with multiplicity. Use the safe-memory matrices

    D_a=[[1,0],[1,0]], D_b=[[0,1],[0,1]], D_e=[[0,1],[0,0]],
    E_e=[[0,0],[0,1]], E_a=E_b=0.

On S×{N,H}, put

    U = Σ_σ A_σ ⊗ D_σ,    V = A_e ⊗ E_e,
    α' = α⊗[1,0],          β' = β⊗[1,1]^t.

Then the number of admissible safe length-m words is `α' U^m β'`; the total admissible number is `α'(U+V)^mβ'`. More strongly,

    Z(t,u)=α'[I−t(U+uV)]^−1β'

weights each occurrence of be or ee by u, **including overlapping occurrences**. The coefficient of u^0 is the safe-word series. An initial e is not erroneously counted as a failure, because the initial memory is N.

For a total word length m, the number with first forbidden pair ending at letter k, 2≤k≤m, is

    α' U^(k−1) V (U+V)^(m−k) β'.

Consequently total minus safe is the sum of these first-hit counts. Keeping separate terminal columns refines this distribution by the final grammar state or marked endpoint. Separate V matrices can distinguish be from ee; letter weights refine it by letter content. If insertion letters have different order costs d_σ, replace their transitions by t^(d_σ) times those transitions rather than silently equating tournament order and word length.

**Proof.** The product state records exactly the admissibility state and whether the previous letter is b/e. Every forbidden pair occurs precisely on a V-transition. Expand the matrix products as paths. In the first-hit expression the first k−1 transitions are safe, the kth is bad, and the remainder arbitrary. These cases are disjoint and exhaustive. The same path expansion proves the multivariate series. ∎

**Boundary.** This proves a language-count theorem only. To obtain observer counts one must additionally prove (i) the map from words to tournaments covers the intended class, (ii) its fibres or canonical representatives are correctly treated, and (iii) observer nonsingularity iff avoidance. A characteristic-polynomial fit proves none of these.

## 4. Cyclic boundary and legitimate necklace counts

For the explicit free test alphabet, define the cyclic-safe language by requiring every adjacent pair, **including last-to-first**, to avoid be and ee. For m≥1 its rooted count is

    c_m = tr(Q^m),   Q=[[1,1,1],[1,1,0],[1,1,0]].

Q has eigenvalues 0,1+√2,1−√2, so

    c_m=(1+√2)^m+(1−√2)^m,
    c_1=2, c_2=6, c_m=2c_(m−1)+c_(m−2) for m≥3.

Length one is important: a and b are accepted, but e is not, because its cyclic adjacent pair is ee. The **formal-language** rotation-necklace count is

    N_m = (1/m) Σ_(d|m) φ(d)c_(m/d).

**Proof.** Each rooted cyclic word is a unique closed walk of length m in its letter adjacency graph. A rotation by r fixes precisely repetitions of a cyclic-safe block of length gcd(m,r). Group the rotations by this gcd and apply Burnside. ∎

No tournament necklace interpretation is asserted. In a general grammar trace of a transfer matrix counts state-marked closed walks; it counts words only if the required cyclic lift is unique, or after correcting its multiplicities.

## 5. Rotation and reversal are not free quotient operations

Linear safe-word membership is **not** invariant under rotation: eab is safe but its rotation abe is not. It is **not** invariant under reversal: eb is safe but be is not. A cyclic-safe word is rotation-invariant by definition, but reversal changes the forbidden pair be to eb; hence the cyclic-safe language is not generally reflection-invariant either (for example eba is cyclic-safe while abe is not).

Thus ordinary Burnside bracelets cannot be applied to this directed language, and necklaces cannot be applied to the linear safe language. If both eab and abe were legal actual construction words, and forbidden-pair iff singularity were proved, their constructed tournaments could not be isomorphic: observer nonsingularity is isomorphism-invariant. This is a **conditional incompatibility**, not an actual tournament counterexample, since the constructor remains missing.

The correct sidecar to retain until a theorem says otherwise is the seed, linear starting cut, endpoint/port markings and order grading. A claimed action on unmarked tournament isomorphism classes must be verified on those objects, not inferred from rearranging symbols.

## Reproduction

    python3 session103_words_necklaces/languages/check.py
    python3 -O session103_words_necklaces/languages/check.py

`check.py` uses separate brute-force words and transfer powers through length 9; checks first-hit partition identities, cyclic closure, and direct rotation orbit sets against Burnside. `receipt.json` records the consequence counts and hostile words. Checks do not evaluate M_2 and do not supply missing insertion definitions. The optimized replay is only a smoke run because Python assertions are disabled; the normal replay is the load-bearing verification.
