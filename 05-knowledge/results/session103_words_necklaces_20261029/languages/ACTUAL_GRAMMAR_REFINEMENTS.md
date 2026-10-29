# Actual insertion grammar: occurrence laws, first-hit factorization and quotient hostiles

**Source recovered:** `d3132eafbfa9c43c6d9f1110566b1e9c3a3286e1`, `05-knowledge/results/operation_closed_synthesis_20261029/{THEORY.md,tournament_probe.py,tournament_results.json}`. Definitions are no longer missing. The earlier `LANGUAGE_BOUNDARIES.md` is a pre-recovery control note; its free-ternary models are not the insertion grammar.

**Inherited, not newly claimed:** the source proves the three-state safe-word transfer, its cubic generating function and minimal recurrence/state counts; cites Shelor's classical bijective construction; and verifies observer nonsingularity iff avoidance through tournament order 12. **New elementary deductions here:** occurrence-count generating functions, first-hit factorization, weighted letter refinement, integral equivalence with the earlier boundary model, and a carefully scoped cyclic grammar count. **FINITE-EXACT:** independent isomorphism hostiles to rotating/reversing insertion words. **PROVED, independently audited:** the all-order equivalence between safe words and native ordered-two-path observer nonsingularity is now established by [avoiding-word triangularity](../triangularity/AVOIDING_WORD_TRIANGULARITY.md) and the [bridge–sink star lemma with persistence](../bridge/BRIDGE_SINK_STAR_LEMMA.md); see the [independent adversarial audit](../audit2/REPORT.md). The language calculations here do not supply those matrix arguments; they now combine with them.

## 1. Exact inherited object

Start at the unique strong four-tournament, type B. Type A has one maximum-outdegree vertex a; type B has two a0→a1. A new vertex v beats every old vertex except its parent p, and is inserted after p in the unique directed cycle. The letters are

|Source|Letter|Parent|Target|
|--|--|--|--|
|A|a|predecessor of a|A|
|A|b|a|B|
|B|c|predecessor of a0|A|
|B|d|a0|B|
|B|e|a1|B|

Every letter adds exactly one vertex. Word length m means order n=m+4; no variable-cost macros occur here. The inherited classical bijection identifies legal words, not their rotation or reversal orbits, with ordinary unique-Hamiltonian-cycle tournament isotypes. Actual cycle rotations are a different operation.

For avoidance of be and ee use states (A,B0,B1), with B1 meaning the previous letter was b or e, and B0 the other type-B state. Initial row α=(0,1,0); terminal column 1=(1,1,1)^t. The safe matrix is

    C = [[1,0,1],[1,1,1],[1,1,0]].

The only forbidden transition is e from B1 to B1. Thus the full legal-word matrix, retaining last-letter memory, is K=C+E33. Indices 33 mean the third state B1. Write D(t)=1−2t−t²+t³. The inherited safe series is (1+t)/D(t); all legal words have series 1/(1−3t+t²).

**All-order observer-count consequence.** By the audited equivalence, the native ordered-two-path observer is nonsingular exactly on these safe words, and then is unimodular. Thus, writing u_n for the nonsingular (equivalently unimodular) count at tournament order n, the classical grammar bijection gives

    Σ_(n≥4) u_n t^(n−4) = (1+t)/(1−2t−t²+t³),
    u_4=1, u_5=3, u_6=7,
    u_n=2u_(n−1)+u_(n−2)−u_(n−3)  for n≥7.

In particular **u_13=2065 is now an exact all-order theorem consequence, not a prediction** (and u_14=4640). Coverage of all unique-Hamiltonian-cycle isotypes uses the separately cited classical completeness theorem; the matrix equivalence itself applies to every finite legal word of the specified seed grammar. The observer convention is ordered nonempty spanning paths, with singleton paths allowed; this does not assert a corresponding theorem for arbitrary tournaments or other observers.

## 2. Exact number of forbidden occurrences, including overlaps

Let z_(m,r) count legal words of length m having exactly r occurrences of be or ee. Occurrences may overlap: eee has two. Then

    Z(t,u) = Σ_(m,r≥0) z_(m,r)t^m u^r
           = α[I−t(C+uE33)]^−1 1
           = (1+t−tu)/(D(t)−u t(1−t)²).                  (1)

Therefore

    Σ_m z_(m,0)t^m = (1+t)/D(t),
    Σ_m z_(m,r)t^m = t² [t(1−t)²]^(r−1)/D(t)^(r+1),   r≥1.  (2)

In particular exactly one forbidden occurrence has series t²/D(t)². This enumerates an explicit combinatorial refinement of the insertion words at every length. By the audited matrix equivalence, r=0 is precisely the nonsingular/unimodular class and r≥1 refines the singular class. The number r is not claimed to equal nullity, rank deficiency, torsion rank or any Smith statistic.

**Proof.** Every legal word gives a unique state walk. Each forbidden pair ends exactly on the extra e-transition at B1, which supplies a factor u. Matrix inversion gives (1). Put B(t)=t(1−t)² and expand `(D−uB)^−1` in powers of u. For r≥1 its numerator is `B^(r−1)((1+t)B−tD)=t²B^(r−1)`, proving (2). The apparently signed numerator in (2) poses no positivity problem: its coefficients count words by (1). ∎

## 3. First forbidden-pair position factors into a safe prefix and a free suffix

Let g_j=[t^j]1/D(t), so g_0,g_1,g_2,...=1,2,5,11,25,... . Let q_j=[t^j]1/(1−3t+t²), so q_j=1,3,8,21,55,... . For a length-m legal word and 2≤k≤m, the number whose **first forbidden pair ends at letter k** is exactly

    g_(k−2) q_(m−k).                                      (3)

Equivalently the two-variable generating series, x marking the prefix through the first violation and y the remaining suffix, is

    x²/[D(x)(1−3y+y²)].                                  (4)

Thus the total rejected series is t²/[D(t)(1−3t+t²)], the difference between all legal and safe words.

**Proof.** A safe prefix of length k−1 must finish in B1 and then append e. Its series by prefix length is `α(I−tC)^−1 e3=t/D(t)`, so the number is g_(k−2). After this e, the underlying legality type is B. Its arbitrary legal suffix count is q_(m−k), independent of the earlier prefix; the flag of having failed no longer restricts legality. Multiplying gives (3), and summing gives (4). ∎

**Proof boundary.** The language argument proves persistence of having seen a forbidden substring, not persistence of a native matrix kernel by itself. The latter is now established separately by the [bridge–sink star argument](../bridge/BRIDGE_SINK_STAR_LEMMA.md) and its [independent audit](../audit2/REPORT.md). Formula (3) consequently also counts singular tournaments by the first forbidden-pair position in their canonical insertion word, using classical grammar completeness.

## 4. Letter-content refinement

Give letters a,b,c,d,e independent commuting weights A,B,Cc,Dd,E. Replace the matrix in (1) by

    M = [[A,0,B],[Cc,Dd,E],[Cc,Dd,uE]].

The weighted series is `α(I−tM)^−1 1`. Explicitly its numerator is

    (1−t(A−Cc))(1+t E(1−u)),

and its denominator is `det(I−tM)`. At u=0 this determinant is

    1−t(A+Dd)+t²(A Dd−B Cc−Dd E)+t³ A Dd E.

This preserves letter-content information discarded by scalar counts. It does not conflate operation-letter counts with path-run counts in M2.

## 5. Why the earlier symmetric boundary realization matched

The pre-recovery note constructed T and boundary vectors

    T=[[1,1,1],[1,1,0],[1,0,0]],  αT=(1,2,0), βT=(1,0,0)^t.

There is an explicit unimodular intertwiner

    S=[[1,0,1],[1,2,0],[1,1,0]],  det S=−1,
    C S=S T,   α S=αT,   S βT=1.

Hence `α C^m 1=αT T^m βT` for all m. This proves integral equivalence of these stationary **linear realizations**, not a word-level bijection. The inverse of S has negative coefficients; a basis transformation is not an operation-preserving identification of construction words or tournament isotypes. The actual five-letter grammar is now known and takes precedence over the guessed walk model.

## 6. Actual insertion-word rotation and reversal fail isomorphism

Two inherited source records at order 7 give a decisive rotation hostile:

|word|canonical cycle-normalized upper-triangle code|det plus|det minus|det M2|
|--|--|--:|--:|--:|
|ecb|100000100001010111111|1|−1|−1|
|cbe|100000100101011111111|64|0|0|

The words are rotations, both legal type B. Their tournaments cannot be isomorphic because nonsingularity differs. Independently of trusting these determinant values, `check_actual.py` decodes their tournament codes and checks every degree-compatible vertex permutation: no isomorphism exists.

Reversal fails already at order 6:

|word|canonical code|det M2|
|--|--|--:|
|de|100001000101111|1|
|ed|100001000100111|1|

The same independent permutation test finds no isomorphism. Thus even safe-word reversal is not a tournament-isomorphism operation. The determinant happens to agree in this pair; determinant equality is not a classifier. A further basic legality hostile is dc (legal from seed B) versus its rotation/reversal cd (illegal). All records and provenance are in `source_rotation_witness.json`.

These facts do not contradict rotation canonicalization of the **unique directed Hamiltonian cycle**. Insertion chronology and cycle-vertex coordinates are different encodings and their operations are not interchangeable.

## 7. A legitimate cyclic-language count, but not a tournament quotient

Define a cyclic grammar word of positive length m as a closed A/B type walk whose cyclic adjacent pairs, including last-to-first, avoid be and ee. This drops the seed-B boundary and replaces it with cyclic closure. Its rooted count is

    c_m=tr(C^m).

Each letter determines its legality source (a,b require A; c,d,e require B); the preceding cyclic letter determines B0 versus B1 when necessary. Thus each cyclic word has exactly one compatible closed productive-state walk. This unique-lift check is what licenses the trace.

The formal-language rotation-necklace count is consequently

    N_m=(1/m) Σ_(d|m) φ(d) tr(C^(m/d)).

Fixed rotations are repetitions of cyclic-safe blocks, so ordinary Burnside applies to this newly defined language. These numbers count **cyclic operation words**, not insertion-word orbits representing tournaments. For example the cyclic one-letter word a is included although a is not a legal first insertion from seed B. The genuine faithful cycle/chord-row encoding being developed separately must retain its own map and isomorphism proof.

## Verification

    python3 session103_words_necklaces/languages/check_actual.py
    python3 -O session103_words_necklaces/languages/check_actual.py

The checker has explicit failure gates, not removable assertions. It enumerates every legal word through m=10, independently checks (1) at u=0,1,2,3, (2) coefficientwise for every relevant r, and (3) at every position; checks cyclic trace/Burnside through m=8; verifies the intertwiner; and performs the two complete permutation-isomorphism hostile checks. It does **not** recompute observer determinants. Those determinant values retain source-provenance status. The new exact enumeration at m=9,10 remains language enumeration, not direct numerical observer evaluation at n=13,14. Its interpretation as exact observer counts at those orders now follows from the separately audited all-order theorem, not from a larger matrix computation.
