# Paired halving: an exact conditional route and the missing reachable-state property

**Session:** 2026-10-19, bounded structural follow-on.  
**Status:** PROVED phase criterion and conditional implications; FINITE-EXACT evidence for the new hypotheses; OPEN all-F paired halving and sharp horizon. This note does not modify the frozen proof packet.

## 1. Result and stopping point

The new even-delay theorem gives a genuine invariant on age pairs, but that invariant does not determine either the paired-halving counts or the sharp dyadic decision horizon. The precise missing information is the nonlinear placement of the value 2 relative to the parity trace of the reachable age vector.

A useful stronger hypothesis is now isolated:

**OPEN even-front annihilator hypothesis.** Let an even number L of binary front inputs act on the initial age vector `(3,3,3,...)`. If at least one input is 1, the resulting vector v has a boundary bit b such that T_b(v) is entirely binary.

This statement would prove the all-F paired-halving claim of HYP-9163. It is strictly more informative than merely preserving even first-defect ages. It has been checked on exact infinite age sequences for every front word through length 14. No proof for arbitrary L is claimed.

## 2. Exact annihilator criterion

Let v be any finite or one-sided infinite word over `{0,1,2}`. Define its cumulative parity before position j by

    h(j) = sum_{i<j} (v_i mod 2) modulo 2.

Call a bit b an annihilator if every entry of T_b(v) is 0 or 1.

**Lemma.** An annihilator exists if and only if h(j) takes a single value on all positions j where v_j=2. If there is at least one 2, the annihilator is unique and equals `1+h(j)` modulo 2. If no 2 occurs, both bits are annihilators.

Equivalently, between each two consecutive occurrences of 2 there must be an even number of entries equal to 1.

**Proof.** Until a dangerous output occurs, the running output is binary and its value immediately before reading v_j is `b+h(j)` modulo 2. An input 0 or 1 preserves binary output. An input 2 preserves binary output exactly when the preceding output is 1, since `|1−2|=1` while `|0−2|=2`. Thus every occurrence of 2 imposes the condition `b+h(j)=1`. These conditions are jointly satisfiable exactly when their h(j) values agree. Induction along the word proves sufficiency, and a binary output necessarily satisfies every condition. ∎

The lemma is valid for arbitrary infinite words: it is a local induction, not a periodicity assumption.

## 3. Why the hypothesis proves paired halving

Let A_r(L) count binary front words of length L whose first dangerous age is r, and let r>0. The Pascal transform is bijective, so these are also the original left-word counts at distance F=L+1.

After the first input 1 the age vector is in `{0,1,2}`. If its first 2 is at age r, exactly one next bit retains that first 2: before it, the two outputs are complementary binary traces. Under the even-front annihilator hypothesis, the other next bit eliminates every 2, not merely the first one. An already binary vector has two binary successors. Under the same annihilator hypothesis, a vector in `{0,1,2}` with first defect at age zero cannot create a positive finite age: its two branches have ages zero and infinity. This is not a claim about arbitrary ternary vectors without that hypothesis. Finally, the all-zero front word has value 3 at age zero, and either extension still has a dangerous age zero.

Therefore, at every even L, each positive-finite-age word has exactly one successor with the same age, and there are no other positive-finite-age successors. This is a literal prefix-extension bijection:

    A_r(L+1) = A_r(L)       for every r>0 and every even L.

In the notation of the main proof packet,

    c6(F) = 2 sum_{r>0} A_r(F−1) 2^(−r).

Consequently `c6(F+1)=c6(F)` for odd F, equivalently

    E(F+1) = E(F)/2,       E(F)=p6(F)−F2^(1−F).

This conditional proof includes exclusion of new positive-delay arrivals from age-zero parents. Checking annihilation only for parents whose delay is already positive would omit that requirement.

## 4. Why the new even-delay invariant is insufficient

Use the pair language C from `GILBREATH_EVEN_DELAY_THEOREM.md`: allowed pairs are `00,01,10,11,20,21,22`, and the two forbidden patterns are `20 (00)* 22` and `22 (00)* 20`.

### A short hostile to annihilation

The infinite word

    00212000...

belongs to C and has first defect at age 2. Its two next inputs have first dangerous ages 2 and 4. The two occurrences of 2 have one intervening 1, so the phase lemma says neither bit can annihilate all defects. The six-entry displayed prefix is the shortest pair-aligned shape with a positive first defect and a later 2 separated by one 1. This is a hostile to an inference from C, not a counterexample to the reachable-front conjecture.

### Adding a correct dyadic period and parity complexity still does not suffice

Consider

    v = (00102120)^infinity.

Its pair cycle is `00,10,21,20`, so it belongs to C, including across the periodic boundary. Its two next inputs have first dangerous ages 6 and 4. Again neither annihilates all 2s.

The parity cycle is `00100100`, with exact binary linear complexity 7 and least period 8. These are exactly the parity constraints for front length L=6. Indeed the actual reachable front word

    alpha = 101100

has age vector `(00100100)^infinity`, which is binary. Thus the legal hostile and the genuine reachable vector have the **same complete infinite parity trace**, the same period, and the same front-length parity sidecar; only their placements of 2 differ. The pair invariant cannot distinguish them.

The checker enumerates all ternary cycles of dyadic periods 2 and 4 and finds no hostile with positive first defect, odd parity complexity, and failure of the phase criterion. Period 8 supplies the displayed example. This is minimality only in that explicitly bounded family, not a general classification theorem.

### The annihilator property is not a two-step invariant on its own

The word `21012` has an annihilator, because its 2s have two intervening 1s. But

    T_0(21012) = 21102,
    T_0^2(21012) = 21002,

and the last word has an odd number of 1s between its 2s. Starting with 21012 followed by zeros, the exact first output is 21102 followed by constant twos, and the next output is 21002 followed by repeating 02. These one-sided sequences also lie in C. Thus induction by two front steps cannot close on C plus the annihilator property alone.

## 5. A second precise route: exact age period would imply the sharp horizon

There is an elementary parity constraint that holds at every front length. Let p_L(z) be the generating series over F_2 of the parity of the age vector after front bits `b_0,...,b_(L−1)`. Prefix summation gives

    p_L(z) = [1 + sum_{s=0}^{L−1} b_s (1+z)^(s+1)] / (1+z)^(L+1).

The numerator equals 1 at z=1. Hence its denominator has exact pole order L+1 there. If

    Q = 2^ceil(log_2(L+1)),

then `(1+z)^Q=1+z^Q`, so the parity sequence is purely periodic with **exact** least period Q. The pole order also proves that any period must be divisible by Q. Equivalently, its binary linear complexity is exactly L+1.

**OPEN exact-age-period hypothesis.** The full integer age vector, not just its parity, is purely periodic with period Q for every binary front word of length L.

Together with the proved even-delay theorem, this hypothesis would give the sharp upper bound

    finite tau <= Q−2,       for L>=1.

The uniform bound, including L=0, is `finite tau <= 2 floor((Q−1)/2)`. At L=0 the initial constant-3 vector has Q=1 and tau=0. For L>=1, Q is even, so the uniform bound becomes Q−2: any first defect lies among ages 0 through Q−1 and cannot have odd age. This proves the proposed upper bound conditionally; it does not prove attainment of the bound at every F.

The exact finite tests below certify pure period Q for every front word through L=14. They do not supply an all-length theorem.

Again, C plus the entire correct parity trace does not force this conclusion. Start with the genuine binary trace `(00100100)^infinity` of the length-6 front 101100, and change only its age-8 value from 0 to 2. The result still belongs to C, because the only nonbinary pair is 20 and there is no 22. Its parity is unchanged, but its first defect is at age 8, beyond the proposed Q−2=6 bound. Moving that isolated 2 to any age divisible by 8 gives arbitrarily late first defects with exactly the same parity trace. None of these altered traces is asserted reachable.

## 6. Exact evidence without finite-preview assumptions

`verify_even_front_annihilator.py` represents an infinite sequence by a finite prefix followed by a repeating cycle. For one T_b application it tracks `(input-cycle position, current running value)`. A repeated pair gives an exact eventual cycle; no guessed age cutoff is used. It then removes redundant period copies and pulls back any periodic prefix. Starting from the constant-3 cycle, this constructs the complete infinite age vector for every finite front word being tested.

The receipt covers:

- 9,840 periodic ternary controls for the general phase criterion;
- all 32,767 front words of lengths 0 through 14;
- the annihilator condition at every even front length in that range, including age-zero parents;
- pure periodicity with the expected dyadic period for every tested reachable vector;
- the short and period-eight hostiles, their exact next-bit first defects, and the actual length-6 source with identical parity;
- bounded hostile minimality among dyadic periods 2, 4 and 8.

Normal and optimized runs agree exactly:

    python verify_even_front_annihilator.py > even_front_annihilator.json
    python -O verify_even_front_annihilator.py > even_front_annihilator_optimized.json
    cmp even_front_annihilator.json even_front_annihilator_optimized.json

## 7. What remains missing

This route has reached a precise stopping point. The pair invariant records enough information to forbid odd first defects, but it forgets nonlinear reachability and multiplicity. The phase criterion identifies exactly what a paired-halving proof would need, while the same-parity hostiles show why neither a local pair argument nor the linear Pascal trace alone is sufficient. A successful continuation must constrain the high bit of the reachable age vector, or give a different direct bijection on front words. No all-F paired-halving, all-age recurrence, or sharp-horizon closure is claimed here.
