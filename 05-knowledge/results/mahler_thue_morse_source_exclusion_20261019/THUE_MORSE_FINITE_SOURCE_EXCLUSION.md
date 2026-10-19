# Every non-erasing binary block image of Thue–Morse misses positive Mahler integers

2026-10-19. **PROVED by the elementary argument below; finite exact controls supplied.** Independently [reviewed](review/INDEPENDENT_REVIEW.md). No external-priority claim. The Mahler 3/2 problem remains open.

## Scope and inheritance

Let M(N)=ceil(3N/2)=(3N+e)/2, with e=N mod2, on positive ordinary integers. The carry word records e at each step. Let t=011010011001... be the Thue–Morse word, the fixed point beginning0 of the substitution sigma(0)=01, sigma(1)=10. Fix arbitrary nonempty binary words U,V, and replace each0 in t by U and each1 by V. Prepend any fixed finite binary word P.

**Theorem.** No positive ordinary integer has carry word P h(t), where h(0)=U and h(1)=V. For each proposed source A there is an explicit finite first-mismatch certificate. Unequal block lengths and commuting blocks are both allowed.

This extends the repository's Mahler carry/repetition viewpoint in THM-4654, [THM-4654](../../../01-canon/theorems/THM-4654-safe-mahler-codes-contract-collatz-logarithmic-memory-and-complexity-approximation.md), and the carry interface [the carry interface](../mahler_carry_interface_20261014i.md) to a uniform morphic template. It neither excludes all automatic sequences nor proves that all safe formal words miss integers. The proof uses a repeated prefix of one particular substitution fixed point; no transcendence theorem is required.

## 1. Carry cylinders and periodic shadows

For a word w=e_0...e_(p-1), define

    C(w)=sum_{j=0}^{p-1} e_j 3^(p-1-j) 2^j.

Following this word gives M^p(x)=(3^p x+C(w))/2^p whenever those are the actual carries. The periodic word w^infinity has the 2-adic source

    y=-C(w)/(3^p-2^p).

The denominator is odd. Its real value lies in [-1,0], because 0<=C(w)<=3^p-2^p. This is a 2-adic source statement, not an ordinary positive orbit.

Any two 2-adic sources have equal first L carry bits exactly when they are congruent modulo2^L. One proof is inductive: equal parity makes their first-image difference 3/2 times their difference, reducing its 2-adic valuation by one. A different parity is exactly the first mismatch. Thus, for a positive integer B, the first carry mismatch with w^infinity occurs at the zero-based position

    v_2((3^p-2^p)B+C(w)).

This integer is positive and nonzero, so the valuation is finite. This also proves uniqueness and existence of finite carry cylinders; the periodic fixed point indeed has the designated first word by its defining congruence and hence repeats it.

## 2. A repeated Thue–Morse prefix with a favorable size ratio

The word01101 is both a prefix of t and a prefix of (011)^infinity. Applying sigma^j and then h preserves both prefix facts. For j>=1 each sigma^j(letter) contains exactly2^(j-1) zeros and the same number of ones. Put

    L=|U|+|V|,  k=2^(j-1)L,
    w=h(sigma^j(011)).

Then |w|=3k, while h(sigma^j(01101)) has length5k and is a common prefix of h(t) and w^infinity. Equal-count substitution is the reason unequal block lengths cause no loss.

## 3. Explicit finite horizon

Let a=|P| and A>=1. If A already mismatches P, stop at that first mismatch. Otherwise write B=M^a(A). Since M(x)+1<=(3/2)(x+1),

    B+1 <= (3/2)^a(A+1).

Choose j>=1, equivalently k in {L,2L,4L,...}, so that

    2^(a+5k) >= 3^(a+3k)(A+1).                         (1)

Such j exists because2^5/3^3=32/27>1. For the periodic shadow word w above, put D=3^(3k)-2^(3k). Then

    0<DB+C(w) <= D(B+1)
                 < 3^(3k)(3/2)^a(A+1)
                 <= 2^(5k).

Therefore its 2-adic valuation is strictly smaller than5k. The first mismatch with w^infinity occurs before5k, while w^infinity and h(t) agree through that position. Consequently A and P h(t) have a first mismatch before total length a+5k. Its exact position is

    a + v_2(DB+C(w)).

This proves the theorem and gives a terminating source-dependent certificate. No search over arbitrary orbit lengths or asymptotic density inference is used.

The certificate is logarithmic in A for fixed blocks and prefix. The source dependence is essential: every finite carry prefix has an integer residue-class realization, so there is no fixed mismatch horizon uniform over all positive sources. This certifies departure from the specified template, not necessarily rejection by the safe follower: the actual integer orbit may continue along a different safe word. The all-depth exclusion is precisely that it cannot retain this template forever.

## 4. A genuinely safe excluded family

Take U=000 and V=100. In every concatenation of these blocks, ones are separated by at least three positions. For every suffix e_0e_1..., its weighted tail satisfies

    sum_{j>=0} e_j (2/3)^(j+1)
      <= (2/3)/(1-(2/3)^3)=18/19<1.

This is the strict Mahler safe-tail condition in the inherited carry model. Thus h(t) is a strictly safe formal word, with a uniform safety margin, yet is not the carry word of any positive ordinary integer. It is an actual formal safety-vs-integer-entry obstruction, rather than a template already rejected by its first unsafe pattern.

The same elementary source theorem applies to any nonempty U,V; other blocks need not be safe. When blocks are loops from the initial follower state0 back to0 and reset along every loop, their concatenations also give strict safe formal paths, by the reset-Buechi equivalence in THM-4072, [THM-4072](../../../01-canon/theorems/THM-4072-mahler-safe-terminal-fibre-product-and-finite-state-obstruction.md), Section1. The finite-source proof itself does not require this extra loop hypothesis.

## 5. Reproduction and separate analytic extension

Run `python3 -B mahler_thue_morse_reset_exclusion.py` and the same command with `-O`. The checker compares certificates with actual positive-integer trajectories, verifies substitution overlaps, and checks selected safe-loop controls. Checks remain active under optimization. The exact outputs are retained beside the script.

The [scalar extension](THUE_MORSE_SCALAR_MAHLER_EXTENSION.md) supplies the stronger rational-versus-transcendental dichotomy, an explicit scalar-theorem hypothesis audit, and a uniform irrationality-exponent bound. Its cited analytic dependency is unnecessary for the elementary finite-source theorem above. Neither theorem excludes arbitrary automatic Mahler carry words.

## 6. Exact audit universe

The standalone checker passes18,802 always-active checks, with byte-identical normal and optimized receipts. It covers all900 ordered nonempty binary block pairs with individual lengths at most4; 13,323 source-specific first-mismatch comparisons; all36 pairs from the six canonical reset loops of lengths1,3,9,12,15,17; and moving positive-integer cylinder representatives at safe-template depths8,16,32,64,128,256. The moving-source controls retain the integer-entry obstruction rather than interpreting finite safety as an infinite candidate.
