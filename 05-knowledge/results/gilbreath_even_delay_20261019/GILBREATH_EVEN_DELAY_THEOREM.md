# A seven-state paired invariant proves that every finite size-six delay is even

**Session:** 2026-10-19, Gilbreath / size-six lane.  
**Status:** PROVED, by an elementary invariant and an independent all-length finite-state closure certificate.  
**Scope:** A single normalized defect 3, binary left sea, and zero right boundary. No prime-specific conclusion and no literature-priority claim.

## 1. Inheritance and the exact claim

Use the definitions and rotated-coordinate identity from [the size-six proof](../gilbreath_alphabet_three_size_six_20261019.md), Sections 4–5. For a binary left word w of length L=F−1, append the normalized defect 3 and then zeros. Let tau(w) be the first number r of appended zeros for which the resulting finite word has a dangerous left edge, with tau=∞ if there is no such r.

The Pascal-front bijection produces a binary word alpha of length L. Its rotated age vector is obtained from v=(3,3,3,...) by successively applying

    T_b(v)_j = |T_b(v)_(j−1) − v_j|,    T_b(v)_(−1)=b,    b in {0,1}.

Its coordinate j is the physical left-edge value at time L+j. All earlier left-edge values depend only on the binary left sea. Thus tau is the first age j at which the final vector has a value at least 2.

**Theorem.** For every F≥1 and every binary left word w of length F−1,

    tau(w) is in {0,2,4,6,...} or is infinite.

In particular, A_r(F)=0 for every odd age r and every F, where A_r(F) counts left words with tau=r.

The proof establishes a stronger statement: after the first front-input 1, every age vector satisfies

    v_(2j+1)=2  implies  v_(2j)=2       for every j≥0.

The previous finite checks at odd ages 1 through 13 are therefore replaced by an all-age theorem. The numerical sharp upper bound on finite tau and the all-age counting recurrence remain separate open questions.

## 2. The pair language

Group the age coordinates into consecutive pairs, beginning at age zero. Let

    S = {00,01,10,11,20,21,22},
    Z = 00,    A = 20,    B = 22.

Let C be the class of finite or one-sided infinite pair words with letters in S and containing neither of the patterns

    A Z* B,        B Z* A.

Here Z* means any number of consecutive copies of the whole pair 00, including zero copies. Other pair letters interrupt a pattern.

Every vector in C has the desired implication v_(2j+1)=2 ⇒ v_(2j)=2, since the pairs 02 and 12 are excluded. That implication alone is not an invariant. For example, the vector 2022 has legal pairs 20,22, but its T_0 image is 2202, which has a forbidden pair 02. The two forbidden pattern families supply the missing information.

## 3. Complete pair-transducer table

If q is the incoming value immediately before an input pair xy, the output pair is

    (|q−x|, ||q−x|−y|),

and the outgoing value is its second coordinate. Since all inputs, outputs, and incoming values after the first front-input 1 lie in {0,1,2}, the following table is complete.

| Input pair | Incoming q=0 | Incoming q=1 | Incoming q=2 |
|---|---|---|---|
| 00 | 00 | 11 | 22 |
| 01 | 01 | 10 | 21 |
| 10 | 11 | 00 | 11 |
| 11 | 10 | 01 | 10 |
| 20 | 22 | 11 | 00 |
| 21 | 21 | 10 | 01 |
| 22 | 20 | 11 | 02 |

**Invariant lemma.** For b=0 or 1, T_b sends C into C.

**Proof, first part: the output uses S.** The only possible output outside S in the table is 02, arising from input B=22 with incoming state 2. Suppose there is a first such output. The incoming state at the very beginning is b≤1. Before the offending pair, the only transition from a state other than 2 into state 2 is input A=20 with incoming state 0, producing output B=22. Once the state is 2, the only transition that keeps it at 2 without already producing a forbidden output 02 is input Z=00. Consequently, the first forbidden output would require the input pattern A Z* B. This is excluded. Therefore every output pair belongs to S.

**Second part: the output avoids A Z* B.** An output A=20 can only come from input B=22 with incoming state 0, and it leaves state 0. With incoming state 0, output Z=00 can only come from input Z=00 and leaves state 0. Finally, with incoming state 0, output B=22 can only come from input A=20. Thus an output pattern A Z* B forces the input pattern B Z* A, which is excluded.

**Third part: the output avoids B Z* A.** An output B=22 leaves state 2, whereas output A=20 requires incoming state 0. Thus B cannot be immediately followed by A. If at least one output Z=00 lies between them, the first such Z, entered in state 2, can only come from input A=20; it leaves state 0. Every remaining output Z then comes from input Z, and the final output A comes from input B. Hence the input contains A Z* B, which is excluded. This proves the lemma. ∎

The binary incoming-state restriction is essential: input B alone is in C, but incoming state 2 maps it to 02.

## 4. The first front-input 1 enters the invariant

Suppose z zeros precede the first 1 in alpha. Before that first 1, the age vector is

    3 (binom(z+j,j) mod 2),       j≥0.

Indeed, on the alphabet {0,3}, T_0 is prefix XOR after division by 3; applying prefix XOR z times to the constant-one sequence gives these binomial coefficients.

Applying T_1 to a {0,3}-valued vector keeps every output in {1,2}: an input 0 leaves the current value unchanged, and an input 3 exchanges 1 and 2. Therefore the vector immediately after the first 1 is exactly

    v_j = 1 + (binom(z+j+1,j) mod 2).                              (1)

This is a prefix sum of the previous binomial coefficients over F_2, using the hockey-stick identity.

**Indexing correction.** The numerator in (1) is z+j+1, not z+j. The case z=0 gives 212121..., making the extra 1 necessary.

For m=z+1, the elementary binary coefficient identity is

    binom(m+j,j) is odd  iff  (m bitwise-AND j)=0.                  (2)

One proof uses the formal power-series identity over F_2

    sum_j binom(m+j,j) x^j = (1+x)^(−m−1).

Alternatively, Lucas's elementary binomial parity criterion gives (2), equivalently that the binary addition of m and j has no carries. The conclusion needed here can also be read directly from Pascal parity: if m is even, the coefficients at j=2k and j=2k+1 agree; if m is odd, the coefficient at every odd j is zero.

It follows that:

- If z is even, then m is odd and every age pair is 11 or 21.
- If z is odd, then m is even and every age pair is 11 or 22.

Both types lie in C: neither contains A=20, so neither forbidden pattern can occur. Every subsequent binary front input preserves C by the invariant lemma.

For completeness, the parity splitting just used needs no named theorem. Write H_m(x)=(1+x)^(−m−1) in F_2[[x]]. When m=2a+1, H_m(x)=(1+x^2)^(−a−1), so all odd coefficients vanish. When m=2a, H_m(x)=(1+x)(1+x^2)^(−a−1), so each adjacent even/odd coefficient pair agrees.

## 5. Concluding the delay theorem

If alpha contains no 1, its first age coordinate remains 3, so tau=0.

Otherwise, after its first 1 all coordinates are in {0,1,2}, and the final vector belongs to C. If its first dangerous value occurred at an odd age 2j+1, that value would be 2. The pair restriction would force value 2 already at age 2j, a contradiction. Thus every finite delay is even. This includes alpha of length zero, which has tau=0. ∎

A useful additional consequence is that, after a 1 has appeared, the first dangerous age never decreases under another front input. The prefix preceding the first 2 is binary, so its image remains binary. Of the two boundary choices b=0,1, exactly one preserves that first dangerous age: their running prefix-XOR values are complementary just before the first 2, and subtracting 2 from 0 gives 2 whereas subtracting it from 1 gives 1. The other branch either has no dangerous age or advances the first dangerous age by at least 2. This last conclusion does not assert that the other branch loses every later 2; that stronger assertion requires additional structure.

## 6. Independent finite certificate and controls

`verify_paired_invariant.py` supplies a second proof of the all-length closure lemma. The pattern-avoidance language has a three-state deterministic monitor:

- State 0: neither partial pattern is active.
- State 1: the suffix is A followed by zero or more Z letters.
- State 2: the suffix is B followed by zero or more Z letters.

A Z leaves the monitor state unchanged. An A rejects state 2 and otherwise moves to state 1; a B rejects state 1 and otherwise moves to state 2. Every other pair in S resets it to 0.

The certificate takes the product of the input monitor, the output monitor, and the transducer state q. It starts at (0,0,0) and (0,0,1), allows exactly the transitions accepted by the input monitor, and checks that no reachable transition is rejected by the output monitor. Only seven product states and 44 allowed edges are reachable. These finite checks prove closure for pair words of arbitrary length, since every such word labels a path in this finite product graph. The complete graph is included in `paired_invariant_certificate.json`.

Additional independent checks are deliberately bounded and are not substituted for the proof:

- Direct first-one construction versus formula (1), for all z=0,...,63 and 128 age coordinates each.
- Original shrinking difference triangles versus the Pascal-front rotated construction, for all 2,047 binary left words of lengths 0 through 10, at 32 age coordinates each.
- Hostile controls: input 2022 maps to 2202 under T_0; input 2220 maps to 2022 and then 2202; incoming boundary state 2 maps legal input pair 22 to illegal pair 02.

The latter two four-cell words show why retaining only the pair alphabet, or retaining only one of the two pattern exclusions, does not yield this closure proof.

Reproduction:

    python verify_paired_invariant.py > paired_invariant_certificate.json
    python -O verify_paired_invariant.py > paired_invariant_certificate_optimized.json
    cmp paired_invariant_certificate.json paired_invariant_certificate_optimized.json

All checks use explicit failures, remain active under optimization, and require only the Python standard library.

## 7. Connection and boundary ledger

**Source:** The exact rotated age sequence for a lone normalized defect 3.  
**Target:** A regular language of age pairs.  
**Map:** Group coordinates (0,1), (2,3), and so forth; retain pair values and whether a suffix is A Z* or B Z*.  
**Preserved predicate:** No first dangerous age is odd; in fact an odd-coordinate 2 always has an even-coordinate 2 immediately before it.  
**Discarded information:** Most of the front word, multiplicities of front words, dyadic periods, and exact later defect ages.  
**Required sidecar:** The two pattern-avoidance monitor states, beyond the local pair alphabet.  
**Cheapest hostile:** 2022, which has legal pairs but maps to 2202.  
**Strongest result proved here:** All finite delays are even at every front length, with an invariant valid on arbitrary one-sided sequences in C.  
**Still open here:** Sharp dyadic horizon, paired-halving at all front lengths, complete all-age recurrence, and a summable tail estimate establishing the repaired asymptotic.
