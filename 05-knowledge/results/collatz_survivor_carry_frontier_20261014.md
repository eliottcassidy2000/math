# Native survivor carries: a stable collision and a minimum-cost frontier

**Update, 2026-10-14:** the strict carry-entropy-gap hypothesis below is now REFUTED by the [polynomial carry-fibre theorem](collatz_polynomial_carry_fibres_20261014.md). Each fibre contains only O(L²) words, so the carry and word exponential rates are equal. The collision identities, exact compilers, finite counts and sourcewise boundaries remain valid.

Date: 2026-10-14. PROVED: the collision identity, its suffix propagation, the projected minimum-cost recursion, weighted-fibre recursion, binary injectivity and entropy bounds. FINITE-EXACT: a complete all-height classification through six odd steps. REFUTED: the proposed strict carry-entropy improvement. OPEN: consequences for sourcewise descent. These are new derivations in this session; external novelty is unclaimed.

## 1. The new signal

The [joint-kernel result](signed_collatz_joint_tower_kernel_20261014.md) turns the number of distinct ternary endpoint carries into an exact rank. This note investigates those carries for actual no-descent prefixes of the positive 3n+1 map.

At six odd steps there are 30 surviving native words but only 24 endpoint carry classes. The reduction 30→24 is the additional collision effect; comparing 24 with all 486 ternary unit classes also includes the much larger, elementary word-count restriction.

The first collision is explicit:

    w=(1,1,1,3),   A=6, B=65;
    v=(1,1,2,1),   A=5, B=73.

Both have endpoint carry 20 modulo 81, and

    Fv(n)−Fw(n)=81(n+1)/64.

This extends to an infinite, suffix-stable substitution family. It explains every carry collision observed through length six. It also yields a constructive compiler that retains only the minimum total valuation needed to reach a carry, when the target is projected reachability. Source measure requires a separate weighted fibre, which cannot be discarded.

## 2. Exact source inequalities, not only slopes

For a positive valuation word w=(a₁,…,aL), write each actual prefix as

    Fi(n)=(3ⁱn+Bi)/2^Ai,   Ai=a₁+…+ai.

Here B₀=0 and Bᵢ₊₁=3Bᵢ+2^Ai. Let Qᵢ=2^Ai and Pᵢ=3ⁱ. “No descent through L” means Fi(n)≥n for every i≤L, with original source n≥3 odd. It does not mean that the orbit is increasing at every individual step.

If Pᵢ>Qᵢ, that prefix inequality holds for every positive source. If Pᵢ<Qᵢ, it instead requires

    n≤Bi/(Qᵢ−Pᵢ).

Thus the exact upper source bound is the minimum of these rational thresholds over the contracting prefixes, or infinity if there are none. The native word sources form one odd class r mod 2^(AL+1), with

    r=(2^AL−BL)(3ᴸ)^(−1) mod 2^(AL+1).

The cylinder is authenticated by replay in the checker. Let n₀ be its least source ≥3. The word supplies actual no-descent sources exactly when n₀ is at most the upper bound. A word with a contracting prefix can in principle supply a finite source family; this possibility must be tested rather than removed by definition.

Call a word *slope-admissible* if 2^Ai<3ⁱ for every prefix. It supplies an infinite native no-descent family. This is the class used by the projected recursions below.

### Why the six-step census covers every positive source

For any positive odd source, U(n)≤(3n+1)/2. Before step i this gives

    nᵢ₋₁≤(3/2)^(i−1)(n+1)−1.

If nᵢ≥n and n≥3, then

    2^aᵢ≤(3nᵢ₋₁+1)/n
          ≤4(3/2)^(i−1)−2/3.

For i≤6 the right side is less than 32. Hence every relevant valuation is at most 4. Exhausting {1,2,3,4}^L and applying the exact native source thresholds therefore covers all odd sources n≥3, not merely sources below a numerical cutoff.

The result is:

| L | Native no-descent words | Distinct carries | All ternary units | Relative odd-source density |
|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 2 | 1/2 |
| 2 | 2 | 2 | 6 | 3/8 |
| 3 | 3 | 3 | 18 | 1/4 |
| 4 | 7 | 6 | 54 | 13/64 |
| 5 | 12 | 10 | 162 | 19/128 |
| 6 | 30 | 24 | 486 | 1/8 |

Every admitted word in this range is slope-admissible. No finite-source-only exception occurs through length six. This equality of the two classes is not asserted at arbitrary length.

## 3. A substitution that survives every common suffix

For m≥0 and b≥3 compare the equal-length words

    w=1^(m+1),b,       v=1^m,2,b−2.

Their total valuations differ by one, Aw=Av+1. Since the valuation-one branch iterated m times gives

    F1^m(n)=−1+(3/2)^m(n+1),

and

    F(2,b−2)(x)−F(1,b)(x)=9(x+1)/2^(b+1),

we obtain the exact affine identity

    Fv(n)−Fw(n)=3^(m+2)(n+1)/2^(m+b+1).

Consequently their constant carries Bw2^(−Aw) and Bv2^(−Av) agree modulo 3^(m+2).

Append any common valuation suffix z of length k and total valuation Az. Its affine slope is 3ᵏ/2^Az, so

    F(v,z)(n)−F(w,z)(n)
      =3^(m+2+k)(n+1)/2^(m+b+1+Az).

The carry collision therefore persists at the full new depth. Algebraically, 2B(v,z)−B(w,z)=3^(m+2+k).

For b=3 and every m≥2, both words are slope-admissible. More generally, whenever the higher-cost word w is slope-admissible, the alternative prefix v is also admissible: its intermediate 2-step is legal for m≥1, and its total valuation is smaller by one. Any common suffix that remains admissible for w remains admissible for v.

This is a real projected-subtree inclusion. It does not authorize swapping the words at the same integer source: their native source cylinders are disjoint. Nor does it automatically yield exponential compression. If a prefix p of length h is prepended, the full-depth collision holds exactly when Bp+2^Ap≡0 mod 3ʰ, equivalently when the prefix carry is −1 mod 3ʰ. Thus arbitrary prepending is not available. In this family the first non-one letter changes, so repeated independent substitutions along one word are not supplied by the identity alone.

### The weight that a quotient must preserve

A native word of total valuation A has relative odd-source mass 2^(−A). The first collision therefore combines masses 1/64 and 1/32 into one carry fibre of mass 3/64. Retaining only the cheaper representative would incorrectly report mass 1/32.

After a common suffix of total valuation C, the correct merged mass is 3·2^(−6−C), not 2^(−5−C). At length six, replacing every fibre by its minimum-cost representative would lose 9/512 of source mass: the correct density is 1/8, while the incorrect value would be 55/512.

## 4. Two exact projected compilers

Let bL=floor(log₂(3ᴸ)), and let E_L(r) be the smallest total valuation among slope-admissible length-L words with endpoint carry r mod 3ᴸ. Initially E₀(0)=0 modulo 1. For every current state and a≥1, form

    r′=(3r+1)2^(−a) mod 3^(L+1),
    A′=E_L(r)+a.

Keep this extension when A′≤b_(L+1), and minimize A′ at each r′.

This computes the reachable carry classes exactly. A higher-cost word with the same carry has exactly the same projected next-carry choices, but less remaining valuation budget. Any future suffix admitted from that word is admitted from the lower-cost representative. This dominance statement concerns projected reachability only; it does not preserve a source identity or the density of its cylinder.

To retain density, use the polynomial

    H_L,r(z)=Σ_{words w in fibre r} z^A(w).

For every coefficient c_A and positive a with A+a≤b_(L+1), add c_A z^(A+a) to H_(L+1),r′ at the same r′ above. Evaluating at z=1 counts words; at z=1/2 it gives their exact source mass. Keeping all costs matters because different costs have different remaining budgets.

Both compilers independently reproduce the six-step exhaustive frontier and all weighted fibre totals.

## 5. Where compression can enter, and the sharpened hypothesis

For valuations restricted to {1,2}, the terminal carry is an injective code of the entire word. Modulo 3 it identifies the last valuation: carry 2 means a=1, and carry 1 means a=2. Then

    r_previous=(2ᵃr−1)/3 mod 3^(L−1)

recovers the preceding carry. Repeat. More generally, injectivity holds for any fixed alphabet containing at most one odd valuation and at most one even valuation. The observed collisions require ambiguity between different valuations of the same parity, here 1 and 3.

Let W_L count slope-admissible words and R_L count their distinct carry classes. Put α=log₂3 and

    β=α^α/(α−1)^(α−1)≈2.8395137305.

The ordinary positive-composition bound gives W_L≤binom(floor(αL),L). Conversely, among compositions of floor(αL) into L positive parts, each cyclic rotation orbit contains a word whose every prefix average is at most the total average: rotate immediately after a maximum of the centered partial sums. It is slope-admissible. Hence

    W_L≥binom(floor(αL)−1,L−1)/L.

Stirling's formula gives lim W_L^(1/L)=β. This is the elementary word-count rate; its comparison with 3 is not a new collision effect.

For a lower bound on carries, binary words in which every prefix has at least as many 1s as 2s have total valuation at most 3i/2 at each prefix. Since 3/2<α, they are slope-admissible. There are binom(L,floor(L/2)) such words, and binary injectivity makes all their carries distinct. Therefore

    2≤liminf R_L^(1/L)≤limsup R_L^(1/L)≤β.

**Historical strict-compression hypothesis, now REFUTED:**

    limsup R_L^(1/L) < β.

The successor now proves equality with β despite the finite collisions: every fibre is polynomially bounded. The observed 30→24 reduction and the suffix-stable family remain a mechanism and finite loss, not an exponential gap. The minimum-cost and weighted-polynomial recursions remain useful for the finer, unresolved distribution of fibre sizes and source weights.

Even a strict gap would concern a ternary shadow. It would not establish that an individual integer descends.

## 6. A mandatory hostile: the coherent Mersenne branch

For every L, the word 1ᴸ survives and has carry −1 mod 3ᴸ. Its least positive native source is 2^(L+1)−1, and its endpoint is 2·3ᴸ−1. These sources really stay above their respective starting values through the prefix.

Thus an arbitrarily thin carry frontier still contains the coherent ternary branch −1 at every level. The finite positive witnesses change with L and tend to infinity. The infinite valuation-one word would require n≡−1 mod 2^(L+1) at every L, whose 2-adic source is −1, not a positive integer. A compactness argument that forgets the original source would retain this counterfeit positive branch.

The new compilers are useful precisely because their loss is explicit: they compute carry reachability and weighted finite source ensembles, not the existence of one positive source following an infinite no-descent path.

## 7. Reproduction

The standard-library checker exhausts only lengths 1–6, using the proved valuation cap 4 and exact rational source thresholds. It independently implements the two projected recursions, verifies binary decoding controls, and checks 48 instances of the general collision identity with common suffixes. All nine non-singleton carry fibres observed across these lengths are explained by that identity.

Normal and optimized Python outputs are byte-identical; 73,020 checks pass. The check count comes mostly from authenticating the native cylinders of the bounded word universe, not from a large integer search.

    python3 04-computation/experiments/collatz_survivor_carry_frontier_20261014.py > 05-knowledge/results/collatz_survivor_carry_frontier_20261014.json
    python3 -O 04-computation/experiments/collatz_survivor_carry_frontier_20261014.py > 05-knowledge/results/collatz_survivor_carry_frontier_20261014.optimized.json
    cmp 05-knowledge/results/collatz_survivor_carry_frontier_20261014.json 05-knowledge/results/collatz_survivor_carry_frontier_20261014.optimized.json

Canonical output-payload SHA-256, excluding its hash field:

    0f6257ef2760ce4e27b1749876c98b3ae733f7759a23b4ef67cd77f150a204ff

No claim is made that the collision family is absent from prior affine-word or sibling-ladder work. The new result here is its explicit use as a weighted carry-frontier compression, connected to the exact rank counter, with a concrete entropy question and a source-preserving boundary.
