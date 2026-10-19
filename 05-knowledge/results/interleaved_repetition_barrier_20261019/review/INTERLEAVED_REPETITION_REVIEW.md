# Independent review: interleaved repetition barrier

Date: 19 October 2026.

## Verdict

**The primary theorem passes this independent proof audit, using the cited Matveev theorem as an external input.** I found no gap in the passage from a lawful repeated expansion factor in an arbitrary earlier context to strict endpoint descent. The proof genuinely allows an unbounded number of expansion/reset episodes.

The initial audit identified one nonfatal scope issue: “every first-crossing bank word of length at most 2^52 descends” and the strict carry statement need **expansion-first** scope. If arbitrary words in `{K,V,F,H}*` are included, `W=F`, `x=1` is a lawful first coefficient crossing with endpoint 1, and its normalized carry is exactly `N/6`. The author has now made the expansion-first restriction explicit in the global definition, standalone lemma, and README, and explicitly excludes this example. I checked the corrected frozen text. **The issue is resolved; no outstanding mathematical objection remains.** The main repeated-factor theorem was unaffected throughout.

The draft’s repeated-expansion theorem, realized-endpoint exclusions, unbounded-episode example, and stated open boundary problem otherwise have appropriate scopes. Neither finite enumeration nor the cited input is being presented as a proof of general Collatz convergence.

## Material audited

- `INTERLEAVED_REPETITION_BARRIER.md`, final frozen SHA-256 `b3ea155d92d0071be52de543651a812e76515b0bda2961a84258866c63ce496d`. The initial reviewed version was `9abeb5538afd8b56c913e826452d99187efffa559a261489645dcbba56809e23`; the only intervening change is the verified scope clarification.
- `verify_interleaved_repetition.py`, SHA-256 `037e0981c9cb70cef4cef3727d6ee99e2e980101d38d2393e1497eff8e1ce829`.
- Companion `POWER_EXPANSION_ZERO_HEAD.md` and its checker, as supporting corollary/provenance.
- The supplied README, manifest, and verification receipts. Corrected README SHA-256: `f2aed4dacbd6172d0f195c44b3b22a253733fe0d05df0a7d47cd7c011dd312f5`.
- The cited theorem was opened directly at the publisher on the review date.

## 1. Guard and fuel identity

Let a chronological parity word have numerator `a=3^h`, denominator `M=2^n`, and carry `b`. Its repeated affine numerator is

`a^k x + b(a^k-M^k)/(a-M)`.

Here `a-M` is odd, including for contracting or all-even words. Multiplying by `a-M` shows that integrality modulo `M^k` is equivalent to

`a^k((a-M)x+b) = 0 mod M^k`,

and `a^k` is invertible modulo `M^k`. Thus the claimed divisibility criterion is exact. Full-word integrality authenticates every prefix parity guard: reduce the full numerator modulo each prefix power of two and cancel the remaining odd power of three. This is not merely an endpoint-only congruence.

For an expanding `U` and positive `y`, the fuel `(a-M)y+b` is a positive integer. Its divisibility therefore gives a genuine real lower bound at least `M^k`. With `b/(a-M)<=C`,

`y+C >= M^k/(a-M) > M^k/a`.

The draft's strict inequality (15) is correct. The argument would fail for arbitrary contracting factors because positivity is essential; those factors are excluded.

The center bound for an arbitrary expansion word follows by induction. Appending a map of slope `s>1` and intercept at most `C(s-1)` sends an accumulated carry at most `C(S-1)` to at most `C(sS-1)`. It does not require equal centers or commuting maps.

## 2. Source floor

Independent reconstruction gives the exact K residue `815464443` and K^4 residue

`588943162011782273059357097228697595`.

The six short initial runs `K^i F`, `K^i H`, `1<=i<=3`, all contract and strictly descend at every source at least the K residue. Since contraction makes the descent inequality stronger as the source increases, checking the stated exact inequality at this lower bound proves the whole source range, not just one residue.

An expansion-first word either begins with at least four K blocks, encounters V before four initial K blocks, or has one of those six short first reset boundaries. V itself begins with K^4. Consequently a hypothetical non-descending first-crossing endpoint must start in the K^4 cylinder. Arbitrarily many later episodes do not change this deduction.

The `F`-only edge case explains the scope qualification in the verdict. It cannot be admitted to this source-floor argument.

## 3. Normalized carry and arbitrary-context transport

For a composition of macroblocks, division of the carry recurrence by the cumulative slope gives exactly

`e_P/S_P = sum_j (e_j/s_j)/S_(j-1)`.

At every macroblock input preceding the first explicit coefficient crossing, `S_(j-1)>=1`. The verified per-letter tariff `b_j/a_j<=n_j/6` therefore gives the prefix debt bound `e_P/S_P<=|P|/6`. A word containing an expansion letter makes the full-word tariff strict, giving `e/S<N/6` and hence `e<N/6` at its crossing.

For the actual source `y=S_P x+e_P` at an interior repetition, `S_P>=1` also permits the added constant C to be transported:

`y+C <= S_P(x+|P|/6+C) <= S_P Z`, where `Z=x+N/6+C`.

Together with the fuel lower bound, this gives `Z>2^r/(S_P s_U)`. Earlier resets are not erased: both their multiplier and their accumulated carry remain in this inequality. This is the key valid extension beyond an initial repeated factor.

K has slope `lambda=3^19/2^30`, with `1<lambda<13/12`; V has slope `lambda^4(3/4)` and length 122. Hence its per-bit rate is at most K's. All omitted expansion multipliers exceed 1 and reset multipliers are below 1, so

`S_P s_U <= lambda^(L_exp/30)`.

With `r>=L_exp/100`, `log 2>1/2`, and `log lambda<1/12`, the resulting strict lower bound is precisely `log Z>L_exp/450`. The position-sensitive version using expansion length in `P U` is valid for the same reason.

## 4. Short-range logarithm certificate

The rational approximation needs only the displayed strict error estimate; there is no unproved continued-fraction assertion. Coprimality and `q>2H0` ensure that `Np-hq` cannot vanish for `1<=N<=H0`, and

`|N log(2)/log(3)-h| > 1/(2q)`.

A hypothetical non-descending source at least r4 has slope defect `u<N/(6r4)<=H0/(6r4)<1/2`. Therefore its positive logarithmic gap is `<H0/(3r4)`, whereas the preceding separation and `log 3>1` give a gap `>1/(2q)`. The exact inequality `2H0q<3r4` supplies the contradiction.

The independent certificate uses a different decomposition: `log 3=log 2+log(3/2)` with 64 rational atanh terms at `z=1/3` and `z=1/5`. It confirms the same `p/q` bracket and exact margins without floating point. The draft's 160-term certificate is also valid.

## 5. Cited all-length input and height ceiling

The publisher's [Theorem 2.1 in Languasco, Luca, Moree and Togbé](https://link.springer.com/article/10.1007/s12188-025-00293-9#Sec2) states the positive-rational product-minus-one bound used in the draft. Setting its two bases to 3 and 2, exponents to `h` and `-N`, and exponent bound to N gives the required estimate for `1-3^h/2^N`. The exponents are nonzero, `h<=N`, and prime factorization excludes equality. Its coefficient is less than `2*30^5*32*2=3110400000<10^10`.

This verifies the specialization, not a new proof of Matveev's theorem. With that external input, `(1-S)x<=e<N/6` gives the stated polynomial-in-N height ceiling. Since `x>=r4>C` and `exp(B)N^B>1`, the bound on `Z=x+N/6+C` and equation (13) follow with ample slack. No hidden assertion about existence of a first crossing is used.

## 6. Episode clock and overlap

Just before the final reset, the cumulative multiplier is at least one. If there are t resets and L expansion bits, the preceding `t-1` reset multipliers are each at most `3/4`, while all expansion multipliers together are at most `(3/2)^L`. Thus

`(4/3)^(t-1) <= (3/2)^L`.

The exact inequality `(3/2)^2<(4/3)^3` implies `t-1<3L/2`. Since each reset has at most 9 bits and `L>=30`, this gives `N<29L/2+9<=15L`. This counts all resets, including consecutive resets; it imposes no episode bound. Internal bits of V remain charged to expansion length consistently.

At `L0=2^48`, `15L0<2^52`, and the strict margin between `L0/450` and `53B+52` is positive. Also `L0>450(B+1)`, making the difference between the linear lower bound and logarithmic upper bound strictly increasing throughout the long range. The two ranges overlap and cover all positive permitted lengths. No asymptotic-only gap remains near the join.

## 7. Nonvacuity and supporting corollary

For the family `K^ell(FK)^(j-1)F`, the multiplier at the j-th reset is `lambda^(ell+j-1)(3/4)^j`. The least crossing index exists because `3lambda/4<1`; preceding reset endpoints are strictly above one by minimality and prime-factor uniqueness. The inequality `lambda^2(3/4)<1` gives `j<=ell`, while any fixed j eventually fails as ell grows. Hence the number of expansion episodes is genuinely unbounded. The initial repeated factor satisfies the density hypothesis for every `ell>=2`.

The companion single-episode theorem also follows directly from the audited main theorem: truncate `U^kR` at its first reset crossing. All expansion bits belong to `U^k`, and `(k-1)|U|>=k|U|/2`, so the repetition test holds. Each later lawful F/H reset is nonincreasing on positive integers, preserving strict descent. Single K or V cases follow from their first-reset contraction and the expansion-first short-range result. Thus the companion's old short-length certificate is not an extra dependency of the primary result or this corollary.

The lower bound on completion length near K-infinity is a necessary condition on actual non-descending endpoint witnesses. It does not exclude the formal 2-adic point from the full realized boundary. The manuscript correctly retains that distinction and leaves arbitrary bank entry open.

## 8. Reproduction and independent controls

The submitted main checker was run both normally and under `python3 -O`. The companion checker was also run in both modes and matches its supplied receipt. Both outputs agree byte-for-byte with the supplied receipt: 255,312 checks, 523 first-crossing leaves, 532 retained depth-eight frontier nodes, and 228,395 literal Terras steps. Its finite forest reaches only two episodes; the separate 1,557-episode family is therefore important and correctly labeled as a different finite control.

`independent_interleaved_controls.py` imports no submitted checker code. It constructs affine carries from the closed sum over odd-bit positions and constructs lawful residues by parity lifting rather than modular inversion. It includes:

- exhaustive fuel equivalence across all binary words of length 1–4, powers 2–4, and all source residues modulo the full repeated length;
- independently reconstructed tariffs, expansion prefixes, source floor, and all six short initial reset cases;
- the alternative exact logarithm certificate and all overlap constants;
- a complete expansion-first macro-prefix forest through depth nine, with actual replay of the least positive source and two additional lifts;
- all expansion-only repeated factors in that forest, including contexts containing earlier resets;
- large interior mixed K/V repetitions after 10 and 100 explicit resets;
- the 1,557-episode family with full chronological replay, strengthening the submitted receipt's coefficient-only large control;
- reset-only equality, wrong-guard, contracting-square equality, and earlier-crossing overshoot probes.

The independent run made **1,225,098 checks**, including **1,156,984 exhaustive fuel/guard residue tests**. Its depth-nine forest contains **1,340 first-crossing leaves**, **1,311 retained frontier nodes**, and **7,492 repeated factors**, **88** of them after earlier explicit resets. The separate interior examples have **11** and **101** expansion episodes. The 172,674-step family replay has **1,557** episodes and agrees with the supplied source hash. The resulting receipt records the exact finite universes and outcomes. The original independent run, the shift-optimized rerun, and the `python3 -O` rerun agree byte-for-byte. The proof audit, rather than the finite controls, supplies the all-length and unbounded-episode quantifiers.

## Resolved wording adjustment

The final primary definition states that a first-crossing bank word starts with a nonempty expansion run over `{K,V}` and explains the excluded `F(1)=1` case. The Section 3 heading and README carry the same restriction. Equation (7) already explains its strictness by the presence of an expansion letter. This resolves the original review finding without changing the principal argument, constants, checker, or receipts. The independent control retains the reset-only equality probe to prevent a future scope regression.
