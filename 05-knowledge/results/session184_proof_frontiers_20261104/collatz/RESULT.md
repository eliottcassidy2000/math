# Exact first descent across two periodic templates

**Status: proved elementary family and obstruction; Collatz remains open.** This extends S179's periodic-shadow parser by actually certifying a template change on infinite arithmetic domains. It does not provide a universal template selector or global rank. No historical novelty is claimed.

## 1. Every slope-first-crossing word is realized by a first-descent progression

For a word `a=(a_0,...,a_(m-1))` of positive valuations, write

    S_j = a_0+...+a_(j-1), P_j=3^j, Q_j=2^S_j,
    B_0=0, B_(j+1)=3 B_j+2^S_j.

Suppose `P_j >= Q_j` for `1<=j<m`, but `P_m<Q_m`.
Set `P=P_m`, `Q=Q_m`, `B=B_m`, and define the odd residue

    R = (Q-B) P^(-1) mod 2Q.

Then **every positive integer n in R modulo 2Q satisfying n>B/(Q-P) has first strict descent at exactly odd step m**, and its descendant is `(P n+B)/Q`.

Proof: the residue equation is precisely the cylinder selecting the entire valuation word and an odd final output. To check intermediate exactness directly, reduce `P n+B` modulo `2^(S_j+1)`. For j<m, the first carry term after prefix j has valuation exactly S_j, and all later terms have higher valuation. Thus `3^j n+B_j` has valuation exactly S_j. The final congruence similarly gives exact valuation S_m. These statements give the claimed odd intermediate states and their native valuation letters. Positivity follows from the positive affine formulas. For j<m, the state `(P_j n+B_j)/Q_j` exceeds n because `P_j>=Q_j` and `B_j>0`. The last inequality is exactly `(Q-P)n>B`. QED.

The arithmetic progression has natural density `1/(2Q)` among all positive integers, or `1/Q` among positive odd integers. The explicit threshold removes finitely many sources and does not change either density. These are first descents, not necessarily descents directly to 1; coverage of all positive sources is not asserted.

## 2. An unbounded, bounded-valuation template switch

Take the expanding template `w=(1,2)`, followed by the contracting template `(2)`. For every integer q>=1 let

    r(q)=min{r>=0 : 3^(2q+r)<2^(3q+2r)}
        = floor(alpha q)+1,
    alpha = log(9/8)/log(4/3).

The word

    W_q = (1,2)^q 2^r(q)

satisfies the theorem's hypotheses. During the expanding part, slopes at even indices are `(9/8)^h>1` and at odd indices are `(3/2)(9/8)^h>1`. During the final two-run, slopes are `(9/8)^q(3/4)^j`; minimality of r gives the first slope crossing at the final step. Therefore an explicit positive-density progression has **first descent exactly at 2q+r(q)**, with every valuation at most two.

The carry has the useful closed form

    B = 3^r * 5(9^q-8^q) + 8^q(4^r-3^r)
      = 5*9^q*3^r + 8^q*4^r - 6*8^q*3^r,
    P=9^q*3^r, Q=8^q*4^r.

This follows from `F(n)=(9n+5)/8`, so `F^q(n)=(9/8)^q(n+5)-5`, followed by `G^r(x)=1+(3/4)^r(x-1)`.

The S179 parser predicts exactly the intended first mismatch: its `(1,2)` center is -5 and the sources here satisfy

    v2(n+5)=3q+1.

Hence precisely 2q template letters match, and the next letter is 2 rather than 1. This is a genuine switch to a different center, from -5 to 1, followed by a proved size descent. It is not a claim that the separate S179 energies decrease across the switch.

Examples: q=1, r=1 gives n=43 and the first descendant below n is 37 at step 3. For q=5, r=3, n=720891 first descends at step 13 to 548047. For q=1000 the compiler gives an exact source with 3,821 binary digits and first descent at step 2410. The machine-readable output records the size exactly; no decimal expansion is needed.

## 3. Sharp resource tradeoff and an obstruction to fixed-margin coverage

If q repetitions of any expanding word of length L and valuation cost A are corrected by a *single* final reset b, a necessary condition for descent, for every positive source, is

    2^(qA+b)>3^(qL+1).

Indeed, otherwise the full affine slope is at least one and its carry is positive. The strict inequality is also sufficient on a tail of its native cylinder by Section 1 whenever all prior prefix slopes are at least one. Therefore a single-reset strategy has a correction valuation growing at least linearly in q. Section 2 replaces that large valuation by an unbounded number of valuation-two steps, whose exact minimum is r(q).

The final slope of W_q is

    rho_q=(9/8)^q(3/4)^r(q),  3/4<rho_q<1.

The number alpha is irrational: a rational alpha=p/s would imply `3^(2s+p)=2^(3s+2p)`, impossible by unique factorization. The fractional parts of q alpha are dense modulo one, so rho_q approaches 1 along a subsequence. Thus actual first-descent words, realized on infinite integer progressions, have arbitrarily small multiplicative contraction margin.

Consequently no finite union of certificate families each restricted by a fixed strict slope bound `rho<=c<1` can capture every W_q. This is a direct obstruction to uniform-margin coverage of all first-descent words, not an obstruction to proving eventual descent by other decompositions.

There is also an exact formal-language statement: the language `{(12)^q 2^r(q): q>=1}` is not context-free. If it were, inverse image under `a -> 12, b -> 2`, intersected with `a^+b^+`, would give a context-free language with Parikh image `{(q,r(q))}`. By [Parikh’s semilinearity theorem](https://doi.org/10.1145/321356.321364) (Rohit J. Parikh, *On Context-Free Languages*, Journal of the ACM 13(4), 570–581, 1966), its image would be a finite union of linear sets. Any nonzero period `(u,v)` in a linear set contained in this graph must satisfy `v=alpha u`, because the discrepancy `r(q)-alpha q` remains in (0,1). Irrationality prohibits every nonzero integer period. Such a semilinear set would be finite, contradiction. This directly realizes the near-critical obstruction behind S179's conservative-grammar discussion rather than assuming the sharp words never occur on positive trajectories.

## 4. General expanding templates

Every word w with `3^L>2^A` has a cyclic rotation for which all prefix slopes are at least one. One proof sets increments `log(3)-a_i log(2)`, chooses a minimum of the cumulative sums before the final endpoint, and rotates there. The total sum is positive, so both the unwrapped and wrapped prefix sums are nonnegative.

For such a rotation, w^q followed by r valuation-two steps satisfies Section 1 when r is the first integer with `3^(qL+r)<2^(qA+2r)`. Thus the construction is not special to (1,2). It handles arbitrary expanding periodic shadows on their exact native congruence domains.

## 5. Tests and remaining obligation

`verify.py` independently replays the native map on 880 sources for all expanding templates of lengths 1 through 5 and letters 1 through 3, rotated as above, at repetition counts 1,2,5,20,100. It checks 135,132 native steps across two-run corrections and single-reset corrections. Long controls additionally test q=1,2,10,100,1000 for the (1,2) family and check the S179 precision identity. Explicit exception-raising checks survive optimized Python; normal and `python -O` outputs are byte-identical.

Sources: [S179 periodic-template rank atlas](../../session179_niche_methods_20261104/collatz/RESULT.md), [S179 conservative grammar analysis](../../session179_niche_methods_20261104/algorithms/RESEARCH.md), and [the older unbounded run-block construction](../../collatz_run_block_family_20261005.md). The older run-block construction is the w=(1) ancestor, not a universal coverage theorem. No forbidden canon file was read or scanned.

**Open:** prove that every positive odd source above 1 has a certified actual descendant strictly below that source (or provide another well-founded global selector). Merely entering a later descent domain is insufficient when its entry value can exceed the original source. These arithmetic families establish genuine unbounded descent certificates and an explicit obstruction to a fixed-margin strategy; they do not close that universal quantifier or exclude nontrivial cycles.
