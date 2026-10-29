# Reflection persistence, rational contact clocks, and a golden height barrier

Repository source pin: `90e704a8f35614d2085778d1b518f03ff0d750ca`. Named packets below are under `05-knowledge/results/` unless a different path is given.

**Status:** PROVED elementary refinements; the contact-clock lemma is CITED classical and reproved here. Finite checks are explicit. No new Lonely Runner lower bound, runner case, or major-conjecture closure is claimed.

## 1. Inheritance: the missing information is not another odd trace

For a nonempty finite set V of positive integer velocities, define

    f_V(t)=min_(v in V) ||v t||,      G_delta={t in R/Z:f_V(t)>=delta}.

The latest `odd_sign_spectra_trijections_20261028/LRC_UNIT_DISTANCE.md`, L1–L2, proves that reflection t↦1−t fixes one component exactly when all speeds are odd. It also proves that for primitive V no new reflection center preserves every labelled phase-distance function. Its mixed-parity control V={1,2} has empty, isolated, and positive-length good sets with the same zero reflection trace. Those results are inherited, not rediscovered here.

The current `minimal_blue_selectors_20261028/THEORY.md` identifies an exact smallest carrier for a specific tournament singleton/complement/reset language. The analogy worth retaining is **operation-relative sufficiency**: a statistic may be exact for one query and blind for another. It does not identify LRC times with tournament frames.

## 2. The entire mixed-parity persistence module doubles

Suppose at least one v in V is even. Work on the explicitly chosen unit-period circle; if minimal-period claims are desired, first divide V by its gcd. For every delta>0, neither0 nor1/2 lies in G_delta. Put

    G_delta^+ = G_delta ∩ (0,1/2).

Then there is a reflection-equivariant homeomorphism

    G_delta = G_delta^+ disjoint union (1−G_delta^+).

This identification commutes with all inclusions G_a⊆G_b for a>=b>0. Hence over Q the full superlevel H_0 persistence module is

    H_0(G_delta;Q) ≅ Q[C2] tensor H_0(G_delta^+;Q),

with every inclusion map equal to the identity on the regular C2 factor tensored with the corresponding half-interval map.

**Proof.** Reflection interchanges the two open semicircles and preserves f_V. The two reflection-fixed points0 and1/2 are excluded by the positive threshold and the even speed. Thus the partition is disjoint at every threshold and compatible with inclusion. Each good set is a finite union of intervals and isolated points, so H_0 is exactly the free vector space on components. The asserted module decomposition follows directly. ∎

In particular, for every a>=b>0, reflection has trace zero not just on H_0(G_a), but also on the image of

    H_0(G_a) → H_0(G_b).

So adding all cross-threshold signed traces cannot recover whether a mixed-parity good set is empty. The ordinary dimensions and birth heights remain informative; the signed trace cancels an entire regular factor. This is an operator-level obstruction to this particular parity strategy, not to LRC itself or to every possible weighted incidence construction.

When all velocities are odd, the same image trace is1 for 0<b<=a<=1/2: the unique reflection-fixed lower component contains1/2 and hence meets every upper good set. Paired components contribute zero. The all-odd witness is already explicit at1/2.

## 3. Exact maxima have short rational contact certificates

Write M=max V and lambda(V)=max_t f_V(t). A maximizer can be chosen with

    t = k/(u+v),       u,v in V, allowing u=v.

Consequently its reduced clock denominator and the reduced denominator of lambda are both at most2M.

This is classical: see Perarnau–Serra, [*The Lonely Runner Conjecture turns60*, §2, equation(4)](https://arxiv.org/html/2409.20160v3#S2), which records the pair-sum clock observation and its longstanding use. The repeated-speed convention here covers a lone active cusp and singleton V explicitly.

**Elementary proof.** Each distance ||vt|| is a continuous sawtooth with slopes ±v, valleys at k/v and peaks at (2k+1)/(2v). The minimum is positive somewhere, since finitely many zeros do not cover the circle. At a positive maximum, either an active sawtooth is at its peak, giving a clock denominator2v, or two active linear branches have opposite slopes. If all active slopes had one sign, moving a little in that sign would increase the minimum. Opposite branches have equations ut−a and b−vt, for integers a,b, so their meeting time is (a+b)/(u+v). At that clock every distance, and therefore their minimum, is an integer divided by u+v. ∎

Thus any true decision lambda(V)>=delta has a rational witness with O(log M) bits in its numerator and denominator. For a rational threshold delta, checking it against every velocity uses exact modular arithmetic and rational comparisons. This does not give a uniform denominator bound independent of the input, nor a fast way to find the witness. Enumerating the pair-sum clocks is pseudo-polynomial in M, not polynomial in its binary length.

This explains the compatibility with the old denominator27 hostile against a universal bound25: short binary certificates can still have unbounded denominator values.

## 4. A quantitative barrier to a literal golden optimum

Let phi=(1+sqrt5)/2 and alpha=2−phi=(3−sqrt5)/2. For every finite integer-velocity instance,

    |lambda(V)−alpha| > 1/(12 M²).

**Proof.** Write lambda=a/q in lowest terms. Section3 gives q<=2M, and 0<lambda<=1/2. The integer a²−3aq+q² is nonzero because alpha is irrational. Therefore

    |lambda²−3lambda+1| >= 1/q².

The other root is beta=(3+sqrt5)/2<3. Factoring the polynomial and using |lambda−beta|<3 gives

    |lambda−alpha| > 1/(3q²) >= 1/(12M²).

This is an elementary quadratic-norm bound, not a new Diophantine approximation theorem. ∎

Golden quantities in packing or Zeckendorf/Fibonacci limits can still be connected through normalized limits, algebraic coordinate changes, or different predicates. What is ruled out is the literal equality of alpha with the optimum of a finite integer-speed instance. More generally, for an algebraic irrational alpha of degree d and a fixed compact interval of rational candidates, its integer minimal polynomial gives a constant times q^(−d) exclusion by the same argument.

The finite contact certificate and irrational limiting geometry therefore live at different scales. The useful data for a transfer are the height M, the chosen normalization, and a proof relating the predicates. Matching the word “golden” supplies none of these automatically.

## 5. Relation to geometry and odd tournament Fourier carriers

The companion geometric Euler repair preserves an entire unit-distance relation while allowing an unconstrained containing region. The present rational certificate preserves existence of a good time but forgets the shape of the good-time set. These are distinct information contracts.

Session90's decorated Walsh correspondence can preserve a full invariant-function space after inverting2; its signed traces must retain the entire operator data to be useful. Here taking only reflection traces intentionally collapses a regular representation factor. One cannot import tournament oddness to recover the lost positive multiplicity without a new map or observation.

This suggests a concrete next question rather than a claimed solution: use rational contact heights and labelled phase owners in a weighted incidence object that distinguishes empty, point, and interval witnesses. Any such proposal must survive V={1,2} at thresholds2/5,1/3,1/4 and the AP13 threshold1/14 control before it can support an LRC argument.

## 6. Exact checks

`check_lrc_persistence.py` uses rational arithmetic only. It compares pair-sum candidates against an independent full sawtooth-cell arrangement for every nonempty subset of {1,...,9} of size at most4. It checks the denominator and quadratic-gap inequalities without a floating approximation to sqrt5; compares all selected cross-threshold image traces; retains the empty/point/interval controls; and checks the six isolated AP13 good times. Run normally and with `python3 -O`; all checks remain active.
