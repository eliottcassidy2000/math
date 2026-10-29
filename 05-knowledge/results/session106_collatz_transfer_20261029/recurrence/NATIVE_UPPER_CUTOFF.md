# Native first-descent survival needs only an upper cutoff

**PROVED:** a sharper carrier for the positive native affine family. **FINITE-EXACT:** the hostile search below. The general two-bound guarded-macro theorem remains useful for other observations.

Fix odd b>=3. At any nonempty native prefix j the affine endpoint is

    y_j=(A_j n+C_j)/B_j,
    A_j=b^j, B_j=2^H_j, C_j>0.

Positivity of C_j follows from C_0=0 and C_(j+1)=b C_j+B_j with B_j>0. Consequently:

* If A_j>=B_j then y_j>n for every positive n. This prefix imposes no survival restriction.
* If A_j<B_j then y_j>=n if and only if

      n <= floor(C_j/(B_j-A_j)).

For a fixed full word u, define K(u)=infinity when none of its prefixes has contracting slope, and otherwise

    K(u)=min_(j:A_j<B_j) floor(C_j/(B_j-A_j)).

Then **exactly** the positive integer sources in the word's exact cylinder with n<=K(u) have avoided first descent at every prefix. The source1 may be excluded separately as the genuine terminal. If the full cylinder is n=r+Mq and the allowed source lower bound is n>=n_min, its survivor parameters are

    max(0,ceil((n_min-r)/M)) <= q <= floor((K(u)-r)/M),

with no upper bound when K=infinity. Thus native survival creates no new lower-bound debt: surviving sources form an initial segment of that cylinder relative to the original source lower bound. With unbounded initial sources, the survivor cell is unbounded **if and only if every prefix slope is at least1**; otherwise it is finite, possibly empty. This equivalence follows because every exact cylinder has infinitely many positive sources and finite K cuts it to a finite set.

The recurrence can retain an absolute source cap K. Each extension updates its inherited A,B,C and replaces K by min(K,floor(C/(B-A))) when B>A, leaving K unchanged otherwise. Exact-cylinder lifting still retains r and the endpoint intercept; the scalar cap does not replace those arithmetic coordinates. It only compresses the historical survival inequalities. Initial finite source bounds can be folded into K from the outset.

## What the arithmetic guard can and cannot be presumed to do

For b=5, word(1,1,5) has K=13 and exact residue13, so its nonterminal survivor cell is the singleton{13}. Word(1,3,3) similarly keeps17. These are genuine finite nonempty cells, not slope-only certificates. The former is the native cycle13->33->83->13.

For b=3, word(2) has the finite survivor{1} if the terminal is not excluded. We searched for a **nonterminal** finite survivor, rather than inferring one exists from a positive real threshold. Exhaustive words of lengths1..8 on letters1..4 (87,380 words) yield none. A separate literal native scan of50,000 odd sources3..100001, examining at most50 steps until their first slope contraction, likewise finds no nonterminal unpaid contracting prefix among the crossings it observes. The receipt records crossings and sources whose crossing exceeds the work horizon.

These are finite negative findings, **not** a proof that every contracting3x+1 word pays every nonterminal source in its exact cylinder. No such general guard implication is used by the theorem or algorithm, and the search does not settle arbitrary word lengths, larger valuations, or larger sources. A positive real cutoff by itself also does not guarantee an integer of the exact cylinder lies below it.

`check_upper_cutoff.py` independently compares the minimum-prefix cutoff with literal prefix trajectories on cylinder representatives, confirms terminal1 and the two5x+1 finite cells, and records both bounded searches. Explicit checks persist under `python -O`.

## Adjacent open conjecture: clocks must not be identified

Rozier and Terracol's [*Paradoxical behavior in Collatz sequences*, §1, Definition1.2 and C1–C3](https://arxiv.org/html/2502.00948v5) describes Terras's open CST conjecture: first descent equals first coefficient below1, for the micro-map odd n→(3n+1)/2, even n→n/2. Their example7→11→17→26→13→20→10→5→8 has final coefficient243/256<1 but endpoint8>7; prior descent to5 prevents a CST counterexample.

**Our deduction:** CST implies no nonterminal finite survivor cells here: micro-step descent during a division run persists to its odd endpoint. No converse/equivalence is claimed; coefficient crossing can occur inside a division run. Our cutoff theorem assumes neither CST nor slope-only payment. The verifier independently checks the micro-clock example, including its carry347/256; it is not an accelerated odd-to-odd counterexample.
