# One of the two E14 obligations closes; the other loses all decorations

## Dependency and scope

This deduction uses the terminal fee certified in `precision_lifts_ramsey_deficits_20261030/COLLATZ.md` C5 at source commit `0578f16726e4e535189fb6074cc2ddb024c95d82`. The inherited E10 certificate was not rerun here. Its consequence for every non-one primitive11 excursion c is

    lambda12 A_c - S11(c) >= gamma,
    gamma = eta-1 = 134277301/48828125,
    eta = 183105426/48828125.

Here lambda_e = 5/8 - 1/(8*5^(e-1)); clocks are the original integer sums of positive valuation letters. A primitive e excursion first returns to the identity affine map modulo 5^e at its endpoint. No abstraction, replacement entry phase, or resampled multiplier is introduced.

## 1. Small-clock exclusion and a new quantized native fee

Set D=8*5^11=390625000 and H=5A/8-S11=j/8. The inherited inequality gives

    H >= gamma + A/D.

Since gamma > 21/8, j>=22. If j=22, then

    A <= D*(22/8-gamma) = 342.

A non-one positive word of length L has 1<=L<A. Its multiplier is 3^L/2^A. A necessary condition for a primitive11 endpoint is 3^L=2^A modulo 5^11. The accompanying finite exact checker excludes this equality for every 1<=L<A<=342 (58,311 pairs). This is a full certificate for the bounded exceptional case, not a sampled all-word inference. Consequently every non-one primitive11 excursion has H>=23/8.

Now lambda14=5/8-1/(25D). Since A/D<=H-gamma,

    lambda14*A-S11
      = H-A/(25D)
      >= (24/25)H+gamma/25
      >= 2+r,
    r = 1062011676/1220703125 = 0.8699999649792.

The threshold 342 and the finite multiplier check eliminate the entire j=22 branch. This strengthens the continuous terminal-fee estimate without claiming an improvement of its abstract certificate.

## 2. Uniform primitive13 bound

Cut a non-one primitive13 excursion w at its actual12 returns. Let a12 count all-one constituents and B count non-one constituents. Cut those B constituents at their actual11 returns. Let a11 count all-one constituents at this inner level and C count non-one primitive11 cores. Necessarily C>=B>=1.

At slope lambda14 an all-one primitive11 block costs31/125 and an all-one primitive12 block costs6/25. Each non-one core costs at least2+r by section1. Every non-one12 constituent adds exactly one level12 charge; w adds exactly one level13 charge. Thus, exactly preserving all original clocks and lower charges,

    delta13(w) >= (2+r)C-B-1+(31/125)a11+(6/25)a12
               >= r+(31/125)a11+(6/25)a12.

This applies to every non-one primitive13 excursion, not merely unpaid candidates. In particular

    g13(-(4,4)) >= r > 4/5.

**The second of the source packet's two E14 inequalities is proved.** Empty tangent classes satisfy the statement automatically.

For delta13(w)<1, C>=2 would give at least2(1+r)-1>1. Hence C=B=1. Moreover r+6/25>1 and r+31/125>1, so a11=a12=0. Therefore every possible E14 failure is a bare non-one primitive11 core which also first returns at levels12 and13 at its endpoint, and whose map is the identity modulo5^14. All four decorated forms of the inherited five-form grammar are excluded.

## 3. Precisely what remains open

For such a bare core c, S13(c)=S11(c)+2, so

    delta13(c)=j/8-2-A/(25D).

We have j>=23 and A<=D*(j/8-gamma). If j>=25 the quantized bound in section1 gives delta13 >= (24/25)*(25/8)+gamma/25-2>1. Thus any strict E14 counterexample must satisfy exactly one of

    j=23, 1<=A<=48,828,467, A=3 mod8;
    j=24, 1<=A<=97,656,592, A=0 mod8.

Both caps are exact integers; the residue follows from j=5A-8S11. All further source constraints remain: native primitive11 status, actual translation, positive-word realizability, and full identity endpoint modulo5^14. Neither scalar branch is proved realizable or impossible here. Thus g13(0)>=1, the full E14 theorem, the all-layer estimate, and Collatz convergence remain open.

## Verification

`check_quantized_native14.py` verifies all fractions, small-clock exclusion, return charge coefficients via exact all-one floor sums, and strict separating inequalities. It also checks the j>=25 exclusion and both residual clock bands. Run normally and with `python3 -O`; identical outputs certify that no mathematical guard relies on `assert`.
