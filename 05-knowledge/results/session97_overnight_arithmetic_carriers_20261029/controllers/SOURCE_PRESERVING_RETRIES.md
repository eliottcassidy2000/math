# Source-preserving retries: repaired rays and an unbounded-clock certificate

**Status: PROVED** for the exact subsets stated; FINITE-EXACT implementations. Global termination of the retry process remains OPEN. Dependencies: the universally terminating cheaper controller and its five-cylinder comparison classifier in this directory.

## 1. A marked controller, not a resampled chain

Keep state `(n0,n,J,D,C)` with

    2^D n = 3^J n0+C.

Here n0 is the fixed original source; n is the current port. Success means `n<n0`, equivalently `(2^D−3^J)n0>C`. A locally descending next transition is insufficient: `27→41→31` has a strict local drop41→31 but still lies above27. Every cumulative update retains exact clocks, source guards and carry. The input and output ports are the actual integers, not independently averaged labels.

## 2. Four repaired nontrivial rays

Let `(r,M)` be one of the four nontrivial expanding first-return cylinders. Write n=r+Mk. Their outputs y have the following exact shifted forms:

|r|M|y+1|
|---:|---:|---|
|125|128|2·81·(k+1)|
|1357|4096|2·729·(3k+1)|
|2269|4096|2·243·(9k+5)|
|2877|4096|2·81·(27k+19)|

If k is even, y≡1 mod4. Its next native valuation a is at least2. For the first-return chart `(A n+C)/2^D`, the combined endpoint is

    z=(3A n+3C+2^D)/2^(D+a).

The worst case a=2 already gives strict descent below **n**, uniformly on the whole source ray:

| Repaired source ray | Positive slope gap 2^(D+2)−3A | Carry 3C+2^D | Odd steps |
|---|---:|---:|---:|
|n≡125 mod256|13|601|5|
|n≡1357 mod8192|1631|50579|8|
|n≡2269 mod8192|1631|38915|8|
|n≡2877 mod8192|1631|31139|8|

In each row, gap·r>carry, and the inequality only strengthens as n increases. This repairs exactly half of each nontrivial expanding cylinder with one extra native step. The division clock is **D+a**, not D+2; a≥2 is an inequality guard, not permission to falsify the actual clock. If the next step initiates a new auxiliary branch, this controller intentionally stops at that intermediate native endpoint; it is a newly specified return rule, not an assertion that a whole second auxiliary tree has already completed.

For the simple expanding cylinder, n≡3 mod16 has first valuation1 and next valuation at least3. The endpoint is at most `(9n+5)/16<n`, giving another repaired ray in two odd steps.

Together these disjoint repaired rays have normalized odd-source density547/4096. This is a finite arithmetic coverage theorem, not an independence claim or a bound on repeated failure probability.

## 3. Compress the entire valuation1 run without losing its clock

For any positive odd n write

    t=v_2(n+1)≥1,     u=(n+1)/2^t odd.

There are exactly t−1 initial valuation1 steps. After these,

    n_run=2·3^(t−1)u−1.

Put c=v_2(3^t u−1)≥1. The next valuation is c+1, and the endpoint is

    z=(3^t u−1)/2^c
     =[3^t n+(3^t−2^t)]/2^(t+c).                  (1)

The compressed transition has odd clock J=t, division clock D=t+c, carry C=3^t−2^t. Its exact source condition on u is

    u≡3^(−t)(1+2^c) mod2^(c+1),                  (2)

together with u odd and n=2^t u−1. Equation(2) retains the exact valuation, not just divisibility. For a preceding chart `(A n0+C0)/2^D0`, composition gives

    A_total=3^t A,
    C_total=3^t C0+2^D0(3^t−2^t),
    D_total=D0+t+c,   J_total=J0+t.

The shifted formulas in the table give t directly from the same original source parameter k: t=1+v_2(k+1), or1+v_2(3k+1), or1+v_2(9k+5), or1+v_2(27k+19). This is a concrete source-preserving arithmetic connection between the four failure ports. The powers81,243,729 are affine derivative factors, not a claimed identification with unrelated tournament symmetry groups.

## 4. A complete hostile family: only t=2,4,8 recover immediately

For n=2^t−1, t≥2, the compressed transition(1) descends below its original input **if and only if**

    t∈{2,4,8}.

The three returns are3→1,15→5,255→205. For every other t≥2 the compressed endpoint is strictly larger. This classifies the first compressed return, not the later fate of those integers.

Proof: u=1 and elementary factorization gives c=v_2(3^t−1)=1 for odd t, and c=2+v_2(t) for even t. For odd t≥3,3^t>2^(t+1), so the slope is greater than1 and the positive carry forces expansion. For even t≥10, `2^c≤4t`, while `(3/2)^10>40` and `(3/2)^t/t` is increasing for t≥10. Thus3^t>2^(t+c), again forcing expansion. The cases t=2,4,6,8 are checked exactly; only6 expands. No empirical extrapolation is used.

This family disproves a fixed upper bound on the number of valuation1 first-controller returns before any local drop: the initial expansion run has length t−1. It also disproves the idea that simply compressing that run guarantees descent below the original source. The unbounded t and exact terminal c are necessary sidecars for this representation.

## 5. Remaining problem and reproducibility

The arithmetic exceptional set for **auxiliary termination** is empty. The first-controller expanding set is exactly five cylinders. Some retry rays are now repaired, and all remaining maximal valuation1 runs have an exact compressed formula. The unresolved global question is whether source-preserving compositions of the remaining rays must eventually cross below their original input; no such claim is made here.

`python3 -B check_source_preserving_retry.py` and its `-O` replay verify the four uniform inequalities,1000 parameter values per repaired ray, the3 mod16 ray through200000, all50001 odd compressed inputs through100001, exact guards/carries/clocks, and the power-family hostile for2≤t≤300. `retry_checks.json` records the outputs. The infinite statements are proved above.
