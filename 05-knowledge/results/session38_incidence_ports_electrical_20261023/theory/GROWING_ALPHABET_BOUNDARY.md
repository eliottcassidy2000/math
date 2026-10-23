# A linear alphabet is necessary for row-multiset class separation

Status: **PROVED elementary counting obstructions.** These concern one tournament's ordinary class signature, even restricted to strong tournaments. They do not classify reconstruction from a full authenticated metagraph inventory, and do not contradict the finite order-eight nine- or twelve-symbol results.

## 1. Model and exact counting bounds

At order n, suppose every outgoing arc receives one of L=L(n)>=1 symbols. Each vertex retains only the multiset of its outgoing symbols, and the tournament signature retains only the multiset of its n vertex rows. No extra data are included. The symbol-assignment rule may be arbitrarily sophisticated; the bound depends only on the output format.

A row is an L-vector of nonnegative counts with sum at most n-1. Loosening this to n, stars and bars gives at most

    R(n,L) = binom(n+L,L)

possible rows. An ordered list of n rows overcounts the outer-multiset signatures. Therefore their number is at most

    U_row(n,L) = binom(n+L,L)^n.                         (1)

There is a stronger bound using the compulsory total mass. Every tournament has exactly

    N = n(n-1)/2

arcs, and every arc contributes exactly one symbol to exactly one outgoing row. An ordered n-row signature is thus an array of nL nonnegative entries with total sum N. Ignoring the individual row-length restrictions only enlarges the set, so

    U_mass(n,L) = binom(N+nL-1,nL-1)                    (2)

is also an upper bound. One may always take the smaller of (1) and (2).

For n>=3, orient one fixed spanning Hamiltonian cycle and choose the other arcs freely. This gives 2^(n(n-1)/2-n) labelled strong tournaments. Each ordinary isomorphism class has at most n! labelled members. Hence the number of strong ordinary classes is at least

    B(n) = 2^(n(n-3)/2) / n!.                          (3)

Whenever either upper bound (1) or (2) is strictly less than (3), two nonisomorphic strong tournaments share the signature. The claim does not require their Hamiltonian path counts to differ.

## 2. First asymptotic bound from the row count

Define the increasing continuous function

    f(a) = (1+a) log_2(1+a) - a log_2 a,    a>0,
    f(0) = 0.

Its derivative is log_2(1+1/a)>0, and f(a) tends to infinity. Thus the equation f(a)=1/2 has a unique positive solution

    a_row approximately 0.10482841666766435.

The standard binomial entropy bound gives

    log_2 binom(n+L,L) <= n f(L/n).

For completeness, this bound follows by evaluating the binomial expansion at probability L/(n+L): any single term is at most the whole sum, which equals one. Thus

    log_2 U_row(n,L) <= n^2 f(L/n).

On the other hand,

    log_2 B(n) = n(n-3)/2 - log_2(n!)
               = n^2/2 - O(n log n).

If these signatures separate all strong ordinary classes at every sufficiently large n, necessarily

    liminf_(n->infinity) L(n)/n >= a_row.              (4)

Indeed, if the liminf were smaller, there would be an infinite subsequence with L/n bounded above by a fixed a<a_row. On that subsequence the upper exponent is at most n^2 f(a), strictly smaller than n^2/2-O(n log n) for large n. This proves more than failure for fixed L: any sublinear alphabet L=o(n) is impossible.

## 3. Sharper asymptotic bound from the fixed arc mass

Bound (2) improves the necessary constant. Put

    t_n = N/n^2 = (1-1/n)/2,    a_n = L/n.

Since binom(N+nL-1,nL-1) <= binom(N+nL,nL), the same entropy bound yields

    log_2 U_mass(n,L)
      <= n^2 [(t_n+a_n) log_2(t_n+a_n)
              - t_n log_2 t_n - a_n log_2 a_n]
       = n^2 t_n f(a_n/t_n).                          (5)

As t_n tends to 1/2, the coefficient tends to (1/2)f(2a) whenever a_n tends to a. The same fixed-gap subsequence argument as above proves the stronger necessary condition

    liminf_(n->infinity) L(n)/n >= a_mass,             (6)

where a_mass is the unique positive solution of

    f(2 a_mass) = 1,
    a_mass approximately 0.14690768667020775.

The fixed-gap argument does not assume that L/n has a limit. If its liminf is less than a_mass, take a subsequence bounded above by a fixed a<a_mass and use monotonicity in a_n and continuity in t_n in (5).

Equations (4) and (6) are necessary, not sufficient. The stronger statement (6) implies (4). Neither establishes a complete signature at linear alphabet size.

## 4. Explicit all-order obstruction: L<=n/16 and n>=64

There is an elementary uniform finite threshold with a completely rational comparison:

    If n>=64 and 1<=L<=n/16, class separation fails
    even when restricted to strong tournaments.       (7)

From the binomial theorem at 1/16,

    binom(n+L,L) <= 17^(n+L) / 16^n.

The hypothesis 16L<=n and the exact integer inequality

    17^17 = 827240261886336764177
           < 1180591620717411303424 = 2^70

give, after raising the row bound to the sixteenth power,

    binom(n+L,L)^16
      <= 17^(16n+16L) / 2^(64n)
      <= 17^(17n) / 2^(64n)
       < 2^(6n).

Thus log_2 U_row(n,L)<3n^2/8=48n^2/128.

Meanwhile n!<=n^n, and for every integer n>=64,

    (n-3)/2 - log_2 n >= 49n/128.                    (8)

At n=64 there is equality. One can verify its persistence without numerical logarithms: (8) is equivalent to

    n^128 <= 2^(15n-192).

For the induction step, (1+1/n)^128 <= (65/64)^128 < 2^15; the last inequality is an exact integer comparison. Hence (3) is at least 2^(49n^2/128), strictly larger than the signature upper bound. This proves (7).

The looser e-based route also works, but requires monotonicity when substituting L<=n/16 into L log_2(e(n+L)/L). The exact binomial argument above avoids hiding that step.

## 5. Reproduction and limitations

Run `python check_growing_alphabet_boundary.py` and then `python -O check_growing_alphabet_boundary.py`. The checker verifies the two exact integer inequalities supporting the all-n threshold, the boundary case n=64, exact information deficits at the maximum allowed integer L for n=64 through 512, stars-and-bars controls, and exact rational brackets around both asymptotic constants. Decimal bisection supplies the displayed approximate constants; the theorems define the constants by their equations, not by numerical approximations.

Default replay compares `growing_alphabet_receipt.json`; `--write` intentionally regenerates it. The proof does not rely on the finite sweep. The checker uses explicit exceptions and no removable assertions.

This is an output-information bound for a single row-multiset signature. An authenticated source/target inventory, endpoint capacities, converse equalities, or a full collection of metagraph observations contains additional information. These bounds alone do not show failure of reconstruction from that larger input. Nor do they prove a Hamiltonian-count distinction, a minimal actual collision order, or impossibility for alphabets whose size grows linearly or faster.
