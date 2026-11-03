# An injective square-root family of Hamiltonian-path counts

Let A(X) be the number of distinct integers at most X that occur as the number of directed Hamiltonian paths in a finite tournament. For every integer r≥1,

    A(4^r) ≥ 2^(r−1),

and all these 2^(r−1) values are realized by strongly connected tournaments of order 2r+1. Consequently, A(X)>√X/4 for all real X≥4. This is a lower bound on the number of attained values, not a positive-density claim or a proof of the full conjecture. It is derived from the one-apex count identity established in the Session 149 packet; no claim of mathematical novelty or literature priority is made.

## Construction and proof

Consider a transitive chain of length 2r and an apex x. Encode by b_i=1 when x beats chain vertex i and 0 otherwise. Choose any binary word b_0...b_(r−1) with b_0=1, followed by r zeros. The first chain vertex is beaten by x, and the last beats x, so the tournament is strong.

Let P be the binary integer represented by the r-bit prefix, and define P_i to be the integer represented by the first i bits (1≤i<r). Set

    C = sum over i=1,...,r−1 of c_i,
    c_i = P_i if b_i=0,  1 if b_i=1 and b_(i−1)=0,  0 otherwise.

The Session 149 one-apex formula gives

    H = 2 + (2^r−1)P + C.                            (*)

Indeed, when the next bit is 0 its inversion contribution is precisely the integer encoded by the preceding bits; when the next bit is 1 it contributes only the 01-adjacency term. The r trailing zeros have preceding-prefix integers P, 2P, ..., 2^(r−1)P, contributing (2^r−1)P. No 01 adjacency occurs in this suffix.

Each c_i≤2^i−1. Thus 0≤C≤sum_(i=1)^(r−1)(2^i−1)=2^r−r−1<2^r−1. Formula (*) can be decoded exactly:

    P = floor((H−2)/(2^r−1)).

Hence the 2^(r−1) prefix choices have distinct counts. Also P≤2^r−1, giving H≤2+(2^r−1)^2+2^r−r−1=4^r−2^r+2−r≤4^r. If r=floor(log_4 X), then 4^r≤X<4^(r+1), so 2^(r−1)>√X/4.

## Reproduction

`python3 check.py` independently builds the tournaments and evaluates their path counts by subset/last-vertex dynamic programming. It checks the formula, injective decoder, strong-family endpoint conditions implicit in the construction, and the claimed bound for all prefixes at 1≤r≤6. The reasoning above is for all r; finite tests only check implementation. See `check.out`. The Session 149 packet at [`FINDINGS.md`](../../session149_tournament_spectrum_20261103/spectrum/FINDINGS.md) supplies the one-apex formula and its proof.
