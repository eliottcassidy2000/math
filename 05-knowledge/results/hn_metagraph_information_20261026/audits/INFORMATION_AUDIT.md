# Independent tournament-response and historical-information audit

2026-10-26. **FINITE-EXACT PASS**, with elementary all-order algebra checked
against the information lane's proof. The independent implementation is
[check_information_independent.py](check_information_independent.py); its
saved output is [information_independent.json](information_independent.json).
Normal and `python3 -O` outputs agree byte for byte.

The computation constructs the entire Hamiltonian-path-count function H by
adding one for every directed Hamiltonian path to all tournaments containing
its arcs. This path-subcube method is independent of a per-tournament subset
dynamic program. Integer fast Walsh transforms are then applied, with exact
rational normalization only at the end. No floating-point arithmetic or
rounding is used.

For each order2..6, all `2^binom(n,2)` labelled tournaments are covered. Every
odd coefficient vanishes, and every pair of incident arcs has degree-two
coefficient magnitude `(n-2)!/2^(n-2)`. The signs reflect the fixed orientation
gauge. The `m-1` sign products along an actual spanning tree of the line graph
of `K_n` have exactly two preimages each, which are global converses. Literal
full flip-response tables are independently grouped through order four; their
fibres are exactly the same converse pairs. Through order six the program
checks the coefficients and all sign-product fibres, without materializing
the quadratic-size response-table matrix. The order-two coefficient assertion
is vacuous because there is only one arc; its response has one converse fibre.

The proof's distinction between information and queries is essential:
`m-1` binary correlations suffice to encode a converse pair, but obtaining
those correlations by a complete Walsh sum can require exponentially many
H evaluations. No efficient query algorithm or A000568 count is inferred from
this information bound.

## Demonstrated errors in the recovered March note

For `P_beta(T)=exp(beta H(T))/Z(beta)`, differentiation of the finite sum gives
`I(beta)=Var_beta(H)` and `I'(0)=E[(H-E H)^3]`. The exact order-five census gives

```
E H = 15/2,                  Var H = 285/16,
I'(0) = 165/32 > 0,
Var(E[H | labelled score vector]) = 1935/112,
sum_{|S|=2} Hhat(S)^2 = 135/8.
```

Thus zero is not a local Fisher maximum, and the two explained-variance
fractions are respectively `129/133` and `18/19`. A sorted score vector gives
the same conditional values here by relabelling invariance. The difference
is conceptual as well as numerical: a score is not a Walsh-degree cutoff.
Higher Fourier terms can contribute to a score-conditional expectation.

The note also confuses connectivity with simple connectedness and higher
persistent-homology vanishing. The independent information-lane check gives
a direct landscape witness: the order-three H=1 cubical sublevel is a six-cycle
with no square faces, so its first Betti number is one. Its Euler-characteristic statement conflates
the filled cube (contractible, Euler characteristic1) with its graph
(`2^m-m*2^(m-1)`). These are type errors; the old sublevel computations would
need separate audits before any topology claim could be recovered. The
historical note is now explicitly bannered, and the mechanism is recorded in
the mistakes ledger. Its elementary converse/odd-Fourier argument survives.

Reproduction from the repository root:

```
python3 05-knowledge/results/hn_metagraph_information_20261026/audits/check_information_independent.py
python3 -O 05-knowledge/results/hn_metagraph_information_20261026/audits/check_information_independent.py
```
