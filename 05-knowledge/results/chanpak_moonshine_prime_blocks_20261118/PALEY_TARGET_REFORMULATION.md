# Exact target selection inside one bounded arithmetic host

**Status:** all-order equivalence and explicit bound, not Chan–Pak coverage. This combines inherited Busch/SCC witness bounds with the classical Paley extension argument. No novelty claim for Paley universality, character sums, or finite decidability.

## 1. A physical witness has logarithmic order

Write H(T) for directed Hamiltonian paths. Strong components of a tournament are linearly ordered, and every Hamiltonian path traverses them consecutively, so their counts multiply. Delete singleton components without changing H. Every remaining strong component has order n_i >= 3 and count h_i >= 3. Busch proves

    h_i >= 5^((n_i-1)/3).

Consequently, if N>1 is realized, it has a witness on at most

    K(N) = floor(log_3 N) + floor(log_5(N^3))

vertices: there are r <= floor(log_3 N) remaining components and their total order is at most r+3 log_5 N. Logs in the definition can be computed exactly with integer powers. N=1 has a singleton witness and we set K(1)=1. For prime N one may instead use K=1+floor(log_5(N^3)), since precisely one nontrivial component remains.

This reduction already appears in the repository's `actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md`; it is not a new decidability result.

## 2. An explicit host containing every required tournament

For K>=1 take q=3^e, with e the least **positive odd** exponent such that

    q > K^2 4^K.

Then q is a prime power, q=3 mod 4, and q<=9K^2 4^K. Let chi be the quadratic character of F_q, extended by chi(0)=0. Orient x -> y precisely when chi(y-x)=1. Because chi(-1)=-1, this is a loopless tournament.

Here are the constants and the distinctness correction in the extension argument. For t>=1 distinct anchors a_1,...,a_t and any signs epsilon_i, let C count x outside the anchors satisfying chi(x-a_i)=epsilon_i for every i. Expand

    F(x)=2^(-t) product_i (1+epsilon_i chi(x-a_i)).

Off the anchors this is exactly the desired indicator. At an anchor it lies between 0 and 1/2, so their total contribution is at most t/2. For a nonempty j-element subset, its product polynomial has j distinct simple roots and is not a constant times a square. The Weil bound is (j-1)sqrt(q); for j=1 its character sum is exactly zero. Since

    sum_{j=1}^t binom(t,j)(j-1) = (t-2)2^(t-1)+1,

we obtain

    C >= q/2^t - ((t-2)/2 + 2^(-t)) sqrt(q) - t/2.

For 1<=t<=K, sqrt(q)>K2^K>=t2^t. The coefficient of sqrt(q) in the subtraction is at most t/2. Thus C>t(sqrt(q)-1)/2>0. For t=0 any field element suffices. Greedily prescribe the orientations of each next vertex to the previous ones. This embeds **every** tournament of order at most K (indeed K+1) as an induced subtournament. The allowed signs incorporate whether an old vertex beats the new one or vice versa. No repeated field coordinates are admitted.

## 3. An exact integer equation, not a residue condition

For N>0 use the host selected above. N is a tournament count if and only if there are 1<=m<=K and distinct field elements x_1,...,x_m with

    sum_{pi in S_m} product_{j=1}^{m-1}
       (1+chi(x_{pi(j+1)}-x_{pi(j)})) = 2^(m-1) N.

The sum and products here take place in **the ordinary integers**, with character values -1,0,1. Each permutation contributes 2^(m-1) precisely when it is a Hamiltonian path of the induced tournament, otherwise zero. Empty products for m=1 are 1. This gives both directions of the equivalence: the finite-field coordinates are a physical tournament, and the witness bound plus universality embeds any existing witness.

For c=1/log(3)+3/log(5), the field size satisfies, for N>=3,

    q <= 9 c^2 (log N)^2 N^(c log 4),

where c log 4 is approximately 3.845919. Thus the ambient host has polynomial size in N, not in the bit length of N. Subset enumeration remains enormous: q^O(K)=N^O(log N), even before path-count evaluation. This does not provide an efficient decoder or imply target solvability. It replaces the sparse one-apex interface with a universal physical arithmetic interface; the remaining theorem would assert solutions to this exact equation for all odd N except 7 and 21.

## 4. Why this changes the question, and what it does not change

`fermat_flip_stern_20261117/README.md` section 3 provides a target-first inverse inside a one-apex family. That interface is powerful but family-restricted. Here there is no such restriction: missed values cannot be blamed on choosing the wrong family. A useful next approach is an exact coefficient/character-sum estimate over subsets, with the target-count indicator retained. Ordinary Paley quasirandomness only supplies adjacency patterns; it does **not** show that the prescribed Hamiltonian-count level is occupied.

This is also a precise limit on the Collatz analogy. Quadratic-character tests live in the cyclic multiplicative group F_q^*, while adjacency uses additive differences. Their interplay packages every bounded tournament in one arithmetic object. It does not construct a decreasing target recursion, and no convergence implication between Collatz and Chan–Pak is asserted. SCC multiplication remains the exact multiplicative renewal law; a universal ambient host does not turn it into additive target control.

## Sources and priority

- Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC 13 (2006), N3: https://doi.org/10.37236/1141 . The uniform lower bound above is weaker than the paper's residue-dependent sharp formulas.
- R. L. Graham and J. H. Spencer, *A Constructive Solution to a Tournament Problem*, CMB 14 (1971), 45–48: https://doi.org/10.4153/CMB-1971-007-1 . Classical quadratic-residue tournament extension via character sums. The mixed-sign prime-power argument and conservative constant are written out here rather than attributed verbatim to that paper's common-dominator formulation.
- C. Even-Zohar, M. Farber and L. Mead, *Ample simplicial complexes*, Appendix A, Theorem A.1, states the finite-field multiplicative-character Weil bound used here: https://arxiv.org/pdf/2012.01483 .
- Repository SCC witness-size priority: `05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md`.
- Repository restricted arithmetic inverse: `05-knowledge/results/fermat_flip_stern_20261117/README.md`, section 3.

## Replay

`python paley_target_replay.py` uses exact polynomial arithmetic over F_(3^7), with an irreducible polynomial independently selected by SymPy. It embeds all 75 labelled tournaments of orders 1 through 4 and checks the integer character equation against independent subset dynamic programming. It also checks all 1,404 signed two-anchor patterns in F_27. These are scoped finite implementation tests, not evidence establishing universal spectrum coverage. `paley_target_replay.json` records explicit matrices and embedded coordinates for the counts encountered.
