# Independent audit of transitive-module template growth

**PASS.** The proof in ../branches/PROOF.md has been checked separately from its author's implementation. No mathematical blocker was found.

The deterministic insertion map is valid at every step and injective because restriction recovers the original path. The C3[TT_a,1,1] count is exactly2^a+1. Every vertex of a nontrivial strong quotient belongs to a directed triangle, so retaining that triangle's relevant module proves H(B)≥2^{a_i}+1 for each module. The nonstrong transitive quotient fails this bound and is correctly excluded.

For a fixed bank of m templates of order≤k, the factors L^k, kL and L bound length vectors, physical marking ranks, and prefix clocks respectively. Named relabelings and observer labels cannot add numerical outputs. A singleton strong suffix contributes only L; arbitrary nontrivial transitive blowups of a singleton quotient are not strong and are outside the class.

For arbitrary quotient templates with at most K modules, the sum over orders j≤K is bounded by K²2^{K(K−1)/2}L^{K+2}. With fixed N0, L_n=floor(log2(3^n(N0+2)−2))=Θ(n). Subtracting the singleton O(n) term and taking logarithms gives n≤K²/2+O(K log n+log K). If K=O(√n), the error is o(n), proving K≥(√2−o(1))√n. If K is larger the lower bound is automatic. For a polynomial menu, its logarithmic contribution is O(log n); if K≤n, n≤K log2 n+O(K+log n), proving Ω(n/log n). Larger K already satisfies that conclusion.

The complexity coordinate is the minimum available number of uniform transitive modules in this specified representation. It is not automatically the order of a prime quotient in a canonical modular decomposition. Non-transitive modules with growing internal structure and arbitrary nonuniform attachments remain outside the theorem.

`verify.py` is an independent forward subset DP, unlike the author's backward recursion. Both normal Python and Python -O pass with explicit, non-removable checks:

-1023 induced-subtournament count comparisons;
-26 labelled strong quotients and400 positive blowups with module lengths1 or2;
-direct permutation checks of the triangle formula for a=1,...,6;
-560 exact summation-bound instances, K=3,...,30 and L=1,...,20;
-branch uniqueness and exact maximum through depth12;
-the nonstrong TT2[TT7,TT7] hostile control H=1.

Receipts: normal.json and optimized.json. Finite checks support the proof audit; they do not substitute for the all-order argument.
