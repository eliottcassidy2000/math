# Pre-registered predictions (opus-2026-10-26-S62), written before the order-7 and order-8 runs

Context: on the marked carrier X_n (2^M states, M = (n-1)(n-2)/2 free pairs), a named guard g_p is the indicator of the
subcube of states whose free pairs along the Hamiltonian word p are oriented as p requires, so a k-fold product
g_(p1) ... g_(pk) is the indicator of a subcube fixing at most k(n - 1) free pairs, and the span of all k-fold products
lies in the Boolean polynomials of degree <= k(n - 1) (PROVED, note section 1). Observed (guard_ranks_paired_cells, n <= 6):
single guards have rank 8, 44, 273 at n = 4, 5, 6; pair products reach the full rank 8, 64, 1024; the paired cells
([x], [tau x]) form a connected bipartite graph with E = 7, 50, 732 edges, V = 8, 24, 112 vertices, cycle rank 0, 27, 621.

Predictions:
1. Order 7: pair products do NOT span the 32768-dimensional presentation-law space (degree bound 12 < 15; proved), and
   the restore depth d(7) (the least k such that every state of X_7 is the common-path cube of k of its own Hamiltonian
   paths, i.e. k paths whose free pairs cover all 15) equals 3, the degree-bound value ceil(15/6).
2. Order 8: on a random sample of 10,000 states of X_8, every state is covered by 3 of its Hamiltonian paths (degree
   bound ceil(21/7) = 3), so the restore depth stays 3.
3. Order 7 paired cells: the bipartite graph of ([x], [tau x]) pairs over the 32768 states is connected (c = 1), so the
   independent terminal responses number 2 x 456 - 1 = 911 and the lost directions are E - 911.
4. The restore-depth sequence follows the degree bound ceil((n - 2)/2) for n = 4..8: 1, 2, 2, 3, 3.
