# Rooted transitive-block extension: exact arithmetic rays

Status: proved by the argument below; independently checked by direct subset DP on all labeled tournaments of orders 2–5, all roots, for each k=1,2,3 (16,212 tests). This is a construction, not a proof that the tournament Hamiltonian-path spectrum has no further gaps. Priority is not established: targeted searches of the session25/26 spectrum packets and tournament theorems found no identical formula, but this is not an exhaustive literature or repository search. Session149 already proves that H=13 is not one-apex; no novelty claimed for that obstruction.

Let T be a tournament with distinguished vertex v, let U=T−v, and write h=H(T), q=H(U), and s for the number of T-Hamiltonian paths starting at v. Let F be the number of ordered pairs (A,B) of directed paths whose vertex sets partition U, both empty paths allowed, such that A is empty or the final vertex of A points to v. As usual the empty induced tournament has one empty path. In particular, A=empty accounts for exactly q pairs.

For k≥1, add k vertices X forming any transitive tournament, orient v→X and X→U, and leave T unchanged. Then

    H(T_k)=h−s+q+(2^k−1)(F−q)
          =(h−s+2q−F)+(F−q)2^k.

The number of Hamiltonian paths of T_k starting at v is exactly q.

Proof. No vertex of U points to X; the only old vertex pointing to X is v. Consequently, in a Hamiltonian path, X can occur only in an initial block and/or a block immediately following v. In each block the order is uniquely fixed by transitivity. Select which X vertices occupy the initial block.

- No X vertex initially: write the path A,v,X,B; this has F possibilities.
- All X vertices initially: the rest is a Hamiltonian path of T which cannot start at v, giving h−s possibilities.
- A proper nonempty subset of X initially: A in A,v,X,B must be nonempty, giving F−q possibilities, independently of the subset. There are 2^k−2 subsets.

Summing yields the formula. A path starting at v must have all X immediately after v and then a Hamiltonian path on U, proving the root count. Equivalently, compute F by summing over S⊆U the product of (i) the number of Hamiltonian paths on S whose last vertex beats v, or 1 for S empty, and (ii) H(T[U−S]). The two path blocks need not concatenate to a T-Hamiltonian path: this is exactly the missing broken-adjacency contribution.

If every u∈U has a directed path to v in T, then T_k is strong: each old vertex reaches v, v reaches X and X reaches all U. In particular, this holds when T is strong. This condition gives strong realizations of the rays.

Example: order the pairs (i,j), i<j, lexicographically and let bit 1 mean i→j. Take n=5, bitmask 41 and v=4. Direct DP gives (h,s,q,F)=(13,4,5,14); the seed is strong. Thus there are strong tournaments with H=5+9·2^k for every k≥1. The finite census in rays.out lists all 39 distinct nonconstant pairs (α,β) found at orders 2–5; this is a finite computation, not a classification at arbitrary order or a completeness proof.

Counterexample to a tempting but false shortcut: take T on {0,1} with 1→0 and v=0. Then (h,s,q,F)=(1,0,1,2), and T_1 is a directed triangle with H=3. It is false that H(T_1)=2h−s or that the v-start count remains s. Deleting the new vertex from the path 0→x→1 creates the non-arc 0→1, hence need not leave a path.

Replay: `python3 test.py` (16,212 checks); `python3 rays.py`. Both use an independent subset-and-last-vertex Hamiltonian-path DP, distinct from the block proof. `test.out` and `rays.out` are frozen stdout. Every small seed is encoded by n, lexicographic edge bitmask, root, and its profile.
