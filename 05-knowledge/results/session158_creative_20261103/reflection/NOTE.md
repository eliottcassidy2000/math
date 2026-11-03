# A reflection-sandwich line and its quadratic response

Status: elementary all-order construction and exact algebraic identity, conditional only on the stated classical odd-cycle formula when that formula is invoked. Finite checks are independent illustrations. No claim of novelty against the entire literature and no assertion about the full Hamilton-path spectrum.

## Construction and theorem

Let A be a nonempty tournament, p=(p_1,...,p_k) one of its directed Hamiltonian paths, A* a disjoint converse copy, and c a new vertex. Form the ordinal sum

    S = A -> {c} -> A*.

Let alpha exchange each vertex with its corresponding converse copy and fix c. Write p* for the corresponding sequence in A*, and use the reflected path P=(p_1,...,p_k,c,p_k*,...,p_1*). Given any tournament T with specified Hamiltonian path q, write K_q(T) for the tournament obtained by reversing every arc outside q.

Then S and its partner Y=K_P(S) are self-converse under alpha; P is a reflection-fixed (blue) path presentation at both endpoints, and Y is strongly connected, since its surviving path P plus the reversed endpoint chord gives a directed Hamiltonian cycle. The first endpoint has H(S)=H(A)^2; the two endpoints are an individual blue line in the converse-merged line multigraph. They need not be different merged nodes.

For the exact second endpoint count set X=K_p(A), a=p_k, h=H(X), s=the number of Hamiltonian paths of X starting at a, and

    r = sum_{U subseteq V(X), a in U} E(X[U],a) H(X-U),

where E(X[U],a) counts directed Hamiltonian paths ending at a and H(empty)=1. Then

    H(Y) = r^2 + 2hr + (h-s)^2.                                   (1)

Let F_t count directed (2t+1)-cycles of Y fixed by reversal and alpha. There is an exact, length-preserving bijection:

    F_t = sum_{U subseteq V(X), a in U, |U|=t} E(X[U],a).           (2)

In particular F=sum_t F_t is the number of all simple directed paths in X ending at a. Redei oddness of every X-U implies F=r (mod 2). Redei oddness of Y and (1), modulo 2, give r=s (mod 2), recovering Session 29's prior portal parity in the restricted case where a is the end of a specified Hamiltonian path in X; Session 29's result applies more broadly. Oddness of h then gives H(Y)=1+2(s mod 2) (mod 4), so F=s (mod 2). Thus an input endpoint-path statistic determines the fixed odd-cycle profile of a new strong self-converse tournament.

For the first endpoint there is also a full polynomial statement: if J_T(z) counts vertex-disjoint collections of directed odd cycles by number of cycles, then J_S(z)=J_A(z)^2. Every cycle is internal to one of the two copies, and cycles in the copies do not conflict. The alpha-fixed collections have polynomial J_A(z^2). At z=2 the first identity agrees with H(S)=H(A)^2. These are exact polynomial identities, not a claim that the strong endpoint has the same polynomial.

## Proof of (1)

Under chord complementation, the two k-vertex blocks become X and X*; every cross-block arc points from X* to X. The central arcs point from a to c, from c to a*, from c to X-{a}, and from X*-{a*} to c. Thus deleting c gives X* -> X, with exactly the two stated portal exceptions.

A Hamiltonian path of Y splits at c. Within each side the X* vertices, if present, must occur before the X vertices, and each block is an internal Hamiltonian path. Split into four cases according to whether the portion of X before c and the portion of X* after c are nonempty.

* Both nonempty: the former must end at a and the latter must start at a*. Summing independently over subsets gives r^2; converse makes the two portal sums equal.
* Only one nonempty: each of the two cases contributes hr. A full path of the other copy is unrestricted, and absence of its portal from the opposite part makes the seam valid.
* Both empty: the full X* path must not end at a*, while the full X path must not start at a. Each count is h-s by converse, giving (h-s)^2.

These cases include empty before/after pieces through H(empty)=1 and are mutually exclusive and exhaustive, proving (1). For n=2k+1>=3 the first and last vertices of P are nonconsecutive, so the endpoint chord is actually reversed. Reflection alpha reverses S by construction and commutes with chord complementation. Thus P closes in Y and all claimed symmetry and strong connectivity hold.

## Exact fixed-cycle correspondence

A reflection-fixed directed odd cycle has alpha-invariant odd support and hence contains c, the unique fixed vertex. Starting it at c, reflection invariance forces its sequence uniquely into the form

    c, v_1,...,v_t, alpha(v_t),...,alpha(v_1), c.

The midpoint arc v_t -> alpha(v_t) and the uniform cross-block orientation X* -> X force v_t into X*. Since an X vertex cannot precede an X* vertex along a directed path, all v_i lie in X*. The only central arc from c into X* is c -> a*, hence v_1=a*. Conversely, every simple directed path in X* beginning at a* produces a fixed odd cycle by this construction. Reversing its corresponding path in X identifies it with a simple path ending at a. Counting the vertex support gives (2) term by term.

## Finite examples and reproducibility

Use vertices 0,...,k-1, force i->i+1, list free pairs (i,j) lexicographically with j>=i+2, and reverse bit q when q is 1. An exact subset DP independently counts both endpoints of every such A through k=5, all 76 presentations in total, checks formula (1), and verifies reflection and the spanning cycle. It also checks k<=3 counts by separate permutation enumeration. Run `python check.py` and `python -O check.py`; both outputs are byte-identical (`check.out`, `check-opt.out`). The separate `cycles.py` directly enumerates fixed odd cycles through k=4 and verifies (2) by length; see `cycles.out`. This is not an exhaustive enumeration of all tournaments of sizes 2k+1.

| A | H(A) | H(S) | (h,s,r) | H(Y) |
|---|---:|---:|---|---:|
| transitive k=1 | 1 | 1 | (1,1,1) | 3 |
| transitive k=2 | 1 | 1 | (1,0,2) | 9 |
| transitive k=3 | 1 | 1 | (3,1,3) | 31 |
| directed triangle k=3 | 3 | 9 | (1,0,4) | 25 |
| k=4, mask 1 | 3 | 9 | (5,1,5) | 91 |
| k=4, mask 4 | 3 | 9 | (5,1,9) | 187 |

The last two examples show why an unmarked scalar gadget law H(Y)=f(H(A)) cannot work for this operation: the portal's r matters even when both H(A), H(X), and s agree. These are blue lines whose endpoints are self-converse and whose second endpoint is strong. They do not furnish the offset operations required in the Session 151 ternary compiler, and they do not resolve the odd-spectrum conjecture.

## Source links

* [Session 151 metagraph structure](https://github.com/eliottcassidy2000/math/blob/3d6628b1804cf80457da16aa52489f016dfc99c4/05-knowledge/results/session151_spectrum_metagraph_20261103/metagraph/STRUCTURE.md) and [path overlap](https://github.com/eliottcassidy2000/math/blob/3d6628b1804cf80457da16aa52489f016dfc99c4/05-knowledge/results/session151_spectrum_metagraph_20261103/bridges/OVERLAP_AND_PATH_WALK.md) establish the line convention, path-fibre, and closed-side caveat.
* [Session 29 converse congruences](https://github.com/eliottcassidy2000/math/blob/3d6628b1804cf80457da16aa52489f016dfc99c4/05-knowledge/results/session29_framed_profiles_20261022/session29-converse-congruences/REPORT.md) supplies the odd-cycle fixed-packing framework, portal parity, and an even-order two-copy seam formula; the present odd-order construction and its four-case identity should be read against those antecedents.
* [Session 157 cycle proof](https://github.com/eliottcassidy2000/math/blob/3d6628b1804cf80457da16aa52489f016dfc99c4/05-knowledge/results/session157_gaps_20261103/cycles/README.md) records the odd-cycle packing formula and established 7, 21 gaps.
