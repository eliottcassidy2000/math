# Lifting the same three-probe wheel through 27 and 81

This is a finite-exact deeper precision test of the [three physical probe operators](README.md), not a new absolute selector. Use factorial-weighted coordinates x_j=(j−1)! f(j), x_0=0. Each row of the physical probes is tridiagonal, and the coefficient of f(j+1) is a multiple of j. Substituting gives the exact integer recurrences

    P: x'_j = j³ x_j + j x_(j+1) − j²(j−1) x_(j−1),
    Q: x'_j = j³ x_j + (j+1) x_(j+1) + j(1−j)(j−1) x_(j−1),
    R: x'_j = (j³−2j) x_j + (2j+1) x_(j+1) + j(1−2j)(j−1) x_(j−1).

Thus modulo m one uses only x_1,…,x_k, where k is the first index with m|k!: x_(k+1) really is zero. This proves a closed finite state model without separately proving that every vector in the state box is a physical polynomial. Start with the transitive two-vertex polynomial f(t)=t². Every reachable vector then comes from a literal probe word. Record also the modulo-five coordinates and order parity, retain only states with x_1 a unit at 3 and 5, and enumerate both outgoing probe choices at each parity. All counts are odd, so CRT gives the modulus 10m.

| m | Full modulus | Reachable states | Largest strongly connected component | Unit residues there |
|---:|---:|---:|---:|---:|
| 3 | 30 | 40 | 36 | 8 of 8 |
| 9 | 90 | 312 | 300 | 24 of 24 |
| 27 | 270 | 3096 | 2988 | 72 of 72 |
| 81 | 810 | 100765 | 97524 | 216 of 216 |

At 27 the seed lies in the largest component. At 81 the original seed is transient (its own return component is a singleton), and one R step (the count-11 tournament) enters the largest component. Thereafter all unit residues can be requested in any sequence by paths within that component. Thus a deeper mod-3 lift needs an actual change of starting model, even though the local projection still covers all units. This is distinct from the earlier theorem for the uniform-context divisor 27. It makes the “widen the admissible class” analogy testable but does not imply all prime powers or any absolute count coverage. The subpower output bound from the main note applies in full.

Replay with `python depth.py | diff -u depth.out -`. The code constructs the reachable graph and independently of any presumed group action obtains its strongly connected components by forward/reverse graph traversals.
