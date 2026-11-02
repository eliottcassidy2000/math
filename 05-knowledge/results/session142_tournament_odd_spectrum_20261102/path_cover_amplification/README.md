# Session142 path-cover amplification

For an n-vertex tournament T, unordered k-path covers satisfy

    k! p_k(T) >= H(T) binom(n-1,k-1).

Combining this incidence bound with the inherited homogeneous run kernel proves

    H(Q[T]) >= (n+2)H(T)

whenever the substitution port of Q is cyclic, with an additional +1 for even n. The context may vary with T. An inactive port instead gives the odd multiplier H(Q-q). Thus cyclic homogeneous-module substitutions cannot support bounded multiplicative radices on an unbounded-order source family. Nonhomogeneous S0/R1 source moves are outside this theorem.

This is a separate bounded note. It does not change the earlier frozen binary-language packet, and it proves no new full-spectrum coverage result. The {7,21}-only gap conjecture remains open.

## Reproduce

    python3 -B verify.py
    python3 -B -O verify.py
    python3 -B verify.py --max-order 3

The selector restricts the source-tournament census; the default is all labelled sources through order5. All pointed contexts through order4 are checked in every mode. Only the Python standard library is used. The checker runs independently of source snapshots and the repository.
