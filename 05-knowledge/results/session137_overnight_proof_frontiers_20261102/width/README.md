# Session137: exact modular contexts and branching resource

New all-order results:

1. The coefficient of the c-path-cover coordinate in every homogeneous
   Hamilton-path context is divisible by `c!(c-1)!`, and this divisor is sharp.
2. Equality in all contexts modulo M is **equivalent** to agreement of
   `p_T(c) mod [M/gcd(M,c!(c-1)!)]`; only finitely many coordinates survive.
   An explicit representative-free substitution kernel computes these states.
3. The exact universal nonempty contextual class counts are5,3,14,11 at
   moduli3,4,8,12, respectively, with width-three witnesses for every class.
4. Arbitrary ordinal joins and unary triangle wraps never give H congruent7
   or11 modulo12. Genuine multi-child triangle branching escapes, giving
   infinite width-three families `11*25^k` and `55*25^k` excluded by the smaller
   grammar.

The old cyclic decoder fails at M8 on L3/L4; the new sharp detector succeeds.
This is an infinite constructor hierarchy, **not** a result that spectral
substitution width is unbounded. The full odd spectrum remains OPEN.

See [PROOF.md](PROOF.md). Reproduce with `python3 -B context_residue.py` and
`python3 -B -O context_residue.py`; outputs `normal.json` and `optimized.json`
match byte-for-byte. No third-party packages or repository imports are used.
