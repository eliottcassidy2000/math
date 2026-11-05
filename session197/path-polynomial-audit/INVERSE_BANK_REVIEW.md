# Independent review: finite-bank inverse search

**PASS.** Reviewed frozen files in `../session197-odd-cycle-path-polynomials/`:

- `INVERSE_BANK.md`: `0ee93289a4ccb875b2f21c083fcc8241ccecd2cdc270fe4e076707b2541a6773`
- `verify_inverse_bank.py`: `1f4df4c4a42c14e6ab152dca08233f520765dd3422b8c865e8b87206d0a4d1f7`
- Previously reviewed dependencies `verify_extension.py`: `3b537a9cbfea66aabe8b1d4cf59c84de45cc41d2ea6f80d053df8c083b1cad76`, and `verify.py`: `10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb`.

Copies were run with normal Python and `-O`. Their identical stdout has SHA-256 `880740e19fc8ae1794fd95dc0112af610e7ef1732e332f45d5bb55feab8d64bc`. Frozen originals were not edited.

## All-order proof

An active root has a crossing triangle. Its substituted graph therefore contains the physical cyclic wrapper D(T), with H(D(T))=H(T)+B_T(2)>=2^|T|. Induced insertion supplies monotonicity. The entire actual input order is bounded, including neutral padding and internal joins. A root in no triangle belongs to no directed cycle and is a singleton strong component; a passive template consequently acts on H by the fixed multiplier H(Q-r)=H(Q).

The witness-order argument strips only a finished outer shell of joins and passive templates. It stops at each active wrapper and keeps its full input. Thus it never equates two same-H profiles across an active continuation. Discarded count-one pieces are removed only from the final scalar product, where that removal is legitimate. Passive nonunit factors have the constructible seed Q=C[K1]. Each retained factor h>1 therefore has a representative of order at most floor(log2 h)+M. The number of factors is at most Omega(P), and the sum of floors is at most floor(log2 P). This proves the stated order bound without assuming seven is absent or deleting any internal padding.

After identity templates are removed, every operation increases degree. Full-profile dynamic programming through the proved order bound is complete by expression induction; count pruning is safe by induced monotonicity. Equal full profiles, not equal H values, are interchangeable. The fixed-bank expression-tree bound is exponential in degree, while the degree cutoff is O(log P), proving polynomial complexity in numeric P. It does not give polynomial complexity in target bit length or an inverse theorem for unrestricted templates.

The 41/43 omissions use exactly the eight stated parity readouts. Their candidate lists are complete. The H=1 candidates must be transitive and fail the power-of-two B(2) test. The remaining H=5 odd-order candidate for 43 would have B(2)=16, hence order at most three, which is impossible. These are grammar omissions at every order, not tournament-spectrum exclusions.

The whole-module obstruction is also correctly scoped: an active context dominates D(T), and a passive context multiplies the odd H(T)>=3. Neither can produce 3H(T)+2 or 3H(T)+4 from the fixed nontransitive input. Nonuniform attachment, splitting the input's neighborhood types, and replacing a different single vertex are not ruled out.

## Independent finite verification

`inverse-bank-audit/independent.py` enumerates the entire bank using integer evaluation tables f(0),...,f(n), forward differences and integer extrapolation. It does not use the source monomial addition, multiplication, shifting, profile storage or generation algorithm. It enumerates both degree orders of joins. At the end, each independently generated full value-table set is compared against a fresh source result, so the test checks every full profile, not merely H values or state counts.

All 12,553 profiles agree, with per-order counts:

`1,1,2,3,7,15,29,61,127,270,498,776,1173,1728,2197,2642,3023`.

Independent trial-division primality reproduces all 1,027 odd primes through 8,191 and the exact same 77 attained primes.

`inverse-bank-audit/walk_check.cpp` separately checks all 77 supplied prime-witness graphs and both explicit 197/211 graphs. It counts length-(n-1) unrestricted directed walks on every induced vertex set, then applies inclusion-exclusion over omitted vertices. A surviving walk visits all n vertices exactly once, so this gives H without either source endpoint-subset dynamic program. All 79 checks pass; 128-bit accumulators avoid cancellation overflow.

No mathematical, scope, or implementation correction is required. The all-order theorem, not the numerical cutoff alone, is essential to interpreting the finite search.
