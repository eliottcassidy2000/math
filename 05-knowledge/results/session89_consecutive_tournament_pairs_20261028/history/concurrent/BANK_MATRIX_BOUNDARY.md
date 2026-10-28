# Nine-prime square-sigma bank = a five-bit rectangle code, not a natural pair invariant

**Status:** elementary all-size rectangle lemma; **FINITE-EXACT independently replayed** instance of the concurrent bank and its affine-action obstruction. Source pin: `b38ec1a1eaf7f4e31a7a38308e2764b5c7650d6d`, path `05-knowledge/results/root_ports_closed_carriers_20261028/BANK_BRIDGE.md`. The exact source is retained as `BANK_BRIDGE.source.md`. This lane independently reconstructs its declared field/probe/prime instance; it does not replay the separate assertion about all eight wedge charts.

## The literal dictionary, not just nine objects on each side

Use old triple rows `(1,4,5),(2,32,34),(8,16,24)` and new triple columns `(1,16,34),(2,4,24)`, in the source's K4 edge-mask coordinates. Their intersections give

```
edge-mask grid          prime grid
[ 1,  4,  5]           [19,  59,  2]
[34,  2, 32]           [ 7,  83, 41]
[16, 24,  8]           [89, 139, 13].
```

For a binary `3×3` matrix M choose exactly those primes whose cells are one, and let `N(M)` be their squarefree product. The independent checker verifies that all nine entries really are distinct primes and reconstructs the square classes of `p+1` using generators `2,3,5,7`.

**Source / target / map:** a subset of this explicitly ordered prime bank maps to its cell-indicator matrix. **Preserved:** the full subset and hence its exact integer product; `sigma(N)=product_(selected p)(p+1)` is square precisely when M lies in the displayed K3,3 cut space. **Lost:** none before a quotient. Forgetting this prime/cell dictionary loses the arithmetic meaning of the matrix.

The source's F64 data are also independently replayed: basis masks `(58,20,1,26,6,11)` in `F_2[t]/(t^6+t+1)`, trace probes `(48,29,25,56)`, and four-bit linear map with columns `(1,2,12,7)`. All 64 trace identities, invertibility of both indicated coordinate maps, norm-one status of all nine ports, and the exact prime square-class dictionary pass. This is a rank-six to rank-four quotient, not a field isomorphism.

## Four rectangles are a complete decision certificate

For an arbitrary `r×s` matrix over F2, the following are equivalent:

1. `M_ij=a_i xor b_j` for row and column potentials;
2. every rectangular four-cell parity is zero;
3. the `(r−1)(s−1)` anchored checks
   `M_00 xor M_i0 xor M_0j xor M_ij=0` hold for `i,j>=1`.

Necessity follows by canceling each potential twice. For sufficiency choose
`a_i=M_i0`, `b_j=M_0j xor M_00`. Then the anchored equation reconstructs every cell. The common shift `(a,b)->(a+1,b+1)` is the only kernel, so the code has dimension `r+s−1`.

For this bank there are exactly four independent checks and dimension five. The 32 accepted subsets among all 512 have weight enumerator

```
1 + 6 z^3 + 9 z^4 + 9 z^5 + 6 z^6 + z^9.
```

Thus first row and first column, counting their intersection once, reconstruct the entire accepted matrix. Rank at most two over F2 is necessary because `M=a 1^T+1 b^T`, but insufficient: the matrix with only `M_00=1` has rank one and fails all four anchored rectangles. The same necessary rank bound holds over Q by `M=a1^T+1b^T−2ab^T`; it is still insufficient.

## The arithmetic predicate does not descend under natural paired relabeling

Use the Session 89 convention `x_ij=M_ij`, `y_ij=M_ji` for `i<j`, and `d_i=M_ii`. Natural simultaneous relabeling of the two tournament cores acts affinely on `x,y`, with the same inversion correction on each.

Start with the zero matrix. Its prime subset is empty and `sigma(1)=1` is square. Transpose the vertex labels 0 and 1 in the natural paired-tournament action. The resulting matrix has exactly `M_01=M_10=1`. Its rectangle syndrome, in order `(i,j)=(1,1),(1,2),(2,1),(2,2)`, is `(0,1,1,0)`. Its selected primes are `{59,7}`, and

```
sigma(59*7) = 60*8 = 480,
```

which is not square. These two matrices encode naturally isomorphic aligned rooted pairs but disagree on the arithmetic predicate. Consequently **no function of the natural pair isomorphism class alone can decide this named-bank square-sigma predicate**.

By contrast the rectangle code is stable under ordinary independent row/column permutations, hence under simultaneous matrix conjugation. This preserves the square/non-square decision for the permuted cell subset, not necessarily the product N itself. The coordinate codec is lossless; the failure occurs only when the wrong quotient is imposed.

The checker exhausts all 512 subsets and all six simultaneous vertex permutations. There are 192 source/permutation pairs for which natural paired relabeling changes the predicate. It finds zero violations for ordinary simultaneous matrix permutations. The exact six-action intersection table is:

| permutation image `(p(0),p(1),p(2))` | accepted codewords whose images remain accepted | zero image accepted? |
|---|---:|---|
| `(0,1,2)` | 32 | yes |
| `(0,2,1)` | 16 | no |
| `(1,0,2)` | 16 | no |
| `(1,2,0)` | 16 | yes |
| `(2,0,1)` | 16 | yes |
| `(2,1,0)` | 0 | no |

Thus **only the identity natural S3 action preserves the full code**. The two three-cycles pass the zero-input test but fail on other codewords; testing the affine offset alone would miss them. Both Python modes pass 3,982 always-active checks with identical receipts, without importing the producer's verifier.

**Boundary:** this explicit finite bank does not classify arbitrary primes with square sigma, produce an infinite number-theoretic family, or settle any major open problem. It is a concrete arithmetic decision/operation mismatch that can guide which information a general pair-matrix representation must retain.
