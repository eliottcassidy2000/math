# Independent review: partial bridge-free C=2/C=3 reduction

**Result: PASS**, with no correction. The packet proves a bounded residual problem, not a complete C=2 classification.

## Bound snapshot

Reviewed files in `session197-lowcycle/`:

- `C2_CUT_REDUCTION.md`: `cf8017dc9d1a97421f9406ae0ee34655d2934a73437c04e253fbdf29d5927cb8`
- `verify_c2_cut_reduction.py`: `d88470cd13c20cf0f1bd06b118d74ca93ddbd02cb32ab84cee6aec2bde6724da`
- `c2_cut_reduction_certificate.json`: `0f1b6593980e4ac9a8feb21716c3ff8bf323c8d2ebf9008595a6213fe19d10c9`

The pinned helper `verify_c1_complete.py` has hash `645ddfc551fd5ce2760f2d1f71318c08f57e1115e54edf9d28d7db3e6ad91d3b`. All files were copied before replay because the verifier writes its certificate beside itself. The original files were not modified.

## All-order argument and primary-source check

The already reviewed separating-vertex product has positive initial and terminal endpoint factors. A factor one exposes the corresponding sole reverse cut arc as an actual strong bridge. Consequently a bridge-free tournament with a cut vertex has at least four Hamiltonian cycles. Both the initial and dual terminal cuts are correct. This proves that bridge-free C=2 and C=3 tournaments are 2-strong. It does not assert sharpness of the bound four.

I opened the specifically cited [official Busch paper](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf), Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC 13 (2006), N3. Corollary 2 on its final page gives the stated cycle lower bound with beta=5^(1/3). The paper's use of 2-connected in its structural discussion is vertex-strong connectivity, not merely the absence of a strong bridge.

The imported inequality is therefore applicable to the residual class. Raising it to the 96th power gives `5^(n-32)<=c^96`. The exact integer comparisons at exponents 41/42 and 65/66 correctly give n<=73 for c=2 and n<=97 for c=3; orders below 32 cause no issue. No floating-point rounding is used. This establishes only the order bound, not the conjectured inventory through those orders.

## Independent controls

Source normal and optimized runs match the frozen receipt byte-for-byte.

The independent checker reuses the reviewer's graph primitives from the C1 review, not the author's enumeration helper. It performs a full directed-cycle DFS and bitset reachability across all 33,867 labelled tournaments through order six.

- All 15,906 low-cycle cut-vertex cases expose a valid endpoint bridge.
- The five-vertex 2-strong histogram is exactly 24 graphs with two cycles. Each is independently checked isomorphic to the cyclic regular five-tournament.
- The six-vertex 2-strong histogram is exactly 720 graphs with four cycles and 1,920 with five cycles. Thus the reported minimum four and absence of small new C=2/C=3 atoms are correct finite statements.
- The seven-vertex example is strong, has nine cycles, has no strong bridge after testing every arc deletion, and has exactly separating vertex 3.
- The hostile six-vertex graph is 2-strong and has exactly five Hamiltonian cycles. Its arc 2 to 5 occurs in none, so the packet correctly rejects the tempting stronger arc-coverage assumption.

The published bound and cut-product proof establish the all-order reductions. The small census remains explicitly finite evidence. No unresolved source, scope, arithmetic, or implementation concern remains in this snapshot.
