# Independent review: square-root fixed-prefix inverse

**PASS.** No mathematical or implementation correction was required. This review concerns the exact final-gadget inverse for a fixed response triple, not selection of a reachable prefix or all-order scalar coverage.

## Frozen source

The reviewed files are in `session197-lowcycle/`:

- `SQRT_FIXED_PREFIX_INVERSE.md`: `084b7c372524f41e99a4087179e0fd30baf35a0de58a0166e3bdcd84c5e407bf`.
- `verify_sqrt_fixed_prefix.py`: `124dd0e43261ddda1f1394d114c564b8eeb40b71431f294b4b1375fadc033aa2`.
- `sqrt_fixed_prefix_certificate.json`: `68ce39a7089b807bbe164072fbebc2e9e62175d130489f8e671ed8e31ada97be`.

The two named physical dependencies were copied unchanged and their hashes were checked by the source verifier. Normal and optimized Python replays both reproduce the frozen receipt byte for byte. Originals were not edited.

## Proof checks

1. Writing `f(a)=g(a+1)`, the append-one and append-zero laws are valid also at `a=0`. In a sequence of appended digits, the coefficient of the original high integer is the reverse-complement low word, not its complement in the original order. Only the first appended zero can retain the high integer's parity. This gives exactly the displayed block identity and its constant term by evaluation at zero.
2. `1<=g(L)<=L` gives `F(L)>=RL` under `0<=e<=h`. Thus every solution lies in `1<=L<=K=floor(N/R)`. For `K>=2`, the selected power `B` satisfies `K<=B^2`, `B<=K`, and `B<2sqrt(K)`. The separate `K=0,1` cases are complete.
3. The decomposition `L=Bx+b+1` is unique. Its allowed high interval lies in `0<=x<=B-1`; parity of `L` is constant on a branch. For consecutive points, `g(x+2)-g(x+1)>=-x`, giving the strictly positive lower bound `RB-hx+h(q-c)`. This does not assume monotonicity of `g`.
4. Binary search followed by equality testing therefore returns every solution and no false solution. There is at most one answer per branch. The implementation covers branches whose final value is below the target, one-element intervals, both parities, and the `L=1` boundary correctly.
5. Tabulation, bit reversal, searches and optional sorting have the claimed `O(sqrt(K) log K)` arithmetic-operation bound and `O(sqrt(K))` integer-word storage. For `K>=1`, `h<=R<=N`; on admissible search arguments the products are bounded by a constant multiple of `N`, so the intermediate bit lengths are `O(log N)`. The claim is not polynomial in input bit length, and the proof says so explicitly. Handling a large input with `K=0` still has the ordinary input-reading/division cost, outside that positive-`K` intermediate bound.

## Separate implementation and graph counters

The independent checker is `sqrt-fixed-prefix/check.cpp`, SHA256 `638e5be2f221e7b3cd78cda89f8344e2281e820420c2faf72550323c19bdae60`.

It imports no source code. Its gadget count assigns each vertex of the transitive deletion to the subsequence before or after the root. A nine-state incidence inventory tracks the last vertex before the root and first vertex after it. This directly counts literal Hamilton words and does not use the recurrence for `g`.

The independent inverse evaluates this path-partition count at the actual integer `L` during binary search. It does not use the block identity or the source's `g` table to evaluate a branch. A second explicit enumeration of every before/after subset checks the small gadget counter. Physical graph counts use inclusion-exclusion over vertex subsets of length-`n-1` walks, rather than the source's endpoint subset DP.

The optimized build and a separate build with library assertions and undefined-behavior sanitization produce identical receipts:

- 4,096 complete literal before/after subset counts;
- 131,838 block identities;
- 710,160 consecutive branch inequalities;
- 140,976 complete inverse comparisons against exhaustive final-`L` scans;
- 224 actual gadget/seam graph counts by walk inclusion-exclusion;
- all three larger exact unique inverses, with the target integers in the proof.

Both independent receipt files have SHA256 `dc7115151fd25d3642560e3f7adef808b19efa796ccab575a2d4af2b07d4c1ba`.

Reproduce from the audit directory with:

```sh
g++ -std=c++17 -O3 sqrt-fixed-prefix/check.cpp -o /tmp/check-sqrt-prefix
/tmp/check-sqrt-prefix
g++ -std=c++17 -O1 -D_GLIBCXX_ASSERTIONS -fsanitize=undefined sqrt-fixed-prefix/check.cpp -o /tmp/check-sqrt-prefix-debug
/tmp/check-sqrt-prefix-debug
```

The finite controls support implementation correctness. The all-target completeness and complexity statements rest on the proved block identity, finite bound and strict branch monotonicity. They do not imply physical inverse closure of arbitrary response triples, nor existence of a suitable prefix for an arbitrary odd target.
