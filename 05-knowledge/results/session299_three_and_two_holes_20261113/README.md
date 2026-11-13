# The two holes and three distinct meanings of multiplication by 3

The [original MathOverflow question](https://mathoverflow.net/questions/232751/the-number-of-hamiltonian-paths-in-a-tournament) and [THM-1370](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md) establish the exclusions; completeness remains **open**. This packet distinguishes the obstruction, its multiplicative reverse, and a physical substitution whose surplus measures return through a root.

## 1. Which steps reverse?

1. **Condensation.** Strong components of a tournament are totally ordered, and paths concatenate, so (H) is the product of their path counts. This is reversible: concatenate any chosen components. It has no additive inverse. In particular 7 could only be carried by a strong component of count 7; 21 would need a component of count 7 or 21.
2. **Insertion.** Every Hamilton path on a deletion extends by inserting the missing vertex: insert before the first vertex it beats, or append if none. The chosen insertion is injective. The reverse deletion is not surjective: some paths become disconnected when the vertex is removed. Even when the full count is unchanged, the available marked continuations can change.
3. **Strong order.** Moon's strong-deletion theorem propagates lower bounds, not exact values. The Moon–Busch minimum is 3, 5, 9, 15, 25 at orders 3 through 7. Larger strong tournaments are above 21. Extremal interval chains give actual witnesses at the floor, but do not fill the interval above it. Compare [THM-4866](../../../01-canon/theorems/THM-4866-single-cycle-floor-bottom-band-uniquely-hamiltonian-one-arc-lorentzian-joins-and-crevasses.md).
4. **Local audit.** At orders 3 and 4 the strong counts are exactly 3 and 5. The strong order-5 and order-6 checks exclude 21; these are finite statements, not a pattern inferred from the two floor gaps. The accompanying independent permutation/subset check reconstructs these small inputs.
5. **Reverse target.** Construct the primitive counts that cannot be factored, rather than trying to reverse a floor inequality. A genuine analogue of the spanning-tree proof would require a family with a tunable exact count: for trees, the missing 2 follows from the block floor, and the constructive converse is the separate identity `tau(C_n)=n`, not the floor inequality. See [THM-4986](../../../01-canon/theorems/THM-4986-the-7-21-omission-proof-has-an-exact-spanning-tree-twin-block-floor-crevasse-lemma-and-the-first-creator-contrast.md).

## 2. Exact multiplicative reverse

Ordinal composition with a directed triangle is a literal three-to-one projection of paths and yields `3H`. Consequently it suffices for full Chan–Pak to realize every positive odd number prime to 3 except 7, and the single value 63 (already known). For if the 3-free part is 7, only exponents 0 and 1 are holes and every higher exponent is an ordinal multiple of 63. This argument does **not** justify division of a particular witness for 63 by 3.

More generally, in the candidate monoid `M={positive odds} minus {7,21}`, the irreducibles are precisely:

* odd primes other than 7;
* `7p` for odd primes `p >= 5`, `p != 7`;
* `49, 63, 343`.

Proof: without a factor 7 any composite splits lawfully. Write a multiple of 7 as `7^e 3^a u` with `gcd(u,21)=1`. For `e>=2`, if any non-7 prime divides it, remove that prime, leaving an allowed multiple of 49. Pure powers of 7 are irreducible only at exponents 2 and 3. For `e=1` and `u=1`, 63 is the first allowed value and higher exponents can peel off 3. For `e=1` and `u>1`, peel off 3 if present, or a prime of `u` if `u` is composite; the sole remaining case is `7p`. The listed candidates have no legal split. Repeated factorization terminates by size. Thus after the already known exceptional atoms are discharged, the unbounded obligations may be taken to be the primes other than 7 and the indicated `7p` family. Every realizable atom requires a strong core; padding cannot replace it. The checker includes literal strong representatives for 35, 49, 63 and 343.

This is not saying that every strong tournament count is an atom. It says exactly what multiplication can do without a new strong witness.

## 3. Local triangle substitution is a different operation

Fix an old root `v` in a tournament `T`; replace it uniformly by a nonempty tournament `C` of order `m`. Let `a_k` count words on the exterior vertices plus `k` indistinguishable root slots, each separated from the next by a nonempty directed excursion. A word is legal when all exterior arcs used and all slot boundary arcs obey the original root orientations. Let `F_k(C)` be the number of ordered spanning path covers into `k` nonempty paths. The inherited uniform substitution identity (THM-3121, [Session 293](../session293_cycle_returns_20261113/control/MODULE_IDEAL.md)) is

    H(T[v←C]) = sum_k a_k F_k(C),      a_1=H(T).

This follows by cutting a path at its maximal C-runs, and the inverse glues the corresponding ordered cover into the slots. In particular

    H(T[v←C3]) = 3H(T)+6(a_2+a_3)
               = H(T[v←TT3])+2H(T).

Here `F(C3)=(3,6,6)` and `F(TT3)=(1,6,6)`: the two three-vertex modules are indistinguishable in every higher-run channel. In particular no uniform triangle substitution produces any 3-free atom. The extra terms are nonnegative; they measure the possibility of leaving the root module and coming back, not an error in an approximate ×3 formula.

### Exact cycle census for two slots

Let `S` be the SCC containing `v` and `P` the product of all other SCC counts. If `v` is a singleton SCC, every `a_k` for `k>1` vanishes. Otherwise, in every two-slot word the excursion stays inside `S` and closes a simple directed cycle `D` through `v` when the slots are identified. For every such `D`, partition the remaining vertices of `S` as ordered sets `U,W`; choose a Hamilton path of `U` ending in a vertex beating `v` and a Hamilton path of `W` starting at a vertex beaten by `v`. Empty factors contribute 1. This gives **all** two-slot words exactly once, after choosing paths in other components independently. In particular, taking `U` to be all remaining in-neighbors and `W` all remaining out-neighbors gives

    a_2 >= P · (# simple directed cycles through v in S).

This is an injection of cycles, not just a count inequality: the excursion recovers the cycle. Every vertex of a strong component lies in a triangle: absent any arc from its out-neighbors to its in-neighbors, no path can leave the former and return to the root. Hence pure ×3 by local triangle substitution holds **iff** the root was a singleton SCC. For `|S|=3` one has `a_2=P,a_3=0`; for `|S|>3` Camion supplies a spanning cycle in addition to the local triangle, so `a_2>=2P`.

The same splitting proves a quantitative bound for any `m>=2`:

    H(T[v←C]) >= H(T)H(C) + (2^m−2) P · (# cycles through v in S).

Indeed every nonempty proper ordered split of C contributes at least one pair of paths, so `F_2>=2^m−2`. Equality in this last inequality holds precisely for a transitive C or C=C3: any other nontransitive C has a proper directed triangle, and splitting it off contributes at least 3 instead of 1. In the displayed overall bound, equality in particular holds for these modules when S=C3. For C3 this makes the smallest possible positive normalized surplus `(H'-3H)/P=6`, attained precisely by a triangle SCC.

### What controls the highest slot?

The exact inherited support bound ([Session 197](../../../session197/prime-atoms/MINIMUM_TERNARY.md), with a simpler fill here) is one plus the matching number of the bipartite arc graph from out-neighbors of v to in-neighbors in S. Each excursion has such a crossing, and distinct excursions give disjoint crossings. Conversely use the matching edges as the excursions, placing all remaining in-neighbors in the prefix and out-neighbors in the suffix. Paths on these two sets exist and are always legal at the slot. This supplies every slot number up to the bound. Thus the underlying control is not the total count H but the number of disjoint return channels at a root.

## 4. A different ternary operation that is genuinely renewable

Replace v by a transitive chain of length m rather than a single triangle. Then `F_k=k! S(m,k)`, by assigning distinct clones to nonempty slots and ordering each slot transitively. If the matching number is 2, the exact formula is

    H_m = a_3 3^m + (a_2−3a_3) 2^m + (a_1−2a_2+3a_3).

A true `H_(m+1)=3H_m+2` clock requires both `a_2=3a_3` and `a_1=3a_3−1`; checking scalar agreement without the matching bound would not certify it. Session 197 gives actual amplitude 60, with the whole ray `60·3^m−1`, and proves minimality in this precise single-root cloning language. At that very same root, replacing by C3 instead of TT3 at m=3 gives `1619+2·179=1977`; the two meanings of 'multiply by three' predict different integers. The third meaning is an arithmetic divisor from a 3-cycle action or an excursion factorial (Sessions 293–295): it constrains observers even when there is no automorphism of the core.

**Open proof target:** synthesize a physically renewable recipe for the multiplicative atoms, or sufficient exact descent with endpoint guards. Multiplication is the easy reverse of condensation. Nothing here reverses the loss of interface data during strong deletion or turns the sparse returning ternary rays into full coverage.

### Matching lower envelope

Write `M_r` for the number of size-r matchings in the out-to-in bipartite graph within S. The constructive support argument is injective on the matching together with a permutation of its r short excursions: these are exactly the length-two blocks between slots and reveal the matched arcs and their order. Use a fixed Hamilton path on each leftover side. Consequently `a_(r+1) >= P·r!·M_r`. All other contributions in that run sector are positive; this is a coefficientwise lower envelope, not an identity replacing the exterior data by a matching polynomial. In particular it distinguishes a single feedback channel, multiple competing channels, and independent ordinal multiplication.

Reproduce the independent finite controls from this folder:

    python -B check.py
    python -O -B check.py

The two modes should match `check.out`. The strong small-order scan uses every labelled orientation, reachability in each direction, and subset path DP; the substitution test separately counts repeated-slot words by literal permutations.
