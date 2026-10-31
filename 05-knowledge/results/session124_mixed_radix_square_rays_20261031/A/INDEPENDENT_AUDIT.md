# Independent audit: basic scale25 compatible prefixes and terminal cycles

## Finding

**PASS.** The frozen basic50 paired certificates and50 direct terminal certificates prove the following, using only B bases42 through1061:

Let `S_l=(85*25^l-1)/2`, and let

`P_l = ((25^t+1)/2,(7*25^t-1)/2) for t=0,...,l`, concatenated.

Then the deterministic new B constructor has common prefix P_l for every input index n>=S_l. For every l>=1 and every N>=S_l, the interval square-sum graph has a Hamiltonian cycle containing the **original unscaled P_l** consecutively. P_l has2l+2 vertices. Every finite prefix of this explicit infinite ray therefore extends to a Hamiltonian cycle in every sufficiently large interval.

Thresholds begin42,1062,26562,664062,16601562. This is a separate prescribed family from the scale49 four-prefix theorem, not a claim that the new result dominates all its retained-path conclusions. Ordinary square-sum Hamiltonicity was already known; no new proof of that previously solved conjecture or historical-priority claim is made here.

## Independent replay

`audit_independent.py` is a self-contained standard-library program, importing no producer verifier, constructor, search code, sibling audit module, or old transition table. It accepts explicit `--pair-dir`, `--cycle-dir`, `--small-base-file`, and `--base-dir` paths.

Normal and optimized runs produce byte-identical receipts with **51,027,430 explicit checks**:

- all50 paired states and all50 terminal states, including complete physical filename/state dispatch;
- all1020 B bases42 through1061, with full vertices, square edges, endpoints and common positional parity;
- all eleven small pairs42 through52 and the four old compressed chunks intersecting53 through1061;
- the frozen old base-manifest hash and full compressed/decoded hashes of every consumed chunk;
-100 full paired expansions, covering all50 states at two successive input depths;
-103 complete streamed cycles, including all50 terminal states at two depths and the depth-four threshold16,601,562;
- eight hostile paired/direct certificate corruptions rejected.

Every one of the103 independent cycle hashes agrees with the producer's separately generated receipt. That comparison is recorded separately; it is not an input to the independent checks. The largest cycle's decimal-lines SHA-256 is `4a3b1fa5c2bc9b9273cddbab8f2b40fad2c282af3eaa7104ecc273b5ead78ca5`.

The depth-four prefix is

`1,3,13,87,313,2187,7813,54687,195313,1367187`.

## Universal affine semantics

For offsets c in[-12,12], lift input position j by `25*v+(-1)^j*c`. Internal edge sums multiply by25. Compatibility of the input B pair ensures that a common label has the same affine representative in either selected member. All25 offsets partition its integer band `[25v-12,25v+12]`; those bands are disjoint and contiguous.

The extra label n+1 has positional parity n+1, by subtracting common odd/even position counts. Its representative belongs below the requested output bound25n+d precisely when

`25+(-1)^(n+1)*c <= d`.

That is the exact member selector checked for every atom. All low singleton labels1 through12 then complete the interval. This proof requires only the weak B invariant and valid positive member lengths; it has no hidden lower bound53. Inputs42 and43 were checked explicitly.

A split replaces the first vertex of one block by the head25-c and leaves the tail beginning75+c. Heads and tails partition the selected block exactly. The auditor checks every required atom exactly once, every forced selection and orientation, every numerical connector square, and the symbolic total length25n+d (or25n+d+1 for the second paired word). The words are complete traversals, not disconnected port matchings.

For paired outputs, it independently tracks each singleton's output parity, each head's output parity, and each tail/whole block's output-minus-input position shift. All common-label profiles agree across the two words. The frozen basic table also satisfies its stronger within-word head/tail shift-linkage convention, which is independently checked. These are universal position identities, not sampled-length approximations.

Every paired word begins with the exact atoms `(single1,single3,forward intact block12)`. Thus it preserves the weak B prefix1,3 and transforms any longer input prefix P into `(1,3),D12(P)`. The numerical connector is3+13=16. Paired endpoints remain8 for even order and10 for odd order.

The direct consumer has one word, the same protected opening, and final endpoint8 for either output parity. It closes to1 by9. It has **no output-pair parity or head/tail shift obligation**, because it is terminal and never recycled as a B pair. Its input still must be compatible for the interval-ownership argument.

## Correct base cutoff and all-size induction

The exact base interval is42 through1061. Define B(n) by returning a base in that range and using the new paired table for every n>=1062. For the recursive branch,

`m=floor((n-12)/25)`, `d=n-25m`,

so42<=m<n and12<=d<=36. Strong induction proves all-size weak B supply without any older transition family.

Since `S_(l+1)=25*S_l+12`, any n>=S_(l+1) has m>=S_l. Applying the protected opening proves that the same fixed constructor has prefix P_(l+1). The terminal table gives the cycle theorem with exactly one additional layer.

The cutoff1062 is essential. Returning old weak bases for n>=1062 could lose the newly required longer prefix and would invalidate the improved depth induction. The independent constructor enforces the new cutoff, and the controls exercise both sides of the second threshold: N26561 and N26562.

## Arithmetic ray identity

Put Q=25^t. The two-coordinate block is `((Q+1)/2,(7Q-1)/2)`. Its internal sum is4Q; its boundary sum with the next block is16Q. All are squares, with roots2*5^t and4*5^t. The blocks are strictly increasing and pairwise disjoint. Their affine recurrence is exactly `(1,3),D12(P)`, so the algebraic ray and the certified prefixes coincide.

A more general arithmetic classification, including Pell-type extreme-offset families, is in the separate `../affine_rays/AFFINE_RAY_LEMMAS.md`. That synthesis is not needed as a certificate dependency here.

## Scope exclusions

This audit concerns only the frozen basic two-prefix, protected-offset12 scale25 catalogs. It makes no Hamiltonian-extension claim for the separately explored four-prefix offset6 or offset8 constructions, nor for relaxed-head candidate catalogs. It does not transfer the scale49 simultaneous-copy counts or482-based seed to these new cycles. No arbitrary prescribed-path extension theorem is inferred.

## Artifacts

- `audit_independent.py`
- `independent.normal.json`, `independent.optimized.json`
- `producer_comparison.json`

The script's base-manifest digest is pinned directly, so its verification requires no prior package's Python code or scale49 transition table. Certificate hashes and the small-base catalog hash are recorded in both receipts.
