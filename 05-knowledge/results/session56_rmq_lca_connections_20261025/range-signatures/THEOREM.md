# Overlap-safe aggregation of fixed-port collision signatures

25 October 2026. Status: **PROVED** elementary algebra and range-query correctness; **FINITE-EXACT** certificates. Classical RMQ/LCA algorithms are **CITED**, not claimed new. No theorem here settles LRC(14) or classifies arbitrary shared-pin networks.

## 1. Inheritance and operation audit

Closest mechanism: Session55 `crt-collision-wiring/THEOREM.md`, fixed row-major cell labels 0,...,8 and all 36 unordered cell pairs. A local array's mask records actual equal pairs, not partition block sizes. CRT over distinct coprime components intersects these relations. The least all-Parker-factor non-Parker product is 35=5*7. Its two local minimal masks are A=0x61089888 and B=0x8a054250 with A&B=0.

Canonical hostile: one F5 component has choices A or B but not their intersection. Duplicating that component as two independent choices fabricates a distinct array. Corrected near miss: bitwise AND of **chosen masks** is idempotent; convolution of **families of alternatives** need not be. Least-used sidecar: the local array multiplicities and the component identities. They distinguish one local choice from two independent choices even when supports agree.

Concept board: fixed-mask meet; alternative-family convolution; dominance compression; weighted histograms; range-query overlap; component provenance. The same minimum-looking syntax does not confer the same algebraic law.

## 2. Family and antichain idempotence: exact characterizations

Let U be a finite interface and let masks be subsets of U. For families F,G define

    F * G = {x intersection y : x in F, y in G}.

Empty families are permitted. This is associative and commutative, with identity {U} and absorbing empty family: expand either parenthesization to all triples and use associativity of intersection.

**Theorem 1 (PROVED).** F*F=F if and only if F is closed under pairwise intersection.

Proof. The diagonal choices x=x give F subset F*F. Therefore equality is equivalent to the reverse inclusion, exactly pairwise meet-closure. The statement also holds for the empty family. Notice that diagonal inclusion alone is not enough.

For a finite family F write Min(F) for its inclusion-minimal elements and define the reduced product A o B=Min(A*B) on antichains.

**Theorem 2 (PROVED).** Min(F*G)=Min(Min(F)*Min(G)). Consequently o is associative, with identity {U}. An antichain A satisfies A o A=A if and only if A is empty or a singleton.

Proof of compression. Every x in F has some x' in Min(F) below it, and likewise y' below y; hence x' intersection y' is below x intersection y. The reduced product is a subfamily of the full product with an element below every full-product element. The two families therefore have exactly the same minimal elements. Applying this identity to triple products proves associativity.

Proof of the characterization. Empty and singleton cases follow immediately. If a,b are distinct elements of A, they are incomparable; c=a intersection b is strictly below both. Choose a minimal element d of A*A below c (finiteness guarantees one). Then d belongs to Min(A*A) but cannot belong to A: otherwise d is strictly below a in the antichain. Thus the reduced product differs from A. In particular it is unnecessary to claim c itself is minimal; another new mask may dominate it.

For a finite nonempty family F, repeated convolution saturates to the nonempty finite-meet closure. If a full pass adds no masks, it is idempotent; otherwise any added mask witnesses the exact invalid independent-choice relaxation. This observation does not authorize replacing a physical component by its closure.

## 3. Weighted alternatives and compatibility

For a histogram h:2^U -> nonnegative integers, define

    (h * g)(z) = sum over x intersection y=z of h(x) g(y).

This counts independent ordered local choices with their multiplicities. It is associative with identity delta_U. The total mass is multiplicative. Thus h*h=h forces total mass W to obey W^2=W, so W is 0 or 1. Conversely the zero histogram and every single unit point mass are idempotent. This characterizes histogram idempotence over nonnegative integers; it is NOT a statement about signed coefficients or probability normalizations.

Dominance compression preserves minimal attainable collision masks and zero attainability, not counts. Never apply the antichain reduction to a counting request and claim its weights remain exact.

Each implemented signature retains component identities, the full histogram, and one list of actual local arrays per output mask. A merger rejects repeated identities. For a CRT interpretation the component moduli must additionally be pairwise coprime, verified before querying. Distinct identities alone do not establish CRT independence.

Independent factors F5 x F7 are legitimate. Independent copies F5 x F5 form a different product ring, not F25 and not Z/25Z. In particular nonzero nilpotents distinguish Z/25Z from a product of fields, while nontrivial idempotents distinguish F5 x F5 from the field F25. No extension-field computation is claimed here. Our exact Z/25Z control is Parker although independent copies of the F5 alternative family admit mask zero.

For shared-pin constraints not known to factor, the primitive combination must instead retain and join compatible assignments on the overlap; the same physical variable must agree, not be selected independently. Associativity of such a richer operation requires its own interface/forgetting proof. This packet proves only the independent CRT case and the explicitly defined abstract monoids.

## 4. Why overlap works, and the exact repair

The ordinary overlapping sparse-table formula covers [l,r] by a left and right interval of length 2^floor(log2(r-l+1)). For a nonempty overlap write the true aggregate as P*Q*R and the returned value as P*Q*Q*R. Associativity plus idempotence removes the duplicated Q. This argument does not require commutativity. If there is no overlap the formula is already valid; if the queried value happens to absorb the extra Q it can also happen to succeed. We do not claim every nonidempotent query fails. Uniform singleton duplication detects the necessity of idempotence for that unmodified formula; special-casing singleton/power-of-two lengths does not repair the length-three hostile below.

The interface operation matters: fixed chosen collision masks with AND are idempotent, as are minimum values. Alternative families, reduced antichains, and weighted histograms generally are not. Meet-closed full families may use overlap for support-only queries, but support idempotence does not preserve physical witness provenance or weights.

**Disjoint sparse table (PROVED algorithm).** At level k partition indices into consecutive blocks of length 2^(k+1), with midpoint at the half-block boundary. Store each left-half suffix and each right-half prefix in the original order. For l<r, k=highest-set-bit(l XOR r). Then l and r lie in opposite halves of the same level-k block. The two stored aggregates cover [l,r] exactly once and in order. Their product is the exact interval result by associativity. For l=r return the leaf unchanged. This proves correctness even for noncommutative semigroups and makes no independence claim beyond that of the supplied operation.

**Segment tree (PROVED algorithm).** Store disjoint dyadic subinterval aggregates. Decompose the query into disjoint canonical nodes; accumulate left nodes from the left and right nodes from the right. Concatenate the two aggregates in that order. Identity padding supplies empty ranges; associativity proves correctness. The code implements queries, not point updates; the standard update extension recomputes the O(log n) ancestor nodes.

These are standard data-structure constructions, with proofs included to audit the application rather than assert algorithmic novelty. The Bender--Farach-Colton paper supplies the classical RMQ/LCA motivation and overlapping RMQ sparse table. This packet does NOT attribute the disjoint construction or the family theorems to that paper.

## 5. Costs, including the interface invoice

For n input values, disjoint-table preprocessing uses O(n log n) semigroup operations and stored aggregate cells, with at most one operation per non-singleton query. A segment tree uses O(n) monoid operations and cells to build and O(log n) operations per query or point-update extension. Integer level selection is O(1) only under the usual word-RAM indexing assumption; one can instead tabulate floor logarithms.

Those are operation counts, NOT constant-time signature operations. Let inputs to a histogram product have a,b supported masks on m interface bits. Straight convolution uses a*b ANDs and a*b count multiply/adds, with at most min(2^m,a*b) output masks. Bit work for masks includes ceil(m/w) machine words; count multiplication/addition depends on coefficient bit lengths. Naive antichain reduction on at most a*b candidates uses O((a*b)^2) subset comparisons, each with its mask-word cost. Local enumeration is additional work, not included in table construction counts. m=36 in the CRT experiment; no claim bounds all reachable signature sizes by a small constant for varying interfaces.

The reference implementation deliberately keeps eager identity tuples and eager realizing-array tuples. Merging intervals of total length L adds O(L) identity checking/copying and up to O(a*b*L) tuple-concatenation work, beyond histogram arithmetic. Storing k witnesses of L local nine-cell arrays costs O(k*L) references (or entries if materialized). Thus O(n log n) aggregate CELLS is not an O(n log n) byte-space claim: eager provenance can reach O(n^2*K) references for a disjoint table with maximum support K. DAG witnesses and interval provenance can reduce copying, but are not implemented. Returning a concrete CRT witness already costs at least its output size. This avoids hiding a large signature/provenance cost behind an RMQ bound.

## 6. A genuine pairwise-coprime length-three hostile

Use the ordered components [2,5,3]. The true quotient ring is Z/30Z. Exhaustive local/CRT/direct checks give 112 magic square-residue arrays, no mask zero, and minimal masks

    0x1000880, 0x80040200.

The overlapping query combines the [2,5] aggregate with the [5,3] aggregate. It independently selects A for the first copy of the F5 component and B for the second; A&B=0, so it fabricates a zero collision mask. Its total count is 784, also wrong. The erroneous abstract product corresponds to four independent factors with total size 150, not the three-factor ring being queried. The identity-guarded merger rejects it before it could be mistaken for a valid CRT result.

Positive controls: fixed-mask AND supports every overlapping query tested; [5,7] is a valid distinct-cell product; all 31 subintervals of five specified pairwise-coprime component sequences agree between disjoint table, segment tree, direct fold, and direct full-modulus histogram enumeration. Their complete supports also agree with the separately reconstructed line-equation enumeration inherited from Session55. Every retained witness reconstructs an admissible CRT array with its claimed exact collision mask.

## 7. Connection ledger and limits

Source: classical overlapping RMQ tables and Session55 labelled CRT collision signatures.
Target: exact interval queries over independent component signatures.
Map: ordered component interval -> histogram/family of intersections of fixed-cell equality masks.
Preserved: full attainable collision support and array multiplicities under full histograms; minimal support and distinctness existence under antichains.
Destroyed by shortcut: component choice correlation when overlap is treated as independence; multiplicity under dominance compression; cell alignment under block-size-only summaries.
Needed sidecar: component identity, pairwise-coprime moduli, fixed labels, and actual realizing arrays.
Cheapest decisive test: [2,5,3] versus duplicated F5, independently checked modulo 30.

Applied method cards: “Correct the object before sharpening the technique” and “Retain local incidence through counting.” No new generic meta-pattern promotion is requested from this single lane. The lesson is an operation/quantifier audit, not a new RMQ asymptotic bound or a universal shared-port gluing theorem.
