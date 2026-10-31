# Starting the positive natural numbers at the square-sum graph G32

Status: definitions and elementary proved consequences, with classical square-sum Hamiltonicity explicitly inherited. This is an enriched representation of ordinary arithmetic, not a new field or changed arithmetic laws.

## Three levels that must not be conflated

Let G_n have original vertex labels1,...,n, with an edge precisely when two distinct labels sum to a square. Distinguish a vertex label v, a graph order n, and a newly assigned numerical rank r. The ports1,3,8 are original vertex labels. They are not themselves the graph objects representing the new arithmetic units.

There are two consistent conventions:

- **Positive-origin convention:** represent r>=0 by G_(31+r). Then the first positive number1 is G32, while zero is the non-Hamiltonian G31. This makes the pre-origin structure concrete. Arithmetic transported to these labelled graph objects is G_a +_31 G_b = G_(a+b-31), G_a *_31 G_b = G_(31+(a-31)(b-31)), for a,b>=31. Zero is G31 and one is G32.
- **Zero-origin convention:** represent r>=0 by G_(32+r). Replace31 by32 in those formulas. Now zero is G32, one is G33, and all represented graph orders have Hamiltonian cycles.

In both cases successor adds exactly the next labelled vertex. The resulting arithmetic is isomorphic to the ordinary nonnegative integers by definition. The operations are not disjoint union, graph join, or port gluing. Under the positive-origin convention, the predicate “has a Hamiltonian cycle” agrees exactly with r>0; this uses the known all-order theorem below.

Simply relabeling vertices v as v-32 would instead change the adjacency equation to x+y+64 being a square. That is a different operation from shifting graph-order indices. We retain the original labels throughout.

## The zero boundary is a genuine structural transition

Gerbicz's classical result gives a Hamiltonian cycle for every G_n with n>=32; the prior Session124 replication supplies its pinned proof dependency and replay procedures. The absence below32 does not follow from downward monotonicity: Hamiltonicity is not inherited by deleting a vertex.

Here is an elementary exclusion for n<=30. Simple cycles require at least3 vertices. For3<=n<=7, vertex1 has degree1; for8<=n<=16, vertex8 has degree1; for n=17, vertex16 has degree1; for18<=n<=30, vertex18 has degree1. At31 every vertex can have degree at least2, but the exhaustive independent certificate in the sibling origin_pairs packet gives exactly two spanning2-factors, with component lengths11+20 and14+17. Neither is connected. At32 the unique Hamiltonian cycle is independently enumerated there.

Thus representing zero by G31 and one by G32 is not representing zero by an empty graph. It places the algebraic zero at a nontrivial, non-Hamiltonian boundary object. Graphs G1,...,G30 can be called a finite prehistory, but they do not by themselves constitute all negative integers or an additive group completion.

## The deterministic arrival rule

For arrival N>=1, its old neighbors are exactly

    {k^2-N : floor(sqrt(N)) < k <= floor(sqrt(2N-1))}.

Consequently its arrival degree is floor(sqrt(2N-1))-floor(sqrt(N)). Every neighbor lies in1,...,N-1; self-loops are excluded. This rule, not a Fibonacci recurrence, specifies the square-sum graph's growth. A numeration system can encode N without replacing this adjacency law.

For example,32 attaches to4,17. Arrival33 attaches to3,16,31. No pair among3,16,31 sums to a square, so no Hamiltonian cycle at32 can be extended by subdividing just one edge with33. The sibling certificate proves the optimal repair deletes two old edges. Graph inclusion, the existence of a new Hamiltonian cycle, and cheap update of a selected old cycle are different predicates.

Original nice-pair certificates add another predicate: paths of orders n and n+1 starting at1, ending at3 for odd length and8 for even length, with identical positional parity for every common label. They first occur at41/42, not32/33. The original vertex label parity is not the positional parity. Under positive-origin ranks r=n-31, order parity is opposite rank parity; under zero-origin ranks r=n-32, it is the same.

In fact **no original nice pair's two closed cycles are related by pure insertion**, at any order. In both paths, labels1 and3 occupy odd positions, because the odd member ends at3. Thus the even closed cycle does not contain edge1–3, while the odd closed cycle does. From even to odd, simple insertion cannot add this old-to-old edge. From odd to even, it would have to subdivide edge1–3, but a fresh common neighbor would make two squares differ by2, which is impossible. If exactly one old cycle edge is deleted on adding a vertex, its two incident new edges force precisely such an insertion. Hence at least two old edges must be deleted for every nice-pair cycle transition. The general pure insertion43→44 therefore belongs to a different witness relation.

## An alternative arithmetic structure: the order semigroup

Keeping ordinary addition of graph orders, rather than shifting it, gives a different exact object:

    S = {0} union {n>=32}.

This is precisely the set of Hamiltonian square-sum orders with an adjoined arithmetic zero. The element0 here denotes a number; it need not assert that the empty graph has a simple Hamiltonian cycle.

S is the ordinary numerical semigroup generated minimally by32,33,...,63. Each listed generator is indecomposable because two positive elements sum to at least64. Every n>=64 is obtained by subtracting32 repeatedly until a generator remains. Its conductor is32, Frobenius number31, and complement in the nonnegative integers has31 elements.

Its Apéry normal form relative to32 is unique: if n=32q+r,0<=r<32, then use q copies of32 when r=0, and (q-1) copies of32 plus32+r otherwise. Its Apéry set is {0,33,...,63}. This is a canonical residue decomposition, not uniqueness of arbitrary sums:66=32+34=33+33. It is also not a constructive operation combining two given square-sum Hamiltonian cycles.

## Where the Zeckendorf connection genuinely enters

Chung–Graham numeration gives a separate canonical address for any nonnegative rank or order. Its weights1,3,8,21,... are even-index Fibonacci numbers; digit constraints and carry rules are proved in the numeration packet. The representation applies equally to r or n, but those are different choices.

For example, encoding the original orders gives31=(1102)_CG and32=(1110)_CG. Encoding positive-origin ranks gives the same graph transition as0 to1. Arithmetic carry complexity is therefore coordinate-dependent; the actual graph change is not explained merely by observing a carry in one chosen encoding.

The useful enriched model retains a labelled graph and relevant realization data above each numeral: a cycle, endpoint ports, positional phase, or a boundary-forest state. Discarding these witnesses still leaves a correct integer address, but loses the information needed to certify future graph constructions. The successor-phase and rooted-forest packets quantify that loss rather than identifying the carriers by numerical resemblance.

## Dependencies and nonclaims

The graph origin and semigroup conclusions use the known theorem “cycles for every n>=32,” credited to Gerbicz, plus the independently replayed finite exclusions. The semiring identities and arrival formula are elementary. The Chung–Graham system and rooted-forest/Fibonacci connection are classical; no historical-priority claim is made.

No proof here identifies numeral normalization with Hamiltonian gluing, makes a new number field, proves arbitrary prescribed-path extension, or proves an all-order coherent nice-pair family.
