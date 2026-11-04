# Vertex insertion is an exact positive directed-cut circuit

Session 176. A proved structural result and an exact synthesis reduction, not a proof of unrestricted tournament-spectrum completeness. No historical-priority claim is made.

## 1. From endpoint information to a circuit

Let T be a tournament on V, n=|V|≥1, and H its number of directed Hamilton paths. Write E(S,u) and B(S,u) for Hamilton-path counts on T[S] ending and beginning at u, respectively; these are zero unless u belongs to S. For distinct u,v define

    K(u,v) = sum over S⊆V with u∈S and v∉S of E(S,u) B(V\S,v).

Thus K(u,v) counts ordered two-path covers, with first path ending at u and second path beginning at v. No edge condition is imposed between those endpoints. In particular K(u,v) is independent of the mutual edge uv. Put

    w_{uv} = [K(u,v)+K(v,u)]/2.

Add a vertex x, specifying b_u=1 for x→u and b_u=0 for u→x. Denote the resulting Hamilton count by F_T(b).

**Positive cut circuit theorem.** Every w_{uv} is a positive integer, and

    F_T(b) = H + 2 sum_{u→v in T} w_{uv} b_u(1−b_v).       (1)

This is an exact arithmetic circuit on Boolean inputs with one positive weighted directed-cut gate per unordered pair. It accounts for every orphan path as well as every legal insertion into an old Hamilton path.

### Proof of the identity

For an old Hamilton path P, legal positions for x number

    b_first + 1−b_last + number of 01 transitions
      = 1 + number of 10 transitions.

Consequently its excess over one is the number of old arcs u→v occurring consecutively with b_u=1,b_v=0. The total number of old paths using that arc is K(u,v): concatenate the two paths of a cover across the existing arc.

Every new path not obtained by legal insertion is an orphan. It contains consecutive v,x,u with u→v in T, b_u=1,b_v=0. Deleting x leaves an ordered two-path cover counted by K(v,u). Conversely each such cover produces one orphan. Summing these two disjoint contributions proves (1), initially with coefficient K(u,v)+K(v,u).

### Integrality

The standard tournament parity theorem states that every Hamilton-path count is odd. Fix all bits other than u,v arbitrarily. The mixed second difference

    F(1,1)−F(1,0)−F(0,1)+F(0,0)

is −[K(u,v)+K(v,u)] by the identity, and is even since its four terms are odd. Thus w_{uv} is integral. This use of tournament parity is explicit; no independent proof of that inherited theorem is claimed.

### Positivity

For n=2, both K values equal one. Otherwise choose any Hamilton path P of T[V\{u,v}], whose existence follows from the standard tournament Hamilton-path theorem. Regard the internal orientation of the consecutive block (u,v) as free. If every vertex of P beats both u and v, append the block. Otherwise let y be the first vertex of P that fails to beat both, and choose z∈{u,v} with z→y. Insert the block immediately before y, ending the block with z. If y is not first, its predecessor beats both endpoints, so the joining edge into the block exists. This gives a Hamilton path using one orientation of the adjacent pair, even if that orientation differs from T's actual uv edge. Removing the mutual edge gives a cover counted by K(u,v) or K(v,u). Their sum is positive. Its evenness then gives w_{uv}≥1. ∎

## 2. Strict submodularity and the complete zero set

For O={u:b_u=1}, equation (1) says F_T(O)=H+2 times the capacity of the directed cut O→V\O in T with capacities w.

For arbitrary A,B⊆V,

    F(A)+F(B)−F(A∪B)−F(A∩B)
      = 2 sum_{u∈A\B, v∈B\A} w_{uv}
      ≥ 2 |A\B| |B\A|.                              (2)

Hence the response is strictly submodular whenever A,B are not nested. This is a quantitative restriction on actual insertion interfaces, not a statement about arbitrary formal profiles.

Moreover F(O)=H if and only if no arc leaves O. A tournament's strongly connected components have a unique transitive order C1→C2→...→Ck. Such a set O must be a suffix union of these components: it cannot split a strong component, since a nontrivial subset of a strongly connected digraph has an outgoing edge, and cannot omit a later component while containing an earlier one. Conversely each suffix union has no outgoing arc. Thus:

**Zero-response classification.** Exactly k+1 orientations of the new vertex preserve H, one at each boundary of the ordered SCC decomposition (including the two outer boundaries). Every other insertion increases H by at least two. In particular a strongly connected base has exactly two count-preserving insertions, the universal source and universal sink.

This identifies precisely which information-flow configurations fail to create a positive feedback cut.

## 3. Exact target synthesis and resource costs

For a fixed base T and odd target m≥H, a one-vertex extension has count m if and only if the positive integer weighted tournament (T,w) has a directed cut of weight (m−H)/2. A witness cut gives all new edges explicitly. This is a complete equivalence for a fixed base, not an assertion that every cut weight or every odd target is possible.

If T is strong and n≥3, its smallest positive insertion increment is twice its minimum nontrivial directed cut. It can be obtained by minimizing over directed s–t minimum cuts for all distinct s,t. Every nontrivial set separates at least one ordered pair, and every such cut is nontrivial. This reduction is polynomial in n and the bit length of the already computed capacities. Constructing those capacities is not claimed polynomial in n: endpoint tables use O(n 2^n) integer cells and O(n² 2^n) arithmetic operations, and all K entries can be accumulated in O(n² 2^n) further operations. Once they exist, each response evaluation takes O(n²) arithmetic operations. Integer bit complexity is additional.

The direct circuit uses binom(n,2) two-literal multiplication gates, weighted addition, and the constant H. It is optimal **in the restricted model of a sum of gates, each depending on at most two input bits**, even allowing unary and constant gates: every unordered pair has a nonzero mixed coefficient −2w_{uv}, and a gate on any other support contributes zero to that coefficient. Thus at least one pair-support gate per pair is necessary. This is not a lower bound for unrestricted arithmetic or Boolean circuits; shared intermediate computations can invalidate that counting argument.

## 4. Quadratic tomography and information loss

There are no higher-order Boolean interactions. All 2^n extension counts are determined by the 1+n+binom(n,2) probes at 0, singletons, and pairs. Explicitly,

    2w_{uv}=F({u})+F({v})−F({u,v})−H,
    F(O)=H+sum_{u∈O}[F({u})−H]−2 sum_{{u,v}⊆O}w_{uv}.  (3)

The degree-two multilinear polynomial is unique, so two labelled base tournaments have identical complete insertion responses if and only if they have the same H, the same undirected pair capacities w, and the same weighted outgoing degree at every labelled vertex. Indeed F({u})−H is twice that weighted outgoing degree. Full insertion probing therefore need not reconstruct the original tournament.

The two opposite labelled directed 3-cycles are a literal information-loss witness. Both have H=3, every pair capacity one, and every weighted outgoing degree one. Their responses in mask order 0,...,7 are

    (3,5,5,5,5,5,5,3).

Thus this endpoint-derived response quotient discards an orientation bit even when the entire insertion truth table is known. One cannot safely replace the underlying tournament by this scalar-response quotient for later arbitrary modifications without proving that those modifications respect the quotient.

### A concrete failure of continuation closure

The opposite triangles can be separated by two successive, identically labelled insertions. Start with T having arcs 0→1→2→0 and with its reversal T*. Add x with outbits (0,0,1), meaning 0→x, 1→x, x→2. Both resulting tournaments have H=5. Next add y with outbits (0,1,0,0) against vertices (0,1,2,x). The resulting counts are 13 for the original cycle and 11 for its reversal.

Their outgoing-mask witnesses on (0,1,2,x,y) are respectively

    (26,12,17,20,2) and (28,9,18,20,2).

The independent `fiber_independent_audit.py` enumerates all 128 pairs of first/second insertion assignments and counts paths by literal permutations, recording this separation. Thus equality of the entire one-vertex response is demonstrably not preserved by an allowed next insertion. This is a specific obstruction to treating that quotient as a compositional state, stronger than merely warning that closure has not been proved.

### Exponentially many indistinguishable labelled bases

Let T be the ordered sum of k directed triangles on fixed labelled blocks, with each triangle independently assigned either cyclic orientation. All 2^k resulting labelled tournaments have the same entire one-vertex insertion response. In particular no number of one-vertex response queries can recover these k hidden orientation bits.

Here is an explicit proof, not an extrapolation from the three-vertex example. Every base has H=3^k. For two vertices in the same triangle, their pair capacity is 3^(k−1): every earlier triangle lies wholly in the first path of the defining cover and every later triangle wholly in the second, leaving the local pair capacity one.

For u in block i and v in block j>i, write t=j−i−1. A cover counted by K(u,v) cannot place any intervening block anywhere. Hence it is zero for t>0; for t=0 its value is 3^(k−2), since the two endpoint-refined triangle counts are each one and every exterior triangle contributes three.

For K(v,u), exterior triangles again contribute three each. The triangle containing u can be split between the two paths in three ways while requiring the second path to start at u: an empty first part contributes one, the unique usable singleton first part contributes one, and the two-vertex first part contributes one. The endpoint-v triangle similarly contributes three. Each intermediate triangle splits freely between two labelled paths, allowing empties, in twelve ways: the two unsplit choices contribute 2·3 and the six proper nonempty subsets each contribute one. Therefore

    K(v,u)=3^(k−2−t) · 3 · 3 · 12^t = 3^k · 4^t,
    w_uv=5·3^(k−2)                    if t=0,
    w_uv=(3^k·4^t)/2                  if t≥1.

These weights are independent of every triangle orientation. Every vertex has exactly one outgoing internal edge of capacity 3^(k−1), and all its other outgoing edges go to later blocks with the displayed invariant capacities. Its weighted outgoing degree is consequently invariant too. The response-equivalence criterion proves the theorem for all k≥1. The cross-block formulas apply only when k≥2. These are distinct labelled tournaments, not 2^k nonisomorphic tournaments: all triangle reversals are isomorphic. The conclusion concerns loss of labelled edge information in the insertion oracle, not a count of unlabeled isomorphism classes.

Let p2(T) denote the number of unordered covers of V by exactly two nonempty directed paths. Every ordered cover is counted by exactly one K(u,v), so

    sum_{unordered pairs} w_{uv}=p2(T),
    average_b F_T(b)=H+p2(T)/2.                       (4)

The average need not be integral. These are optional consistency checks, not a claim that the full path-cover profile is determined by H.

## 5. An explicit streaming benchmark

For the transitive base 1→2→...→n, with u<v and d=v−u,

    w_{uv}=1                    if d=1,
    w_{uv}=2^(d−2)              if d≥2.               (5)

To prove this, K(u,v)=1 for adjacent u,v and zero otherwise: a cover that ends its first increasing path at u and starts its second at v cannot place a vertex strictly between them. Conversely K(v,u)=2^(d−1): all vertices before u are forced into the first path, all vertices after v into the second, and each of the d−1 intermediate vertices can choose either path. This proves (5), including the adjacent case. Substitution into (1) recovers an alternative positive-cut expression for Session 175's binary feedback streaming constructor. The exponentially growing long-range gate capacities make explicit why a small number of nonuniform feedback bits does not behave like a free additive arithmetic circuit.

## 6. Relation to the spectrum frontier

The inputs here are actual tournaments and actual endpoint counts, not unverified numerical profiles. The mechanism strengthens Session 175's orphan accounting by combining its two terms into an exact positive integer cut circuit, classifies its entire minimum set, and supplies a complete fixed-base synthesis problem and a quadratic response quotient.

The global exclusions 7 and 21 remain exclusions. Unrestricted completeness for all other positive odd integers remains open; the strong-prime and strong-7p construction lanes from Session 175 are not closed here. Sequential insertion does not evade this obligation: the capacities must be recomputed after each new vertex, and the new pair capacities are not free parameters. A useful next construction target is a family of *realized* weighted tournament kernels whose directed-cut spectra contain provably overlapping intervals.

## 7. Reproduction and provenance

Run from any directory:

    python3 check.py > certificate.json
    python3 -O check.py > certificate.optimized.json
    cmp certificate.json certificate.optimized.json
    python3 fiber_check.py > fiber_certificate.json
    python3 -O fiber_check.py > fiber_certificate.optimized.json
    cmp fiber_certificate.json fiber_certificate.optimized.json

The standard-library script enumerates every labelled base tournament through order five and every orientation of one added vertex. It checks the cut identity against endpoint dynamic programming, positive-integral capacities, SCC zero-response counts, tomography, reconstruction of every truth table, and the average identity. For bases through order four, every extension is independently counted by literal permutations. Transitive kernels are checked through order ten. A separate checker tests all independent triangle orientations for one through four ordered blocks (up to twelve vertices), including the explicit exponential-fiber kernel formulas. Failures raise explicit exceptions, so optimized Python does not disable checks.

Inherited source: [Session 175 index](https://github.com/eliottcassidy2000/math/blob/897db58335ce4f7f5c332a238caf572e46537806/00-navigation/S175.md), especially its construction/README.md §5 and interface/PROOF.md §§2–3 under 05-knowledge/results/session175_tournament_spectrum_20261104/. The present proofs derive the circuit theorem from those endpoint/orphan definitions. This work inspected only targeted prior packets; it did not read or scan the repository's mistakes ledger.
