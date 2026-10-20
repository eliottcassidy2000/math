# Closed half-edges, Hamiltonian bypasses, and the nonstrong boundary

Status: PROVED deductions from the cited tournament bypass theorem and inherited tiling definitions; the order-three to order-five classification is FINITE-EXACT with direct exhaustive certificates. No all-order connectivity or complete arrangement theorem is claimed.

## 1. Definitions and provenance

Let n>=3. A tiling is a tournament with one distinguished directed Hamiltonian path, normalized to a fixed vertex order. Complement tau reverses every nonpath arc. Reflection sigma reverses path positions and represents converse after relabelling. A literal line is a two-element tau orbit; it is blue when sigma fixes its two tilings, and black otherwise. A merged node is an isomorphism class identified with its converse class. Loops and repeated lines remain literal instances, each with two half-edges.

H(T) counts directed Hamiltonian vertex words. HC(T) counts directed Hamiltonian cycles up to rotation. The ordinary class fibre is H(T)/|Aut(T)|. All statements below count these exact objects, not labelled cycles modulo a different equivalence.

[THM-4729, closed/open tilings](../../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md) introduced the explicit closure inventory in the concurrent source at commit466113f5c06. Its nonstrong independent-set and first-cut results independently overlap the [strong-core packet](../merged_metagraph_purity_and_strong_core_20261020/STRONG_CORE_COUNT_RECURSIONS.md). This note credits that closure coordinate and adds the bypass classification, directed/cut consequences, and corrected leading asymptotic.

Call a path closed if its last vertex beats its first; otherwise call it open. An open Hamiltonian path together with the shortcut from its first vertex to its last is exactly a Hamiltonian bypass, obtained by reversing one edge of a directed Hamiltonian cycle. This is a predicate of the path presentation, not merely the node.

## 2. Closure inventory and all-shore constraints

### Theorem 1: exact inventory

For an ordinary tournament class,

    c(T)=n HC(T)/|Aut(T)|,
    o(T)=(H(T)-n HC(T))/|Aut(T)|.

A merged node containing two converse classes has twice each count. For every literal line, one half-edge is closed and the other open. In particular the total closed and open tiling masses both equal half the tiling cube; the same equality holds separately in blue and black.

Proof. Closing Hamiltonian paths are rooted Hamiltonian cycles, so there are n HC(T) of them. An automorphism fixing a Hamiltonian word fixes every vertex, so the action is free. Dividing its two invariant path subsets by |Aut(T)| gives the formulas. Converse preserves both subsets and their counts. Complement flips the nonpath endpoint arc, proving the line assertion. Reflection permutes the nonpath coordinates, so it commutes with complement. Thus complement preserves the sigma-fixed tilings and their complement, proving the assertion within each colour. This is THM-4729's closure mechanism. QED.

### Theorem 2: directed and cut constraints

Orient each line from its open to its closed end. A loop is a directed loop. At a node u the indegree is c_u and the outdegree is o_u, with each loop counted once in each. For a shore S of nodes and any one of the total, blue, or black line systems, write C(S), O(S) for its closed/open half-edge masses, I(S) for internal line instances including loops, and d(S) for crossing line instances. Then

    incoming(S)-outgoing(S)=C(S)-O(S),
    I(S)<=min(C(S),O(S)),
    d(S)=C(S)+O(S)-2 I(S)>=|C(S)-O(S)|,
    d(S)=C(S)-O(S) modulo2.

Proof. Each internal line consumes one closed and one open half-edge, including a loop. Each incoming crossing line contributes one closed half-edge, while each outgoing crossing line contributes one open half-edge. Thus C=I+incoming and O=I+outgoing. Subtracting and adding give every assertion. QED.

This is a bipartition of half-edges, NOT of merged nodes: one strong node generally carries both sorts. The constraints need not determine the pairing. For a nonstrong shore C=0, so d=O and I=0, recovering the nonstrong independent periphery. Camion's theorem gives c_u>0 exactly at strong nodes, hence every line's closed endpoint belongs to the strong core.

## 3. Fully closed classes are completely classified

### Cited input, with exact scope

Abdelhamid Benhocine and A. Pawel Wojda, *On the existence of specified cycles in a tournament*, Journal of Graph Theory 7(4) (1983), 469–473, [DOI and primary publisher abstract](https://doi.org/10.1002/jgt.3190070413): every tournament of order n>=7 contains D(n,p) for each 2<=p<=n, where D(n,p) is obtained from a directed n-cycle by reversing p-1 consecutive edges. Only p=2 is used here, giving a Hamiltonian bypass. This is the already source-audited published input used in the [independent rigid-strong proof](../merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md#3-rigid-strong-tournaments-are-completely-classified). We cite it, not a broader digraph exception classification, and do not claim to reprove it.

### Theorem 3: the only fully closed classes at n>=3

Every Hamiltonian path closes if and only if the tournament is C3 or the nonregular five-vertex tournament

    Q5 = a -> C3 -> b -> a.

Equivalently, these are exactly the classes with H(T)=n HC(T). In Q5, a beats all three cycle vertices, those three beat b, and b beats a.

Proof. For n>=7, the cited bypass theorem supplies an open Hamiltonian path, so equality fails. For every even n, Redei's parity theorem makes H(T) odd, whereas n HC(T) is even; equality again fails. The only remaining orders are3 and5. Direct exhaustive enumeration over all 2^binom(n,2) labelled tournaments and all n! vertex words, followed by full-permutation isomorphism canonicalization, yields the following exact classification:

- n=3: one fully closed class, C3, with H=3, HC=1, |Aut|=3 and scores(1,1,1).
- n=5: one fully closed class, Q5, with H=15, HC=3, |Aut|=3 and scores(1,2,2,2,3).

The portable verifier includes these full finite universes, explicitly checks existence of Hamiltonian paths, computes the closing-path count directly, and canonicalizes only the surviving classes. An independent implementation checks the same cases. Thus the finite base classification plus the all-order exclusion proves the theorem. QED.

The elementary shape of the positive Q5 is also transparent: starting at a, a closing Hamiltonian path must proceed through the cycle vertices in cyclic order and then b; there are three choices. These are all three Hamiltonian cycles through a. Their five cuts give 15 closing paths, equal to the direct full path count. C3 has its three cyclic words. The regular cyclic C5 instead has H=15, HC=2 and five open paths, so it is not this exception. No use is made of calling C5 a Paley tournament.

### Corollary 4: no sinks and directed cycles at every n>=6

For every strong T of order n>=6,

    H(T)>n HC(T)>0.

Thus its merged node has positive incoming and outgoing degree in the full open-to-closed orientation. All nonstrong nodes have indegree zero. Every outgoing line from a strong node ends at another strong node, so the induced directed strong core has positive outdegree at every node. It follows that every strong node reaches a directed cycle, possibly a loop.

Proof. Theorem3 excludes fully closed classes at n>=6; Camion supplies a cycle in every strong tournament. Every line's closed endpoint is strong. Repeatedly following outgoing edges in the finite strong core eventually repeats a vertex. QED.

There is at least one strong tournament at every n>=3: orient a Hamiltonian cycle and choose the other arcs arbitrarily. Thus the core is nonempty and the cycle assertion is nonvacuous. Positive indegree in the whole metagraph does not assert positive indegree inside the strong induced subgraph. Nor does the cycle conclusion imply strong connectivity, underlying connectivity, or any all-order blue-graph connectivity. This is the endpoint-closing orientation, not the cyclic-triangle-gradient orientation.

## 4. Correct leading asymptotic of nonstrong tilings

Set T_0=1 and T_n=2^binom(n-1,2) for n>=1. Let a(n) count strong fixed-path tilings, including a(1)=1, and NS(n)=T_n-a(n). The unique first strong component of a path presentation gives the exact recurrence

    NS(n)=sum_(k=1)^(n-1) a(k) T_(n-k).

The cross-component arcs are forced; the first block is strong and the remaining block is arbitrary. Both our strong-core packet and concurrent THM-4729 prove this bijection.

### Theorem 5: factor2 and an explicit uniform error bound

For every n>=4,

    2-(n-2)*2^(-(n-3)) <= NS(n)/T_(n-1)
                               <= 2+(n-3)*2^(-(n-3)).

Consequently

    NS(n) ~ 2*2^binom(n-2,2),
    NS(n)=2*T_(n-1)*(1+O(n*2^-n)).

Proof. First, the endpoint terms k=1 and k=n-1 are exactly T_(n-1) and a(n-1). For 1<=k<=n-1 an elementary exponent identity gives

    T_k*T_(n-k)/T_(n-1)=2^(-(k-1)*(n-k-1))<=1.

Applying this with n-1 in place of n shows

    NS(n-1)<=sum_(k=1)^(n-2) T_k*T_(n-1-k)
           <=(n-2)*T_(n-2),
    a(n-1)/T_(n-1)>=1-(n-2)*2^(-(n-3)).

Dropping the nonnegative middle terms proves the lower bound. For 2<=k<=n-2,

    (k-1)*(n-k-1)>=(n-3),
    a(k)*T_(n-k)/T_(n-1)<=2^(-(n-3)).

There are n-3 such middle terms. Since a(n-1)<=T_(n-1), adding their bound to the two endpoint terms proves the upper bound. Both errors tend to zero, establishing the asymptotic. QED.

The endpoint contributions are asymptotically equal. Thus the coefficient1 suggested in HYP-9345(f) is false; the fact that k=1 is the largest individual displayed term does not make k=n-1 negligible. The recurrence gives NS(n)/T_(n-1)=1.972900390625 at n=8, approximately 1.9998779894098946 at 16, and approximately 1.999999999992724 at 40. These are checks, not the proof.

## 5. Boundaries and reproduction

This resolves the closure question in HYP-9345(a) and the leading asymptotic in (f), not all remaining questions there. All-order connectivity, blue self-line uniqueness/near-regularity, self-line count formulas, full metagraph determination, and higher asymptotic terms are not proved here. The separate group-theoretic ML statement is not proved by removing its necessity for rigid classification.

Run `python3 verify.py` from this directory. Its JSON output must match `verification.json`; running under `python3 -O` gives identical output. It requires only Python's standard library and no repository paths or remote access. See `VALIDATION.md` for independent audit and provenance.
