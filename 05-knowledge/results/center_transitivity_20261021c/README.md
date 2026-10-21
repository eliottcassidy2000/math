# Full vertex orbits, regular subgroups, and the order-25 exclusion

**Current follow-up:** [ORDER35.md](ORDER35.md) subsequently excludes 35 and
lowers the center-degree bound to 28. The next two paragraphs and list (1)
retain the explicitly order-25-only intermediate result; they are not the
current combined unresolved list. [VT15.md](VT15.md) separately excludes the
vertex-transitive part of the still-open order-15 lane.

October 21, 2026. **PROVED:** a deficit-six center of order 25 or 35 must be
vertex-transitive. **FINITE-EXACT, with a complete necessary universe:** no
order-25 tournament is such a center. This packet alone reduces the inherited
ten unresolved orders to the following nine:

\[
\boxed{12,14,15,18,20,21,24,28,35}. \tag{1}
\]

The order-35 coefficient remains **OPEN**, with vertex-transitivity a necessary
condition. This packet does not supply its center or metagraph wiring. It does
not incorporate a concurrent order-12 exclusion; a later synthesis may combine
that separately authenticated result with (1).

## 1. Inheritance and the information that must be retained

The input is a literal strong self-converse tournament \(T\). Write
\(A=\operatorname{Aut}(T)\), let \(H\) count its directed Hamilton words, and
put \(h=H/|A|\). A Hamilton word is *blue* when its positional reversal is an
anti-automorphism. Let \(b\) count blue \(A\)-orbits and
\(\Delta=h-b\). Thus \(\Delta\) is the number of black path orbits, not the
number of vertices or an unlabelled graph count. For self-converse sources the
converse/reversal operation pairs the black orbits, so \(\Delta\) is even.

Closest proved mechanism: the score-run obstruction in
[the incoming odd-shore proof, §8](../odd_shore_deficit_centers_20261021/session22-centers/PROOF.md).
Its half-supported prefix lemma, proper-shore deletion argument, and exact
odd-shore balance are in §§1–3. The current combined
[degree-35 proof](../deficit_six_centers_20261021b/COMBINED_CUTOFF.md)
retains the ten orders
\(12,14,15,18,20,21,24,25,28,35\), and asserts only regularity at 25 and 35.
The general closing-cycle input is
[cycle-cut obstruction, §5](../paired_black_leaf_reduction_20261020/CYCLE_CUT_NEAR_PURITY_OBSTRUCTION.md),
with its explicit closing words in the
[prefix classification](../pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md).

Corrected near miss: the nonautomorphic successor is an **existential choice
of Hamilton cycle**. It is not a property of every cycle of the source.
Canonical hostile: regularity alone does not imply vertex-transitivity; §6
supplies an exact seven-vertex example. The least-used sidecar here is the full
vertex-orbit partition of \(A\), together with the action of an
anti-isomorphism on that partition.

The live concepts are rotation shores, full vertex orbits, anti-isomorphism
normalization, regular prime-square subgroups, and full ordered adjacency.
Their connection contract is:

| Source | Map and target | Preserved predicate | Lost information / sidecar |
|---|---|---|---|
| A marked Hamilton cycle of \(T\) | Label each vertex by its full \(A\)-orbit | Automorphic words have the same orbit-label word | Labels alone lose edges; retain the actual four-cycle and suffix |
| A transitive order-25 action | Choose a regular subgroup and its group coordinates | Every source becomes a Cayley tournament in the enumerated universe | Group labels and connection mask must be retained; repetitions are not new classes |
| A Hamilton word on its source | Complete ordered adjacency key | Equality of keys is **iff** equality of full \(A\)-orbits | The key is not a metagraph adjacency rule; retain source and word for ports |

The present advance uses actual permutations and source words. A count of
seven black orbits is sufficient to exclude \(\Delta=6\), but supplies no
unspecified black-port connections.

## 2. Orbit runs upgrade regularity to vertex-transitivity

Fix a Hamilton cycle with successor \(R\), and put

\[
m=|A\cap\langle R\rangle|,\qquad q=n/m.
\]

The inherited cut lemma gives exactly \(q\) path orbits among the cycle cuts,
at most \(\gcd(2,q)\) of them blue. The \(q\) orbits of \(R^q\) have odd
size \(m\), induce regular strong tournaments when \(m\ge3\), and are called
rotation shores. Every full \(A\)-orbit is a union of shores.

**Orbit-run lemma.** Suppose \(T\) is strong and self-converse,
\(q\ge5\), and \(m\ge5\). If a full vertex orbit consists of at most four
shores, then

\[
\Delta\ge q-\gcd(2,q)+4. \tag{2}
\]

**Proof.** First let the vertex orbit \(O\) consist of one shore. Choose a
directed four-cycle \(S\) in that shore, using the inherited strong-tournament
pancyclicity input. Against any other odd shore, the exact number of arcs out
of \(S\) is \(4k\). It cannot give exactly two exits at every vertex of the
other shore, because that would say \(4k=2m\). Thus some terminal supplies at
least three exits or three entries. Since \(S\) is a proper subset of a
shore, its complement is strong. The inherited prefix construction therefore
supplies at least four black path orbits, counting their converse mates. Each
constructed word has a consecutive run of four vertices in \(O\). A marked
cycle cut has no adjacent vertices in \(O\).

Now let \(O\) contain \(c\) shores, \(2\le c\le4\). Choose \(S\) in one
shore and choose the terminal in a **different shore of the same orbit**.
The same balance calculation supplies three entries or exits. The constructed
word has a run of five vertices in \(O\), including that terminal. A cut word
has the repeating shore-label pattern of the cycle, which contains just
\(c\le4\) positions from \(O\) per period and at least one other position.
It has no such five-vertex run.

Both arguments survive the entry dual and the converse mates. Indeed an
anti-isomorphism \(\sigma\) normalizes \(A\):
\(\sigma A\sigma^{-1}=A\). Hence it sends a full vertex orbit to another
full vertex orbit of the same cardinality. The latter also contains exactly
\(c\) rotation shores, whether or not \(\sigma\) maps each individual
shore to a shore. Reversal preserves the existence of the consecutive run.
Each of the four resulting black orbits is therefore disjoint from every
marked-cycle cut orbit. Membership in a fixed full \(A\)-orbit is preserved
by every automorphism, so this is an actual orbit distinction. Adding the at
least \(q-\gcd(2,q)\) black cut orbits proves (2). ∎

Consequently, when \(\Delta\le6\), \(5\le q\le8\), and \(m\ge5\),
every full vertex orbit has at least five shores. There cannot be two of
them. **The source is vertex-transitive.** This does not use the false
implication “regular implies vertex-transitive.”

The inherited necessary parameter table has only
\((q,m)=(5,5)\) at order 25, and \((5,7),(7,5)\) at order 35.
The actual-cycle choice and the exact cut index are part of that table's
proof. Thus both orders satisfy the new transitivity condition.

## 3. Why the order-25 Cayley universe is complete

We prove a group lemma for any prime \(p\).

**Lemma.** Every faithful transitive permutation group of degree \(p^2\)
contains a regular subgroup isomorphic to \(C_{p^2}\) or \(C_p\times C_p\).

Let \(G_v\) be a point stabilizer. Choose a Sylow \(p\)-subgroup \(P\) of
\(G\) containing a Sylow \(p\)-subgroup of \(G_v\). Comparing the exact
\(p\)-parts of \(|G|=p^2|G_v|\) gives
\([P:P\cap G_v]=p^2\). Thus \(P\) is transitive.

Choose a central element \(z\in P\) of order \(p\). It is semiregular:
if it fixed one point, centrality and transitivity would make it fix every
point, contradicting faithfulness. Its orbits are \(p\) blocks of size
\(p\). The induced transitive \(p\)-group on these blocks has order \(p\),
because the \(p\)-part of \(p!\) is \(p\). Choose \(g\in P\) acting on
the blocks as a \(p\)-cycle. Its permutation order is \(p\) or \(p^2\).
If it is \(p^2\), it is a regular \(p^2\)-cycle. Otherwise \(z\) and
\(g\) are commuting independent order-\(p\) elements. The group
\(\langle z,g\rangle\) is regular: \(g\) moves between the blocks and
\(z\) acts freely inside each one. ∎

Apply the lemma to the full automorphism group of a vertex-transitive
25-vertex tournament. Choosing an origin identifies its vertices with one
of the two abelian groups. The source is then

\[
x\longrightarrow y\quad\Longleftrightarrow\quad y-x\in S,
\qquad S\mathbin{\dot\cup}(-S)=G\setminus\{0\}. \tag{3}
\]

Each group has twelve inverse pairs, so all its connection sets are encoded
by twelve bits. There are exactly \(2\cdot2^{12}=8192\) frames in this
complete necessary universe. It is **not** a list of nonisomorphic sources:
the two presentations already share sixteen literal labelled sources, and
additional isomorphism duplications are harmless.

Inversion is an anti-automorphism of every source in (3). No claim that an
arbitrary transitive degree-35 group contains a regular subgroup is made.
The prime-square argument does not apply at 35.

## 4. A seven-word certificate suffices

For an actual Hamilton word \(W=(w_0,\ldots,w_{24})\), retain all 300 bits

\[
K_T(W)_{ij}=\mathbf1[w_i\to w_j],\qquad 0\le i<j<25. \tag{4}
\]

Two words have equal keys iff the unique map taking one ordered word to the
other is an automorphism of the entire source. This is an equivalence, not a
sufficient score invariant. A word is blue iff its positional reversal is
anti on **every** pair; the consumer checks this directly.

For every frame, [the literal certificate bank](certificates.txt) contains
seven Hamilton words starting at zero. Their keys are pairwise distinct and
all seven words are black. Thus every frame has \(\Delta\ge7\), and since
the sources are self-converse it has \(\Delta\ge8\). In particular none
has deficit six. Combined with §§2–3, this excludes every order-25 center.

All 57,344 words are replayed by [check.py](check.py). It reconstructs both
group laws and every connection mask independently of [produce.cpp](produce.cpp),
checks exact permutations, all path edges, all blue tests, and all seven full
keys. Production verification does not search for paths, enumerate
automorphisms, infer a full Hamilton count, or use the discovery counters.

The deterministic capped DFS producer merely found these finite witnesses.
Its observed total was 643,594 nodes; the worst frame used 8,377 nodes and 32
completed words. The cap was ten million nodes per frame, with exhaustion an
error rather than an exclusion. No cap was exhausted. These are discovery
statistics, not completeness hypotheses.

## 5. Resulting frontier and what it does not establish

This packet proves the necessary transitivity restriction and excludes order
25. Applied to the frozen ten-order source table, it gives exactly (1).
The existing degree bound remains 35, and the existing recurrence order 18
and safe onset 48 are unchanged. There is no new order-35 coefficient, no
full center classification, and no metagraph reconstruction from \(h,b\) or
\(\Delta\) alone.

The most concrete next obligation is now a **vertex-transitive** order-35
source with an actual marked Hamilton cycle of type \((q,m)=(5,7)\) or
\((7,5)\). A cyclic Cayley search would cover only a stated subfamily until
a separate group-theoretic completeness argument is supplied. The
seven-word consumer can be reused on any correctly supplied source; its
current completeness theorem concerns order 25 only.

## 6. Controls and reproduction

Independent small controls enumerate all Hamilton words and all permutations
of \(C_3\), regular \(R_5\), and Paley7. They obtain respectively
\((H,|A|,h,b)=(3,3,1,1),(15,5,3,3),(189,21,9,3)\). Literal permutation
orbits independently agree with the complete ordered keys. Thus pure-blue
and deficit-six positive controls both remain present.

The regular nontransitive hostile has rows
\((14,84,56,50,97,67,13)\). Every vertex has outdegree three, but the number
of cyclic triangles in its outneighbourhood is
\((1,1,0,1,0,1,0)\), a nonconstant automorphism invariant. Twelve interface
hostiles reject bool/float aliases, wrong groups or masks, changed sources,
duplicate/negative vertex labels, duplicate keys, and wrong containers.

From the repository root:

```text
python -B 05-knowledge/results/center_transitivity_20261021c/check.py
python -B -O 05-knowledge/results/center_transitivity_20261021c/check.py --verify
```

Normal and optimized verification agree at **163,904 body checks**
(the verifier additionally compares the saved receipt). The receipt LF SHA256
is `1cd59dbfdc5ff5c296fc69aa7f5327a41991e237ba2c291d991697e2792ad0c8`.
The literal bank LF SHA256 is
`135cbc88608f3e6233cc82012984fda4fbd5b6a2b66ac7e1461191850965c151`.
The bank has 8,192 frames and 8,176 distinct labelled adjacency matrices;
these numbers are not isomorphism counts.

During development the first saved-receipt comparison compared Python tuples
with JSON lists. It was repaired before promotion by comparing canonical
serialized text, with LF-normalized reading. No mathematical payload or
certificate changed. The producer is optional for verification and its
executable is not included.

Independent full proof/API review: tiling_fold, October 21, 2026, **PASS**.
The parent independently replayed the normal consumer and reproduced all
163,904 checks and the saved receipt. The author separately replayed optimized
verification with identical output. The subsequent
[order-35 proof and complete certificate](ORDER35.md) extend this packet;
the nine-order list above is expressly the result of the order-25 step alone.
