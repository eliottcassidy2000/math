# A new marked source closes the deficit-nine line; exact rooted port-order extension

**Status: PROVED scoped infinite constructions and exact finite certificates; same-agent proof review with separate adjacency and edge-walk consumers. Universal raw tournament spectrum completeness remains OPEN.**

The inherited source is packet f at `b9aa14f16edd`: closed subcubic square-sum carriers, a complete edge compiler, fixed-affine fan rigidity, and the bounded-loss theorem. The decisive change here is to enlarge the **parent tournament**, keeping the same marked-state interface. The previously closed singleton/D alphabet is not a universal class of parents.

## Inheritance / concept board

| Lane | Closest mechanism | Hostile or missing information | Action / result |
|---|---|---|---|
| Anchor: tournament spectrum | Native state (F,q,g,h), positive cut/interior weights | Fixed singleton/D sources have positive dyadic intercept; bounded one-apex losses cannot supply deficit9 indefinitely | Search actual arbitrary parents for q a power of two and negative intercept |
| Anchor: square-sum ports | Primitive disjoint path switches | Coarse codes tie switches; an unoriented cycle has no rooted port order | Normalize at the marked edge, then retain region and bit ownership |
| Niche: exact source census | Every tournament has a Hamiltonian observer order | Labelled enumeration repeats classes but must retain every marked cut | Fix the observer path and enumerate all remaining arcs and cuts |
| Wildcard: weighted queries | Hamilton cycles of the refined carrier form a Boolean cube | Separately feasible edge and order requests can conflict on one bit | Joint compiler, exact Hamilton polynomial, and signed-cost optimizer |

The method cards used are **Close semantic reading classes before judging a theorem**, **Verify the consequence, not the model's own assumptions**, and **Use redundant paths as detectors**. The finite census certifies sources; the recurrence, not a long list of examples, proves each infinite family.

## 1. The actual source that closes deficit9

Use vertices0,...,6 for a parent tournament A. Its outgoing-neighbour bit masks are

    [12,121,2,100,45,69,21].

Add a marked vertex x=7 with outgoing set S={0,1,6}, mask67. Every other parent vertex points to x. The exact marked state is

    (F,q,g,h)=(449,512,54,117).

Here h counts Hamilton paths of A, F counts those of A+x, g counts A-paths whose last vertex belongs to S, and q counts ordered two-path covers(P,Q) of A with Q nonempty and its first vertex in S. P may be empty. These are raw labelled path counts; no automorphism division occurs.

Append a vertex z dominated by every existing vertex other than x, with z->x. This is the native zero operation:

    (F,q,g,h) -> (F+q+g,2q,0,h).                 (1)

For r>=1 successive zero vertices, (1) gives

    H_r=2^r q+F+g-q=2^(r+9)-9.                 (2)

The restriction r>=1 is essential: the unextended marked source has F=449, whereas2^9-9=503; the initial g=54 has not yet been paid into F and reset. The resulting tournament has8+r vertices. Therefore for every k>=10 it has exactly2^k-9 Hamilton paths on **k-1 vertices**.

**Complete positive-line theorem.** For every integer k>=5, the number2^k-9 is realized as a raw tournament Hamilton-path count. Among positive values on this line, the only excluded exponent is k=4, which gives the already proved forbidden value7.

The finite boundary k=5,...,9 is given by actual adjacency witnesses in `deficit_nine_family.json`. At k=8 the target247 is an ordinal join of witnesses for13 and19; Hamilton paths of an ordinal join concatenate, so counts multiply. The remaining boundaries come from the inherited exact inverse. They are independently recounted, rather than inferred from an old coverage headline.

### A six-arc description and a guarded scalar return

There is a particularly short description with no transfer-state notation. On the ordered vertices1,...,n, n>=9, orient every edge forward except

    7->1, 5->2, 6->3, n->1, n->5, n->6.

This tournament has exactly **2^(n+1)-9** Hamilton paths. The parent relabelling is(1,4,3,5,6,0,2) in the original zero-based labels; all tail vertices follow it, and the marked vertex is last. `deficit_nine_family.py` independently checks this direct six-arc adjacency against the source construction and recounts it.

From n=9 onward, its actual counts obey H_(n+1)=2H_n+9, with H_9=1015. This scalar return is lawful because the realized state lies on the invariant relation q=F+9, g=0, h=117. It is not a universal scalar action on arbitrary source words. More generally every certificate(4), after its first zero step, realizes the invariant relation q=F+c and hence the guarded return F->2F+c. Source entry is the mathematical content; once it is certified, repetition is automatic.

### A direct proof of the native recurrence

A Hamilton path of the enlarged tournament either ends in z or has z immediately before x. In the first case delete z: the old path must not end in x, giving F-(h-g) possibilities. In the second case remove z and x: the remaining ordered parent paths(P,Q) contribute q when Q is nonempty, and h when Q is empty. These cases total F+q+g.

In an ordered parent two-path cover counted by q after adjoining z, z is last in either P or Q, since all old parent vertices dominate it. Deleting z gives one old cover in each of two disjoint ways; the marked out-neighbour set excludes z. Thus q doubles. Every Hamilton path of the new parent ends in z, so its marked-end count is zero, while h is unchanged. This proves(1) with the physical vertex clock retained.

### Why this is a lawful escape from the previous obstructions

The new parent has minimum feedback-vertex deletion number3, checked over every deletion subset. Its h=117 is outside the singleton/D grammar's powers-of-five invariant. It is not a transitive chain with one apex, and its periodic tail is constant, not mixed. None of packet f's obstructions applies to this source. Those results remain true in their stated domains.

The complete labelled parent scan through order6 finds no dyadic-slope source with deficit9. A seven-vertex parent does exist. Thus **seven is the minimum parent order for this particular pure-zero-tail source condition**, not a claim about the minimum tournament order for individual counts2^k-9. The witness census and its absence boundary have different scopes.

### Exact constant defect: three feedback vertices and six feedback arcs

For every r>=1, the family in(2) has feedback-vertex number exactly3 and feedback-arc number exactly6.

Deleting S={0,1,6} leaves the parent in transitive order(4,3,5,2), followed by all zero-tail vertices, followed by x. Conversely the original parent is induced and deleting any two of its vertices leaves a directed triangle. `feedback_profile.json` gives all21 two-deletion triangle witnesses. This proves the vertex claim at every r.

For the arc claim, the r tail vertices are a transitive clone block: all old parent vertices dominate them, and they dominate x. In any linear order, arrange the clone labels in their internal transitive order, without changing external costs. Then move all clones into the cheapest of their external slots. Their common external adjacency makes this cost at most the original sum, so an optimum has consecutive clones.

For each order of the eight old vertices and each of its nine slots, the backwards-arc cost is B+m*r, where m counts old parent vertices after the block plus the indicator that x lies before it. Complete enumeration gives the minimum intercepts at m=0,...,8:

    [6,5,5,5,5,5,5,5,7].

Thus the minimum is6 for every r>=1. An independent subset-DP check gives parent minimum3, marked-core minimum5 and the extended minimum6 at r=1,...,7. The finite slot enumeration plus clone compression proves the all-r claim; it is not extrapolated from the seven extended examples.

This is a fixed-defect, fixed-core construction outside the one-apex grammar. The previous bounded-loss obstruction and the present six-arc family are compatible: the number of feedback vertices is part of the model, not an interchangeable statistic.

## 2. A general cut-energy certificate for a whole infinite family

For any actual nonempty parent A and cut S, let b_uv=(W_uv+W_vu)/2 be the nonnegative integral symmetric ordered-two-path-cover weights. Put

    L(S)=sum_{unordered {u,v} subset S} b_uv,
    C(S)=sum_{u in S,v outside S,u->v} b_uv.

The inherited exact identities are

    F=h+2C,      q=2C+g+2L.

Consequently

    F+g-q=h-2L.                                  (3)

A finite actual cut satisfying

    q=2^a,       2L-h=c>0                         (4)

is a certificate for **every** count2^k-c with k>=a+1, by adjoining k-a zero vertices. The parent and marked set are part of the certificate. Merely solving the numerical equations in nonnegative synthetic coordinates does not establish an actual source.

`power_sources.cpp` computes literal subset endpoint tables and cover weights for actual parents. `audit_power_sources.py` independently recovers h,g,F and q by ordinary path counting on A, A+x and its first zero extension, without importing the cover-weight code. The identity q=H(first extension)-F-g gives an independent audit of every reported source. The deficit9 tail is directly counted through k=17.

### Complete fixed-path enumeration

Every tournament admits a Hamiltonian path: insert a new vertex before the first path vertex it beats, or at the end if it beats none. The preceding path vertex, if any, beats the inserted vertex. Relabel that path0->1->...->n-1. Therefore fixing those n-1 arcs and enumerating all other arcs, together with **all** marked subsets, covers every parent/cut up to relabelling. It is a presentation cover, not an isomorphism-class census.

Through order8 this scans2,131,020 fixed-path presentations and539,001,930 nonempty cuts. It finds548 distinct positive deficits with a dyadic q. The admitted set is complete for parent order<=8; each retained witness has a separate adjacency audit. The all-labelled scan through6 and the fixed-path scan through7 agree on their overlapping deficit sets. No division by Hamilton-path count or automorphism order is used.

In particular every odd c in1..115 has a retained certificate with a<=11 and parent order+1-a<=1. Hence, for every k>=12, **all58 odd integers in[2^k-115,2^k-1] are Hamilton-path counts**, with each witness having at most k+1 vertices. This is a uniform near-dyadic window, not coverage of every odd integer. The first missing positive odd deficit in the order-eight pure-tail source bank is117; this does not mean that2^k-117 is globally forbidden.

### A wider uniform window and its exact threshold obligation

The deterministic larger-parent sample and three targeted nine-vertex certificates extend the contiguous deficit bank to every odd c<=235. The targeted cuts for117,169,223 have states respectively

    (7755,8192,320,1363),
    (3821,4096,106,469),
    (7745,8192,224,1153).

Their complete adjacency data and independent audits are saved. Combining sources by the smallest available slope yields118 certificates, all with a<=13 and parent order+1-a<=1. Consequently, for **every k>=14**, every odd count in[2^k-235,2^k-1] is realized on at most k+1 vertices. `window_certificate.py` reconstructs and counts all118 actual tournaments at both common exponents14 and15, in addition to the source-level audits.

The first absent deficit in the merged bank is237. It was not found in the declared50,000-parent order-nine sample, the1,998 one/two-arc neighbours of the three recorded parents, or their1,536 one-vertex insertions. These are separate bounded searches, not nonexistence results for other parents or for the raw count line.

The source gate has a useful exact limitation. Since F+g>0, every certificate(4) has c<q. If it is used with r>=1, then

    N=2^(a+r)-c > 2^(a+r-1).

Thus its exponent k=a+r is necessarily the first power of two above N. One cannot increase k arbitrarily to bypass an expensive source. For a target N between2^(k-1) and2^k, the sufficient source obligation is precisely an actual certificate for c=2^k-N with dyadic q<=2^(k-1). The existence of a source at a larger q does not discharge that target.

A cheap re-marking probe also records its boundary: among eight natural uniform re-markings of each of twelve sources, only one retains a dyadic q. This96-case calculation rules out no general re-marking scheme; it shows why copying the old cut, complement or scalar count does not automatically preserve the needed source interface. Higher cover response remains available as an additional coordinate.

## 3. Exact rooted port orders in the square-sum carrier

See [PORT-ORDERS.md](PORT-ORDERS.md) for the full theorem and proofs. Each primitive switch permutes one contiguous interval of the reference cycle, with fixed ordered ends and unchanged label set. These intervals are disjoint. A coarse source bit may control several intervals, and those ties must be retained.

Normalize the cycle at its invariant edge1--3, or1--8 for the protected native carrier. Required rooted orders of selected labels are then constant region-order checks followed by at most one binary decision per region. Shared-bit contradictions remain visible. In a fully refined carrier the compiler is complete for **all Hamiltonian cycles of that subgraph**.

For any k distinct marked vertices, at least2^max(b-floor(k/2),0) cycles preserve their reference rooted order. One marker cannot detect a permutation of its region, so each paid bit requires two markers in a disjoint region. This preserves more freedom than prescribing complete incident neighbourhoods. It does not extend arbitrary port orders in the full square-sum graph.

Joint edge and order requirements compile to the same shared-bit constraints. Weighted cycle sums factor, and arbitrary signed additive edge costs can be optimized independently per free bit. These are exact consequences for the constructed graph, not an algorithm for unrestricted Hamiltonicity.

## Open targets and reproduction

The new deficit9 theorem supersedes the previous packet's open task for that particular line. **Universal completeness of all positive odd raw-H counts except7,21 remains OPEN.** The sharper question is now quantitative source entry: how large a parent and dyadic slope are needed for each deficit, and do those costs permit a full interval rather than fixed-width windows? A universal source-existence theorem alone would still need an exponent threshold matching the requested integer.

Run the scripts listed in `VERIFICATION.json`. The complete order-eight fixed-path census is the long run; its final witness bank is independently auditable without rerunning the discovery search. `MANIFEST.json` pins inherited inputs and packet bytes. Normal/optimized consumers, independent adjacency counts, actual cycle walking and complete small order-query controls are recorded separately. No subagent audit or external novelty claim is made.
