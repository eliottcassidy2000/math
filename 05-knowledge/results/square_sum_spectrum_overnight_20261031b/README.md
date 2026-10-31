# Recursive forests, finite scalar charts, and a reversible four-vertex block

**Status: PROVED scoped composition/decoder/obstruction results; FINITE-EXACT restricted tournament spectrum; unrestricted spectrum completeness OPEN.** Same-agent mathematical review, with independent literal path counts and a separate C++ word-tree traversal. No external novelty claim is made.

## Inheritance and concept board

Two incoming Session124 checkpoints changed the work rather than merely adding citations. Commit `146717912b` supplies an explicit infinite square-sum ray preserved by reusable compatible-pair transitions. Commit `727bfdf464` lowers its base boundary to42 and independently certifies every allowed odd **unrestricted** tournament count through999999. That second checkpoint also supplies a positive four-coordinate block action and proves why repeated triangle blocks lose non-3-divisible outputs. The present results use these mechanisms and keep their scopes.

| Lane | Inherited mechanism | Hostile / repaired near miss | New use |
|---|---|---|---|
| Anchor | Exact positive marked-block action on (F,q,g,h) | Repeated triangle blocks have a uniform three-adic loss | The strong four-vertex block has H=5,E=32 and an exactly decodable cut history |
| Anchor | Fixed-source scalar-controller rigidity | Switching finitely many terminal words might appear to restore amplification | Every closed return is neutral; a finite append-only scalar atlas remains too small |
| Niche | New compatible B recursion and its threshold42 extension | The new ray and the old native path conflict on small vertices if simply superimposed | Lift2303 disjoint copies into intact interiors while preserving the native path on singleton ports |
| Wildcard | A scalar arithmetic port can retain a word | The same (h,q) loses generator order when singleton and four-vertex blocks mix | Keep the full marked state and integer clock; an explicit seven-vertex collision measures the loss |
| Information | Four labelled cut choices give16 subsets | Sixteen valid residues need not form a linear four-bit code | Five fixed linear probes are necessary, but four adaptive linear probes suffice |
| Incoming reroute | Session124's canonical weighted-cut law | A nonnegative cover coordinate can still conceal negative interior weight | Positive cut/interior coordinates give a new exact inverse guard |

This uses the META-PATTERNS card **Canonicalize mechanics; retain the predicate**. The source labels, allowed operations and retained endpoint observables determine which compression is valid. No graph-class identification follows from a shared number.

## 1. One native path and 2303 simultaneous ray prefixes

Keep the previously proved native path

    R=(1,8,17,19,6,10,15,21,4,12,24).

Let P_0=(1,3) and, in one-based order, put

    P_(l+1)=(1,3,6,19) concatenated with (49p_j+(-1)^j4)_j.
    S_l=(85*49^l-1)/2.

These are the incoming [threshold42 ray and compatible-pair construction](../session124_square_prefix_threshold42_20261031/README.md), with |P_l|=4l+2. Define

    U=[-1200,1200] minus ([-24,24] union [74,122]),  |U|=2303,
    R_(l,t)=(2401p_j+(-1)^(j+1)t)_j,  t in U.

**Theorem.** For every l>=0 and every N>=S_(l+2), one square-sum Hamiltonian cycle on[1,N] simultaneously contains:

- the original native path R;
- the original seed482–3362–6242;
- all2303 pairwise vertex-disjoint paths R_(l,t), each of length4l+2 vertices.

The thresholds are sufficient, not claimed optimal. The first two are102042 and5000082. Every finite prefix of each of2303 fixed infinite rays can thus be retained simultaneously with R, for every sufficiently large initial interval. This does not assert an infinite Hamiltonian cycle or arbitrary prescribed-forest embedding. The older exact threshold72 for R alone remains unchanged.

### Proof with both boundaries retained

Use the new42-based B constructor, whose base cutoff is2081 and which recurses from2082. Its certified common prefix is P_l whenever its index m is at least S_l. Reusing the old cutoff2620 would invalidate this invariant; the constructor here uses the incoming repaired boundary.

Feed this B pair to the old protected Aplus/Aminus transition. Those tables require only the weak B contract, which the new pair satisfies. They split a1–3 prefix only at offsets0 or−2. Every other offset c gives an intact copy

    (49p_j+(-1)^j c)_j.

There are47 uniformly safe offsets. All these copied labels are at least25. In the final lift, the only internal cuts are at the original A vertices2 and7; neither can cut one of the retained copies. The [native final392-state table](../square_sum_spectrum_overnight_20261031a/README.md) still joins exactly the same typed segments and singleton ports, preserving R and the482-based seed.

The second affine lift changes the copy parameter into 49u plus a signed offset in[-24,24], where u=-c. The sign is the common positional parity of the copy in the A pair. As the second offset varies, either sign gives the entire interval[49u-24,49u+24]. Omitting c=0,-2 is exactly omitting u=0,2, giving U above. Hence every required path has its original displayed labels, independent of N and of which input member or direction was selected.

Distinct p_j lie in disjoint width2401 bands centered at2401p_j. Within a band, different t give different labels. Thus the2303 paths are pairwise disjoint and avoid the singleton range[1,24]. The t=961 ray begins3362–6242, so the old seed attaches at its endpoint without consuming a third incidence. At l>=1 it continues through15367, as expected from the inherited four-vertex seed extension.

Two successive quotient steps require

    m=floor((N-1200)/2401) >= S_l.

This follows from N>=2401S_l+1200=S_(l+2). The old A and final port identities depend on parity and the weak pair contract, not on the former numerical base53. At m=42 the blocks, prefixes and cut segments remain nonempty; the actual controls explicitly test this lowered boundary. This proves the all-size statement using the complete finite tables.

`joint_rays.py` replays all98 new paired states, all11 new small bases, the intact-offset condition in all196 old A states, and all392 native final states. It then independently stream-checks complete cycles at102042,102043 and5000082, including every vertex, square edge and prescribed forest edge. Normal and optimized results match:10646608 checks, with11526 distinct prescribed edges in the largest control. The incoming independent audit remains the source of the new recursion theorem; the infinite conclusion here is the intact-segment proof, not an extrapolation from these three cycles.

### Late incoming corollary: the same interface accepts mixed-radix rays

Commit `4338480c37` arrived during close-out. Its [mixed-radix theorem](../session124_mixed_radix_square_rays_20261031/TERNARY_SCHEDULE_THEOREM.md) supplies three operators with(radix, offset, fixed prefix)

    A=(25,12,(1,3)), B=(25,8,(1,3,6,19)), C=(49,4,(1,3,6,19)).

Every finite code w gives a common weak-B prefix P_w from index S_w=(85R_w-1)/2 onward, where R_w is the product of its radices. Its proof explicitly guarantees the same per-label parity agreement required above, including independently colored split heads. Therefore the preceding intact-segment argument applies verbatim with P_l replaced by P_w: the native path, seed and2303 lifted copies coexist for every

    N >= 2401S_w+1200 = (85*2401R_w-1)/2.

This is a composition corollary of the incoming audited weak-pair contract and our audited outer interfaces. The present `joint_rays.py` controls cover the all-C cases only; no new mixed-code full-vector audit is claimed here. For pure B at depth l the new rays have4l+2 vertices with the25^l threshold factor, but their actual labels differ from the old C rays. This is not a threshold improvement for an unchanged prescribed ray. All2303 copies in this corollary use the same code w; independent code choices per offset need an additional ownership proof.

## 2. A finite scalar chart atlas still cannot generate the spectrum

An **append-only scalar chart** consists of a fixed source alphabet Sigma and a fixed terminal word A. Its objects are all words w in Sigma*, observed by H(wA). A directed controller edge from(Sigma,A) to(Sigma',A') appends a fixed word B and promises

    H(w B A')=a H(w A)+b  for every w in Sigma*.

The target domain is part of the promise. For these pure alphabet domains it forces Sigma to be contained in Sigma' and B to belong to Sigma'*. Finite atlases have finitely many charts and edges and start from finitely many actual source words. Controllers with source-dependent guards, re-marking, word reversal/reencoding or vector-valued observations are outside this definition.

**Theorem.** If every source alphabet has at least two letters, every closed controller walk acts identically on H, and the entire finite atlas attains only finitely many counts from its finite seeds. Allowing singleton alphabets increases this to at most O(log X) distinct counts up to X. In particular no such fixed atlas proves the unrestricted odd spectrum conjecture.

### Closed returns and cycle removal

On a directed cycle the alphabets must all coincide, by inclusion. Concatenating its inserted blocks gives B in Sigma*, and the composite edge has the form H(wBA)=aH(wA)+b. For two-letter and four-letter alphabets, the [previous packet](../square_sum_spectrum_overnight_20261031a/README.md) proves a=1,b=0 whenever this holds universally. The three-letter cases follow from the same certificates: each contains{0,3} or{1,2}, whose zero-mass prefix identity has a nonzero constant value for every suffix, forcing a=1. If3 is present, the 3^m odd-cycle injection forces unbounded insertion debt except for sink padding. For{0,1,2}, any inserted letter misses a marked out-bit; choose the constant source word1^m or2^m that beats that marked vertex. The same injection gives debt at least2^m for every nonempty insertion. Thus closed returns are neutral in all alphabets of size at least two.

Erase a closed portion of a controller walk. The actual prefix changes, but the count and chart at the return do not. Every later controller is universal on that entire chart, so all later counts remain unchanged. Repeating this removal reduces any walk to a simple chart path with the same final count. There are finitely many such paths and seeds. This proves the first assertion without imposing a word-depth bound.

For singleton alphabets, an unbounded portion of a walk produces just a constant-letter run a^m. Alphabet inclusion prevents changing from one singleton alphabet to another, or returning to a singleton after entering a larger alphabet. Erase all cycles in the larger-alphabet portion. The final spectrum is consequently contained in a finite bank of fixed-context one-run templates p a^m q, plus finitely many counts.

The variable run is a transitive module: every outside vertex has the same relation to all its vertices. If its representative lies on no directed cycle, insertion or erasure of the whole module leaves H unchanged by the odd-cycle-collection formula. Otherwise it lies on a directed triangle: if there were no such triangle, every in-neighbor would beat every out-neighbor, preventing a return from an out-neighbor to the vertex. Replacing that triangle's module vertex by any odd-sized increasing subset of the m clones gives2^(m-1) distinct odd cycles. Their singleton collections alone give H>=2^m+1 for m>=1. Thus at most floor(log2 X)+2 distinct values from each template lie below X. This proves the logarithmic bound. The general bounded-run obstruction was already present in [run_density.md](../feedback_word_hamilton_spectrum_20261020/run_density.md); the new step is deriving the finite one-run normal form from scalar chart returns.

The singleton boundary is real: for every w in{1}*,

    H(w110)=2H(w10)-1.

This one-chart loop gives3,5,9,17,... . It shows why the theorem distinguishes unary pumping from genuinely branching source languages.

### Complete small atlas control

`finite_terminal_atlas.py` enumerates every terminal of length at most3 and every appended{0,1}-block of length at most2:85 charts and2015 exact affine edges. Four actual prefix rows form a closed observable source space **including the affine constant coordinate**. The unreduced36-state reachable space is larger; confusing it with this observable dimension failed the first internal test and is not used.

All957 edges within strongly connected chart components have consistent affine gauges, so every cycle has identity composite. Among these,450 append a nonempty block; all such blocks contain only1s. The full closure from the85 empty-prefix chart seeds has265 chart/count pairs and19 counts, the largest59. An independent vertex-subset algorithm checks8060 literal edge instances. This census illustrates the theorem; it does not supply its all-size proof. `clone_run_controls.py` independently checks100 module templates and400 literal path counts.

## 3. The strong four-vertex block retains its entire cut history

The incoming [million-prefix proof](../session124_million_path_spectrum_20261031/PROOF.md) gives the following exact positive block action. For a parent A and marked external vertex x, retain

    F=H(A+x), h=H(A),
    q=number of ordered path covers(P,Q) of A with Q nonempty and x->first(Q),
    g=h-t, where t counts parent paths whose last vertex beats x.

For a marked right ordinal block B, its profile(H,E,p,q_B,s,t_B,F_B) acts by

    F'=H F+p q+(H-s)g+(F_B-H)h,
    q'=E q+q_B h,       g'=(H-t_B)h,       h'=Hh.

The definitions, disjoint path decomposition and positivity proof are inherited unchanged. The four coordinates are counts of actual path/cover objects, not arbitrary nonnegative lattice points.

Take the strong four-vertex tournament D with out-neighbor masks(8,1,3,6), in vertex order0,1,2,3. Literal enumeration gives H(D)=5 and E(D)=32. The q_B contributions of the four vertices are respectively

    (5,6,9,7).

Thus a cut mask M, where x beats precisely those vertices, has digit d(M) equal to the corresponding subset sum. All16 digits are distinct:

    Ddigits={0,5,6,7,9,11,12,13,14,15,16,18,20,21,22,27}.

There is a small exact K4 connection: the six two-vertex cut sums are exactly11,12,13,14,15,16. Singles, pairs, triples and the full set occupy disjoint numerical ranges, while the six pair sums are distinct. This explains subset uniqueness. Complementing the cut sends d to27-d, since E-H=27. These are path-cover weights on actual marked cuts, not an identification with ordinary graph-cycle data.

**Decoder theorem.** Fix an actual initial state with odd h_0. Any sequence of these D blocks has a unique cut word determined by its terminal pair(h,q). If h=5^r h_0, set h_previous=h/5. The last digit is

    d = q * h_previous^(-1) modulo32,

chosen in Ddigits, and then

    q_previous=(q-d*h_previous)/32.

Distinct admissible residues give a unique cut mask. Repeating recovers every cut and eventually the initial q. Conversely this procedure, including divisibility, nonnegativity and the final seed check, decides membership in this pure-block orbit. Replaying the recovered word checks the other coordinates F,g. It is therefore a genuine actual-word membership criterion for this restricted semigroup, unlike membership in the old17-dimensional observation lattice.

Equivalently,

    q_r=32^r q_0+h_0 sum_(j=1)^r d_j 32^(r-j)5^(j-1).

There are exactly16^r distinct residues q_r modulo32^r for the fixed seed and length. The five-adic height and two-adic clock are both retained; multiplying by the odd source height must be inverted before reading the terminal digit. `diamond_profiles.py` checks all4369 cut words through length3 and all16 literal profiles.

### One adaptive probe saves a bit

View the normalized digit as a five-bit vector. For every nonzero v in F_2^5, the packet gives two admitted digits whose XOR is v. Any fixed linear map to four bits has a nonzero kernel direction, hence identifies an admitted pair. Five nonadaptive linear bit probes are therefore necessary and sufficient.

Nevertheless four **adaptive** linear probes suffice. The recorded depth-four decision tree first reads the two low bits, then chooses the remaining parity queries according to those answers. All16 leaves decode correctly, meeting the information lower bound log2(16)=4. This is a precise5-versus4 separation on the normalized digit oracle; it is not a claim that normalization of an unknown raw q is free, or that four fixed linear probes suffice. The31 collision witnesses and complete decision tree are in `diamond_profiles.json`.

### Mixing generators exposes the missing coordinate

Let S_b denote a singleton with x beating it iff b=1. Its actions are

    S_0:(F,q,g,h)->(F+q+g,2q,0,h),
    S_1:(F,q,g,h)->(F,2q+h,h,h).

From the same source(1,0,0,1), the two five-vertex additions

    D_mask2 ; S_0   produce (17,12,0,5),
    S_0 ; D_mask9   produce (11,12,2,5).

Both resulting tournaments have seven vertices. They agree on h,q and the actual vertex clock but disagree on F,g. Full adjacency masks and independent direct counts are retained. Thus the pure-block decoder cannot silently decode arbitrary mixed words. This is a useful missing-coordinate witness, not a reason to abandon the richer state.

Repeated mask2 also gives a concrete escape from the triangle valuation obstruction. From that same seed, its outputs satisfy, for r>=1,

    F_r=(100*32^r+(216r+359)*5^r)/675.

For r>=3, F_r=3*2^(r-2) modulo5, so repeated use does not force divisibility by5. A different, source-specific filter survives: F_r=(r-1)*2^r modulo3. No primality conclusion follows. The first values9,181,5041,156477 illustrate why avoiding one divisibility obstruction is insufficient.

## 4. A compact new candidate grammar and its exact current boundary

Let G be the class obtained by adjoining one marked vertex to an ordinal sum of singleton blocks and copies of D, with arbitrary cuts. Starting with one parent vertex and either cut gives the two initial states(1,0,0,1),(1,1,1,1). These generate its entire count spectrum: a parent starting with D can first be given an extra global source, which changes no Hamiltonian-path count.

**FINITE-EXACT over all word lengths:**

    spectrum(G) intersect [1,100000]
      = {positive odd values <=100000} minus {7,21,43,93,173,247}.

The primary exact-state traversal visits2095168 states. A separately implemented C++ traversal enumerates literal block profiles and visits63403536 word-tree nodes **without state deduplication**. Both obtain49994 counts and the identical100001-byte support bitmap hash. This result is about the restricted class; all six omitted values are not being declared globally forbidden. In particular the four extra omissions are known globally attainable.

Here is the all-word pruning proof. F and q never decrease, 1<=h<=F, and0<=g<=h. If F exceeds the ceiling it cannot return. If q exceeds the ceiling, any subsequent block with p>0 overshoots F at once. The only p=0 generators are S_1 and D_mask15; they are actual ordinal products and multiply F by1 or5. The consumer records **all powers-of-five tails** before stopping a q-overflow branch. Among states with F,q bounded, h and g are bounded too, so exact-state deduplication terminates. In the independent word tree, every positive q at least doubles at each step; at q=0 only exact singleton self-loops are removed, while a zero-cut D increases h by5. Thus that traversal also terminates without an arbitrary depth cutoff. The singleton-only pruning rule would have lost product tails here; the repaired rule explicitly preserves them.

**OPEN strengthened target:** Is spectrum(G) exactly the positive odds excluding those six values? Equivalently, does this specific grammar attain every odd value at least249, together with its verified smaller values? The current computation is evidence only. If proved, the known witnesses for43 and173, and the ordinal products93=3*31 and247=13*19, would complete the unrestricted spectrum apart from7,21. This is a concrete sufficient target with one primitive four-vertex block and a positive four-coordinate state.

The incoming unrestricted coverage through999999 remains the stronger **finite unrestricted** result and is credited to Session124. We did not replace it with this smaller bound. The point of the new computation is the much more specific grammar, exact decoder and four explicit grammar obstructions. Further large censuses should serve that proof target rather than merely repeat known unrestricted coverage.

### A complete inverse with a decreasing integer rank

There is now an exact inverse membership procedure for the mixed grammar, not just the pure-block decoder. Every reachable h is5^r. Its possible g values are{0,h} when h=1, or{0,h/5,2h/5,3h/5,4h/5,h} otherwise, because the last generator overwrites g.

For h>1, a proposed last D block has previous height h/5. Its cut digit is forced by q modulo32 after division by that odd height. The displayed forward formulas then determine previous q and, for each of at most six possible previous g values, previous F. Require all divisions to be integral and all state inequalities. A proposed last S_0 requires g=0 and even q, giving previous q=q/2 and previous F=F-previous_q-previous_g. A proposed last S_1 requires g=h and q-h nonnegative even, giving previous q=(q-h)/2 and unchanged F.

After deleting exact self-loops, **every inverse step strictly decreases the positive integer F+q+h**. A D step decreases h by5 and q by at least a factor32; an S_1 step decreases q; an S_0 step decreases q unless q=0, in which case any nontrivial predecessor decreases F. The recursion therefore terminates.

At h=1 there have been no D blocks. Membership is exactly

    g=q mod2,   F=f(q),
    f(0)=1, f(2k)=f(k)+k+(k mod2), f(2k+1)=f(k).

The binary expansion of q reconstructs a singleton word from the fixed source(1,0,0,1). Thus the inverse stops at an actual source, not an arbitrary positive lattice point. A successful inverse yields the complete ordered word and its vertex clock.

**Finite target-entry bound.** Remove trailing p=0 generators; they contribute only a power-of-five factor and sink padding. If the remaining core is not the source, its last generator has p>0. For S_0, q_new<=2F_new. For every D cut with p>0, the exact16-profile table verifies p>=5 and5q_B<=32F_B. Since F_previous>=h_previous,

    F_new >= p q_previous+F_B h_previous,
    5q_new=160q_previous+5q_B h_previous <=32F_new.

Consequently a target F is in this grammar iff, for some trailing factor5^k, one of the finitely many canonical states at F/5^k has

    h a power of5, h<=F/5^k, g<h in the stated six-value set,
    0<=q<=floor(32(F/5^k)/5),

and is accepted by the well-founded inverse. The source case F/5^k=1 is included. This is a complete decision criterion for the restricted grammar. **It does not prove that any candidate is accepted for every allowed target.** That uniform entry obligation is the remaining conjecture in a concrete form.

`diamond_inverse.py` independently exhausts the six grammar omissions and constructs249,251,1001. It also constructs the important incoming target103061=7*14723, whose smaller factors cannot supply an ordinal-product repair:

    103061 = H(S_1 S_1 S_0 S_1 S_1 S_1 S_1 S_1 D_mask10 S_1 S_0 S_0 S_0),

from the fixed source, with final state(103061,114424,0,5) on18 vertices. This is a new witness **inside this grammar**, not a claim that the count was previously unattained globally; Session124 already supplies a16-vertex witness using a different block. A separate literal-adjacency DP checks all four positive inverse witnesses, with no transfer or inverse-code imports. The inverse retains sources, generator order and integer vertex clocks.

## 5. Incoming weighted cuts reveal a positive interior coordinate

Two further incoming checkpoints, `80e780001b` and `f5076a605f`, were read during the session. Their fixed-ear, bounded-width and triangle-packing obstructions limit specific carriers, not the raw spectrum. The [canonical weighted-cut law](../session124_tournament_cut_response_20261031/README.md) supplies a useful new coordinate for the actual four-state block action.

For a **nonempty** parent A with marked out-neighbor set S, let b_uv be the incoming positive symmetric cover weights. Define

    C=sum_(u in S,v outside S,u->v) b_uv,
    L=sum_({u,v} subset S) b_uv.

Then the existing marked state satisfies the exact universal identities

    C=(F-h)/2,    L=(q-F+h-g)/2,    q=2C+g+2L.

To prove the second identity, let s_v,t_v count parent paths starting/ending at v, and W_uv count nonempty ordered two-path covers ending at u and starting at v. The coefficient of v in q is

    q_v=s_v+sum_u W_uv=t_v+2 sum_(v->u) b_vu,

by the incoming endpoint-divergence identity. Sum over S. Each internal pair contributes once to the weighted outdegree sum; the rest contributes C, while sum_(v in S)t_v=g. This proves q=g+2C+2L. Thus the previously hidden coordinate is a literal nonnegative integer, not a formal sign correction.

Every actual state consequently satisfies

    q>=F-h+g,    q-F+h-g even.                         (interior guard)

The source-specific inverse now enforces this necessary guard at every predecessor and at entry. It is not sufficient for reachability; the inverse still must recover the actual source. For target103061, this reduces the number of memoized candidates from5849860 before the guard to409690 after it, while recovering the identical18-vertex word. Both optimized and normal executions agree. This is a concrete improvement obtained by merging incoming mathematics, not a universal complexity estimate.

### Positive transfer in the new coordinates

Use the state(C,L,g,h). For any nonempty right block with profile(H,E,p,q_B,s,t,F_B), put

    alpha=(p+H-s)/2, beta=(E-H-p+s)/2,
    a=(F_B-H)/2, ell=(q_B-F_B+t)/2.

All four coefficients are nonnegative integers. Here a is its weighted outgoing cut, ell its internal weight, alpha=(H-s)+weighted_indegree(V-S), and beta=s+weighted_indegree(S). Also p<=E-H and alpha+beta=E/2. Substituting into the inherited block action gives

    C'=(H+p)C+pL+alpha*g+a*h,
    L'=(E-H-p)C+(E-p)L+beta*g+ell*h,
    g'=(H-t)h,    h'=Hh.

This is an exact positive integral conjugate, not a loss of coordinates. For singleton additions it becomes especially simple:

    S_0:(C,L,g,h)->(2C+L+g,L,0,h),
    S_1:(C,L,g,h)->(C,C+2L+g,h,h).

Both double C+L and add the old g, while the new bit selects which of the two weights changes and sets the next g. This gives a precise two-channel interpretation of the binary clock. It neither asserts that arbitrary positive coordinates are realizable nor pays any unrelated Collatz debt.

The independent literal-cover consumer exhausts all1099 labelled nonempty parents of order at most5, all33866 cuts and actual child-DP counts, and541856 comparisons with the diamond block action. It imports no incoming cut-response implementation. The algebra above proves the universal statement; the exact controls audit its conventions, positivity and source interpretation.

The next entry problem can therefore be phrased using genuinely nonnegative cut and interior weights. A theorem covering all C at some actual height h, rather than a scalar controller on unrestricted words, would address the remaining spectrum gap. Universal coverage of those states is still OPEN.

## Reproduction and limits

Run from this directory, retaining the sibling packages named above:

```bash
python3 joint_rays.py
python3 finite_terminal_atlas.py
python3 clone_run_controls.py
python3 diamond_profiles.py
python3 diamond_spine_scan.py --limit 100000
python3 diamond_inverse.py
python3 diamond_inverse.py 103061
python3 audit_inverse_witnesses.py
python3 cut_interior_coordinates.py
c++ -std=c++17 -O2 diamond_tree_audit.cpp -o /tmp/diamond_tree_audit
/tmp/diamond_tree_audit 100000 /tmp/diamond_tree_support.bin
c++ -std=c++17 -O2 ../session124_million_path_spectrum_20261031/audit/adjacency_dp.cpp -o /tmp/diamond_adjacency_audit
/tmp/diamond_adjacency_audit inverse_witnesses/H249.txt inverse_witnesses/H251.txt inverse_witnesses/H1001.txt inverse_witnesses/H103061.txt
```

Python normal/optimized receipts agree. The two tournament traversals agree by exact support hash, not merely support size. The full joint forest control checks its own actual cycle. Small mixed-block controls count full adjacency paths and path covers independently of the inherited block action. `MANIFEST.json` freezes inputs, source and receipts. No claim of universal Collatz, LRC, plane coloring or Jacobian progress follows from these scoped constructions.
