# A shared-C3 global observer, with exact sandwich resource profiles

24 October 2026. **PROVED:** the all-size labelled construction, normalization and capacity statements below. **INHERITED COROLLARY:** the sandwich square/diagonal profile and looplessness specialize earlier center and nonstrong-source laws. **FINITE-EXACT:** the complete controls stated in `normal.json`. No new classification of all metagraph neighbors, cross-order harmonic transport, or general global-observer theorem is claimed.

## 1. Inheritance before extension

All archival inputs are pinned at repository commit `c09b642a629b7d12b529b88c2697abaee78e73e7`; original paths and byte hashes are in `SOURCE_PINS.json`.

- **Closest proved mechanism:** THM-2454, `01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md`, §2, gives palindromic-stack orbit and blue multiplicity. `deficiency_three_20261020/centers.md` §2 already states `h(A -> C -> A^op)=h(A)^2 h(C)` and `b(A -> C -> A^op)=h(A)b(C)`.
- **Existing quotient mechanism:** `split_word_metagraph_spectra_20261021/session23-slack-spectrum/README.md` §4 explicitly identifies an earlier path-state square with coordinate-swap reflection. Its hostile controls show that incidence counts do not determine distinct neighboring classes.
- **Existing loop-free mechanism:** THM-4729, `01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md`, L3, proves that every nonstrong source has only strong complement partners and no loops. We reproduce its elementary closing-arc argument in our gauge rather than claim it as new.
- **Canonical hostile:** `session47_resilient_native_wiring_20261024/global-observer/THEOREM.md` exhibits coherent local partner hosts that cannot share one observer: two required outgoing observer arcs at one owner. Its directed-block compiler is inherited, not contradicted.
- **Corrected near miss:** two different full phase lifts of one normalized resource are not two available resources. Moreover, distinct path orbits of an arbitrary R need not survive as distinct local F/reflection resources of `R -> C3`. Our positive construction uses the authenticated native prefix keys, not that invalid inference.
- **Least-used sidecar:** the actual anti-isomorphism of the shared C3. It fixes the source orientation on the overlap and forces the right native phase to be `1-k`, not `k`.

Live board: path-orbit square / converse involution / normalized unit capacity / contiguous observer overlap / marked endpoint gluing / two integer work counts. Anchor: a positive shared global observer. Niche: exact direct-metagraph profiles. Wildcard: literal Paley11/1729 on a 25-owner, 24-path-arc sandwich, without unrelated numerical identifications.

## 2. Labelled agreement paths and normalization

For a tournament A on n labelled owners and a directed Hamiltonian path
`P=(p_0,...,p_(n-1))`, define `K_P(A)` by keeping the n-1 consecutive path arcs and reversing every other arc. The directed agreement graph of A and `K_P(A)` is exactly P. Conversely, if two tournaments on the SAME labelled owner set have a directed agreement graph that is one Hamiltonian path P, the second is exactly `K_P(A)`. The path is unique. Thus the map

    (A,P) -> (A,K_P(A))

is a bijection onto labelled tournament pairs whose agreement graph is a directed Hamiltonian path. It is equivariant under common owner relabelling. Also `K_P(K_P(A))=A` and

    K_(reverse P)(A^op) = K_P(A)^op.

For n=1 the sole empty edge set corresponds to the sole path; all later sandwich applications have n>=5. Forgetting the labels or retaining only the two tournament isomorphism classes loses the agreement graph and does not supply a gluing witness.

Normalize along P by naming p_i with position i, and record the nonconsecutive upper-triangle bits

    w_A(P)_(i,j)=A(p_i,p_j),     0<=i<j<n, j-i>=2.

 Consecutive arcs are all 1, so this word determines the full normalized tournament. Equality `w_A(P)=w_A(Q)` holds iff the positionwise map `p_i -> q_i` is an automorphism of A. Consequently normalized words are exactly Hamiltonian-path orbits under Aut(A); the action is free because an automorphism fixing an ordered spanning path fixes every vertex.

Let F complement every word bit, and let rho permute coordinates by

    rho(w)_(i,j)=w_(n-1-j,n-1-i).

Then `w_(K_P(A))(P)=F(w_A(P))`, and

    w_(A^op)(reverse P)=rho(w_A(P)).

The latter formula has NO extra bit complement: the reversal of positional order and tournament converse cancel. F and rho commute and are involutions. A normalized unit resource is the orbit under `<F,rho>`, not one of its full labelled lifts. A blue word satisfies rho(w)=w; a black word does not.

## 3. The complete sandwich profile: an inherited law with explicit resource map

Let R be ANY nonempty finite tournament, set `h=H(R)`, `g=|Aut(R)|`, `a=h/g`, and form fresh labelled copies

    S = R_L -> C3 -> R_R^op.

No strongness or self-converse hypothesis on R is required. Put q=|R|, n=2q+3, label the middle vertices c_0,c_1,c_2 with `c_i -> c_(i+1)`, and identify labels in the two copies of R by the indicated converse map.

Every Hamiltonian path of an ordinal sum visits its summands in order and consecutively: after crossing to a later summand there is no arc back. Aut(S) preserves the unique ordered strong-component list, hence each summand and all its internal components. Therefore

    H(S)=3h^2,     |Aut(S)|=3g^2.

This remains valid if R itself is decomposable or has C3 strong components: an order-preserving automorphism of the finite SCC chain fixes each position, so equal adjacent SCC types cannot be permuted.

Let X be the a path orbits of R. A path of S has the unique form

    P + (c_k,c_(k+1),c_(k+2)) + reverse(Q),          (1)

with P,Q Hamiltonian paths of the respective R copies. Independent automorphisms in the two halves and the middle cycle show that its normalized word is determined exactly by `([P],[Q]) in X x X`. Conversely equality of normalized words gives an S-automorphism, whose restrictions identify these two path orbits. Thus this is a bijection, not merely a numerical product.

An anti-automorphism beta of S swaps the two labelled R copies and sends `c_i` to `c_(-i)`. Applying beta to the reversed path in (1) gives

    Q + (c_(1-k),c_(2-k),c_(-k)) + reverse(P).

Since the three middle phases are in one automorphism orbit, rho acts on `X x X` by exactly `(x,y) -> (y,x)`. The blue locus is the diagonal, of size a, and the black locus has size `a(a-1)`.

Every complement partner is strong. In (1) the first vertex lies in R_L and the last in R_R^op, so their source arc points first-to-last. They are not consecutive because n>=5. K reverses that arc while retaining the whole forward path, producing a directed spanning cycle. Hence every target is strong, while S is nonstrong. There is no F loop, either before or after converse merging. In particular a source word never lies in the F-image of the source fibre.

It follows that the F/rho resource orbit meets this source fibre in `{(x,y),(y,x)}`. For x=y it has two total words; otherwise four. Distinct unordered pairs give disjoint resource orbits. Consequently the exact loop-free incidence/resource profile is

    normalized source words                  a^2
    blue incidence / blue resources          a
    black weighted incidence                 a(a-1)
    fundamental paired-black resources       a(a-1)/2
    all normalized F/rho resources           a(a+1)/2.

The black weighted incidence counts two source words per fundamental paired-black resource. These are NOT counts of distinct neighboring isomorphism classes. The source is a point of symmetry (pure blue) iff a=1. This is the C3-centered specialization of the inherited center law; the new use here is the fully retained resource map and the construction below.

The nonempty-R condition is needed for this particular strong-target proof. With empty wings S=C3 and its complement along a Hamiltonian path is transitive, not strong.

## 4. Native two-cell construction with one actual global observer

Use the authenticated rigid eight-owner source T from the archived `CUBES.md`, canonical code29558785, with

    P0=(3,0,6,7,1,5,4,2),
    P1=(3,1,6,7,0,5,4,2),
    g=(4,0,2,3,1,5,6,7).

The literal g is an isomorphism `K_P0(T) -> K_P1(T)` and fixes the actual first and last owners3 and2. Let `T_r=T -> ... -> T` have r>=1 blocks, q=8r, and let `P_b` be the concatenation selected by a binary key b in `{0,1}^r`. Pins at six common clocks in every block give exactly these `m=2^r` paths: the only unpinned owners are0 and1, at the two available clocks, and both assignments are precisely P0 and P1. The seed has H(T)=409 and Aut(T)=1, verified independently below. These m pinned paths are a selected family, not all `409^r` prefix paths.

Construct

    S_r = T_r^L -> C3 -> (T_r^R)^op.

The left local cell is `L=T_r^L -> C3`; the right local cell is `D=C3 -> (T_r^R)^op`. They overlap in the SAME three actual C3 owners, not three copies identified after the fact. The right cell is a relabelled converse of L via the map that sends each prefix owner to its right counterpart and sends native cycle owner c_j to shared owner c_(-j). It is this anti-isomorphism that makes the two induced sources agree on the overlap.

For any keys a,b and phase k define

    O_(a,b,k)=P_a^L + (c_k,c_(k+1),c_(k+2)) + reverse(P_b^R).   (2)

This is a Hamiltonian path of S_r. Its left and right restrictions are contiguous segments with exactly the same three-owner path on the overlap. Therefore

    K_(O_(a,b,k))(S_r)|L = K_(O|L)(L),
    K_(O_(a,b,k))(S_r)|D = K_(O|D)(D).                         (3)

Indeed a pair internal to either contiguous segment is globally consecutive iff locally consecutive; this checks every positive and negative adjacency flag, not only the local subsequence order. The two local partner tournaments consequently agree on their shared C3, and (3) is an actual single global-observer realization of both.

In the original native L frame, the right local record is the converse of key b at phase `1-k`: reversing `(c_(1-k),c_(2-k),c_(-k))` and applying `c_j -> c_(-j)` gives the shared phase-k path. Thus its normalized word is rho(w_b), with no invented independent phase.

### Exact local capacity and its global black image

The ordered SCC size profiles of L and its converse are `(8,...,8,3)` and `(3,8,...,8)`, so they are not isomorphic. Every L complement partner is strong, while L is nonstrong. Source words from distinct keys are distinct because the prefix is rigid. Together these three separations prove that the m four-word sets

    {w_b, rho(w_b), F(w_b), F(rho(w_b))}

are pairwise disjoint. This is the authenticated local resource reservoir.

The left record in (2) consumes local resource a; the right consumes local resource b. Hence the two slots satisfy capacity one EXACTLY when a!=b. If a=b, the global observer and both local partners still exist, but the two local slots double-use one normalized resource. Different phases cannot repair that collision: they are kernel lifts of the same resource.

Globally, (2) has normalized word indexed by ordered pair (a,b), and the global resource is indexed by the UNORDERED pair `{a,b}` from §3. Thus every capacity-valid two-local-slot assembly is off-diagonal and globally black. The global resource quotient forgets which local slot holds which key; swapping a and b keeps the global resource but swaps the marked local assignment. There are exactly

    3m(m-1) capacity-valid labelled global records,
    m(m-1) capacity-valid normalized global source words,
    binom(m,2) capacity-valid global black resources.

There are two pin scopes. For a FIXED ordered pair (a,b), pin both private prefixes completely; exactly the three shared phases remain, and only those phases move in the route below. For the UNION over keys, impose instead the six common seed clocks in each left block and their reflected clocks in each right block, with C3 unpinned. Ordinal factorization forces every block to occupy its fixed segment. Each block has only the two unpinned owners0 and1 and exactly the two legal assignments P0,P1 (reversed in the right half); C3 has exactly three paths. Hence the union is also a complete fibre for these weaker common pins, not for the fixed-pair full-prefix pins. No phase route below releases a fixed prefix or moves between different key pairs.

The full selected weaker-common-pin sandwich menu, including local-capacity-invalid diagonals, has 3m^2 full records, m^2 source words, m blue resources and `binom(m,2)` black resources. The complete UNPINNED sandwich has prefix orbit count `409^r` and source fibre `409^(2r)` instead. The selected m must never replace the complete a in a whole-node profile.

This local-resource injection is specific to the native family, or to another family with a separately authenticated injection. It is not proved merely by choosing different path orbits of arbitrary R.

An explicit FINITE-EXACT hostile is `R=C3 -> U`, where U has adjacency rows `0001,1000,1100,0110`. It is the strong four-owner tournament; R has `(H,Aut,a)=(15,3,5)`. The distinct R-path orbits represented by `P=(0,1,2,3,6,5,4)` and `Q=(0,1,2,5,4,3,6)` become reflection-paired after appending the common C3 path `(7,8,9)`. Thus in `L=C3 -> U -> C3` their different local source words consume the SAME F/rho resource. This authenticates the precise failed generalization at order10.

### All phases preserve both local resources and one global resource

The literal permutation `c_j -> c_(j+delta)` fixing all other owners is an automorphism of S_r. It sends (2) at k to (2) at k+delta. Equivariance of K makes it an actual global partner isomorphism, restricting to actual source and partner isomorphisms on both local cells. Every one of the three phases is available; their normalized words and normalized resources do not change.

The two local units lie in the common order-(8r+3) reservoir; the global unit lies in the separate order-(16r+3) carrier. This is not a claim of three units in one capacity pool. Each local resource and the global resource are continuously held, with no release or later reacquisition. The three lifts supply no extra capacity. Any target-class-dependent per-slot character and either the ordinary or converse-merged CLASS-incidence vector is constant separately on each local cell and on the global record. No cancellation between cells is required for these particular class-incidence differences, which are each zero.

This does NOT authenticate arbitrary owner-sensitive current rows, target-side attachments, or a harmonic bank on a different-order metagraph. Ordinary local source classes on opposite cells are converses and need not be equal; their own phase-invariance is the claim. The local carrier order is 8r+3, while the global carrier order is 16r+3. Equal endpoint classes within a carrier do not transport numerical harmonics between these orders.

### Two integer work counters, not one

For any nontrivial phase change, each of the three middle owners changes its indexed clock in the global path. Therefore global owner-clock work is exactly3. The left local record includes those three occurrences at local clocks `q,q+1,q+2`, and the right includes them at clocks `0,1,2`. The sum of per-local-slot owner-clock occurrence work is exactly6. These are two different declared counters. Collapsing the common physical owner set does not make six equal to three; it changes which occurrences the counter counts. Both remain attached to the full records.

The accepted phase graph is the three-state complete graph, with local-work edge cost6 and global-work edge cost3. It is a genuinely positive counterpart to the earlier two-outgoing-arc global-observer hostile, not a route around that hostile under its unchanged constraints.

## 5. A common target: many resources are not many neighbors

The endpoint-fixed seed maps compose blockwise across ordinal sums. A change from key0 to any a uses g in each changed left block. A change from key0 to any b uses the SAME g on each changed right-copy block, now between converse partners. Both actual endpoints are preserved. A central phase change uses its literal C3 rotation. Between different blocks all partner arcs reverse except the last-to-first observer junctions; the marked endpoint conditions preserve exactly those exceptions.

Consequently an explicit product permutation carries `K_(O_(0,0,0))(S_r)` to `K_(O_(a,b,k))(S_r)` for EVERY selected a,b,k. All m blue resources and all `binom(m,2)` black resources in the selected pinned sandwich menu reach ONE neighboring isomorphism class. That target is self-converse because a diagonal presentation is blue and F commutes with rho. It is strong by §3, hence different from the source.

This is a common-target theorem for the selected pinned menu only. It does not assert that the entire `409^(2r)`-word source fibre has one neighbor. It gives an explicit all-size demonstration that even the positive global-observer construction must count resources separately from distinct neighbors.

## 6. Complete small profiles and Paley11 control

The checker computes Hamiltonian paths by DFS, independently recomputes H by subset DP, enumerates full automorphism groups by adjacency-pruned backtracking (no assumed circulant group), and obtains every path orbit by its normalized word. Complete R-orbit representatives then give every S source word without enumerating all tournaments of order |S|.

- R=C3: `(H,Aut,a)=(3,3,1)`; S has order9, one blue incidence, zero black incidence, one resource.
- R=cyclic C5 with forward differences{1,2}: `(15,5,3)`; S has order13, blue3, black weighted6, resources6.
- R=Paley7 with forward differences{1,2,4}: `(189,21,9)`; S has order17, blue9, black weighted72, resources45.
- R=Paley11 with quadratic-residue differences{1,3,4,5,9}: `(95095,55,1729)`. The 1729 distinct normalized R words are computed exhaustively. The theorem then gives its order25 sandwich 2,989,441 source words, 1,729 blue incidences, 2,987,712 black weighted incidences and 1,495,585 resources. Its Hamiltonian observer has24 path arcs.

The sequence1,3,9 does not continue by powers of3: the next listed a is1729, not81. These values and that failure were already recorded in THM-2454's historical §4; they are independently reproduced here, not newly discovered numerology. The 1729 and24 occurrence is the stated concrete construction, not a Leech-lattice, string-theoretic or modular identity. The large Paley11 S fibre is derived from the proved bijection, not claimed to have been enumerated pair by pair.

The independent finite controls also cover every one of the75 labelled nonempty tournaments R of orders1..4. For each, direct whole-S path and automorphism enumeration checks the factorization and exhausts the normalized S fibre. Native tests r=1..4 exhaust the selected pinned menus, all three phases, both local partners, actual source anti-isomorphism, common-target maps, both resource counters, converse phase relation and work3/work6. All checks use explicit failures, so normal and `-O` runs have identical force.

## 7. Scope and remaining boundary

The map preserves actual owner labels, contiguous local observer edges and nonedges, normalized local capacity, selected global resource identity, the full kernel phase, true endpoint transports and both declared integer work counters. It forgets the ordered key assignment only when taking the global resource quotient, a loss recorded explicitly above. It preserves class-incidence values within each existing carrier, not arbitrary owner-sensitive electrical constraints or a cross-order harmonic bank.

No arbitrary overlapping-cell family is claimed to admit such a sandwich embedding. The directed-block compiler and its negative flags remain the correct general gate. Here the overlap is one actual contiguous C3 segment, with a specified anti-isomorphism and off-diagonal local-resource condition; those are the positive construction's exact hypotheses.
