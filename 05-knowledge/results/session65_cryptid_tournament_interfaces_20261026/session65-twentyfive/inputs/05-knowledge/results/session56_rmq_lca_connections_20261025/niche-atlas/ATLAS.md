# Session56: niche connections with RMQ/LCA and current gluing work

## Source pin and status

Historical repository pin: `4f1760fd1eff346767c032524c932d5e317401d2`, inspected through `session50-publication/clean-repo`. Paths below are relative to that commit unless explicitly described as local Session55 packets. Read AGENTS, START-HERE and active guardrails; targeted theorem files, not chronological logs, supplied the inheritance. No literature novelty claim. `PROVED (inherited)` describes the cited packet's scoped status, not a fresh full audit. New elementary deductions are proved here. `OPEN application` means the proposed use has not been wired into the actual research implementation. No all-size native tournament theorem, LRC closure, integer magic-square existence result, or Collatz termination claim follows.

The closest current mechanism is Session55's exact CRT meet convolution on named equality partitions. The canonical hostile is the aligned-versus-crossed three-block pair A,A versus A,B: identical block sizes, different global distinctness. The corrected near miss is that existence-antichain compression preserves feasibility but destroys exact counts. The least-used sidecars recovered here are the actual address/scale map, target-conditioned continuation state, and higher-order interaction data.

Live concepts: ordered extrema; named-port equality; hierarchical versus crossing partitions; address/valuation scale; continuation kernels; higher-arity defects.

## 1. Stable Cartesian trees preserve a witness, not a numerical minimum

**Source:** `01-canon/theorems/THM-2253-online-dyadic-contrast-tournament-extractor.md`, §§1–3, retains the actual first terminal dyadic node; `05-knowledge/results/square_sum_ports_20261023/README.md` retains named exposed endpoints rather than only their count.

**Target/map:** For any finite ordered array a, recursively choose the least index minimizing `(a_i,i)` as root, with the left/right subarrays as children. Map interval [l,r] to the LCA of vertices l,r. This is its leftmost minimum position: if the recursive root lies outside the interval, recurse into the unique child containing both endpoints; if inside, it is both the minimum and their LCA. This proves the identity for every finite array, including ties.

**Preserved by the Cartesian/LCA map alone:** the chosen argmin owner index, with the declared leftmost tie rule. The minimum value is recovered only by looking up that index in the original array or retained key sidecar. **Lost without sidecars:** numerical values, all other tied owners, and actual chronological admissibility of the terminal event. Store the original array, owner sets if needed, and the admissibility predicate.

**Hostile:** [1,1] has two minimum owners; [1,2] has the same leftmost-minimum Cartesian shape. Tree shape alone does not preserve tie multiplicity. Likewise [1,1] and [5,5] have identical Cartesian trees and argmin answers but different minimum values. **Cheapest test:** all ternary arrays of length at most seven and all intervals, plus [1,1]/[1,2] and [1,1]/[5,5]. **Status:** PROVED elementary RMQ/LCA identity; FINITE-EXACT test; OPEN application to current packets. Nothing says a dyadic terminal-event query already is this RMQ.

## 2. CRT collision meet is an idempotent range-query algebra

**Source:** local `session55-crt-collision-wiring/THEOREM.md` §§2–3; local `session55-weighted-crt-counts/THEOREM.md` §§1–3.

**Target/map:** Given a FIXED sequence of selected local arrays, store their 36-bit equality masks C_i. A factor interval [l,r] maps to AND of C_l,...,C_r. Associativity composes adjacent blocks. Idempotence allows two overlapping power-of-two blocks covering the interval: repeated factors do not change the mask. This is the exact algebraic reason the overlapping-block RMQ strategy transfers to mask intersections.

**Preserved:** surviving collisions of that selected tuple, including its zero-mask distinctness predicate. **Lost:** how many tuples yield the mask and which array choices were available. Restore weighted histograms for counting and witness pointers for reconstruction. This is not an algorithm for selecting good local arrays merely by inspecting one tuple.

**Hostile:** existence-neutral modulus 3 multiplies counts by two in the Session55 counting theorem. An overlapping-block product of zeta counts repeats factors and therefore overcounts; idempotence is absent. Even sum(1,2,3)=6 differs from sum(1,2)+sum(2,3)=8. **Cheapest test:** exhaustive small-mask intervals and a duplicated count factor. **Status:** PROVED elementary range identity; OPEN implementation; not a complexity claim about the full Bell(9) search.

## 3. Partition lattices are not secretly LCA trees

**Source:** `01-canon/theorems/THM-2330-partition-lattice-gordian-lift-spectrum.md`, status and partition-coarsening lift theorem; Session55 weighted theorem §§1–2.

**Target/map proposed:** encode the refinement order of equality partitions as ancestry so that meets/joins become LCA-like queries. **Preserved if restricted:** a truly hierarchical family of clusters with verified nesting can be stored as a rooted tree. **Destroyed by general tree replacement:** crossing compatibility states.

**Hostile/proof of obstruction:** on labels {0,1,2}, the discrete partition is below both 01|2 and 02|1 in refinement order, but these two partitions are incomparable. Under an order embedding into ancestry with coarser states as ancestors, they would be incomparable ancestors of one node, impossible in a rooted tree. Reversing the order does not rescue the whole lattice: the top partition then supplies the same obstruction. Thus there is no faithful ancestry encoding of the entire partition order.

**Sidecar:** retain the lattice/DAG, or prove that the actual admissible subset is a hierarchy before using LCA. **Cheapest test:** run this three-label obstruction and a laminarity test on the proposed cluster family. **Status:** PROVED no-go for faithful whole-order encoding; OPEN restricted-family opportunity. This does not forbid using a tree as a computational expression tree for lattice operations.

## 4. Persistent versions keep branch provenance; they do not supply compatible ports

**Source:** `05-knowledge/results/session25_spectrum_wiring_20261021/session25-ports/PORT_SIGNATURES.md`, §§1–3; Session55 weighted theorem §§2–3.

**Target/map:** represent each arithmetic or boundary-gluing branch by an immutable version root. Each node stores the full named-port state and a proof/witness pointer; updates create a new path of nodes rather than mutating an ancestor. Equal states can be shared only after the continuation semantics are fixed.

**Preserved:** prior branch state, rollback provenance, and exact named compatibility when stored in the payload. **Lost by the version tree alone:** all mathematical contents, multiplicities, and attachment names. The LCA of two version roots says where their edit histories diverged, not whether their arithmetic witnesses can glue.

**Hostile:** on three named boundary vertices take G=edge12 and H=edge13, anonymous-isomorphic fragments. Attach K=edge12. Their k-coloring counts become k²(k−1) and k(k−1)². A persistent tree containing only anonymous shapes still returns the wrong compatibility result. **Sidecar:** fixed boundary labeling and the full partition-indexed extension vector. **Cheapest test:** k=2,...,5 on this three-vertex hostile; branch once and verify the old payload remains unchanged. **Status:** inherited PROVED signature law; OPEN engineering design, no implemented persistent data structure or benchmark claimed.

## 5. Two-adic ultrametrics become LCA only in the correct digit order and scale

**Source:** `01-canon/theorems/THM-3503-rule30-odometer-ultrametric-regrading-and-orbit-closure-dimensions.md`, Theorem 2.1, equations (9)–(10): ν2(iota(t)−iota(s))=v_m when m=ν2(t−s).

**Target/map:** for a bounded quotient, insert phase digits into a trie LOW BIT FIRST. For distinct phases x,y, their common-prefix length is ν2(x−y), since equality modulo 2^m is exactly equality of the first m low bits. The LCA depth gives phase distance; the source theorem's regrading m↦v_m then gives the ambient Rule30 state distance 2^(−v_m).

**Preserved:** finite quotient phase-prefix depth and, with the innovation sidecar, the inherited ambient ultrametric. **Lost:** physical innovation depth if only trie depth is kept; topological conjugacy alone does not preserve dimension or metric scale.

**Hostile:** 1 and 5 share two low bits but not a leading-zero high-bit prefix of that length in a three-bit representation. A sparse innovation sequence can preserve the phase tree while radically changing ambient scale. **Cheapest test:** all phase pairs below 64, then check actual innovation data rather than assuming v_m=m. **Status:** inherited PROVED regrading; elementary PROVED low-bit-prefix lemma; OPEN data-structure use. No new Rule30 dimension lower bound.

## 6. Fibonacci/Berggren addresses warn against treating tree depth as arithmetic size

**Source:** `01-canon/theorems/THM-3382-fibonacci-ray-dual-index-harmonic-bifurcation-and-ternary-heap-addresses.md`, status and connection contract; `THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md`.

**Target/map:** retain the same ray parameter r but compare Fibonacci-time address 3r+4 with canonical ternary-heap address (33·9^r−9)/8. An LCA or trie can organize common prefixes of branch words, but any density, harmonic mass, or height budget must use the declared arithmetic embedding.

**Preserved:** branch ancestry and ray membership. **Lost:** counting density, harmonic class, primitive-content phase, and affine owner if only depth is kept. **Sidecar:** exact injection into N, branch label and normalization/owner current.

**Hostile:** the source proves divergent harmonic mass and density 1/3 for {3r+4}, but convergent harmonic mass and density zero for the heap addresses of the SAME ray objects. Depth also identifies two distinct odd-depth rays. **Cheapest test:** print paired addresses and the first split between equal-depth branches; use the explicit formula to prove geometric versus linear growth, not a finite sum as evidence of convergence. **Status:** inherited PROVED; OPEN use as a guardrail in indexed searches. Shared Fibonacci words do not prove a Collatz or tournament bridge.

## 7. Zeckendorf carry is a missing update coordinate, not a tree-address decoration

**Source:** `05-knowledge/results/duck_zeckendorf_20260925.md`, equations Z8–Z11 and §7.

**Target/map:** use canonical Fibonacci supports for integer labels; store the source's two-coordinate charge Q and carry δ(a,b). The exact update is Q(a)+Q(b)−Q(a+b)=(−2δ,δ), and the color update includes the bit δ mod2. A prefix/sum data structure must aggregate this lawful operation or retain full labels, not XOR colors blindly.

**Preserved:** sum color with carry; full charge retains integer n=A+2B. **Lost:** integer magnitude and squarehood in the color-only quotient; support overlap is not itself a complete carry classifier.

**Hostile:** distinct integers 1+4 and 2+3 both sum to 5, but their naive color XORs differ because 1 and 4 share a Fibonacci atom. Repeated 1+1 also differs from the singleton representation of 2. **Sidecar:** δ or the exact charge/label. **Cheapest test:** those two decompositions and the source's carry cocycle before designing a persistent update. **Status:** inherited PROVED arithmetic decoder; OPEN data-structure application. No square-sum Hamiltonicity follows from the color.

## 8. Collatz binary lifting must retain the affine itinerary and its admissibility

**Source:** `01-canon/theorems/THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md`, Setting and statements 1–2. The source explicitly separates elementary equivariance from FINITE-EXACT bijectivity through period 22.

**Target/map:** on a fixed finite parity itinerary of the standard map, each step is x↦3x+1 or x↦x/2. Compose affine maps by (a,b) after (c,d) ↦ (ac,ad+b). Doubling blocks therefore computes exact rational affine formulas for the chosen itinerary. Store the residue cylinder on which that itinerary is actually followed; its intersection under composition is essential.

**Preserved:** finite orbit values on that cylinder and the step count. **Lost:** admissibility if one keeps only the affine formula; actual time if one replaces a run of even steps by an accelerated edge. **Hostile:** odd inputs 1 and 5 agree modulo 4 but ν2(3x+1) is 2 versus 4, so coarse residue agreement does not fix accelerated edge length. The source's even-period golden reader also collapses 0,−1,−2 to zero, warning that an address quotient can lose orbit identity.

**Cheapest test:** finite exact compositions for inputs 1..100 and 12 standard steps; explicit valuation pair 1/5. **Status:** PROVED elementary composition; inherited bounded golden dictionary only; OPEN fast-query implementation. No all-integer termination conclusion and no universal odd-period bijection promoted from bounded evidence.

## 9. Tropical range composition requires disjoint blocks, unlike scalar RMQ

**Source:** `01-canon/theorems/THM-2248-higher-interaction-defect-complex-and-tropical-trace-spectrum.md`, §§4–6; `THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder.md`, equation (7).

**Target/map:** a finite state transition kernel K is composed using (K⊗L)ij=min_k(Kik+Lkj). A segment tree can store interval kernels using associative composition. The elementary proof expands both bracketings into min over the same pair of intermediate states. Source THM-2897 instead minimizes sums of primitive cap costs over decompositions; its F_p<h remains only a sufficient no-cover gate.

**Preserved:** finite-state optimal costs when the state boundary and kernels are exact. **Lost:** minimizing witnesses and multiplicities under scalar minimum; the original complete infinite continuation kernel if one truncates its state space.

**Hostile:** the 2×2 all-ones kernel squares to the all-twos kernel. Overlapping blocks therefore count repeated transitions and are not legal merely because scalar RMQ overlaps. **Sidecar:** disjoint interval decomposition, state labels, argmin witness, and any truncation certificate. **Cheapest test:** all triples of 2×2 Boolean-entry matrices for associativity, plus all-ones nonidempotence. **Status:** PROVED elementary algebra; OPEN application; no automatic finite-state reduction for Gordian or tournament continuations.

## 10. Higher-arity interaction obstructs pairwise-tournament compression

**Source:** `01-canon/theorems/THM-2248-higher-interaction-defect-complex-and-tropical-trace-spectrum.md`, §3; `THM-2330-partition-lattice-gordian-lift-spectrum.md`, status and fixed-arity obstruction.

**Target/map:** compress a packet to pairwise interaction or an orientation, then build a tree of pairwise joins. This can be sound only if the summary contains every future composite-context interaction, not merely original vertex pairs.

**Preserved:** tested binary interactions. **Lost:** irreducible higher-order savings. **Hostile:** in Z^r give generators ±e_i cost 1 and ±(e_1+...+e_r) cost r−1. Every proper subset of the basis has additive length |S|, but the full sum has length r−1. For r=3 every original pair is additive and the triple is not. More generally the cost of a size-s subset is min over k in Z of |k|(r−1)+s|1−k|+(r−s)|k|, exactly s for s<r and r−1 for s=r.

**Sidecar:** higher-defect complex or full target-conditioned continuation kernels. **Cheapest test:** r=3 first, then r=3..8. **Status:** inherited PROVED abstract metric-monoid hostile; FINITE-EXACT controls. It is NOT an exhibited knot symbiont of higher arity, and it does not refute a specific native tournament theorem without identifying its carrier.

## Priorities and stopping rules

1. Implement card 2 only for fixed selected masks, and card 9 only for exact finite kernels. Stop if the target query needs multiplicities or unrepresented boundary states.
2. Pursue card 5 if the actual task has phase-prefix queries; retain the innovation regrading explicitly.
3. Before proposing an LCA reduction of collision signatures, run card 3. If the family crosses, keep the partition DAG rather than inventing a hierarchy.
4. Use cards 4,7,8 as state-design tests: branch history, carries, and itinerary cylinders are separate coordinates.
5. Keep cards 6 and 10 as hostile controls whenever a tree or tournament quotient claims to preserve global arithmetic consequences.

## Reproduction and scope

`python3 check_connections.py` and `python3 -O check_connections.py` produce identical JSON. Declared finite universes are in that output. Tests independently reproduce the elementary controls, not every inherited theorem. No code imports another lane's implementation. Initial writes in the shared workspace failed because its filesystem was full; the deliverable was rebuilt under `/tmp/session56-niche-atlas`, with the failed empty workspace files not treated as certificates. No remote writes or changes to existing packets were performed.

