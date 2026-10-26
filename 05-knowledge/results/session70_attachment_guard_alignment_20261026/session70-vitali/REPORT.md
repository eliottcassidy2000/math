# Native re-marking guards realize the full Boolean moment hierarchy

2026-10-26. Session70 recovery/extension lane. **PROVED:** A1–A4 below from explicit elementary constructions. **FINITE-EXACT:** the independent checker covers every constructed alignment and every simultaneous-subset guard for `k=1,...,7`; normal and `-O` outputs agree. **INHERITED:** the paired-rooted collision and its C4 construction, Session69 P3 recursive response closure, THM-840's general rank criterion, and the abstract parity-law obstruction. No literature priority, universal prime-tournament classification, or classical nonmeasurability is claimed.

## Immediate results

1. The inherited `2^k` alignment fibre at `n=7k+2` lives in a **rigid strong source tournament for every k**. Growing alignment memory is not caused by nontrivial source automorphisms or failure to choose an equivariant path.
2. A single legal re-marking guard can test any conjunction of the k hidden bits. The complete guard matrix is the Boolean zeta matrix, with determinant one. Its rank is exactly `2^k`, over every field.
3. Retaining all named guard moments through block order t has rank `d_t=sum_{j=0}^t binom(k,j)`. Its exact missing linear sidecar rank for all simultaneous-subset guard responses is `2^k-d_t`. On probability laws, full alignment-law dimension is `2^k-1`; this is not a lower bound on nonlinear encodings of individual presentations.
4. Uniform even/odd alignment-parity laws have the same every-proper-block marginal, yet the actual k-block re-marking guard has probabilities `0` and `2^(1-k)`. Thus every fixed locality order fails at a larger native tournament scale. Full named moments repair the loss by finite Boolean inversion.

The new part is **native realization and exact size of the moment obstruction inside tournament re-marking**. Abstract parity distributions and linear-kernel rank were already proved in the repo; their previous uninstantiated hypotheses must not be presented as new theorems.

## 1. Source pins, inheritance, and bounded recovery

The native construction and historical bytes were read at git commit
`fb872272125510c317b8114dd3f697b433c595f9`.
`SOURCES.json` records exact paths, git blob IDs, SHA-256 hashes, and local snapshots. The current worktree is sparse; `git show PIN:path` was used for missing physical paths, not inferred absence.

The closest mechanism is C4–C5 of
`05-knowledge/results/strong_atom_continuation_vitali_20261026/README.md`.
C4 already proves the `2^k` paired-rooted alignment fibre and k-bit deterministic lower bound. C5 already defines the coarsest guarded continuation quotient. Session69 P3 in
`05-knowledge/results/session69_prime_ports_power_aliquot_20261026/session69-tournament/REPORT.md`
already closes the **ordinary uncoloured all-arity ordered path-cover response recursion**. The current lane does not reopen that result. The unresolved compression is of aligned individual choices, coloured continuation, or economically represented responses under the required operation alphabet.

The six-concept board was: observation fibre; named versus averaged event; finite moment order; relative frame; equivariant selector; recursive attachment. The new native guard conjunction changes each as follows: one fibre already contains a Boolean cube; averaging needs its whole joint law; correlation order grows with block count; endpoint frames do not identify that cube; source rigidity rules out blaming symmetry; uniform substitution physically realizes every conjunction while P3 still governs uncoloured aggregate incidence.

Additional targeted recovery beyond the immediately preceding HISTORY ledger:

- **THM-840**, `01-canon/theorems/THM-840-hamming-five-continuation-congruence-boundary.md`, §1a: exact minimum linear sidecar rank `dim N(ker R)`. This is stronger and more useful here than merely citing its set-theoretic operation-kernel criterion.
- **THM-2237**, `01-canon/theorems/THM-2237-truncated-boolean-moment-interval-and-parity-top-atom-majorants.md`, §3: both parity laws have equal proper Boolean marginals. That result warns that native realization must be established separately. A2 supplies that realization.
- **THM-561**, `01-canon/theorems/THM-561-miss-zeta-factorial-atom-inversion.md`: rank-summed zeta coordinates recover the count law, not named configurations. The finite named Boolean inverse below is the required stronger carrier.
- `05-knowledge/results/hamiltonian_observer_restrictions_20261023/compiler/README.md`, §6, and `history/AUDIT_COMPILER.md`: no equivariant selector is an **unpointed-input** theorem. Retaining the reference path or an external label order already supplies a gauge. The audit corrected this exact scope error before publication.
- **THM-796**, `01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md`, §8: the node deletion kernel first fails at n=5; at n=7 there are 1,206 non-lumpable nonzero blocks among 1,312. Direct two-step deletion differs from the normalized Markov product in 1,778 entries and 270/272 rows, max error 1/2. These are inherited finite counts, not rerun here.
- `05-knowledge/results/three_sorted_metagraph_continuation_minimization_codex_S11.out`: the older operation alphabet is top/bottom deletion plus reflection, not current arbitrary guarded re-marking. At n=7, labelled continuation of the bare boundary has 16,359 cells on 16,384 lines; joint Xi has 16,382, while retained coherent phase gives all 16,384. These counts must not be transferred to another alphabet. Its residual relative holonomy has order five; our seed mismatch below has order two, so even the holonomy order is carrier-specific.

`agents/processor.py --status --peek --limit 8` could not read the coordination inbox because `.machine-id` is absent. No registration was invented; no mathematical step depends on that inbox.

### Incoming fresh re-marking work: different operation menu

During this lane, incoming main `7d5c27c1f53ef90f2da1efad4793b90f53a75cd5` added `05-knowledge/results/fresh_remarking_quotient_and_divisor_arrangements_20261026.md` and **HYP-9416**, `05-knowledge/hypotheses/HYP-9416-fresh-re-marking-two-port-quotient-joint-classes-are-minimal-and-line-responses-classify-strong-atoms.md`. Both were read and separately pinned in SOURCES.json. Their all-order assertions remain OPEN; their through-seven tables are explicitly not independently audited by the incoming packet and are not rerun here.

V1 keeps the multiset of complemented successors over all s-to-t paths; V2 allows all paths and resets ports. These operations deliberately forget which named permutation produced a successor. Our rank bound instead retains each individual named guard expectation for arbitrary laws on aligned states. Thus finite fresh-chain compression does not contradict this result, and this result does not refute that compression. Aggregating operator names is a change of experiment, not an exact repair for all the original named experiments. Also, C4 already applies to a **single global marked Hamilton path** through k blocks: the incoming phrase “beyond one run / multi-run alignments” must not narrow its actual scope.

## 2. Exact native construction

Write a tournament by outgoing-neighbour bit rows. Take

```
X=(2,12,41,17,39,75,31),
Y=(2,20,9,51,37,71,31),
p=(0,1,2,3,5,6,4).
```

Both have the standard marked path `0->1->...->6`. Let tau reverse every arc except these standard path edges. The maps, interpreted as old vertex to new vertex, are

```
s=(0,1,4,2,3,5,6): X -> Y,
t=(0,2,4,1,5,3,6): tau(X) -> tau(Y).
```

They fix both endpoints. Both presentations are black. The permutation p is not a directed path in X (`3->5` is absent), and is one in Y. These are the inherited C3 seed, not a new minimality claim.

For `b=(b_1,...,b_k) in {0,1}^k`, let `S_b` be the ordered sum of k blocks, block i being X or Y according to `b_i`. Add u,v with

```
S_b -> u -> v -> S_b.
```

Call this tournament `T_b`, with vertices in block order then u,v. Every standard path is legal. The outer directed triangle makes T_b strong. Extending s and t independently inside each chosen block gives the same rooted source type and the same rooted tau-target type for every b; target inter-block arcs are reversed except for the fixed endpoint hinges. Source maps and target maps preserve those hinge endpoints. This is precisely C4's fibre at `n=7k+2`.

All presentations are black: the standard first vertex has degree `7k-5` and last vertex v has degree `7k`. If converse plus positional reversal fixed the tournament, these two degrees would sum to `n-1=7k+1`; instead their sum is `14k-5`, unequal for every positive integer k.

### A1. The growing fibre has rigid sources

**PROVED.** Every source T_b has trivial automorphism group.

The degree sequence of X, in vertex order, is `(1,2,3,2,4,4,5)`. Thus vertices 0,2,6 are fixed by every automorphism. The degree-two vertices 1,3 are distinguished by their opposite orientations to fixed vertex 0. The degree-four vertices 4,5 are distinguished by their opposite orientations to fixed vertex 2. Hence X is rigid; so is the isomorphic Y. X is strong because it contains `0->1->...->6->0`.

In T_b, u has degree one, v has degree 7k, and every block vertex has degree between 2 and `7k-1`. Thus u and v are uniquely identified and fixed. The remaining induced tournament S_b has k strong components, ordered linearly. Every automorphism preserves that order and therefore every component, whose rigidity finishes the proof.

For additional seed control, tau(X) is also rigid: its degree sequence is `(6,4,3,4,2,2,0)`; 0,2,6 are fixed, and each repeated-degree pair is distinguished by its orientation to vertex 2. The two unique seed isomorphisms s,t differ by

```
t^-1 s=(1 3)(4 5).
```

This permutation is not an automorphism of X. It is a mismatch between two identifications, not residual source symmetry. Consequently the unbounded alignment obstruction and the finite no-equivariant-selector obstruction are logically distinct.

## 3. A2: every Boolean conjunction is an actual one-step re-marking guard

For `A subset [k]`, define the named permutation `P_A` as follows. Read block i in order p if `i in A`, and in its standard order otherwise; keep blocks in their original order and finish with u,v. This is a permutation of all n vertices, independently of b. Let

```
g_A(b)=1 iff P_A is a directed Hamiltonian path of T_b.
```

**PROVED, all k and A:**

```
g_A(b)=product_{i in A} b_i = 1_{A subset supp(b)}.       (A2)
```

Within a selected block the p-guard succeeds exactly for Y. Within an unselected block the standard path always succeeds. Every cross-block boundary goes forward in an ordered sum, regardless of the chosen internal endpoint, and all final wrapper transitions are legal. Thus precisely the selected block conditions remain. The empty A is the always-legal standard path.

These are **simultaneous one-step** re-markings, not a tacit multiplication of success probabilities of a sequential process. Re-marking changes coordinates, so treating several sequential guards as independent without transport would be invalid. The literal P_A construction avoids that issue.

Let `Z[A,B]=1_{A subset B}` index guards by A and states by `B=supp(b)`. In numerical subset-mask order, Z is upper triangular with diagonal one. Hence `det Z=1` and Z is invertible over every field. In particular the native guard observables span every function on this entire `2^k`-state fibre.

**Deterministic scope:** the k singleton guards already recover b, so the minimal number of extra binary labels on this fixed coarse fibre is k, as inherited C4 proved. The number `2^k` below is a dimension of a linear function/distribution representation, not a need for `2^k` deterministic memory bits.

## 4. A3: exact linear memory and positive full-law carrier

Take the vector space V of signed weights on the `2^k` alignments over a field F. Define

```
N mu = (sum_B Z[A,B] mu(B))_{A subset [k]},
R_t mu = (N mu)_A for |A|<=t,       0<=t<=k.
```

Thus R_t retains every **named** guard expectation up to t selected blocks, including total mass at A empty. All subsets of rows of an invertible matrix are independent. Therefore

```
rank R_t = d_t = sum_{j=0}^t binom(k,j),
dim N(ker R_t) = 2^k-d_t.                                (A3)
```

By the elementary rank criterion of THM-840 §1a, the minimum rank of any **linear** sidecar S permitting exact recovery of all N from `(R_t,S)` is exactly `2^k-d_t`. This can also be proved directly: N is injective, so `(R_t,S)` must be injective; rank subadditivity gives the lower bound, and appending the missing zeta rows attains it.

For probability laws over the reals, mass is already one. The normalized fibre has affine dimension `2^k-1`; after R_t it has generic remaining affine dimension `2^k-d_t`. Strictly positive laws admit small perturbations in every kernel direction, so positivity does not erase this dimension bound. The statement is not about arbitrary nonlinear, discontinuous, or unbounded-precision encodings.

The complete **positive minimal law carrier on this family** is equivalently the alignment distribution `mu(B)` or all named moments `m(A)=Nmu(A)`. The decoder is finite Boolean inversion:

```
mu(B)=sum_{A superset B} (-1)^(|A|-|B|) m(A).             (A4)
```

An arbitrary vector m is the moment vector of a probability law iff `m(empty)=1` and all inverse masses in (A4) are nonnegative. This exact feasibility guard prevents treating arbitrary moment-compatible formal records as native laws. Here every B has already been physically realized by T_b, so these conditions are also sufficient for a mixture on the constructed family.

If every guard is to be recovered as a linear functional of a d-dimensional state representation `v_b`, then `Z=L V` gives `d>=rank Z=2^k`. Allowing an affine constant instead needs at least `2^k-1` nonconstant coordinates. This is a representation bound, not an algorithmic complexity or entropy theorem without a probability model.

## 5. A4: proper-locality statistics miss a native global guard

Let mu_0 and mu_1 be uniform probability laws on bitwords with even and odd total parity. For every proper subset A,

```
E_mu0 g_A = E_mu1 g_A = 2^(-|A|).
```

More strongly, every proper-block marginal is the same uniform law. But

```
E_mu_(k mod2) g_[k] = 2^(1-k),
E_mu_(1-k mod2) g_[k] = 0.                              (A5)
```

The proof is direct: prescribing a proper set of coordinates leaves at least one free bit, and exactly half the completions have each parity. The all-one word has only one parity. THM-2237 proves this abstract cube fact; A2 proves that its missing top event is an actual tournament re-marking guard, not an invented observer.

Thus for every fixed t, choose k=t+1. Two mixtures on one paired-rooted strong-source fibre agree under every test depending on at most t named blocks, but disagree on a lawful one-step full-block guard. In particular all block-preserving re-marking permutations that change at most t blocks are indistinguishable in law. This does **not** assert blindness of all arbitrary re-marking permutations: a global permutation may access all k bits and is precisely the allowed distinguisher.

The gap `2^(1-k)` tends to zero; the theorem concerns exact observability, not a uniform approximation-hardness bound. Scaling an observable is not automatically a physical probability operation.

## 6. Finite section versus genuine nonmeasurable selector

Three distinct mechanisms must remain separate.

1. **Finite observation fibre.** Every subset of a finite tournament space is measurable for the discrete sigma-algebra. A target varying on a lambda or endpoint fibre fails to be a function of that observation. It is not thereby non-Lebesgue-measurable. Current THM-169 still has its repaired finite-evidence status, not its former all-order proved claim.
2. **Finite equivariant section.** An equivariant surjection `pi:E->B` of G-sets has a G-equivariant section iff each base orbit has a point b whose stabilizer G_b fixes some point of `pi^-1(b)`. Necessity follows from equivariance; sufficiency follows by choosing such an e at an orbit representative and defining `s(gb)=ge`, which is well-defined because G_b fixes e. For the forgetful map `(T,P)->T`, automorphisms act freely on full ordered spanning paths. Hence a nontrivial Aut(T) forbids an unpointed equivariant single-path selection, but a retained reference path or label gauge changes the input and can select. Our A1 sources have trivial Aut, yet still need growing alignment memory: the two obstructions cannot be identified.
3. **Infinite measure-theoretic selection.** A dense rational-translation transversal, or the separately proved Haar-nonmeasurable E0 parity-tail transversal, needs an infinite invariant measure and a countable group/tail action. Those hypotheses are absent here. The lattice quotient C/Z[i] has a measurable half-open square section. The historical `the-vitali-atom.md` correction banner explicitly repairs that conflation.

The rigorous common pattern is **forgetting a coordinate that a declared target or continuation still observes**. Nonmeasurability, equivariance, deterministic state size, and linear moment rank are four different predicates with different assumptions.

## 7. Connection contract and surviving frontier

| Source -> target | Actual map | Preserved | Lost on compression | Sidecar / decisive control |
|---|---|---|---|---|
| C4 alignment b -> marked tournament T_b | ordered sum and C3 wrapper | both rooted endpoint types, black colour, strongness; now rigidity | named relative alignment | singleton p-guards read all k bits |
| Native guarded re-marking -> Boolean events | P_A -> g_A | actual legality equals subset containment | averaging/truncation discards correlations | complete Boolean zeta matrix, determinant1 |
| Alignment law -> order-t responses | mu -> R_t mu | every t-local named expectation | `2^k-d_t` linear directions | missing moments, exact inverse and nonnegative masses |
| Abstract parity law -> tournament-law hostile | push law through b -> T_b | proper marginal equality and actual global failure gap | no physical information is invented | explicit P_[k] reads top event |
| Unpointed path fibre -> source | forget P | tournament class | one chosen physical path | stabilizer fixed-point criterion, not Vitali measure theory |
| Ordinary response recursion -> individual continuation | forget individual alignment inside aggregate coefficients | Session69 P3 incidence weights | named state/guard outcome | admissible simultaneous source-target alignment, inherited C1 |

For the full operation menu, the inherited minimal set-theoretic object remains **context-equivalence under every named guarded continuation**, and the aligned source-target forest is a sufficient closed carrier. This report does not claim to solve economical representation on arbitrary primes or prove all-order guard reconstruction of every tournament; the parallel guards lane addresses that question. What is closed here is a full native Boolean subfamily and its exact deterministic/linear/stochastic cost.

## 8. Replay and evidence

```
python3 -B check_alignment_moments.py > check.normal.json
python3 -B -O check_alignment_moments.py > check.optimized.json
cmp check.normal.json check.optimized.json
```

No producer module is imported. The checker builds directed edge sets, tests the supplied seed maps, enumerates every 7! seed automorphism twice, builds all native sources for `k=1,...,7`, checks strongness, wrapper degree separation, black colour, and both rooted identifications, and evaluates every one of the `sum_k 4^k=21,844` native guards. It independently checks GF(2) ranks, every truncation rank, every proper full marginal, direct matrix moments against a fast superset transform, and inversion. The determinant-one proof supplies the every-field statement; GF(2) computation is an independent control, not its only justification.

Positive controls: full named inverse and complete guard distinguishers. Hostiles: top-deleted parity responses collide; a global guard repairs them; the directed triangle has three paths and no rotation-fixed one, but a retained reference is a lawful selection. Boundary k=0 is the one-state constant-guard case, with no parity pair.

Hashes and current check totals appear in the final `CHECKS.md`/`SHA256SUMS`. Mathematical files are outside the shared worktree; parent owns integration and publication.
