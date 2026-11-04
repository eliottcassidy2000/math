# Session 176: information that a tournament circuit must retain

**Five independently reviewed packets establish exact interfaces, constructive edge control, and resource limitations.** They do not settle the Hamilton-path spectrum conjecture or prove a new unrestricted Boolean-circuit lower bound. The organizing question is operational: which information must survive a proposed composition, and what does retaining it actually cost?

## 1. The inherited target

[Session175](https://github.com/eliottcassidy2000/math/blob/897db58335ce4f7f5c332a238caf572e46537806/00-navigation/S175.md) distinguishes three different quantities: tournament isomorphism counts, ordinary Hamilton vertex-word counts H, and path-orbit counts. This session concerns H and interfaces that compute it.

Ordinary H=7 and H=21 are globally impossible. The conjecture is that every other positive odd integer occurs. The inherited reduction leaves two infinite construction obligations: strong witnesses for every odd prime other than 7, and for every 7p with p an odd prime other than 3. Known carry witnesses 49, 63, 343 complete that multiplicative reduction, not those missing constructions. The inherited computational prefix is strictly below 80,557; this session neither extends that certified prefix nor claims 80,557 is missing.

The progress here is not another unexplained numerical search. It identifies exact state spaces and supplies a strong unit edge gate that persists at arbitrarily large orders, while proving why several tempting simplifications cannot deliver completeness.

## 2. Information levels depend on the allowed observer

These results form an observer-relative hierarchy, not a claim that every listed statistic refines every other one.

**A scalar H is insufficient even for uniform substitution.** The [information packet](../information/PROOF.md) gives two five-vertex modules with H=15 but different ordered path-cover profiles q. The same two-vertex cyclic wrapper sends them to strong tournaments with counts 123 and 135. The lost coordinate is therefore observable by an actual tournament, not merely by a formal linear functional.

**The complete q profile is exactly sufficient for all uniform one-hole contexts.** If exterior vertices each see the whole module uniformly, every Hamilton word splits into module runs and an exterior skeleton:

    H(C[T]) = sum_k c_C(k) q_T(k).

The explicit alternating transitive probes P_r have 2r exterior vertices and triangular diagonal r!, recovering all q coordinates. At least 2k−2 exterior vertices are needed to expose the formal coordinate q(k), and these probes attain that bound. Equality of q is equivalent to equality under every such context and is preserved by actual nested tournament substitutions.

This is a compositional theorem, with a precise boundary: a nonuniform attachment can distinguish even the two orientations of a labelled two-vertex module, despite identical q. The linear dimension required on actual n-vertex profiles is at least floor(n/3)+1, witnessed by ordered triangle/singleton sums. This is not an n-bit lower bound, a claim that formal rank n is realized, or an obstruction to packed integer encodings.

**All one-vertex nonuniform responses are richer in endpoints but not closed under continuation.** The [insertion packet](../circuits/PROOF.md) proves

    F_T(b) = H(T) + 2 sum_(u→v) w_uv b_u(1−b_v),

with positive integer capacities w_uv=(K(u,v)+K(v,u))/2. It includes legal insertions and orphan paths. The capacities, weighted outgoing degrees and H determine the entire response table; quadratic tomography recovers them. Exactly one insertion at each SCC boundary preserves H, giving k+1 zero-increment orientations for k SCCs. All others increase H by at least two.

Yet opposite labelled triangles have identical full one-insertion responses and separate after the same two successive attachments, with final counts 13 and 11. Ordered sums of k independently reversed triangles give 2^k distinct labelled bases with identical entire one-insertion response. All those bases are isomorphic when labels are forgotten. This is exponential labelled information loss, not exponential diversity of isomorphism classes.

The uniform-context theorem and the one-insertion theorem concern different operation sets. Uniform-context equivalence is closed under its declared substitutions; equality for one nonuniform insertion is explicitly not closed under a second. Neither should silently replace the other.

## 3. Constant arithmetic width is not constant information

The [shared-gates benchmark](../gates/NOTE.md) studies a transitive backbone with one arbitrary feedback vertex. Its Hamilton count has a nonzero mixed derivative for every pair of input bits, but a three-register streaming recurrence uses at most 6n binary arithmetic operations. Every ordered cut with at least two bits on each side has bilinear interaction rank exactly two; end cuts have rank one.

Nevertheless, after k input bits, two suffix probes reconstruct the prefix's entire binary integer. Exact deterministic streaming therefore requires at least 2^k states, or k bits, at that cut. Sending the prefix attains the corresponding one-way communication bound. Dense pairwise dependence, bilinear rank, register count and state bits are genuinely different resources.

[Independent bit-level verification](../gates/INDEPENDENT_BIT_AUDIT.md) supplies width n+1 and an explicit O(n²) Boolean implementation of the binary output H. Arithmetic gates with growing integer inputs are not constant-size Boolean gates. Similarly, the insertion packet's binom(n,2) lower bound applies only to a sum of pair-supported gates, not circuits with unrestricted shared intermediates. No general B2 lower bound is inferred.

## 4. A positive construction, and the exact reason it does not yet scale

The [spectrum packet](../spectrum/PROOF.md) constructs strong tournaments L_r,U_r differing by one literal edge, with

    H(L_r)=6·2^r−1,    H(U_r)=6·2^r+1.

This gives a two-state unit-response interval at every order r+4. It is stronger than an isolated finite cube: the unit edge persists at unbounded counts. It is not growing interval coverage, since the centers double and these intervals do not overlap across orders.

Replacing the transitive module by an arbitrary tournament A gives responses 11H(A)+6q_A(2) and 13H(A)+6q_A(2). The edge increment is exactly 2H(A), hence remains two precisely when A is transitive. Recursive substitution of a nontransitive gadget amplifies the increment rather than furnishing an independent binary instruction.

The same packet gives a complete counting circuit for a transitive backbone with k arbitrarily interconnected, nonuniform feedback vertices. It retains empty runs and endpoint legality, automatically counting all orphans. Its explicit state bound is k!·3^(k+1), with O((n+1)k!(k+1)3^(k+1)) integer additions. This is an exact realized construction language, not a minimal automaton, polynomial-time unrestricted algorithm, or existence proof for chosen output values.

## 5. A quantitative lower bound on one constructor language

The [resource packet](../resource/PROOF.md) sharpens Session175's finite-clone obstruction. For prime outputs, discard every quotient with either zero or at least two nontrivial SCCs. Each remaining quotient contributes exactly its unique strong core; singleton padding is irrelevant. Deduplicate those cores up to isomorphism and reversal.

If the resulting K cores have orders q_j≤d and L=floor(log_2(N−1)), positive transitive inflation, reversal and ordered sums supply at most

    sum_j L^(q_j) <= K L^d

distinct prime counts through N. Covering all required primes with K≤N^(alpha+o(1)), alpha<1, consequently requires

    liminf d log log N / log N >= 1−alpha.

An elementary central-binomial argument suffices; no prime-number-theorem assumption is needed. This is a library-diversity versus quotient-order tradeoff for the specified model. It does not bound arbitrary feedback number, arbitrary module substitutions, Boolean gates, or global tournament attainment. The persistent unit gate is compatible with this limitation: fixed-template inflation can have useful infinite families without covering the prime lane.

## 6. Concrete next mathematical steps

1. Grow realized strong response components and prove their odd intervals overlap. The unit edge supplies an exact local control, but its doubling gap and 2H(A) substitution law must be overcome explicitly.
2. Minimize the growing-feedback circuit only after specifying an operation-closed response equivalence. Compare future continuation responses, not just present H or one-step insertion tables. Any savings must include integer bit costs.
3. Seek genuinely growing interface families or richer nontransitive modules. The resource theorem quantifies one obstruction to fixed-core cloning; it is not a reason to abandon constructive search outside that class.
4. For Boolean lower-bound work, retain the original-coordinate gate ledger, shared descendants and protected outputs. The earlier [S118](https://github.com/eliottcassidy2000/math/blob/main/00-navigation/S118.md) deletion accounting and [S123 protected forest repair](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session123_realizers_and_binary_capacity_20261031/circuits/RESOURCE_COST_AND_BOUNDARY_KERNEL.md) remain separate inputs. A new accumulation theorem must prevent recharging the same shared resource; these tournament identities do not supply it by analogy.

## 7. Verification and provenance

Each packet supplies proofs, finite exact controls, matching normal/optimized receipts and an independent review: [gates](../gates/AUDIT.md), [insertion circuits](../circuits/AUDIT.md), [insertion fibers](../circuits/FIBER_INDEPENDENT_AUDIT.md), [uniform information](../information/REVIEWER_AUDIT.md), [spectrum constructions](../spectrum/REVIEWER_AUDIT.md), and [constructor resources](../resource/INDEPENDENT_AUDIT.md). Their finite universes are stated locally. During spectrum review the empty-backbone aggregation bound was corrected from n to n+1; the final statement above uses the corrected bound.

The core inherited mathematical source is Session175's endpoint/orphan interface, path-cover constructions and clone obstruction, linked in the packets. The contextual viewpoint also continues [Session171's response framework](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session171_knots_information_20261103/framework/NOTE.md). Existing classical inputs, including tournament Hamilton-path existence and parity, remain attributed inputs. No literature-priority claim, formal-proof-assistant certification, hosted-CI pass, or resolution of the unrestricted spectrum is asserted here.
