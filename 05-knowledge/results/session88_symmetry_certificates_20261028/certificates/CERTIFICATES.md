# Session 88: exact state is not the same obligation as an exact certificate

**Status:** the elementary statements and proofs below are **PROVED** within their explicit assumptions; scripts are **FINITE-EXACT** controls, not all-order experimental proof. The group-action mechanisms are classical. The arbitrary-set Gram rigidity used in the interval companion is **CLASSICAL/CITED: Fulkerson–Gross (1965), Theorem 2.1**, with a self-contained independently rederived proof, not a new theorem. The main recovered implementation opportunity is the interval-pattern maximum-spanning-tree formula in [INTERVAL_PATTERNS.md](INTERVAL_PATTERNS.md), which replaces the displayed exponential inversion used in an older repository result and recovers all realized patterns through a positive-prefix trie. No literature novelty claim is made.

## Inheritance, scope, and concept board

Read-only source: `session87-four-vertices/verified-repo`, supplied at `7a0bab9bfa7bb87c5a10c7ba239592d90cfafd51`; additional materialized canon was read in `session87-four-vertices/repo`. The latter is a supplied source collection, not represented here as an independently refreshed current GitHub checkout. No source files or shared mistakes ledger were changed.

* Closest tournament mechanism: Session87 `05-knowledge/results/session87_four_vertices_enumeration_20261028/EXACT_GENERATION.md`, §§2–3: pointed states are aligned `(parent, attachment)/Aut(parent)`; a parent minimum and an independently chosen mask minimum are invalid. Its `check_generation.py` retains every attaining parent map.
* Hostile: the same file §6, the decorated-C3 81/27 pair; agreement of local stabilizers or subset-H data does not identify attachment actions.
* Closest graceful mechanism: **THM-2789**, `01-canon/theorems/THM-2789-interval-gram-tomography-and-graceful-gap-tail-quadratic-detector.md`: Gram determines the incidence-column multiset but Boolean Möbius reconstruction is exponential. **THM-2793**, `01-canon/theorems/THM-2793-oriented-ramp-reference-and-rooted-graceful-gap-reconstruction.md`: one *vector* of ramp moments and a root-zero anchor restore the missing order/side.
* Classical primary-source recovery: Fulkerson–Gross, [*Incidence matrices and interval graphs*](https://msp.org/pjm/1965/15-3/pjm-v15-n3-p11-p.pdf) (1965), Theorem 2.1, p. 837. Their point-row incidence convention is transposed here; their Gram-equality/row-permutation result is exactly the arbitrary-set realization rigidity used in [INTERVAL_PATTERNS.md](INTERVAL_PATTERNS.md), §6. Primary theorem checked 2026-10-28. This is a recovered classical mechanism and specialized algorithmic application, not a priority claim.
* Corrected near miss: **THM-2761**, `01-canon/theorems/THM-2761-graph-edge-sum-discriminant-codegree-factorization-and-graceful-sign-gauge.md`: oriented differences `(-1,1)` are distinct but their absolute differences coincide.
* Square-class source: corrected **THM-4840**, `01-canon/theorems/THM-4840-tutte-camion-dichotomy-for-arithmetic-square-class-banks-hamming-15-is-the-equal-parity-degree-space-of-k6-rank-6-prime-banks-unit-cayley-tournaments-c3-t5-p7-and-fano-chromatic-numbers.md`, including its correction banner. No uncorrected parity-extension or prescribed-TU-signing claim is inherited.

Concept board: faithful operation state; named root-switch blocks; action closure versus reachability; stabilizer/transporter reuse; one specified interval-column decision; arithmetic/graceful local guards. The anchor is exact tournament classification reuse; the niche is reusable certificates; the wildcard is the old graceful Gram algorithm.

## 1. What the information lower bound actually says

Let `X` be a state space, `T_a` its allowed operations, and `q_b` its demanded outputs. Treat legality of partial operations as an observable too, or totalize operations with a failure symbol. A representation `R` sufficient for all finite future operation words must distinguish states whenever some future output does:

    R(x)=R(y)  implies  q_b(T_w x)=q_b(T_w y) for every legal w,b.

If the right-hand observations separate every pair of states, `R` must be injective. Without that separating hypothesis, faithful literal-state reconstruction does **not** follow. This is an elementary indistinguishability condition, not a universal impossibility theorem for compression.

For one predicate `P(x)`, a verifier can instead consume `(x,c)` and accept a witness `c`. It need not reconstruct `x` from `c`, simulate future operations, or decide unrelated predicates. A successful witness can be much shorter than a complete operation carrier. But four distinct costs remain:

1. the input read by the verifier;
2. the witness length;
3. witness discovery;
4. witness verification, including trusted/preverified setup.

Calling a certificate short does not make its discovery cheap. A verifier using an enormous complete parent bank has not replaced that bank with a small object. An existential positive certificate does not automatically furnish short negative certificates or exact counts.

## 2. A named, aligned root-switch row is a certificate for one whole block

### Setup and exact assumptions

Let `U` be a tournament on `n>=1` named vertices. A **verified pointed bank** `B` has exactly one representative of every pointed isomorphism class under consideration. Every entry `j` consists of a parent tournament `T_j` on `0,...,n-2` and a mask `S_j`, with child `E_j` obtained by adding root `* = n-1` whose out-neighborhood is `S_j`.

For the full order-n universe, a complete bank is obtained from a complete pairwise nonisomorphic order-(n−1) parent bank and exactly one mask per orbit of the **full actual parent automorphism action**. That setup is assumed verified here; it is not established merely by checking one row.

For each vertex `v` of `U`, supply:

* a **named bank ID** `j(v)`;
* a bijection `f_v:V(U)->V(E_j)` sending `v` to `*`;
* equivalently, its restriction from `U-v` to `T_j`, with attachment mask transported in the **same** frame.

The verifier checks permutation validity, root image, and every directed adjacency. This subsumes deletion-arc and mask checks. It performs no canonical search.

### Row theorem

If all `n` root witnesses pass, then

    Q_U(j) = #{v : j(v)=j}

is exactly the root-switch row, including every nonzero entry and its multiplicity. Its support is precisely the pointed bank entries whose underlying ordinary tournament is isomorphic to `U`. Therefore one verified row identifies its entire ordinary equivalence block.

**Proof.** Each witness identifies `(U,v)` with the indicated pointed bank object. Uniqueness of the pointed bank makes this its unique ID. Conversely any pointed bank object with underlying tournament isomorphic to `U` is the image of `(U,v)` for the preimage of its distinguished root under an unpointed isomorphism. Every `v` was included. This proves both inclusions, and counting vertices proves multiplicities. ∎

Consequently rows are constant on each ordinary block, and one can safely fill all named entries in the support after the first verified row. This conclusion is about the actual dictionary, not an anonymous multiset of row weights.

### Size and verifier cost

For bank size `b`, explicit maps take `O(n^2 log n+n log b)` bits, apart from the bank and child input. Verification takes `O(n^3)` adjacency comparisons plus dictionary lookup/permutation checks. Bit-level costs depend on adjacency storage and ID encoding. Storing only one map per root is enough: retaining **all** isomorphisms is optional stronger information, not necessary for this theorem.

If discovery already knows a deletion-to-parent map and an action word moving its transported mask to the bank mask, compose them. Verification can check the resulting single map instead of rerunning either canonicalization. A word certificate is also valid, but its explicit length/action costs must be counted.

### Hostiles and the exact completeness boundary

* **Missing roots:** valid maps for a proper subset of vertices certify only those memberships and lower bounds on their multiplicities. For a transitive order-three tournament all three roots have different pointed IDs; omitting one leaves one block member undiscovered.
* **Duplicate root records:** `n` records do not suffice unless their source roots are exactly the `n` distinct named vertices.
* **Separate minima:** an isomorphism of deletion parents and an unrelated convenient mask permutation need not assemble into any pointed isomorphism. Check the same map on both.
* **Anonymous shape:** many nonisomorphic vertex-transitive tournaments have anonymous row `[n]`; the named ID contains the distinguishing information. The blockfill lane independently tests the cyclic order-seven connection sets `{1,2,3}` and `{1,2,4}`: both shapes are `[7]`, while their actual group orders are `7` and `21`.
* **Incomplete bank:** a successful local row proves the class of this child has all its pointed types represented; it does not prove every other ordinary class exists somewhere in the bank.
* **Wrong ID equality:** a hash is an indexing aid, not a substitute for exact pointed-state equality. A duplicate representation under two IDs can falsely split a block; merging nonisomorphic objects under one ID invalidates the bank contract.

A collection of partial membership edges can still be useful: every edge safely unions states known to belong to one class. It yields the exact partition only if each true block is connected by the supplied edges and all bank states are covered. Otherwise it gives a refinement and hence an **upper** bound on the number of classes, not exact class counting. Each verified same-class witness by itself gives a **lower** bound on that block's size. These directions must not be confused.

## 3. The same row yields the full child automorphism group

Choose a root `r`. Let `j=j(r)` and `f_r:U->E_j` be its certificate. In the known full parent action, let

    H = Aut(T_j)_{S_j}.

Extend every element of `H` by fixing the new root, and conjugate by `f_r` back to `U`. This is **exactly** `Aut(U)_r`, not merely a subgroup.

For each `v` with `j(v)=j`, set

    a_v = f_v^{-1} o f_r.

Then `a_v` is an automorphism of `U` carrying `r` to `v`. The following is the usual stabilizer–transversal argument:

    Aut(U) = < Aut(U)_r, a_v : j(v)=j >,
    |Aut(U)| = Q_U(j) |H|.

**Proof.** The row identifies the full orbit of `r`. Given any automorphism `a`, put `v=a(r)`. Then `a_v^{-1}a` fixes `r`, so is in the stated stabilizer. The generators are themselves automorphisms. The order formula is orbit–stabilizer. ∎

No child-sized canonicalization or enumeration of all child automorphisms is logically needed. The parent mask-stabilizer computation itself may be nontrivial; it is not free if only parent generators, rather than a full action list or certified stabilizer chain, are supplied.

### Adaptive first-row construction

This is more than caching rows after they are found:

1. Lift the known original parent-mask stabilizer to a subgroup `K<=Aut(U)`.
2. Retain the already known original rooted ID/map.
3. Propagate that map along the `K`-orbit of its root. If `a(u)=v`, the map at `v` is `f_u o a^{-1}`.
4. Classify one root not yet covered by known `K`-orbits, retaining a transport.
5. If its ID equals an earlier root's ID, compose the two maps to get a new actual child automorphism; enlarge `K` and propagate again. Lifting its parent-mask stabilizer is another safe source of generators.
6. Stop only when all roots have assigned, verified maps.

Even a proper subgroup safely avoids classification calls inside its own orbits. It never licenses merging different orbits just because they have equal sizes or weak invariants. At termination the full row is exact, and the theorem above certifies the full group from a fixed root's stabilizer plus its reached transporters.

The sibling blockfill implementation tests this adaptive construction against complete exact rows through order eight. This file states the theorem; its independent small checker audits group reconstruction for every ordinary class through order five. For asymmetric children no orbit savings need occur. Fewer canonical calls do not alone imply a runtime gain after orbit and permutation bookkeeping.

## 4. Exact orbit certificates: reachability plus closure

Let `X` have exact decidable equality, and let `S` be a finite list of specified bijections of `X` generating `G`. The action maps and their semantics are fixed, not guessed from observed samples. The ambient `X` and `G` need not be finite.

A certificate for a finite orbit of `x` consists of:

* distinct literal states `B=(b_0=x,...,b_{M-1})`;
* a rooted transport tree: for every `i>0`, indices `p(i)<i` and `s(i) in S` with `s(i)(b_p)=b_i`;
* a closure table `c(i,s)` with `s(b_i)=b_c(i,s)` for every `i,s`.

**Theorem.** These checks imply `B=Gx`.

**Proof.** Tree induction gives `B⊆Gx`. Closure gives stability under positive generator words. Since every `s` is injective and `B` is finite, `s(B)⊆B` implies `s(B)=B`, so `B` is also stable under inverses. Thus `Gx⊆B`. ∎

For monoid actions, the same tree and closure proof gives the exact *forward reachable set*, without the inverse conclusion. For an infinite candidate set, finite-set inverse closure is unavailable; include inverse generators or prove inverse stability explicitly.

A closure table is optional if the verifier recomputes every successor and searches the literal-state dictionary. The table makes all destinations explicit and prevents treating a filtered traversal's missing transitions as completed work. The official [GAP `orb` documentation](https://gap-packages.github.io/orb/doc/chap3_mj.html) describes Schreier trees and incremental generators; importantly, it warns that grade-restricted enumeration can report closure despite omitting part of the full orbit. Our contract requires literal closure under every declared generator, with no hidden filter.

### Honest costs and negative answers

For state length at most `ell` bits and `t=|S|`, explicit certificate size is

    O(M ell + M log(M(t+1)+1) + Mt log(M+1)).

Verification uses `M-1+Mt` action/equality checks, plus validation and state distinctness. Let `C_act` and `C_eq` denote their actual costs; table checks cost `O((M+Mt)(C_act+C_eq))`, with distinctness/dictionary overhead separately. A deterministic sorted dictionary adds comparison costs. Hash tables must compare full states after hash matches. A tree shares prefixes; expanding one full word for every orbit state can need quadratic total word length even when the tree has `M-1` edges.

The verifier is output/certificate-sensitive. `M` may be exponential in the original instance length; **this is not a blanket polynomial certificate theorem**.

* Tree without closure: only a reached subset/lower bound.
* Closure without root reachability: may include unrelated orbits.
* Closure after a noninjective encoding: an exact quotient orbit, not necessarily a source orbit. For example S4 acts transitively on its six two-subsets, while the cardinality encoding has a one-point orbit.
* A target found on a valid path is a positive reachability witness. A target absent from a partial search is no negative witness. An exact closed orbit excluding it is a negative witness, potentially huge.
* The certificate proves the orbit under the **supplied generators**. Claiming this is the orbit under a larger intended group requires a separate generation/completeness proof.

## 5. Arithmetic and graceful certificates obey different guards

### Square-class arithmetic

Fix distinct primes `p_1,...,p_m`, and include **all** prime factors of their `p_i+1` in a valuation-row bank. Let `A_li=v_l(p_i+1) mod 2`. For a subset vector `x`, multiplicativity of divisor sum on squarefree products gives

    sigma(product_{i:x_i=1} p_i) is square  iff  Ax=0.

This preserves a specific multiplicative square predicate, not integer magnitudes, addition, prime labels, or sigma on repeated prime powers. A positive certificate can instead give an integer `z` satisfying `z^2=product(p_i+1)`; direct big-integer equality verifies it without reconstructing the factors. A nonzero syndrome row is a negative certificate once the complete, correct valuation setup is trusted. If a claimed odd-valuation obstruction uses a new base `q`, its primality must be checked: `q=4`, `N=4` is a hostile to treating an odd exponent of any composite base as a nonsquare certificate.

A column permutation `P` with a verified invertible row map `L` satisfying `AP=LA` preserves `ker A` and therefore square-selection decisions. It need not preserve the actual selected product or square root. This is a lawful action on the **code**, not automatically an arithmetic symmetry of addition or of prime values. The corrected THM-4840's equal-degree-parity K6 description applies to exactly this typed code; its extra parity coordinate is edge-count parity, not common vertex-degree parity.

### Square-sum paths

For `Q_n` on `1,...,n`, an affirmative Hamiltonian-path certificate is a permutation `(v_1,...,v_n)` and, optionally, square roots `s_i` such that `v_i+v_{i+1}=s_i^2`. Checking distinctness, coverage, and `n-1` local equalities proves the decision. There is no requirement to carry the whole graph or every possible continuation. The checker verifies the inherited THM-4505 order-15 path.

A failed local guard refutes that *candidate*, not every path. The degree obstruction used in **THM-4505**, `01-canon/theorems/THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md`, is a separate negative certificate: for `n>1`, an isolated vertex or more than two degree-one vertices prevents a Hamiltonian path. This does not promise a short obstruction for every non-Hamiltonian graph, nor transfer to a restricted move component without coverage.

### Graceful labels and missing edges

For a fixed graph with `m` edges, a certificate of gracefulness is an injection `f:V->{0,...,m}` whose **absolute** edge differences are exactly `1,...,m`. Verification is local plus a difference inventory. Distinct signed differences alone fail on the path labelled `(0,1,2)`, as THM-2761 records.

**THM-2786**, `01-canon/theorems/THM-2786-binary-golomb-universal-edge-difference-separation-and-graceful-compression-boundary.md`, identifies an important over-solving trap: separating all vertex-pair differences requires span at least `binom(n,2)`, whereas a graceful tree asks only for its `n-1` edge differences. A compressed successful witness is allowed to collide on nonedges. That is a concrete reason decision certificates can be smaller than universal recoverable representations.

For one graceful gap-tail query, THM-2789's interval Gram is sufficient for the *existence/count of an unmarked column pattern*. [INTERVAL_PATTERNS.md](INTERVAL_PATTERNS.md) now computes that exact count in `O(m^2)` arithmetic operations; a partial-pattern version also recovers every realized column pattern through a positive-prefix trie in `O(m^4)` arithmetic operations, never enumerating the absent patterns. It still does not identify the gap address or check which side contains the proposed attachment vertex. THM-2793's ramp/root sidecars restore those data when needed. A certificate must match the requested decision, not silently upgrade an unmarked existence statement to a marked extension.

## Reproduction and scope

    python3 -B check_certificates.py > receipt.json
    python3 -B -O check_certificates.py > receipt-optimized.json
    cmp receipt.json receipt-optimized.json

The checker uses explicit exceptions, not Python assertions. It checks all ordinary tournaments through five by unrestricted permutations; exact orbit positive/hostile controls; the seven-prime square-class bank's 128 subsets; the order-15 square-sum path; the signed/absolute graceful hostile; and the interval formula on the universe detailed in the receipt. No trusted legacy code is imported by this checker. These finite controls supplement the proofs, not establish unbounded claims.
