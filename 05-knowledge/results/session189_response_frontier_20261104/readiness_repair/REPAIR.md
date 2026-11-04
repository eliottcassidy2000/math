# Repair of the incoming marked cycle-count correspondence

The incoming PROOFS.md, Corollary 2.1, at commit `a92e0b49befb6ba1ac3fd4f42cca5130db76033b` asserts a bijection of marked strong tournaments that is not onto as written. The elementary cycle-count identity preceding it is correct. The corrected marked statement below works for every positive digit. A canonical-choice argument additionally repairs the unmarked correspondence for digits one and two, without asserting that the original marked statement was true.

## 1. Definitions and exact image

For a tournament T with a named sole-exit vertex a, let r be its unique out-neighbor. Put U=T-a. Define Phi(T,a) by reversing every arc of U and retaining the sole-exit attachment a->r, with every other old vertex pointing to a. All physical names are retained. Phi is an involution on all such marked tournaments.

The elementary deletion identities are

    C(U+a_exit)=s_r(U),
    sigma(U+a_exit,a)=e_r(U),
    C(Phi(T,a))=e_r(T-a)=sigma(T,a).

Here cycles are counted modulo cyclic rotation, and the fixed marked vertex identifies their unique starting position.

**Corrected marked bijection.** For every sigma>0, Phi gives a relabeling-compatible bijection between:

1. Marked strong tournaments (T,a) with a a sole exit of readiness digit sigma.
2. Marked tournaments (Y,a) with a a sole exit, C(Y)=sigma, and Y-a strong.

C(Y)>0 already implies Y is strong. The extra deleted-core condition is necessary; it is not implied by strongness of Y alone.

**Proof.** First suppose T is strong and e_r(U)>0. Strongness of T and the fact that a's only exit is r imply that r reaches every vertex of U. In the linearly ordered strong-component decomposition of U, r therefore lies in the first component. A Hamilton path of U ending at r puts r in the last component. Consequently U is strong.

Thus the image Y has strong deleted core U^op and exactly sigma Hamiltonian cycles. Conversely, if Y-a is strong, its converse is strong. Attaching a with sole exit r preserves strongness: a reaches all old vertices through r, and every old vertex reaches a through an old vertex other than r. The resulting T has readiness digit s_r(Y-a)=C(Y)=sigma. Applying Phi twice restores every named arc. This proves both directions. QED.

## 2. Concrete failure of the original marked statement

Take outgoing rows

    Y=(2,4,9,3),  a=1,  r=2.

Y is strong, a's only out-neighbor is r, and Y has exactly one Hamiltonian cycle. Yet

    Phi(Y,1)=(14,4,0,6)

is transitive, hence not strong. Since Phi is an involution, this marked Y cannot be the image of any marked strong positive-ready tournament. Its deleted core is not strong.

This does not disprove the corresponding unmarked counting formula. The same Y has another sole-exit vertex, 0, whose deletion is strong. That distinction is essential.

## 3. Canonical unmarked repair for digits one and two

Call a sole-exit vertex admissible if deleting it leaves a strong tournament.

**Lemma.** Let Y have order n>=4, at least one sole-exit vertex, and exactly sigma Hamiltonian cycles, where sigma is one or two. Then Y has exactly one admissible sole exit.

**Uniqueness.** If two vertices a,b each have outdegree one, one points to the other, say a->b. Deleting b makes a a sink, so b cannot be admissible. Thus two admissible exits cannot coexist.

**Existence.** Choose any sole exit a->r. If Y-a is strong, a is admissible. Otherwise write its ordered strong components as C1->...->Ck, k>=2. Because a Hamiltonian cycle of Y deletes to a Hamilton path starting at r, r lies in C1 and

    C(Y)=s_r(C1) product_(j=2)^k H(Cj).

Every H(Cj) is a positive odd integer. Since C(Y) is one or two, every factor H(Cj), j>=2, must be one. A tournament with Hamilton count one is transitive, so a strong component with count one is a singleton. In particular Ck={b}.

The vertex b loses to every other vertex of Y-a and beats a. Hence its unique exit is a. In any Hamiltonian cycle of Y, the final portion is z,b,a,r. Because n>=4, z is not r, and therefore z->a. Delete b and replace z,b,a by z,a. This leaves a Hamiltonian cycle on Y-b, proving that Y-b is strong. Thus b is admissible. QED.

Positive-ready strong tournaments likewise have a unique admissible sole exit: the proof in Section 1 shows any positive-ready exit is admissible, and the uniqueness argument applies. Conversely an admissible exit in a strong tournament is positive-ready.

It follows that choosing the unique admissible exit is canonical and commutes with every relabeling. The corrected marked bijection therefore descends to a genuine bijection of unmarked isomorphism classes between:

- Strong sigma-ready tournaments.
- Strong tournaments having a sole exit and exactly sigma Hamiltonian cycles.

This repairs the unmarked correspondence for sigma=1,2 and n>=4. It does not replace the false original assertion about all possible markings.

## 4. Boundaries and what this repair does not certify

Order three is exceptional: the directed triangle has one Hamiltonian cycle and sole exits, but deleting any vertex leaves a non-strong pair. There is no admissible exit.

The unmarked repair is not asserted for arbitrary cycle counts. For example,

    Y=(14,20,24,18,1)

is strong with exactly three Hamiltonian cycles and a sole-exit vertex, but no admissible exit. Its deletion leaves a singleton dominating a directed triangle. The corrected marked statement with the deleted-core condition still applies.

Applying a cited external enumeration requires checking that it counts exactly the unmarked class above, with the same order range and convention. The repair does not independently audit Garey's count, the digit-two enumeration, an all-order heredity theorem, the pruning code, or the reported finite census completeness. The insertion/block lemma is a necessary local restriction on extensions, not by itself a proof that every larger ready class has an eligible smaller parent.

The adjacent verifier independently checks the four-vertex counterexample, its canonical admissible marking and inverse, the order-three boundary, the three-cycle boundary, and the incoming order-12 digit-two witness. Its normal and optimized receipts agree.
