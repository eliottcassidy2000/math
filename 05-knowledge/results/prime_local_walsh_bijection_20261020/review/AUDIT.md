# Independent review: prime-local Walsh correspondence

2026-10-20 UTC. **PASS** for the three notes identified below. No mathematical correction is required. The stronger addendum supersedes the initial note's narrower open-problem boundary.

## Reviewed revisions

- `PRIME_VALUATION_THEOREM.md`: `e7498c1aabdec73f391cb41520727550133a0fb4996fd038353797cf339f71b3`
- `SYLOW_ACTION_ADDENDUM.md`: `17431560b037d5dc12b027d8a55a99edbac7ac5b5772243f6f0ec30d49a9b6ec`
- `COKERNEL_ADDENDUM.md`: `f35dad137b7b41ab348f7142a7511f8f6b8039470a1ff62acfec34c508658ded`

These files are in the sibling `session10-walsh/` directory. The audit covers the proofs, declared construction, and independent finite controls through n=5. It does not claim a fresh n=6 enumeration, external literature priority, or a formal proof-assistant verification. The separate converse/parity refinement is outside this audit.

## 1. Weighted pairing and Walsh basis

The argument is valid over R=Z[1/2]. Inverse entries need only be p-integral for odd p, not integral over Z. The inverse formula gives the row-to-column valuation inequality on p-unit entries; a nonzero determinant modulo p supplies a support perfect matching. Equality of the total valuations forces equality on each matched edge. Every support perfect matching works. Equality is not asserted for every support edge.

The signed Walsh invariant-basis construction is sound over R: a sign-odd stabilizer forces an invariant coefficient to be its negative, hence zero; a sign-even stabilizer permits a consistently transported signed orbit sum. No averaging by n! is needed. Orthogonality of the original Boolean characters gives the stated weighted identities and dyadic invertibility.

Sorting by valuation gives square invertible diagonal layers. This follows from the already proved equality of valuation multiplicities, not from an unstated assumption about row and column block sizes.

## 2. Cofactor construction

At every stage the retained matrix A is an actual row-column minor of the original matrix. The identity (A Q)[r,r]=1 guarantees some column j with A[r,j] and Q[j,r] both nonzero. The latter coefficient certifies that deleting row r and column j leaves an invertible minor.

The stated update

    Q_new[k,l] = Q[k,l] - Q[k,r] Q[j,l] / Q[j,r]

is correct with inverse indices ordered columns-by-rows. This is the inverse of the retained original minor, not a replacement matching matrix. Thus no fill-in entry is mistakenly admitted as an original support edge. Column sign changes preserve both nonvanishing tests at every stage. The stated cubic field-operation bound is sound.

## 3. Independent finite controls and joint-prime obstruction

`independent_small_check.py` imports no target checker. It independently reconstructs complete labelled tournament and graph orbits for n=0,...,5, using explicit directed edges and permutation transport, and checks the signed Walsh matrices, canonical catalogues, automorphism orders, and determinants against the frozen data.

For Q and F_3,F_5,F_7 it recomputes a fresh inverse at every deletion step and compares it with the proposed downdate, as well as the target matching and cofactor transcript. All 84 deletion steps agree. Normal and optimized runs give byte-identical receipts; the code uses explicit exceptions rather than removable assertions.

The joint odd-part multisets agree through n=4 and differ at n=5. The elementary independent obstruction is valid: tournament automorphism groups have odd order; an order divisible by 15 inside S_5 would have order 15, hence be cyclic by the Sylow counts, but S_5 has no element of order 15. The sign-even empty graph has automorphism group S_5.

There is also a two-row Hall certificate in the exact n=5 matrix. The tournament masks 5 and 76 have automorphism orders 3 and 5. Each has just one column with coefficient coprime to 15: graph mask 0, the empty graph. Their two-element row set therefore has only one neighbor. This checks the failure of a simultaneous unit-support matching directly, despite invertibility modulo 15.

## 4. Sylow-action addendum

The crucial sign cancellation is correct. If Y is a labelled copy of a sign-even graph and h fixes T, the transport formulas multiply by the same orientation sign in epsilon(hY) and phi_(hY)(T). Their product is unchanged. Therefore the signed summand is constant on every orbit of any subgroup of Aut(T).

For a p-group P, all non-singleton orbit sizes are divisible by p. A p-unit signed sum consequently forces at least one P-fixed graph copy. Taking P to be a Sylow p-subgroup of Aut(T), the matched valuation equality makes the common P also Sylow in the output graph automorphism group. Relabelling then preserves the S_n-conjugacy class of the Sylow action, not merely its abstract group or order.

The predicate that Aut(object) contains an S_n-conjugate of a specified p-group H is determined by this Sylow-action class: every p-subgroup lies in a Sylow subgroup, and Sylow subgroups are conjugate within the automorphism group. Hence one matching preserves all such predicates, and their finite Boolean combinations, for the chosen odd prime. This proves precisely the prime-power portion of HYP-9340(a), whose definition was checked in the repository.

The refined action-block argument is also valid. Within each already square order layer, distinct action classes have no connecting entries; full rank forces each rectangular class block to be square and invertible. No labelled equivariance, fusion data, or mixed-prime fixed-point claim is inferred.

## 5. Cokernel addendum and publication boundary

The equality M(D_b R^N)=D_a R^N follows directly from M D_b M^t=c D_a and the invertibility of c and M. Thus [x] maps to [Mx] by a well-defined isomorphism of the two cokernels. Identifying R/dR with the cyclic group of order odd(d) is valid, with powers of two inverted modulo the odd modulus. The n=5 decomposition (Z/3)^4 plus Z/5 versus (Z/3)^3 plus Z/15 is correct. This is a direct-sum isomorphism, not a bijection of the displayed cyclic summands.

When assembling the publication, identify the initial theorem's fixed-non-Sylow-p-group limitations as superseded by the Sylow addendum. The inherited HYP-9340(b) suggestion that a single bijection might preserve every odd-H symmetry predicate simultaneously is refuted by the n=5 joint-prime obstruction; equality for each predicate separately never implied such a simultaneous coupling. The conjecture for each mixed-prime odd-order H separately remains outside these proofs.

No repository checkout or Git data were edited by this review.

## Final Sylow revision reconciliation

The final Sylow addendum appends direct finite action controls without changing the preceding proof. Its verifier and frozen JSON/transcript were inspected: the nontrivial p-unit edge counts are 1,4,13,97 at n=3,4,5,6, with 23 matched transporter witnesses in total. Normal and optimized transcripts agree. The bounded whole-automorphism-group shortcut is explicitly checked and correctly restricted to this range. These added n=6 controls were inspected, not independently rerun by this reviewer.
