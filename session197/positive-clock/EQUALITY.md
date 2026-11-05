# Equality classification for the marked self-converse positive-unit clock

This is a separate extension arising during independent review of the peer's sharp amplitude72 theorem. It is not a change to that frozen packet.

## Theorem

Let a rooted tournament (T,r) admit an anti-involution fixing r. Its transitive root clones have H_m=72*3^m+1, for every positive m or eventually, if and only if T is an ordered join

    transitive k vertices -> S -> transitive k vertices,

for some k>=0, where the rooted core S has outgoing rows

    [414,260,8,50,294,327,415,318,12], root4.

The core is rigid, has a unique root with this clone identity, and its unique anti-involution is reversal. Consequently, for each odd n>=9 there is exactly one unlabelled rooted source and one unlabelled unmarked source with this property. There are n! labelled rooted pairs and n! labelled unmarked tournaments; if the root label is prescribed in advance, there are (n-1)! labelled sources. There are none at other orders.

## Proof

The established root-run argument turns an eventual identity into the exact one. Its first two values are217 and649. The root cannot be a singleton SCC. The Hamilton-path product of all components outside its strong component divides gcd(217,649)=1, hence those components are singletons. Removing them preserves all clone counts.

The anti-involution fixes the root component setwise and reverses the ordered SCC list. Thus the strong root core has odd order, and the singleton prefix and suffix have the same length k.

The strong core has H=217. Busch's exact Theorem1 gives minimum225 at order11. The safe inequality 5^(n-1)<=H^3 and 5^11>217^3 exclude order at least12. Therefore the odd strong core has order at most9. This step uses the exact order11 bound: the weaker uniform order11 estimate alone is not sufficient at H=217.

The complete physical fixed-root anti-involutive universe through order9 is independently regenerated in `check_positive_clock.cpp`. It enumerates every tournament on O=N+(r), every symmetric O-to-tau(O) crossing matrix, and forces the remaining arcs. Direct root-run word counts and exact matching rank select precisely24 presentations with profile(217,216,72), all at order9. A degree-filtered exact isomorphism backtracking check maps every one to the displayed core with its root preserved, and finds exactly one such map. Thus there is only one rooted strong equality core.

The same exact checker enumerates all automorphisms of the core, with and without fixing the root; both groups are trivial. It also checks all nine possible clone roots, with rectangular cross-matching tests where needed. Only root4 has matching rank2 and the required profile. This proves root uniqueness without assuming an arbitrary candidate root is balanced.

The converse construction is immediate: ordered singleton padding contributes one to every Hamilton-path product. Extend reversal of S by pairing opposite singleton positions. This gives an anti-involution fixing the same root and preserves the clone identity.

An automorphism of a padded source preserves its ordered SCCs, fixes each singleton, and restricts to an automorphism of S. Hence every padded source is rigid. There is exactly one allowable root, giving the stated labelled counts. Any two anti-isomorphisms differ by an automorphism, so the displayed anti-involution is unique as well.

## Additional independent mass check

Choose one vertex from each reversal pair {0,8},{1,7},{2,6},{3,5}, order the chosen four vertices, insert fixed root4, and complete by the reversed anti-image. There are384 candidate words; direct adjacency checks retain11. Their full list is in `equality-controls.json`.

This agrees with the peer's11 normalized amplitude72 presentations. Since the equality core is rigid, the anti-involution centralizer is trivial, and the fixed-word/normalization fibre has no further quotient. The24 physical O/I presentations likewise equal4! relabellings of the O side. These are checks, not assumptions in the equality proof.

## Verification and dependencies

`check_positive_clock.cpp` and `independent-opt.txt` provide the complete independent bank, run counts, exact isomorphism/automorphism tests, all-root check, and physical clone controls. Optimized and debug outputs agree. `equality-controls.json` records the literal fixed-word enumeration and exact order-bound integers.

The peer's frozen `session197-orbit-wildcard/positive-clock/PROOF.md` supplies the already established sharp positive-unit theorem and run notation. The root-run/matching proof is `session197-prime-atoms/MINIMUM_TERNARY.md`. Busch's source is *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC13(2006)N3, Theorem1, DOI10.37236/1141, https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf.
