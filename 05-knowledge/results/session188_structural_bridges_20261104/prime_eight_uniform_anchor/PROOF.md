# An endpoint-balanced prime-eight quotient has a universal uniform anchor

## 1. Statement and exact scope

Fix the eight labelled quotient vertices with outgoing bit rows

    Q=(238,4,240,118,195,146,34,74).

Let A_0,...,A_7 be arbitrary tournaments on fixed disjoint physical name sets, with every |A_i| at least two. Put T=Q[A_i] and U=Q^op[A_i]. Target child interiors may instead be any tournaments with the same complete ordered path-cover profiles as their respective source children.

Attach ports x->y. The single query z=(a,c)=(1,0), on quotient masks, means every name of A_0 points to x, x points to every other old name, and y points to every old name. The exact full augmented profile contrast satisfies

    q_h(T+ports)-q_h(U+ports) = 4 if h=1, and 0 otherwise, modulo 8.   (1)

Here q_h counts ordered covers by h nonempty directed paths. This holds for arbitrary unequal child sizes and interiors, not just pair or triangle children. It also holds for the converse-complement query Jz=(255,254).

Consequently the same physical root-incidence lift of z separates T from U in every fixed ambient tournament containing the old core as a module. Its Hamilton contrast remains 4 modulo eight. There is no restriction on the total core order, and no invalid generic modulo-eight lift for odd cores is being used.

This is a theorem for one specified prime quotient. It is not a classification of all prime-eight quotients, a theorem for arbitrary prime quotients, or unrestricted reconstruction. Singleton children are not covered by the theorem.

The independently verified packet `even_endpoint_obstruction` proves that Q is prime and that its Hamilton start and end vectors have the same parity mask 55. Thus the simpler endpoint-parity criterion cannot prove this result. Its failure did not imply a response collision.

## 2. Complete child-profile domain

The all-order five-type normalization from `five_profile_types/PROOF.md`, Section 1, says that every nontrivial child has, after division by its odd Hamilton count modulo eight,

    q1=1, q2=2+4s, q3=4s+2t, qk=0 for k>=4,             (2)

where s,t are bits. The four possibilities (s,t)=(0,0),(1,0),(1,1),(0,1) are respectively the pair, even-order-at-least-four, odd-cycle-even and odd-cycle-odd types. All actual nontrivial children are covered. Their Hamilton counts have not been replaced by one as integers: the product of these eight odd counts will be restored after evaluating the normalized functional.

For a physical block-uniform query, the ordered-cover run expansion is multilinear in child profiles, with exactly one q-coordinate factor for each child. An ordered cover is encoded by concatenating its ordered paths, marking each intervening cut, and then replacing every maximal same-child run within each path by its child label. Each port occurs exactly once. Internal runs of a child are assigned to their ordered occurrences in precisely q_k(child) ways. Cuts prevent runs from merging, including when the same child lies on both sides of a cut.

The skeleton determines all quotient and port adjacencies within paths. Across a cut there is no adjacency requirement. This construction is a bijection, including empty prefix/suffix positions and ports that are endpoints of their cover paths.

Modulo eight, any child with four or more runs contributes zero. A repeated child, with two or three runs, contributes an even factor. Thus every term with three or more repeated children also vanishes. It is therefore sufficient, and exact, to retain run-count vectors r in {1,2,3}^8 with at most two coordinates larger than one. There are

    1 + 8*2 + binom(8,2)*4 = 129                         (3)

such vectors. No low-degree truncation of the two-port Boolean response has been made. The truncation concerns child-cover residues only; each retained skeleton represents complete physical ordered covers.

## 3. Complete finite symbolic certificate

For h=1,2,3, let K_h(Q;a,c;r) be the number of cut-marked run skeletons with r_i occurrences of child i, both ports once, and h nonempty ordered paths. Record K modulo eight. The exact recurrence in `covers.cpp` has state

    (r, port-used mask, final symbol, number of cuts),

with r initially zero except for an initial old symbol, or zero when starting with a port. Appending an old symbol increases its coordinate; appending a port sets its unused bit. Without a cut, the relevant quotient/port edge must exist, and consecutive identical child labels are forbidden. With a cut, every unused next port or permitted next child is allowed, including the same child. Cuts are inserted only between existing nonempty segments. The final state has every child present, both ports used and exactly h-1 cuts.

Counts above three occurrences and states with three repeated children are discarded as justified in Section 2. No discarded state can return to the retained domain. All transitions strictly increase the word length, so the increasing base-four state/port-mask order is acyclic.

The verifier retains every one of the 129 source and target kernel entries for each of the three cover coordinates and each of the two queries. It expands, coefficientwise over the integer polynomial ring reduced modulo eight,

    P_h(s,t)=sum_r [K_h(Q;z;r)-K_h(Q^op;z;r)]
                         product_i q_(r_i)(s_i,t_i).     (4)

The complete symbolic result is

    P_1 = 4,   P_2 = 0,   P_3 = 0                       (5)

for both z and Jz. Every nonconstant coefficient is zero modulo eight. This is a coefficient identity, not sampling of the 4^8 bit assignments and not an assumed relation among actual profile coordinates. The executable checks all coefficients, including zeros, with exact integer arithmetic.

Restoring the product of the eight odd Hamilton counts preserves the value 4 and the zero values. For h>=4, h! divides every ordered-cover count and eight divides h!, so all further contrast coordinates vanish. This proves (1).

The finite quotient-kernel computation is part of the proof. Its domain is explicitly bounded and exhaustively verified. The child-order conclusion is unbounded because (2) and the run bijection are all-order statements.

## 4. Arbitrary fixed ambient lift

Let K be a fixed tournament with a distinguished vertex replaced by the old core. Give both ports the original core's incidence to every exterior vertex. They join the core in a new module. The full ambient Hamilton count has the occupied-run expansion

    H(K[B])=sum_(h>=1) c_h q_h(B),                      (6)

where the nonnegative integers c_h depend only on K and its other fixed children, and c_1=H(K) is odd when all exterior vertices are displayed individually. Equivalently, c_1 is the odd Hamilton count of the contracted ambient tournament.

Applying (1) in (6) leaves exactly 4c_1=4 modulo eight. All other cover-coordinate differences are zero. Thus the proof works in both even and odd total core order. It relies on the complete profile certificate, rather than assuming that a nonzero odd-core Hamilton contrast always lifts.

## 5. Original-module root families

Suppose any family of modules of the original core is simultaneously reversed, with the whole core selected an odd number of times and all selected sets contained in the core. Since Q is prime, every proper original module lies within one child. The inherited hereditary original-module theorem gives coordinatewise nondecrease of every target child cover profile relative to its source profile. The root operation converses the quotient; it does not change these cover inequalities because whole-child converse preserves all q-coordinates.

Let the baseline be quotient converse with source child profiles. The exact whole-augmented-profile duality for the child-uniform query J(a,c)=(1-c,1-a) gives opposite baseline integer contrasts for z and Jz, in every fixed ambient. One way to see this is to reverse all paths of the augmented core and interchange the two ports, then restore child interiors using equality of their q-profiles. Exterior port bits remain the fixed root incidences. The inherited `monotone_dual_queries/PROOF.md` states this precisely.

By (1) both baseline contrasts are nonzero. Choose the member of {z,Jz} with negative source-minus-baseline contrast. Nonnegative run substitution, first through the children and then through the ambient, makes the actual family target at least as large as that baseline. The same chosen query therefore still has strictly negative contrast.

Hence the fixed two-query battery {z,Jz}, physically expanded and lifted as above, detects every original-module family containing this root oddly. For a fixed source and ambient, the sign-selected one of these queries works for every such family. This is source-relative detection, not pairwise injectivity between arbitrary orbit outputs or changes outside the fixed core.

## 6. Exact controls and reproduction

The proof certificate consists of the complete two-query, three-coordinate, 129-vector kernels and their exact polynomial contraction. Independent literal controls include 53 small child models and all their first three augmented cover counts, larger nontrivial and mixed-parity child graphs, fixed ambient lifts, eight root-family controls, and retained signed physical cubic and quartic batteries. The batteries keep every actual query value and every sign; shared physical labels are not collapsed or discarded.

Run `python3 verify.py > normal.json` and `python3 -O verify.py > optimized.json`, then compare. All checks remain active under optimization. The installed C++17 compiler is used in temporary directories and nothing is installed. Physical Hamilton counts use unsigned 64-bit integers only up to augmented order 20, where 20!<2^64. Complete first-three-cover counts are computed only up to order 12, where 12!*binom(11,2)<2^64. The kernel arithmetic is exactly modulo eight, independently of these literal integer controls.

Frozen source allowlist: `PROOF.md`, `README.md`, `verify.py`, `covers.cpp`, `physical.cpp`, `normal.json`, `optimized.json`, `SHA256SUMS`. Exploratory files in the separate `prime_eight_profiles` directory are not part of this packet.
