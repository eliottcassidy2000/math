# Independent bounded-slack / prime-amplitude audit

2026-11-05. **PASS for the stated all-order normal form, six coefficient
rows, both amplitude-seven exclusions, and the explicit H13 sign control.**

Frozen source hashes:
- PROOF.md: `eb3fd77753f216e36e1ab8ce26c9cbd9185cfe54ddd96acc47503d2c44863f31`
- verify.py: `dc500c0f04b24a3d35e67ac16346fd55868d79be33bd60051da2d2faed22a2fb`
- Both receipts: `43602fc71baa06269de0ff00bf5e142a6977df42b228c420e9c51365e0a35cc0`

## Independent checks

`audit.py` independently builds forward-extension path tables, start tables
on the converse, and the complete boundary-cover matrix by all two-part
vertex partitions. It checks every one of the 5,405 rooted tournaments
through order five against literal permutations and cuts, including the
one-vertex case and empty deletion graph. It also checks 1,240 larger
random rooted cases through order nine. Every profile, transport identity,
boundary-slack identity, quantitative inequality, and normalized profile
matches. The H3 and H5 strong-component isomorphism classes were independently
recovered from the entire labelled universe through order five.

The exact unit-source table was rebuilt from the deletion graph's
boundary-cover matrix, rather than by the source's full rooted tournament
profile enumeration. All six rows agree:

- d=1: beta +1 has bases {1}; beta -1 has {3,5}.
- d=3: beta +1 has {6}; beta -1 has {9,15}.
- d=5: beta +1 has {11,12,16}; beta -1 has {15,39}.

Each base is multiplied by arbitrary nonnegative powers of two. Independent
representative profiles are saved in `AUDIT.json`.

`padding.py` adds twelve explicit adversarial checks with multiple
nontrivial SCCs, singleton components between them, deleted suffix singletons
on both sides of a later nontrivial SCC, empty root-outneighbor sets, and
randomly interleaved vertex labels. All profile normalizations agree.

An independently ordered clone construction gives the H13 control counts
13,25,49,97,193,385. The first four are additionally checked by literal
permutations. The source verifier passes normally and with `-O`; outputs
are byte-identical and equal to the frozen receipts.

## Proof audit

**Transport and positivity.** Removing the initial root of the second path
gives `U=sum_O R_x(S)` with its singleton contribution d. The boundary
identity follows by separating R-covers whose concatenation through the
root is illegal: covers with singleton second path contribute d-e, and
the remaining illegal covers give exactly the ordered pairs counted by f.
There is no missing empty-path term; H(empty)=1 is explicit in this packet.

The inherited two-endpoint lower bound can also be proved directly. Given
a Hamilton path P on S minus {x,y}, find its first vertex not beating both
x and y. If there is one, choose an endpoint q beating that vertex, and
place the other endpoint p immediately before q at that position. Earlier
vertices beat both endpoints, so all external edges of the inserted block
are legal. If no such vertex exists, append the block. Cut between p and q;
the internal orientation is immaterial. This produces one of the two
specified boundary covers and recovers P on deleting x,y. Choosing the
endpoint deterministically gives the required injection. Thus
`B_xy+B_yx >= H(S-{x,y}) >=1`, including S={x,y}, and
`f>=binom(|O|,2)` is valid.

**SCC normalization.** The SCC ordering used is the tournament's abstract
ordered join, not an interval ordering of vertex labels. Components before
the first root outneighbor can contribute only to the first boundary path;
components after the last outneighbor can contribute only to the second.
Consequently removed prefix singletons are neutral, while every removed
suffix singleton multiplies U by two and leaves d and f unchanged. If a
suffix singleton exists, the final deletion SCC is disjoint from O, hence
e=d. This remains true when removed suffix singletons occur before or
between later nontrivial SCCs. Moving them to a terminal K-chain is an
equality of profiles, not a claimed graph isomorphism.

The interior-singleton construction is valid: choose an outneighbor in the
first occupied SCC as the start of the second path, and one in the last
occupied SCC as the end of the first. Entire endpoint SCCs and all other
non-singleton components can be assigned as stated. Each designated
interior singleton independently chooses a path, yielding `2^u` distinct
covers. If the occupied SCC is unique, u=0. Otherwise f>0, so no logarithm
of zero is taken. At most d vertices lie in nontrivial deletion SCCs,
using `H(C)>=|C|` for strong tournaments and the product of their counts.
Together these facts prove the claimed all-order core bound.

**Exact finite reduction.** At deletion count 1 there is no nontrivial SCC;
at deletion counts 3 and 5 exactly one nontrivial SCC exists. Its order is
at most its path count, so the exhaustive labelled enumeration through
order five suffices to establish the only possible strong component types.
Every normalized graph is then a fixed strong component, when present,
inserted among at most `k(f)+ell(f)` singleton components, plus arbitrary
root incidences. The enumeration covers all positions and all incidences.
The source's reduction test removes only non-normalized candidates.
Consequently the finite table plus the universal theorem is an all-order
classification, not an empirical cutoff extrapolation.

**Unit scope and prime exclusions.** Fixing d bounds f for unit offsets
only in the singleton-partner/K source class. The document explicitly
excludes general varying partners from that inference. For odd prime
amplitudes, the complete inherited star synthesis and its physical
inequalities independently force the partner reduction: U=1 would imply
the profile (1,0,2,1), a transitive partner, and a power-of-two L, which
cannot be the odd prime. Thus U=p and L=H(Y)=H(Y-c)=1.

The new endpoint insertion injection is valid: in a path ending at r,
remove r and insert it before the first outneighbor. All preceding vertices
beat r, and the new path does not end at r. Deleting r recovers the input,
so `e<=h-e`. Oddness gives `h<=2(U-f)-1`, and therefore
`3d<=2p+2beta-1`. At p=7 this leaves d=1,3 for beta=-1 and d=1,3,5 for
beta=+1. No corresponding table row or its dyadic multiples contains 7.
Both all-order cloning obstructions follow without assuming a global
exclusion of the Hamiltonian-path count seven.

**Sign and scope.** The scalar return is `H -> 2H-beta`. The forbidden
H13 family has beta=-1 and counts 13,27,55,...; the actual H13 control has
beta=+1 and counts 13,25,49,... . The H15 obstruction is for the returning
minus step. Neither excludes a single accidental clone step, a higher-run
family, an unrelated graph operation, or the counts H13/H15 themselves.
No global {7,21} spectrum theorem is claimed.

No mathematical or verifier blocker was found. All source edits remained
with the producer; this auditor wrote only its own scratch artifacts.
