# Independent Round2 audit: zero-defect and sharp growth

**Verdict: PASS.** The all-order zero-defect classification, sharp singleton floors and equality interfaces, composite boundary, and two-input transfer are correct. No correction is required. The exact low-growth classification is a real strengthening of the finite-only observation in THM-4876, not a proof of all-target coverage. The reversed-extreme-arc examples and the outdegree/endpoint inequalities remain inherited.

Inputs are `../session231-physical-descent/ROUND2.md`, `round2_verify.py`, and `ROUND2_CHECK.json`. This audit also inspected the precise current statements in THM-4873 `end-count-lemma-and-c2-digit-heredity-the-strict-ternary-parent-dies-at-the-grandchild.md`, part (ii)(D), and THM-4876 `deletion-and-flip-increments-of-the-hamiltonian-path-count-an-increment-of-two-forces-at-most-three-paths-the-increment-ladder-thm-013-for-all-orders-and-what-extension-spectra-forbid-for-chan-pak.md`, statement 7. THM-4876 explicitly labels the relevant old growth/equality observation FINITE-EXACT through order 9. Targeted prior-source searches did not locate an earlier all-order version of this precise classification; this is not an exhaustive literature-priority claim.

## Proof audit

### Sole exit and zero-defect SCC factorization

Let a have sole outneighbor w in B=A−a. Every path ending at a is a B Hamilton path not ending at w, hence e=d−q. Every Hamilton path not ending at a uses a→w and uniquely corresponds to an ordered B cover whose second path starts at w; the first path may be empty. Thus h=e+R. Root-starting covers either have the singleton path (a), counted by d, or use a→w, counted by R. Consequently u=d+R and δ=q, with no boundary exceptions.

The inherited δ≥r(r−1) implies that a nonsink with δ=0 or 1 has exactly one exit. For δ=0, q=0. A vertex has positive endpoint count exactly when it lies in the final SCC: paths follow SCC order, and Camion provides each endpoint in a nontrivial final SCC; a singleton final SCC is immediate. Therefore B=P→C with P ending in the w SCC and C nonempty, uniquely.

The w-starting cover bin is already nonempty inside P. The other bin may be empty there. Each C two-bin cover extends these bins, including the possibility of an empty bin in C, and each full cover restricts uniquely. Hence R=zQ and (h,e,u)=(pc+zQ,pc,pc+zQ). This also proves the converse construction, because w is then outside B's final SCC. There is no implicit assumption that a numerical (p,c,z,Q) is physical.

Since z≥p and Q≥2c, h≥3e. All variables are positive. Equality requires z=p and Q=2c. If C has at least two vertices, a nontrivial ordered bipartition and Hamilton paths in its parts give a cover with both bins nonempty, in addition to the 2c covers with an empty bin. Thus Q=2c iff C is a singleton. This proves both the inequality and its exact equality conditions.

For e=1, pc=1 gives p=c=1. A tournament with H=1 is transitive: a triangle has three Hamilton paths and canonical insertion of remaining vertices is injective; every nontransitive tournament contains a triangle. Hence P,C are transitive and w is P's sink, giving z=1. With t=|C|≥1, Q=2^t, so h=u=2^t+1. Remove P−w, which is a universal transitive prefix. The remaining order w,C,a is transitive except for a→w, and a is the former last vertex. This proves the asserted full topology, not just the scalar values.

### Old and improved singleton floors

With b=(h−1)/2−e≥0, L=2h+δ−e=(3h+1)/2+δ+b. All quantities δ,b are nonnegative integers. Old equality means δ=b=0, whence (h−1)/2=e≤h/3, so h≤3. The h=1 zero-defect nonsink alternative is excluded because the zero-defect formula gives e=pc≥1 and h≥3e. At h=3 equality gives e=1; the preceding e=1 classification gives exactly a transitive prefix followed by a rooted triangle. Conversely all such graphs attain L=5.

For h≥5 old equality is impossible, so δ+b≥1 and L≥(3h+3)/2. Improved equality has precisely two cases:

- δ=0,b=1: e=(h−3)/2≤h/3 gives h≤9. Since e=pc is odd, h=7 is impossible on this face without using any global theorem excluding scalar 7. Thus h=5 or 9. At h=5 the e=1 classification gives the reversed-extreme-arc four-vertex core rooted at the former last vertex. At h=9, e=3 and h=3e force C singleton and H(P)=u_w(P)=3, with w in P's final SCC. This is exactly the retained physical interface condition.
- δ=1,b=0: q=1, e=d−1, h=2d−1, and h=d−1+R imply R=d. Thus B rooted at w has zero defect and endpoint 1. If w were a sink then q=d=1 and h=1, excluded. The e=1 zero-face classification gives d=2^t+1, t≥1, hence (h,e,u)=(2^(t+1)+1,2^t,2^(t+1)+2). Adding a immediately before the former last vertex gives exactly the reversed-extreme-arc tournament rooted at its penultimate vertex, with its universal transitive prefix unchanged.

Converses follow from the displayed constructions. Therefore the equality input scalars are precisely 2^k+1, k≥2. Every equality output is 3·2^(k−1)+3≥9, hence composite. For a prime L>5 with h≥5, integer strictness gives L≥(3h+5)/2 and h≤floor((2L−5)/3). No restriction on all prime representability follows.

The wording “7 is impossible” in the producer is dispensable: its immediately following parity proof already supplies the needed exclusion locally, avoiding any dependence on a global H=7 impossibility theorem.

### General two-input transfer

Using u=h+δ,

    N−k(2h+δ−e)=uv−ku+e(k−f)=u(v−k)+e(k−f)≥0.

For h≥5 at a nonsink port, this yields 2N≥(3h+3)k. Equality requires equality in both nonnegative stages, so the first input lies in the complete singleton equality list. Every such input has e>0 and u>0, forcing v=k and f=k. Physicality then makes the second port a sink. Conversely a listed first input and any sink second input do attain equality. Since N=kL and L is composite, all equality targets are composite. In particular two nonsink ports cannot attain equality. This transfer concerns the actual rooted interface, not just h,k.

## Hostile topology and an optional explicit refinement

The exceptional state (9,3,9) cannot be omitted or forced into a transitive-prefix-plus-one-strong-core family. A concrete hostile is C3→C3 rooted in the second triangle. Its state is (9,3,9) and L=15. Here P=C3→{w}, C is a singleton, and w is a sink of P with H(P)=u_w(P)=3. Removing the SCC prefix changes the state by factor 3; that prefix is not transitive. The producer correctly keeps this case.

The physical condition at h=9 can optionally be expanded further. Every H=3 tournament is a single C3 with arbitrary singleton SCCs before and after it: Hamilton counts multiply across SCCs, and a nontrivial strong component of count 3 has at most three vertices by Camion. Consequently:

- If w is a sink in P, then the full source is TT_j→C3→TT_l→C3, rooted in the last triangle, j,l≥0.
- If w is nonsink in P, the old-floor classification gives P=TT_j→C3 rooted in that triangle. Appending C singleton and adding a with sole exit w gives a fixed five-vertex strong core preceded by TT_j, with root a.

These describe the h=9 zero-defect case without introducing any unproved scalar realizability assumption. The producer's interface formulation is already complete and need not be replaced by this refinement.

## Independent computation

`round2_check.py` imports only this audit's own Round1 helper. It never imports either producer checker or `inherited.py`. Its state calculator enumerates all simple directed paths by DFS, accumulating induced-subset totals and root-start counts. On all 1,099 graphs through order 5 these states agree with this audit's independent literal permutation/split counts. SCC decomposition uses independently computed reachability sets. Rooted graph isomorphism checks verify the claimed normalized topologies directly, rather than merely accepting matching scalar states.

Exhaustive root-zero labelled universe through order 6:

- 33,867 graphs; 5,405 sole-exit cases; 1,598 zero-defect nonsinks.
- 152 old-floor equality topology checks; 748 improved-floor equality topology/interface checks.
- 254 distinct physical triples, giving 64,516 ordered transfer tests.
- 60,198 sharpened-transfer tests; all 35 transfer equality pairs have sink second ports and composite targets.
- Literal reversed-extreme-arc family parameters k=2,…,9, yielding h=5,…,513.

Observed improved equality states are exactly (5,1,5), (5,2,6), (9,3,9), (9,4,10), (17,8,18) in this finite universe. These finite controls do not establish the all-order quantifiers; those follow from the proof above.

Reproduce with:

    python3 -I -B round2_check.py > ROUND2_CHECK.json
    python3 -I -B -O round2_check.py > ROUND2_CHECK.optimized.json
    cmp ROUND2_CHECK.json ROUND2_CHECK.optimized.json

Checks use explicit exceptions, not assertions. `AUDITED_INPUTS.json` preserves Round1's original snapshot. `ROUND2_AUDITED_INPUTS.json` records the later Round1 header/status and Round2 snapshot rather than silently rewriting earlier provenance.
