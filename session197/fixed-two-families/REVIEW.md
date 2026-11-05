# Independent review: three fixed-two-port families

**PASS.** Reviewed the all-order physical construction, response identities, arithmetic formulas, converse symmetry, and finite inverse regions. The conclusion is exact membership and constructive witnesses for the three stated families, not completeness for all strong C2 tournaments, all feedback-two tournaments, or the full Hamilton-path spectrum.

## Frozen source binding

Source directory: `session197/fixed-two-families/`.

- PROOF.md: `c48c6deca5f5a0fd5b8336ec11a297dd8bf8d97aac83bca726390044bb22d545`
- core.py: `dc32dc1ef3647b694015cb4ff1ef99415c162bfbaf58d3b545175b09f0378042`
- families.py: `482b2e553b85ecb299cdab77893a68d65c9f92feb4ffcf8ef499fffe00a3b9e2`
- verify.py: `f7d3ebc7aee5cb3d597d918ac8fa771378c719428f1ea04188d6f0d661f7a58e`
- Both source receipts: `5f418f3dd20b555d774514c88d3317ca8eac2f2bde84a1851e20cb62aa2ab5fd`

The complete source manifest was checked after review. One wording clarification was requested and accepted before this freeze: the P suffix emits the single mask 2 z_i per digit; it does not repeat each digit twice. No mathematical or code correction was needed.

## All-order proof checks

1. The port convention is physical: start with a->b, append a vertex beating precisely the selected old ports, and forget only its newborn exposure. All old nonports point to every later newborn. Deleting the two fixed ports therefore leaves a transitive tournament at every length.
2. The one-apex profile (g(u),u mod 2,u,1) follows from the actual retained-root O/Z moves, including the singleton and leading-bit boundary. U(u)'s two start-one ports and all displayed end/cover responses are consistent with the graph.
3. N's gate retains the FIRST clone, whereas M's relabeled construction retains the LAST clone. The distinction is essential. For M, independently splitting the reserved path gives contributions d0, h(U), D, u, and 4u-1, exactly the displayed R0. The unique forbidden placement in the last contribution is accounted for.
4. P's two intermediate births form V->C3. Deleting the designated old port gives precisely the transitive entry-last/exit-first three-vertex cell. The displayed Dcell matrix and subsequent J move give all four profile coordinates. The argument remains valid for u=1.
5. The O/Z suffix formula is a genuine retained-root iteration. The v=1 boundary preserves the old endpoint count, whereas each nonempty encoded suffix ends in Z and sets it to d. Start count two is preserved because a root-start path cannot leave a new sink or a new vertex whose only old out-neighbor is the already visited root. Cutting at the final Z gives exactly two Hamilton cycles and hence strongness. P's empty suffix already has two cycles.
6. Reversing every arc, swapping ports, and reversing the nonport order sends masks 0<->3 and fixes 1,2. It reverses the three central words 2,03,023 to themselves after that mask change, and exchanges the two binary parameter words. This is an actual graph map, not merely symmetry of the formulas.
7. The inverse bounds are globally valid. Substituting g(x)>=1+2(x mod2) gives the N lower bound with extra nonnegative terms 2ut+2vp+pt; the M bound has extra p+t+4pt+3ut+3vp. The P bound follows from g>=1 and p,t<=1. These imply the advertised finite regions after using parameter symmetry. Equality at powers of two is correct for N and M.

## Independent computation

Copied source normal and Python -O replays both reproduced the exact frozen receipt. Original files were not modified.

A new checker `fixed-two-families/independent.py` constructs every arc directly from the stated mask words. It imports neither the source core nor the source family implementation. It enumerates literal directed path words of every vertex subset, retaining their starting vertices. Thus Hamilton counts, ending/starting counts, deletion counts, Hamilton cycles, and ordered reserved two-path covers are obtained by word enumeration, not by the source's endpoint subset dynamic program.

The independent checks pass:

- 64 directly built one-apex graphs, providing g values independently of its recurrence;
- 386 family graphs with u,v through 12, including exact Hamilton counts, C=2, transitive port deletion, and every edge of the converse maps;
- 94 complete gate profiles with u through 32;
- 176 retained-root suffix checks, including v=1 endpoint preservation;
- 573 finite inverse comparisons against an oversized parameter rectangle.

The independent normal and -O receipts agree byte for byte.

- Independent checker: `763614cf113bc8717c2e105c3baed6555b911bbe0844dbde4f50593aad253fe6`
- Independent receipt: `25fc7f23702aac520ba6a9ca403da7ddb4d2b5b7999b85720bdc030cf9406db7`

These computations corroborate, rather than replace, the all-order proof. No unresolved issue remains.
