# Independent audit: strong scale25 four-prefix extension

## Result

**PASS.** The frozen50 paired and50 terminal certificates, together with bases42 through1061, prove the all-size unscaled-prefix cycle theorem for

`P_0=(1,3)`,

`P_(l+1)=(1,3,6,19), (25*p_j+(-1)^j*8)_j`.

For `S_l=(85*25^l-1)/2`, the deterministic paired constructor has common prefix P_l for every n>=S_l. Every N>=S_l, l>=1, has a square-sum Hamiltonian cycle containing original P_l consecutively. Its length is4l+2 vertices.

This statement concerns the offset8 four-prefix family. It does not assert a cycle extension for the separate offset6 prefix(1,3,6,30), nor for arbitrary numerical square-sum rays.

## Independent evidence

`audit_independent.py` is a standalone standard-library auditor. It imports no producer code or old transition table. The required paths are `--pair-dir`, `--cycle-dir`, `--small-base-file`, and `--base-dir`. Its optional `--partial` mode explicitly declines to claim completeness; the final evidence uses full mode.

Normal and optimized full runs match byte-for-byte with **51,028,631 explicit checks**:

- all50 paired and50 terminal states, complete physical catalog and filename dispatch;
- all1020 bases42 through1061;
-100 full paired expansions, covering all50 states at two depths;
-103 complete streamed cycles, through the depth-four threshold16,601,562;
- eight hostile word-certificate corruptions;
- an additional mutation flipping exactly one head-color entry, rejected specifically by the common-head parity check;
- a concrete independently expanded Hamiltonian-pair counterexample to necessity of the old within-word head/tail linkage.

All100 certificate hashes match `FROZEN_SOURCE_MANIFEST.json`, and all103 independent cycle hashes agree with the producer's separate run. These comparisons are recorded in `producer_comparison.json`, not used to manufacture the independent result.

At N16,601,562, the independently checked cycle begins with the full18-vertex P4:

`1,3,6,19,17,83,142,483,417,2083,3542,12083,10417,52083,88542,302083,260417,1302083`.

Its decimal-lines SHA-256 is `2bbe959be4176097b648e42ac56a7647697e68cc3474019631c46e7beca10231`.

## Universal semantics and the essential relaxation

The affine coverage proof is unchanged: the25 offsets partition each source vertex's band `[25v-12,25v+12]`, input compatibility fixes its representative, and the forced member choice is `25+(-1)^(n+1)c<=d`. Singleton1 through12 complete the interval.

The exact protected opening is

`single1, single3, single6, head6, forward intact block8`.

The fourth label is25-6=19, which is the head of the split offset6 block—not an additional low singleton. The following intact block starts25-8=17. The connector sums are4,9,25,36, and the entire input prefix survives in that block. The split offset6 is tracked as part of exact ownership, so the high literal19 is neither omitted from nor duplicated in its source block.

Each head has its own actual output-position parity. Each tail/whole block has a constant output-minus-input position shift. The independent verifier compares every head color across the two words and every tail/whole shift across the two words. It does **not** require a head to have its tail's shift within one word. This is the necessary/sufficient compatibility criterion under the disjoint source-band decomposition; see `FREE_HEAD_LEMMA.md` for the proof.

The distinction is substantive. At input parity0, residue18, the head at offset6 has output parity0, while its tail shift is0. The old rule would require head parity1, but the actual pair on1068/1069 is compatible and square-sum Hamiltonian. Thus the relaxed certificate system is not merely a relabeling of the earlier stronger condition.

All structural obligations remain independently checked: atom identities and arities, exact coverage, input choices, orientations, symbolic output lengths, every connector square, protected opening and required endpoints. Paired outputs end at8 for even order or10 for odd order. The terminal consumer ends at8 and closes to1; as a single final cycle, it has no output-pair parity obligation.

## All-size construction

Use certified bases only at42 through1061. For every larger n, use

`m=floor((n-12)/25)`, `d=n-25m`,

and the new paired state `(m mod2,d)`. Then42<=m<n and12<=d<=36, proving termination and all-size weak B supply. Since `S_(l+1)=25*S_l+12`, the same constructor retains P_l whenever n>=S_l. Applying the matching terminal state gives the cycle theorem.

The cutoff1062 is mandatory: retaining old weak bases above1061 would not establish the longer-prefix induction. The full replay checks the boundary cases and every state at both initial and recursive depths.

## Arithmetic simplicity, not monotonicity

Writing Q=25^t, the four-coordinate block is

`((2Q+1)/3,(10Q-1)/3,(17Q+1)/3,(58Q-1)/3)`.

P_l consists of the first l complete blocks followed by the first two entries of the next. Its root multipliers are2,3,5,6, scaled by5^t. Lifted vertices have residues17 or8 modulo25, while the newly prepended literals1,3,6,19 have neither residue. The bounded offset makes the lift injective on any simple positive input, so induction proves simplicity at every depth.

The ray is not globally increasing: every block boundary descends by `(2-8Q)/3`, starting19→17. Simplicity is supplied by residue exclusion, not by an incorrect monotonicity assumption.

## Publication scope

The frozen basic two-prefix scale25 audit is unchanged. The present packet is a separate stronger prescribed family using a justified semantic relaxation. Neither packet requires old scale49 transition tables. Mixed scale25/scale49 constructions may invoke both audited implications, but their additional parameter formulas and scope should be stated separately.
