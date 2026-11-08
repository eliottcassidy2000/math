# Independent adversarial audit: graft inverse and scalar contexts

Verdict: PASS. No mathematical correction is required within the packet's stated scope. No publication or producer-file edits were performed.

## Proof review

- The inherited all-order closure proof was read in full through its matrix theorem, including all three cases for the rooted two-path cover and the exceptional singleton fourth row. The new determinant expansion has the correct cofactor sign and gives det(K)=s e(qQ-Hb).
- The endpoint-access theorem is valid for tournaments of every nonempty order. For a longest root-starting path and an outside vertex beaten by a path vertex, inserting it after the last such path vertex extends the path. Thus, if maximal, every outside vertex beats every path vertex. There is no outgoing edge from the path's vertex set, contradicting root reachability. Reverse orientations for the end condition. Having both endpoint-access conditions at one root is equivalent to strong connectivity, including the singleton convention.
- For strong nonsingleton gadgets H,q,s,e are positive. When Delta=0, row 4 of the transformed matrix is q/H times row 1; row 2 is independent by its zero last entry and positive first entry. Column 2 supplies one further independent dimension. This proves rank exactly three, not just determinant zero. Delta nonzero proves rank four. The argument does not make an unjustified exact-rank claim about nonstrong gadgets.
- The minimum order-four assertion has a complete finite argument: exceptional order one, nonstrong order two, and the sole strong order-three isomorphism class. Its Delta is zero. The displayed order-four example attains full rank. This does not infer minimality from a truncated census.
- Both elimination equations and every sign in the inverse check. Since inversion holds on a rational basis, this is genuinely an inverse on the entire rational vector space. The source split difference U-d is nonnegative, vanishes for sinks, and is positive for any nonsink via the path (a,v) and a Hamilton path of the complement, including an empty complement. The claimed deletion-ratio inequality consequently has the correct direction.
- The observation matrix has determinant four. Its rows are precisely the empty, singleton, twice-singleton, and edge contexts. Physical source states have determinant minus four, so the response rank is four. Every one-layer first row lies in the stated two-dimensional span, including singleton gadgets. Singleton and edge first rows are independent, and the original-count row raises the rank to three. Depth two is necessary and sufficient in the stated rational-linear sense.
- Equality of states is equivalent to equality of all scalar graft-context responses: closure gives the forward direction and four explicit contexts give the converse. This says neither that tournament isomorphism is recovered nor that nonlinear encodings require four dimensions.

## Prior work and novelty boundary

Read Session225 PROOF.md and, using git show HEAD:path in session226-work, adaptive_state_induction_20261103/BLOCK_SWITCH.md. B2 explicitly proves rational linear dimension four for its R/N/J1 language with a response determinant of -32. The new packet correctly acknowledges this and limits its observer contribution to the graft language and the depth distinction. Its determinant factorization, inverse, rank/deletion-ratio analysis are developments of the inherited general matrix, not a claim to have first introduced the four-state closure. No literature-priority assertion is warranted or made.

## Executed checks

All commands completed successfully:

- Producer verify.py, normal isolated mode and optimized isolated mode. Both outputs are byte-identical to producer CHECK.json.
- inherited_verify.py, normal and optimized isolated modes. Outputs agree byte-for-byte. The inherited script itself agrees byte-for-byte with Session225's verify.py.
- independent_check.py, normal and optimized isolated modes. Outputs agree byte-for-byte. This program imports no producer module. It generates tournaments as edge sets, computes path counts via independent recursive completion and endpoint permutation counts, and computes determinants by the Leibniz formula rather than producer elimination.

Independent coverage: all 5,405 labeled rooted presentations at orders one through five, including every root. Checked root-start/reachability equivalence, reversed endpoint equivalence, strongness equivalence, split/sink boundary, determinant factorization, strong singular rank-three minors, and minimal full-rank order four. Checked 10,304 inverse images of rational coordinate basis vectors, hence every tested invertible matrix on all rational inputs. Checked 1,140 literal grafts using every rooted gadget through order four on all four physical spanning sources, including singleton gadgets and sources.

The producer additionally checks all 1,099 root-zero source contexts literally and all nine gadget statistics against permutation/split enumeration. Independent checks are corroboration, not the basis of all-order validity.

## Audited input hashes (SHA-256)

- PROOF.md: ebb93392b7230ecae282371cb486478acf623d40916da6a025a7101f2c32125d
- verify.py: 783937df54a60e991bdf6438ee1ccb21512fcc8bb18fb8a999f33125532135f1
- inherited_verify.py: 9e834801a3dcaeb74a22a9fd7787d2bdf8c6c1b49e594f7177f0031bdfe088d5

Receipts: replay.json, replay.optimized.json, inherited.json, inherited.optimized.json, independent.json, independent.optimized.json. Independent executable: independent_check.py.
