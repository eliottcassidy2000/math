# Independent review of the full opposed-binary sparsity theorem

Verdict: **PASS**, bound to:

- `session197-opposed-binary/PROOF.md`: `8eb0b6345476320f00a75304927e3e192a12347ff618f9fab425d4101b305e64`
- `session197-opposed-binary/verify.py`: `1375f06162bc73be425c8e4057a4c9dde5c4d2ed903f5d7dfd123abc2a966938`
- `session197-opposed-binary/receipt.json`: `99da9628e529af69c4de1ec84b32fad394b19ed881618e43e09216c5768b891b`
- `session197-opposed-binary/bounded_image.json`: `72500619242a9043c1a34eeed8d31acdf439b6fa3b2ccd9eb688211f4fe25c22`

## Proof review

Both outer-run bounds use genuine induced triangle-clone tournaments. For an outer 2-run the two apexes give the cycle; for an outer 1-run the adjacent opposite-letter vertex and the specified apex give it. The extra vertices are distinct from the run, including the short two-run cases. The elementary clone-path count and induced-path injection are valid.

For c>=2 the two boundary blocks are disjoint. The four anchors enforce precisely the four nontrivial incidences for the three slots in apex order a,b: last(S0)=2, first(S1)=last(S1)=1, and first(S2)=2. All other boundary vertices are placed on the safe sides of these anchors. Arbitrary assignments of middle vertices remain legal and are recovered from their positions relative to the apexes, proving injectivity. Cases c=0 and c=1 are correctly separated rather than using overlapping anchors.

The presentation inventory is exhaustive and unambiguous: maximal first/last runs determine a,b,sigma,tau and the middle word. Its four contributions sum exactly to L+8L^2*2^R. The all-1 family contributes only the discarded scalar one; all-2 words obey the separate clone bound. Thus weak padding does not create an order-cutoff loophole.

The exponent alpha=log_3(2) and its comparison with 2/3 are correct. The Euler product is indexed by distinct integer generator values and only overcounts repeated factorizations. It gives the stated product-support bound and the prime-omission conclusion. No odd-cycle identity or external strong-order theorem is needed.

The finite H19 exclusion is genuinely all-order because every relevant word is included by the two proved bounds. The displayed aligned-letter escape is a strong graph: its vertex order 0,1,2,3,4,5,0 itself is a directed Hamiltonian cycle, and the independent count is 19.

## Independent replay

I reran the complete frozen verifier under normal Python and `python -O`; outputs agreed and all four frozen digests matched. The replay covers 511 binary words against independent graph DP, 988 induced boundary clone graphs, 8,940 explicit distinct free-slot paths on 396 words, two independent transfers on all 16,392 candidate words through501, and graph-DP checks for all97 positive witnesses. The exact output has96 strong values and139 product values through501. The independently bounded ceiling19 inventory has516 words and no19.

The inherited general slot mechanism and source boundaries are explicitly separated from the new quantitative consequence. No correction or unresolved issue was found in this final snapshot. The theorem applies to opposed incidences only, not the full four-letter two-apex class or the unrestricted tournament spectrum.
