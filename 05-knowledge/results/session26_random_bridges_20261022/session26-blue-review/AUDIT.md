# Independent audit: one-apex colored wiring

Audited source: /workspace/scratch/07d5532f2fce/session26-blue/RESULT.md
SHA-256: 168abb3d2ea09304b7e000c234c824e968959acec20b2155e09133ba413c661a
Verdict: PASS for mathematical statements and proofs as scoped, with two nonmathematical editorial observations below.

## Independent computation

check.py and check_n8.py enumerate literally every tournament with the fixed directed Hamiltonian path 0→1→…→n−1, assigning all non-path arc bits independently. The stratum test intersects all cyclic triangles, rather than generating deletion-apex words. Ordinary class identification subsequently minimizes the apex word over all deletion-transitive vertices, using the uniquely degree-sorted transitive complement. Converse is applied directly to the full adjacency matrix. Blue is tested by a[i,j]=a[n−1−j,n−1−i]. Adjacency is tested by directly toggling every non-path arc.

n : ordinary classes, self-converse classes, tilings, blue tilings
2 : 1, 1, 1, 1
3 : 2, 2, 2, 2
4 : 4, 2, 8, 4
5 : 9, 5, 43, 9
6 : 21, 3, 209, 9
7 : 48, 10, 926, 34
8 : 106, 4, 3898, 18

All four formulas match every computed order. Only-blue classes are exactly the transitive class and the centrally padded C3 when n is odd. Direct internal directed tiling incidences at n=3,…,8 are 2,8,30,52,48,24 respectively. Thus the n≥9 independence qualification is not cosmetic.

## Proof checks

- Consecutive disjoint triples prove the feedback sum inequality on every edge. Isomorphism and converse invariance justify the quotient statement. The general d threshold, neighbor lower bound, loop threshold, and transitive sharpness all follow with the stated rounding.
- Alternating six windows at n=8 have packing number two in each color, giving separate bounds two and two. This does improve the universal sum bound two.
- Cyclic-triangle intersection identifies the possible deletion apices. In the two-apex strong case, the opposite cuts must reach the extreme ends or a nonempty prefix/suffix separates the tournament; the C3(1,1,n−2) exception and its two rooted words follow. Odd automorphism group order, together with the at-most-two apices, gives rigidity for strong orders ≥4. C3 and arbitrary outside padding are treated separately and correctly.
- The strong class and self-converse counts, root-forgetting multiplicities, and prefix/suffix summations yield (2)–(3), including n=2 and odd n=3 exceptions.
- The strong rooted Hamiltonian mass decomposition and subtraction of one exceptional double-apex path count yield (4). The C3 contribution must be H/Aut=1, and the source uses this correctly. Outside padding leaves Hamiltonian counts and automorphisms unchanged.
- The rigid self-converse anti-isomorphism is unique. For a unique apex it fixes that apex and reverses the base. A reversal-pair transversal with maximum i>r has 2^(i−r−1) choices; the all-left transversal supplies 1−b_r. Hence b=J+1−b_r exactly. Summing gives the stated odd strong blue mass, and the double-apex subset-path proof supplies the even case.
- A symmetric path through the unique fixed apex must have it central, whereas the apex-first Hamiltonian path is a nonsymmetric witness. The double-apex witness and C3 padding establish precisely the within-stratum points-of-symmetry statement. Nothing classifies outside nodes.
- Reflection commutes with κ. Every nonfixed reflection pair contributes one fundamental black incidence; reflection-fixed tilings contribute blue incidences. At n≥9 no selected endpoint is adjacent to another selected endpoint, so these incidences are boundary edges, with no further division by two. This argument does not silently apply to small-order internal edges or loops; the document explicitly excludes that inference.
- Both rational OGFs expand from the correct strong-component generating series with two independent padding sequences. Initial a_0=a_1=v_0=v_1=1 and asymptotic constants are correct.
- n=9: A=227, C=19, merged=123, V=15989, B=131, black=(15989−131)/2=7929.
- n=10: A=475, C=5, merged=240, V=64751, B=35, black=(64751−35)/2=32358.

## Editorial observations

1. The relative companion classification link did not resolve at audit time because the companion was still being written. The source repeats enough of the mechanism for this independent audit.
2. For a fully standalone notation declaration, explicitly state b_i=1 iff apex→base_i and 1≤j<i≤m in the Hamiltonian-path sum.

No claim was checked or inferred concerning unrestricted tournament class enumeration, full metagraph wiring, or global points of symmetry. The supplied theorem packet appropriately makes none.

## Final editorial-delta verification

The bit convention and exact summation bounds were added, with no mathematical change. Final source SHA-256: b979dd680b70baea1ee9aedc0eae7b2475a71055ffb1ac50acdb98b652ffe43f.

PASS is rebound to this final source. The companion exists at /workspace/shared/session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md; publishing the specified sibling session26-blue and session26-tournament directories resolves the relative link. Both editorial observations are resolved for that publication layout.

All enumeration scripts and JSON receipts are frozen by SHA256SUMS in this audit directory.
