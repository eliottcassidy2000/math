# Coprime mobile module of order five

Status: exact construction and transfer identities verified; independent audit requested. This does not prove Chan–Pak completeness.

## Module and binary compatibility

Let B have adjacency masks [28,5,16,6,10], with root pointing only to vertex0. Literal path enumeration gives H(B)=13, Q_B(2)=104 and per-start two-bin weights [26,18,14,18,15]. Thus its permitted-start weight is26=13·2, and Q/H=8. It replaces the three bits010 while preserving normalized U/d. The complementary root neighborhood replaces101. In original labeling [18,8,3,5,14] these subsets are {4} and {0,1,2,3}. There are no other nontrivial integer-multiple-of13 subsets in this representative.

Write S0,S1 for the inherited singleton transfer matrices, A=S1S0, M=matrix(B,1), and G=M−13S0S1S0. Exact local enumeration gives

M=[0 0 0 0 0 0;4 13 0 0 0 0;11 0 11 0 11 0;22 77 26 90 22 77;2 0 2 0 2 0;0 14 0 14 4 27].

G has rank3. Let F denote the inherited triangle defect T1−3A. For every nonempty binary word w, with singleton transfer W, all four products FWF, FWG, GWF, GWG vanish. This is stronger than a finite spacer check: form K=image(F)+image(G) and L=S0K+S1K. Exact rational rank calculations prove dim(L)=2, FL=GL=0, and S0L,S1L⊆L. Induction on word length proves the statement. The verifier includes this finite invariant-subspace certificate.

## An unbounded physical fiber with deletion count13

For r≥2 and 0≤j≤r−2, replace bits2j through2j+2 of (01)^r by B with the prescribed root neighborhood. All between-block arcs point forward. The resulting rooted tournament has2r+3 vertices and

D=13(4^r+6r−1)/9 + 4^(j+1)−10 + 4·1_(j=0),
E=0, d=13, U=13(4^r+2)/3.

Proof of D: with r=j+k+2 its defect term is ell A^k S1 G A^j v. The inherited identity A^4−6A^3+9A²−4A=0 supplies the order-four recurrence in each variable. The verifier checks its4×4 initial grid equals4^(j+1)−10+4·1_(j=0), whose expression satisfies the same recurrence. The baseline is13 ell A^r v. The sidecar follows directly from binary compatibility and the untouched terminal1.

For1≤j≤r−3, remove the initial singleton0 and final singleton1 to obtain a strong tournament on2r+1 vertices, with the same D and

E=d=13, U=13(2^(2r−1)+1)/3.

Both boundary singletons and the new boundary blocks remain untouched under this restriction. The new first singleton1 is beaten by the root and new last singleton0 beats the root, proving strongness through the ordered blocks. The B module itself is strongly connected, as can be checked directly from its arcs. For fixed r these r−3 scalar values are distinct. At r7 they are23731,23779,23971,24739 and their gcd is1.

## Mixed defects: an escape and its boundary

Any separated collection of triangle01 and B010 replacements in a binary word has exact additive expansion

D=(product H_i)D_base + sum_i (product_(h≠i) H_h) delta_i,

provided at least one ordinary singleton separates every two modules. Terms with multiple defects vanish by the invariant-space identity. With t triangle modules and u B modules, t,u≥1, every scalar is divisible by3^(t−1)13^(u−1). Thus the H13 module removes the universal restriction to powers of3, but repeated use of finitely many module types still has a product/lcm divisibility obstruction. One triangle and one H13 module need not force either divisor. This is not a prime-coverage theorem or a prime-capable exponential-fiber theorem.

## Exhaustive finite module bank

search.py examines every labeled tournament on n≤6 vertices, using subset dynamic programming for Hamilton paths and ordered two-bin counts. For each binary-compatible tournament Q/H=2^k it exhausts all root-outneighbor subsets and records representatives by (H,Q/H,sorted alpha,available intermediate m values). This grouping is a statistics grouping, not a claim of isomorphism classification.

Compatible (H,Q/H) and labeled counts:
- n1: (1,2):1.
- n2: (1,4):2.
- n3: (1,8):6; (3,4):2.
- n4: (1,16):24; (3,8):16.
- n5: (1,32):120; (3,16):120; (13,8):120; (15,8):24.
- n6: (1,64):720; (3,32):960; (9,16):80; (13,16):1440; (15,16):288.

In particular no H11 binary-compatible module occurs through order6. The H15 order5 representative is regular, alpha_v=21 for all five vertices, so it has no intermediate alpha(S)=15m. Order6 contributes no new H beyond1,3,9,13,15. No nonexistence beyond order6 is asserted.

Artifacts: search.py, bank1.json through bank6.json, verify.py, VERIFY.json. The independent module check in verify.py uses literal permutations, separate from search.py's dynamic programming. The transfer implementation is inherited from ../session228-fiber-switch/transfer.py.
