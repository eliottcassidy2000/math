# Additional parity structure, without an impossibility claim

Status: independently audited theorem and literal counterexample. The baseline universal-pair reciprocity theorem is in the separate frozen carriers packet.

## Two-adic refinement

For a tournament U, its unrooted ordered path-cover polynomial B_U(k) has expansion

B_U(k)=sum_pi 2^(c_nontrivial(pi)) k^(c_total(pi)),

where pi ranges over permutations whose nontrivial cycles are directed odd cycles of U. Singleton cycles are allowed. One derivation uses the squarefree path series f and the tournament matrix identity: f=det(I+A^T X)/(det(I-A X)det(I-X)). Taking the formal logarithm cancels even cycles, leaving sum_i x_i plus twice the sum of directed odd-cycle monomials. Exponentiating k log(f) gives the formula. All operations are in a nilpotent squarefree ring over the rationals, so the logarithm is finite and coefficient extraction is exact.

For even order n>2, B_U(2) is divisible by8. The identity permutation contributes2^n; every nonidentity term has an even number of cycles, at least2, and at least one nontrivial cycle, so its power of2 is at least3.

A hypothetical universal pair has opposite source-order parity and common C=B_U(2), b=B_(U-root)(2). Its even-order source and its even-order deletion both have order>2: the previously proved positive odd d divisible3 excludes deletion order<=2, and positive even starts exclude the exceptional tiny sources. Hence C,b are divisible by8. Combining the already proved divisibility by12 gives

C,b divisible by24.

The common Z descendant has C_Z=2C+2T-2L, also divisible by8. Hence T congruent L modulo4. The baseline L congruent6 modulo12 and T divisible6 sharpen this to

L,T congruent6 modulo12.

These restrictions do not force an infinite divisibility tower. The congruence class C,b=0 modulo8 and L,T=2 modulo4 is closed under O/Z/J (use the inherited formulas). A physical modular source already realizes the combined necessary congruences, with profile

(123,4,33,630,123,33,2760,672,1422).

It is not itself a paired source. The tempting stronger statement that every strong even-order tournament has B(2) divisible16 is false: the strong six-vertex tournament (48,5,9,3,14,30) has H=33 and B(2)=360.

## Equal global moments do not force order parity

There are explicit tournaments of opposite parity with identical H and B(2): orders17 and18, both H=2025 and B(2)=5971968. Their complete adjacency rows are in GLOBAL_PRODUCT_CERTIFICATE.json, independently evaluated by verify_global_product.cpp.

The construction uses ordinal sums, for which H and B(2) multiply. Let A6 have (H,B2)=(45,432), B9=(75,3456), C9=(27,1728). Then

U=A6 ordinal-sum A6 ordinal-sum TT5,
V=B9 ordinal-sum C9.

H(U)=45^2=75*27=H(V), and B2(U)=432^2*32=3456*1728=B2(V).

The finite complete unrooted bank through order9 contains no direct opposite-parity equal (H,B2) pair; the product witness shows that this finite observation must not be extrapolated. It also rules out an attempted proof of universal-pair nonexistence based solely on those two global moments. This particular pair does not supply bounded-start sources or matching rooted coordinates.
