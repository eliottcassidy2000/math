# The arithmetic matroid of the magic branch arrangement

**Status:** PROVED by elementary Smith-normal-form inclusion–exclusion plus an explicitly certified finite table. CITED topology is separately identified. This is a continuation of `RESULT.md`; that frozen base packet is unchanged. No novelty, rational square-point, or general F1-descent claim.

## 1. A universal nonzero-entry theorem for finite abelian groups

Let A be any finite abelian group, written additively, let N=|A|, and put h_d=|A[d]|, where A[d]={a:da=0}. Count labelled3×3 group-magic arrays, with no distinctness assumption, whose nine entries are all nonzero. The normal form in `RESULT.md` identifies their entries, in a fixed permutation, with

    E+i x+j z,  i,j∈{−1,0,1}.

**Theorem B1.** Their number is exactly

    N³−9N²+(20+8h_2)N
      +1−10h_2−4h_3−4h_4−3h_2².                 (1)

Thus only2-,3-,4-torsion can affect this statistic. Equality of the group order alone does not determine it: C4 gives13 whereas C2×C2 gives21; C9 gives224 whereas C3×C3 gives200. These count arrays with all entries nonzero, NOT injective arrays or square-valued arrays.

### Proof and explicit finite certificate

Take the nine integer row vectors

    V={(1,i,j):i,j∈{−1,0,1}}⊂Z³.

For a subset S of these rows with rank r and positive Smith invariant factors d_1,...,d_r, the simultaneous zero set in A³ has size

    N^(3−r) ∏_i h_(d_i).                         (2)

Indeed the integral row and column operations used in Smith reduction induce automorphisms on powers of every abelian group, and a diagonal equation d_i a_i=0 has exactly h_(d_i) solutions. Inclusion–exclusion now counts the complement. The complete512-subset signed table is

| Nonzero Smith factors | Sum of(−1)^|S| |
|---|---:|
| empty |1|
|1|−9|
|1,1|20|
|1,2|8|
|1,1,1|1|
|1,1,2|−10|
|1,1,3|−4|
|1,1,4|−4|
|1,2,2|−3|

Substitution of h_1=1 gives(1). The certificate includes every subset's mask and Smith factors, and the full subset-size/type histogram. Every Smith result is independently recomputed from gcds of all1×1,2×2,3×3 minors. A separate parent implementation uses integer minors only and matches the table. This is a finite exact proof component with an explicit complete universe, not an interpolation from small group orders. ∎

## 2. Toric and additive finite-field counts are different specializations

### 2.1 Generic multiplicative branch complement

In the multiplicative magic torus, let parameters be(e,x,z)∈(F_q^×)³ and delete the nine hypertori

    e x^i z^j =1,  i,j∈{−1,0,1}.

Call the resulting open set U_×. Set n=q−1 and g_d=gcd(n,d). Since F_q^× is cyclic of order n, Theorem B1 gives

    #U_×(F_q)=n³−9n²+(20+8g_2)n
             +1−10g_2−4g_3−4g_4−3g_2².          (3)

It is a characteristic quasipolynomial in n with period12. This counts points of the BASE complement; it must not be multiplied by512 to obtain rational points of the ordinary-square cover.

For characteristic at least5, q lies in one of the four unit classes modulo12, giving

| q mod12 | #U_×(F_q) |
|---:|---|
|1|q³−12q²+57q−105|
|5 or7|q³−12q²+57q−97|
|11|q³−12q²+57q−89|

For example q=3,4,5,7,11 give1,2,13,57,417. Formula(3), not the restricted residue table, covers characteristics2 and3.

### 2.2 Additive special-fiber branch complement

Let U_+ be the set of(E,x,z)∈F_q³ for which all nine E+i x+j z are nonzero. Here A=(F_q,+), so h_d=q if the characteristic divides d, and h_d=1 otherwise. Thus

    #U_+(F_q)=
      q³−4q²+6q−3                         if char=2;
      q³−9q²+24q−16=(q−1)(q−4)²          if char=3;
      q³−9q²+28q−20                      if char≥5.       (4)

The9-coordinate branch arrangement is an arithmetic object, not a single matroid unchanged in all characteristics.

### 2.3 Exactly where the characteristic changes occur

The84 triple determinants of V have absolute-value histogram

    0:8, 1:32, 2:32, 3:4, 4:8.

Consequently its matroid is unchanged in all characteristics at least5. In characteristic3 precisely four previously independent triples become dependent; the nine points(1,i,j) form the affine plane AG(2,3), with12 three-point lines rather than the8 straight real3×3 grid lines. These are exactly Session58's four omitted wraparound lines, whose ordinary magic sum defects were±3x and±3z.

In characteristic2 the nine normals collapse to four distinct normals, with multiplicities1,2,2,4. Any three distinct ones are independent. This gives the U_(3,4) rank structure after removing repeats. Neither a rank count alone nor the abstract three-torus records these changes: the indexed integral characters are the necessary sidecar.

## 3. Discriminant and precise branch meaning in the flat family

For the finite-free square-root algebra from `RESULT.md`,

    B=R[r_1,...,r_9]/(r_i²−z_i),

the squarefree monomials r^ε form an R-basis. Trace pairing entries are zero unless the two exponent vectors coincide; the ε-diagonal entry is512∏z_i^ε_i. Therefore

    disc(B/R)=2^4608 (z_1⋯z_9)^256                 (5)

in this specified basis. There are512 basis vectors, each z_i occurs in256 diagonal factors, and512^512=2^4608. This proves the identity over the universal integer coefficient ring and hence after any base change.

Over Z[1/2,t], the cover is étale exactly where every z_i is invertible. At t=0, z_i are the nine additive magic linear forms. For t invertible, the variable change w_i=1+tz_i identifies z_i=0 with w_i=1, precisely the nine hypertori in §2.1. Thus(3) and(4) count the two base opens over which the corresponding ordinary-square covers are étale in odd characteristic. In characteristic2 the same base complements still make sense, but the covers are nowhere étale.

The source's relative differential criterion is2r_i invertible for every i. At a point where some z_i vanishes, its root also vanishes and that factor is ramified. Ordinary square membership is NOT preserved by the coordinate change w_i=1+tz_i. A finite-free degree512 cover need not have512 rational points above a rational base point.

## 4. No split-torus stratification of these opens over finite fields

A finite stratification by split tori over F_p would imply

    #U(F_(p^r))=Σ_d a_d(p^r−1)^d, a_d∈Z_≥0,

for every r. The additive polynomials(4), rewritten in n=q−1, are

    n³−n²+n                  in characteristic2;
    n³−6n²+9n                in characteristic3;
    n³−6n²+13n               in characteristic≥5.

Each has a negative coefficient. Equality at infinitely many extension sizes forces polynomial equality, impossible for the nonnegative torus-count coefficients.

For U_×, restrict to even extension degrees; then all relevant prime-to-characteristic2-,3-,4-torsion splits. The fixed polynomials in n are

    n³−9n²+28n−28            in characteristic2;
    n³−9n²+36n−51            in characteristic3;
    n³−9n²+36n−59            in characteristic≥5.

The coefficient−9 of n² again rules out a finite stratification by split tori over any finite base field. For a base F_(p^a), use its even extensions, so the argument is unchanged.

This is a scoped obstruction to the specified split-torus stratification of these OPEN COMPLEMENTS. It does not contradict the actual monoidal descent of their ambient torus, does not obstruct blue presentations, and does not assert no possible broader F1 model. If p≥5 and p≡5,7,11 mod12, the odd/even extension counts also demonstrate that U_× has no single counting polynomial over F_p.

## 5. Actual complex topology, with an explicit external theorem

The preceding counts and constructions require no topology. The following is a **CITED application**, not inferred merely from numerical resemblance.

Ye Liu, Tan Nhat Tran and Masahiko Yoshinaga, *G-Tutte polynomials and abelian Lie group arrangements*, arXiv:1707.04551v3, §§2.1–2.2 and5.2, provide the standard common framework for arrangements defined by integer lists over different coefficient groups. Primary full text inspected at https://arxiv.org/html/1707.04551v3 . Formula(2.1) is the hyperplane-complement Poincaré specialization; formula(2.2) is the complex toric-complement specialization. Corollary5.6 identifies the fully split quasipolynomial constituent with the arithmetic characteristic polynomial. The paper's general framework and classical specializations are NOT new claims of this session.

For the present list, the ordinary and arithmetic characteristic polynomials are respectively

    χ_lin(u)=u³−9u²+28u−20,
    χ_ar(u)=u³−9u²+36u−59.

With ambient rank3, the cited formulas become

    P_(U_+(C))(s)=−s³χ_lin(−1/s)
                =1+9s+28s²+20s³,

    P_(U_×(C))(s)=−s³χ_ar(−(1+s)/s)
                =1+12s+57s²+105s³.

Thus their Betti-number lists are(1,9,28,20) and(1,12,57,105), and their Euler characteristics are0 and−59. This is a concrete topological change across the nonproper family, not a contradiction of flatness. No smooth-proper invariance theorem applies to these open families.

## 6. Verification and connection contract

Run `python3 check_branch_arrangements.py` and `python3 -O check_branch_arrangements.py` (SymPy required). Both use explicit exceptions. They compare two Smith methods on all512 subsets; modular elimination in characteristics2,3,5,7; all cyclic group orders1≤n≤120 by direct forbidden-center counting, with full triple enumeration through n=12; and16 additive finite-field cases up to q=81, implemented as elementary-abelian vector groups, never as Z/qZ for nonprime q. Four small root algebras also receive direct trace-determinant checks. The ordinary and optimized certificates/stdout agree byte for byte. `branch-certificate.json` gives every subset and result.

Connection: integral character list → graph/magic torus and additive tangent arrangement; preserved data are the integer relation vectors. Lost data under tangent/rank passage are torsion components of character kernels. The missing sidecar is the full Smith invariant list. The decisive tests are the2-,3-,4-torsion formula, the four characteristic3 wraparound dependencies, and the inequivalent cyclic/elementary-abelian counts at group orders4 and9. Arithmetic squarehood and rational sections remain independent gates.
