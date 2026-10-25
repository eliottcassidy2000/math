# Direct polynomial verification of Theorem 3

This supplements repair.md and closes its boundary-specialization caveat. It does not change the correction to the reduction preceding (5.5).

Source: Ajai Choudhry, arXiv:2508.07806v1, equation (5.1), printed p.11, https://arxiv.org/pdf/2508.07806v1. The eight input polynomials were independently transcribed from the printed formula and compared with its HTML rendering. The certificate does not rely on substituting into the rational (5.4) construction.

Put u=u1, v=u2, x=v1, y=v2, and abbreviate
V=x²+y², W=x²−y², H=x⁴−6x²y²+y⁴,
E=u²−v², F=u²+v², T=x²y².

Here is a polynomial square-root matrix R for precisely the row/column ordering printed in (5.1). Each numbered line is a row, with semicolons separating its four entries:

1. HWF; 2uvVH; 2W(4TE+uvV²); 2W(4TE−uvV²).
2. HVE; 2xyFH; 2xy(V²E+4uvW²); 2xy(V²E−4uvW²).
3. W(V²E−16uvT); 4xy(W²E−uvV²); 4xyFWV; 4xyWV(E−2uv).
4. W(V²E+16uvT); 4xy(W²E+uvV²); 4xyWV(E+2uv); 4xyFWV.

For every i,j, direct expansion gives a_i+b_j=R_ij² as an identity in Z[u,v,x,y]. The executable certificate checks these sixteen displayed identities exactly. Separately it factors each a_i+b_j, verifies a square integer leading coefficient and even multiplicities for all polynomial factors, constructs another polynomial root, and verifies the resulting identity by expansion.

Consequently, the square-valued conclusion of printed Theorem 3 holds for all integer parameter choices, including every boundary where the earlier rational parametrization is undefined. There are no denominator restrictions in this check. Repeated entries or zero values can occur for some parameter choices, so the claim proved is all sixteen indexed sums are squares; this does not assert four distinct entries in both sets for every specialization.

Files: polynomial_certificate.py is executable; polynomial_results.json contains the verified factor-derived roots and scope. All computations use exact arithmetic. No full-paper reproduction is included.
