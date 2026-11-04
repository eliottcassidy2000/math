# Independent review of the cyclic-child disjoint kernel

**PASS, Sections1--3 only.** The exact disjoint-quartic kernel, the arbitrary-module same-port zero identity, and restoration of singly marked endpoint-balanced exterior children are correct. The alternating-hub positivity and cube corollaries in Sections4--6 have a separate reviewer; this audit does not substitute for that review.

The pinned proof is2c625e1a7c7cff166813d7f7243220b160fe04c3f32418e75bb07f267d9fb1be. No mathematical correction was needed.

## Algebraic review

1. In the squarefree exterior path algebra, every excursion from the module back to itself has the same series tau. Eliminating the exterior therefore replaces the module adjacency by B_A+tau J, including its diagonal. A diagonal excursion would repeat its module endpoint and is killed by its squared vertex variable. Repeated exterior support is likewise killed. The resolvent and the full-path Schur expressions therefore count exactly the allowed paths, without a lost factorial or an independent-marginal assumption.
2. Conversing the module transposes its effective path matrix, exchanges row and column sums and leaves H fixed. For p->q->r->p, direct enumeration of every internal word of length1,2,3 gives the displayed M skew and all three subsequent C/R identities. The zero and factor2 are exact integers, not merely parity statements. The tau dependence in the M skew is essential and is retained.
3. Substituting the complete eight-term disjoint-quartic normal form and extracting all three module names gives exactly the four terms of the claimed exterior kernel. The independent formal checker treats alpha,beta,tau,E_i,E_j,S_i,S_j,G_ij,G_ji,ell_i,ell_j,rho_i,rho_j as algebraically independent commuting variables and only imposes squarefreeness on the three module names. Its resulting polynomial has ten monomials and agrees identically with the claimed kernel. Thus this is a formal identity, not a finite collection of exterior examples.
4. For arbitrary module size, the same-port normal form factors through C_p R_q+C_q R_p. Its exterior factor can be written

       alpha(rho_i S_j+rho_j S_i)
       +beta(ell_j E_i+ell_i E_j)
       +2 alpha beta H(rho_i ell_j+rho_j ell_i).

   Both this factor and the C/R factor are invariant under transpose. This proves the arbitrary-module same-port zero claim without extrapolating from the triangle calculation.
5. In the quotient-converse application, each exterior marked child contains exactly one physical mark. Its contribution to each ordered-cover word is therefore a single specified endpoint profile, with no missing multiplicity: run slots remain distinguished. Whole converse swaps E and S; endpoint balance permits restoration. Every unmarked child uses its ordinary converse-invariant ordered-cover profile. No child contains both exterior marks in this theorem. That restriction correctly avoids requiring a joint marked profile.

## Independent checks

`formal.py` enumerates all three-name internal words and performs polynomial arithmetic in a mixed ring: three nilpotent module variables and thirteen unrestricted exterior variables. It checks all four displayed internal identities, the complete exterior kernel, and the same-port zero before extracting module support. It imports no producer code and does not use a symbolic-algebra package.

`literal.cpp` separately counts Hamilton paths by full physical subset DP, without run compression or the kernel formula. It retains every query value for four complete sixteen-query interpolation records in three orientations: source, local triangle converse, and quotient converse. This is192 literal augmented-tournament counts on11 and17 physical names. The disjoint quartics are60 for C3[C3,C3,C3] and-109584 for the regular five-quotient with C3 children. Both corresponding same-port quartics vanish. Local and root contrasts agree exactly, independently checking the balanced-restoration step.

All arithmetic is integral. The largest literal tournament has17 names, so its count is at most17!, and the signed sixteen-term interpolation is safely within signed64-bit range. No response coefficient is silently dropped: every finite difference is computed from all sixteen actual query values.

## Reproduction and scope

Run `python audit.py > normal.json` and `python -O audit.py > optimized.json`, then compare. Both modes recompile and rerun the literal checker in a temporary directory and execute every formal check with explicit exceptions. Publish REVIEW.md, audit.py, formal.py, literal.cpp, normal.json, optimized.json and manifest.sha256. The standalone exploratory executable and literal.json beside them are not needed for publication.

This review does not assert that every prime quotient has a nonzero kernel, that all larger children can replace C3, or that an arbitrary strong outer ambient preserves the exact integer sign. Those boundaries remain explicit in the producer proof.
