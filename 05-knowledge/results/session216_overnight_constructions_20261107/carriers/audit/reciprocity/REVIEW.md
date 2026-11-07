# Independent audit: rooted reciprocity and universal paired sources

Status: APPROVED in the stated necessary-condition scope. The proof does not establish existence or nonexistence of every universal pair, numerical tail coverage, or an obstruction to arbitrary paired-source constructions.

## Interface and inherited dependency

The inherited S203 synchronization criterion equates all-word O/Z/J adjusted-base synchronization with equality of the eight nonstart coordinates (h,e,L,R,d,C,b,T) and start difference s1-s2=2, for nontrivial positive-start sources in the two staggered residue classes. Its proof and the rooted-cover closure are frozen under dependencies/ with canonical provenance. This justifies using the eight-coordinate premise as the universal-pair interface, not merely a sufficient special subclass.

Here C is the ordinary ordered two-bin cover count, b the corresponding root-deletion count, and T the ordered three-bin cover count whose first nonempty path starts at the root. Ordinary nondistinguished bins may be empty. The interpretation F(3)=T is essential. The manuscript refers to these inherited coordinates; readers should retain those definitions or link the dependency explicitly.

## Polynomial integrality and reciprocity

Grouping other nonempty path components as an unordered set gives integer coefficients in the falling-factorial basis (k-1)_j. Since each such falling factorial lies in Z[k], F,G lie in Z[k], not just the ring of integer-valued polynomials. This distinction is load-bearing for the divisibility claims.

The squarefree-ring proof is sound. AX is nilpotent; M=I-AX is invertible. With X diagonal, X^2=0 and M(I-X)=M-X. The tournament relation yields I+A^T X=M(I-X)+1 x^T. The rank-one inverse denominator is exactly P because x^T(I+X)=x^T. Multiplication by x_r kills the extra diagonal factor, giving E_r(-x)=-S_r/P. Reversing paths identifies the unrooted generating functions for A and A^T, yielding P(-x)=1/P. All inverses exist in the indicated squarefree ring because their constants are1. Repeated vertices vanish, and the outer x_r excludes repeated visits to the starting/ending root.

Taking top-degree coefficients in S_r P^(k-1) and E_r P^(k-1) proves F(1-k)=(-1)^(n-1)G(k). Negative exponents are legitimate unit powers; their coefficients agree with the polynomial extension. The singleton case has F=G=1 and satisfies the identity.

## Necessary conditions

For same-order parity, D=F1-F2 has roots -1,0,2,3. Division by their monic linear factors takes place in Z[k]. Evaluation at1 is therefore a multiple of4, ruling out start difference2.

For opposite parity, D=(k-2)(k-3)Q with integer-coefficient Q. Its values at0 and-1 yield 2 sigma e=6Q(0) and 2 sigma L=12Q(-1), hence3|e and6|L. Comparing F(1) and F(-1) modulo2 forces both starts even. Together with the width-eight staggered labels, exactly(4,2),(10,8),(16,14),(22,20) remain. This excludes singleton inputs before the nontrivial O/Z/J transfer formulas are used, so the singleton Z exception is not being ignored.

The inherited closure preserves equal nonstart coordinates and starts after every common letter, while adding one vertex to both sources preserves their relative order parity. Applying the same local divisibility statement to descendants Z,J,ZJ,JJ gives all eight multiples of3. The chronological formulas e_Z=d, e_J=h-e, L_J=C-L, L_Z=L+b, e_ZJ=h+R-e-d, L_JJ=C+2T-L are correct.

Since D(1)=2, Q(1)=1. Thus Q(-1) is odd and L=6 modulo12. Applying this to Z,J,JJ gives12|b,12|C,6|T. The usual tournament parity theorem makes h,d positive odd multiples of3. These deductions are necessary conditions only; modular compatibility is not a physical realization theorem.

## Independent exact controls

independent.py does not use the producer's SymPy or ordered-bin convolution construction. It enumerates unordered vertex-set partitions, counts directed paths literally inside each block, singles out the root-containing block, and forms the integer falling-factorial polynomial. It checks all5405 rooted labeled tournaments through order5, all coefficients of the reciprocity identity, and every relevant profile evaluation. Same-parity collision controls cover5340 repeated keys.

It separately computes inverse unrooted path series by squarefree subset convolution and checks P(-x)=P^-1 and S_r/P=-E_r(-x) coefficientwise, with168,146 rooted subset checks. Normal and optimized receipts are byte-identical and use explicit exceptions rather than assertions. Producer verification also replays exactly in both modes, with its1099 root-zero labeled controls and1034 collisions. The producer verifier requires SymPy; the independent verifier uses only Python's standard library.

## Finite search replay and scope

The root-incidence generator reproduces its exact receipt:1,823,702 orientations,63,546 bounded-start retained presentations,809 necessary-modular candidates. candidates.cpp independently reconstructs each of those809 literal graphs with the new root at a different vertex index, and computes all nine profile entries by direct path/cover DP. All agree, normally and in a debug build. Exactly425 candidates pass L=6 mod12,12|C,b and6|T; their rows match ROOT_STRONG_MODULAR_CANDIDATES.txt exactly.

The deep search reproduces809 sources, depth7,2,653,520 visited presentations,939,817 keys,0 hits. Its key retains all eight nonstart coordinates plus s+2(s mod3). A collision between its two residue slots is exactly the required start-gap-two pair, not merely scalar adjusted-base equality. Every nontrivial gate preserves the positive bounded start and residue. The bank is supplied, not asserted exhaustive. In particular sources that enter the modular class only after other prefixes remain outside this search.

Source order is at most9 and seven appended vertices give at most16 vertices. Hamilton paths are bounded by16!, ordinary two-bin covers by17*16!, and three-bin covers by binomial(18,2)*16!, all within signed64-bit range. No arbitrary-order arithmetic safety is asserted.

The historical C3 and circulant5 external-family searches also reproduce exact receipts:3060 and15120 graphs respectively, maximum order12. These finite failures are now subsumed for their same-order universal-pair goal by the all-order parity obstruction; they are not needed for the theorem.

## Reproduction and frozen files

    python independent.py > independent-replay.json
    cmp independent-replay.json INDEPENDENT.json
    g++ -O3 -std=c++17 candidates.cpp -o /tmp/session216-candidates-audit
    /tmp/session216-candidates-audit frozen/BANK.txt frozen/ROOT_MODULAR_CANDIDATES.txt
    cd frozen && sha256sum -c ../FROZEN_SHA256SUMS

frozen/ preserves the reviewed producer manuscript, scripts and receipts. FROZEN_SHA256SUMS pins exact bytes; PROVENANCE.json pins the inherited synchronization and closure dependencies. Replay receipts are retained separately. No source files were changed and no publishing was performed by this auditor. Parent remains sole publisher.
