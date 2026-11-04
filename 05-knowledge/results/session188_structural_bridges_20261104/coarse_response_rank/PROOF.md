# Exact rank ten for the complete block-uniform cyclic contrast

**All-order exact linear reduction.** Let A, B, C be arbitrary nonempty tournaments on three fixed, disjoint sets of physical names. Compare

    T=C3[A,B,C],    U=C3^op[A,B,C],

with every child interior unchanged. Ports x,y have fixed edge x->y. For three-bit masks a,c, every name in a selected child points to the indicated port, and that port points to every name in an unselected child. Write Delta(a,c)=F_T(a,c)-F_U(a,c).

All 64 contrasts are universal rational linear combinations of these ten actual physical block queries, in this order:

    (0,1), (0,2), (1,0), (1,1), (1,2),
    (1,3), (1,4), (2,0), (2,1), (2,2).

Here bits 1,2,4 mean A,B,C. The coefficients have denominator at most two and absolute value at most two. The complete 64 by 10 doubled integer coefficient matrix is retained in the receipt. Hence all 64 block-uniform contrasts vanish if and only if these ten vanish.

The universal vector space of these contrast functions has dimension exactly ten, already witnessed by transitive children of total order at most six. This is a rank statement about linear reconstruction of all coarse contrast data. It is **not** a lower bound of ten queries for separation, adaptive identification, or a nonlinear decision procedure. Nor does it prove that the ten contrasts cannot vanish simultaneously for some children. It makes no assertion about arbitrary nonuniform named queries or ambient lifts.

## 1. Exact physical-word functional

We use the complete block-uniform word functional developed in the companion `coarse_mod8_boundary/PROOF.md`, Section 2, and the maximal-run argument in `covered_response/PROOF.md`. The formula is restated and implemented independently here so reproduction needs no sibling Python imports.

Let q_i(D) count ordered covers of child D by i nonempty directed paths. Put q_0(D)=0. Define

    R=1-XYZ,
    N=[X,XY,XYZ; XYZ,Y,YZ; XZ,XYZ,Z].

The (i,j) entry of N/R is the generating function of all nonempty cyclic quotient run words starting in child i and ending in child j; each occurrence of a child contributes its variable. Reversing the cyclic quotient transposes N. This works because a directed cyclic word has a unique continuation and additional full turns multiply its weight by XYZ.

For a,c in {0,1}^3, put o=1-a and p=1-c. Define

    L_a=R+1^t N a,     R_o=R+o^t N 1,     M_(o,c)=o^t N c,
    W_(a,c)(N)=L_a (R+M_(o,c)) R_p + L_c M_(p,a) R_o.

The first term represents words with x before y. Their intervening old segment may be empty because x->y. The second term represents y before x and requires a nonempty old segment because y does not point to x. Prefixes and suffixes may be empty. Thus W/R^3 includes every port position and all repeated-child run patterns.

Apply the linear functional

    L(X^rY^sZ^t)=q_r(A)q_s(B)q_t(C).

Each run word is filled by one ordered path cover in each child, in occurrence order. This is a bijection with physical Hamilton paths in the augmented tournament. Therefore

    Delta(a,c)=L(D_(a,c)/R^3),
    D_(a,c)=W_(a,c)(N)-W_(a,c)(N^t).

Every child profile is finite, so the formal series is well-defined under L. The common denominator expands as

    1/R^3 = sum_(k>=0) binom(k+2,2) X^kY^kZ^k.

The statement includes every physical contribution at the queried masks, including cubic and quartic terms. It does not substitute aggregate endpoint skew for the actual two-port response.

## 2. Complete finite polynomial identity

The standard-library checker forms all 64 integer numerator polynomials D_(a,c). Their combined support has 48 monomials. Exact Fraction Gaussian elimination finds the ten listed basis polynomials. For every mask pair, the checker obtains integers b_(a,c),j and checks the full coefficientwise identity

    2 D_(a,c) = sum_(j=1)^10 b_(a,c),j D_j.

All b have absolute value at most four. The receipt lists the complete doubled matrix, including the ten identity rows, in lexicographic query order `(a,c)` with a outermost and c innermost. Applying L after division by the common denominator proves the all-order reconstruction identity, independent of the child sizes or any inequalities between profiles. This is a finite exact algebraic certificate for a universal identity, not extrapolation from a tournament census.

## 3. Sharp linear dimension on actual tournaments

To prove the bound ten cannot be lowered for universal linear reconstruction, use the following ordered triples of transitive child sizes:

    (1,1,1), (1,1,2), (1,2,1), (2,1,1), (1,1,3),
    (1,2,2), (1,3,1), (2,1,2), (1,2,3), (2,1,3).

Form the ten by ten matrix whose rows are these examples and whose columns are the ten basis contrasts, all on their actual physical names. Its exact determinant is

    -3057647616.

The nonzero determinant proves independence as functions of actual tournaments, not merely of formal unrestricted profile variables. The checker obtains every entry by literal Hamilton-path dynamic programming on the augmented tournament, without the run-word formula, and also checks agreement with that formula.

The result concerns unrestricted nonempty children. It does not assert rank ten on every fixed size pattern, the equal-profile subfamily, or the pair-free subfamily.

## 4. Reproducible physical controls and scope

Both exact responses at all 64 masks are retained for each of the ten rank witnesses, plus two nine-old-name controls: C3[C3,C3,C3] and C3[TT3,TT3,C3]. This gives 1,536 independent literal physical counts. The latter retains the known individual zero contrast without confusing it with vanishing of the complete coarse response. All 64 reconstructed contrasts agree with literal counts in every control.

Run:

    python3 verify.py > normal.json
    python3 -O verify.py > optimized.json
    cmp normal.json optimized.json

Only the Python standard library is required. Explicit errors keep every check active under optimization. The four-file packet is `PROOF.md`, `verify.py`, `normal.json`, `optimized.json`.

The theorem provides an exact ten-coordinate description of restricted coarse contrast. A zero ten-vector would establish only equality of all 64 block-uniform queries. Full named-response equality remains a strictly separate question, and no such collision is claimed here.
