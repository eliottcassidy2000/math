# Four noncirculant prime-nine quotients with arbitrary children

**Status: finite exact certificates, independent review pending.** This packet concerns four explicit row tuples only. It does not assume any incoming order-nine census is complete and does not raise the global highest-changed-node frontier beyond nine.

## 1. Statement and four fixed specimens

Use quotient names0 through8 and outgoing row masks. Set

    Q0=(192,289,99,71,15,153,306,94,189),
    Q1=(320,193,83,71,43,15,432,61,190),
    Q2=(448,209,195,71,45,15,432,312,62),
    Q3=(384,177,147,135,105,77,15,368,126).

For each displayed Q, let A_i,B_i be arbitrary nonempty children on exactly the same physical sets, with q_h(A_i)=q_h(B_i) for every i,h. Matching is childwise; different children may have different profiles and orders. A child-uniform actual two-port query separates Q[A_i] and Q^op[B_i] in every fixed ambient.

The physical ports obey x->y. Bits a_i,c_i mean all names in child i point to x,y. Both ports inherit the original core root's incidence at each exterior name. Source and target are never independently relabeled.

Every proper nonsingleton subset is explicitly checked for the module property; none is a module in these four quotients. Their Hamilton counts are2251,2239,3027,2479, so the four are pairwise nonisomorphic. Their sorted score sequences are respectively

    (2,3,4,4,4,4,4,5,6),
    (2,3,4,4,4,4,4,5,6),
    (3,4,4,4,4,4,4,4,5),
    (2,4,4,4,4,4,4,4,6).

They are genuinely outside the regular circulant-nine family. These observations classify only the displayed four specimens.

## 2. Singleton-position reduction and two nonanchor residuals

All four quotients have endpoint parity p=0 and all-nonsingleton second-endpoint residue sigma=0. Thus those two uniform tests alone miss them. The independently reviewed all-order singleton-sensitive anchor theorem applies without a quotient-order restriction.

Let v_j be one exactly when child j is nonsingleton, and let V be its physical child-label mask. Exact endpoint and repeated-run primitives give the first and second anchor contrasts through the inherited affine formulas

    Delta(0,e_i)=d_i(V) modulo4,
    Delta(e_i,0)=2t_i(V)+2d_i(V) modulo8.

Every nonzero first anchor has complete contrast(2,0,...) modulo four. Every nonzero second anchor has complete contrast(4,0,...) modulo eight; p_i=0 removes the odd-total cycle correction. These gaps survive every fixed ambient.

The complete512-pattern evaluations leave precisely

    Q0: V=254,511;
    Q1: V=511;
    Q2: V=511;
    Q3: V=54,511.

For the two singleton-containing residuals, the same actual query(a,c)=(1,2) gives a gap two modulo four. Specifically, if K0 is the bare two-port Hamilton count and K_j counts words with child j twice, the complete modulo-four reduction gives

    Delta H = Delta K0 +2 sum_(j in V) Delta K_j modulo4.

At V=254 for Q0 and V=54 for Q3 its value is2. The common odd child Hamilton normalization fixes this residue. Triple-run contributions vanish by the return-segment involution, and all later terms vanish modulo four. The complete augmented-profile contrast is consequently(2,0,...) modulo four and survives every ambient. Actual singleton/pair physical graphs and cyclic-root ambient lifts are retained independently.

Thus every assignment containing a singleton is covered, leaving only V=511.

## 3. Small cycle-safe certificates cover every nontrivial profile

For a nontrivial child the complete normalized profile modulo eight is one of

    P=(1,2,0), E=(1,6,4), O=(1,6,6), C=(1,2,2),

with higher coordinates zero. Set b_i=1 for E,O and c_i=1 for O,C. In particular c_i is actual child-order parity. For an actual block query, let K_i,L_i,K_ij count the repeated-run patterns with child i twice, child i three times, and children i,j twice, respectively. All physical port incidences and orders remain present. Then

    Delta H / product_i H(A_i)
      = Delta K0 +2 sum_i Delta K_i +4 sum_(i<j) Delta K_ij
        +4 sum_i b_i Delta K_i +2 sum_i c_i Delta L_i modulo8.   (1)

The inherited return-segment argument makes L_i even. Each retained baseline contrast is divisible by four. Therefore the common odd normalization fixes the contrast and its division by four is an exact affine bit in eighteen variables.

A query is called safe only when the two exact bare quotient-extension Hamilton-cycle counts have the same parity. This supplies the full-profile check for odd total child order: the inherited formulas give Delta q2=Delta q3=2 Delta H-4 Delta Z modulo8, and safe queries have Delta Z even. Every safe nonzero gap has complete profile(4,0,...) modulo8, for even or odd total order, and survives every fixed ambient.

The finite certificates use these safe actual-query batteries:

    Q0: (0,7),(1,3),(272,34),(66,66),(144,272),(258,10),(24,64);
    Q1: (36,32),(6,260),(260,10),(10,260),(144,40);
    Q2: (3,132);
    Q3: (3,1),(272,34),(272,129),(33,17).

For each row, the displayed affine contrast bits sum to one. Thus no assignment of the eighteen profile bits can annihilate the whole battery. Every listed query is checked for cycle safety separately; local nonvanishing is not mistaken for arbitrary-ambient separation.

The Q2 certificate is especially simple: at(a,c)=(3,132), the baseline contrast is4 modulo8 and all eighteen profile derivative bits vanish. Both bare cycle counts are even. This one fixed query therefore separates every all-nonsingleton assignment for Q2, irrespective of all child sizes and orientations.

`linear_certificates.json` retains each exact XOR identity, with its indices into the complete selected-query bank. The verification checks all2^18 profile assignments for each row. This finite algebraic certificate, together with the complete profile reduction, gives unbounded child-size scope. The surrounding206-query systems and their ranks are retained but are not required to extrapolate beyond the listed batteries.

## 4. Source-relative original-module corollary and limits

For a fixed source core and ambient, use the appropriate singleton-pattern battery or the listed nontrivial-profile battery and add its core converse-complement queries J(a,c)=(511-c,511-a), with exterior bits unchanged. A baseline-negative member exists by complete augmented-profile duality. If an XOR family of modules of the original source is supported inside the core and selects the root oddly, all proper selected original modules lie inside prime children. Hereditary profile monotonicity and nonnegative run coefficients prevent their profile increases from canceling the negative gap.

This is a restricted source-relative root-family consequence. It is not all-pairs named-response injectivity, does not authorize exterior operations or operations on newly created modules, and is not a complete order-nine classification. The quotient list is exactly Q0 throughQ3 above.

## 5. Reproducibility and physical retention

The self-contained packet retains all four exact9-by9 joint endpoint matrices, all684 anchor primitive coordinates, every512-pattern anchor bank, and both orientations' exact55 run channels plus one rooted cycle count at each of206 actual queries. No quadratic response truncation or endpoint-product approximation is used.

Every one of the55 channels on both sides of all seventeen certifying safe queries is independently reconstructed from actual transitive-child physical graphs by exact interpolation:1,870 graphs of order at most thirteen. The two singleton residuals and their cyclic-root ambient lifts add eight literal sides of order at most twenty. Original child sets, augmented rows, port masks and exact Hamilton counts are retained.

The complete named Boolean response banks for all four bare quotients and their converses retain every coefficient and every basis value, including all cubic and quartic coefficients and all zeros. Their polynomials reconstruct the full named Boolean cubes; this is separate from the selected-query affine certificate for large substitutions.

Run `python3 verify.py` or `python3 -O verify.py`. Only Python standard-library tools and the installed C++17 compiler are required. Checks use explicit exceptions and survive optimized Python. All literal Hamilton counts fit unsigned64-bit integers because the maximum augmented order is twenty. Independent review is required before promotion.
