# Independent audit: cyclic three-block completion cost

**PASS, no correction.** The bivariate ordinary-class generating function, exact cost coefficients, linear-scale limiting law and mean1/9 are valid for the specified C3 transitive-blow-up family. They do not preserve or prescribe Hamilton-path counts.

## Frozen binding

- C3_COMPLETION_LAW.md: `4b20e89a7944302a25c921ea687975cacd169de16e95f565b376c513c24ef773`.
- verify_c3_cost.cpp: `126ed0d7f66519a781ee458ac1e4967abb36277a239c5eca9dc07c5d86d04d56`.
- Both supplied receipts and both copied source replays: `bc4ea08f7df50bc542cba9d51c3754470495f761e192562af5eb7095ed48f8bd`.

The source remains unchanged in session197-orbit-wildcard/multiblock-orbits/. Its main completion theorem is now independently PASS, report session197-multiblock-orbits-audit/REPORT.md SHA64a181aae0ae4f660fdc325547535884b8bf7e912cd1cc37ad3b4eaf285cf37d, bound to main proof934f76db…. This audit is distinct from that main review.

## Exact class and cost bijection

The quotient automorphism group consists of the three cyclic rotations. A sorted triple of three distinct weights therefore gives two ordinary classes, interchanged by converse. A repeated-weight triple gives one. All-equal triples also give exactly one, so their nontrivial stabilizer does not introduce a missing factor.

The main theorem gives D equal to the minimum pairwise absolute difference. In sorted order x<=y<=z this is min(y−x,z−y). The D0 class is uniquely specified by its repeated value a and third value b; even a=b has one representation. Hence its series is z³/[(1−z)(1−z²)].

For d>=1, imposing D>=d is exactly the substitution y=x+d+u, z=x+2d+u+v with x>=1 and u,v>=0. It contributes two classes and total order3x+3d+2u+v. Subtraction of the d+1 tail cancels the (1−z³) denominator and proves the exact-D coefficients. Summing in d yields the claimed rational bivariate series. The coefficient floors, ordinary Burnside specialization, maximum floor((n−3)/3) and its attainment are all correct, including n3,4,5.

## Limit distribution and moment

For fixed interior x in(0,1/3), d=ceil(nx) leaves residual order n−3−3d=n(1−3x)+O(1). The positive coin series with weights1,2,3 has coefficient N²/12+O(N): its pole at1 has order3 and all other root-of-unity poles have smaller order. The tail numerator is twice this coefficient, while the ordinary class denominator is n²/6+O(n). Their ratio is exactly asymptotic to (1−3x)².

At zero the survival probability is identically1, and at1/3 it is identically0 because the maximum cost is strictly below n/3. Thus the stated continuous limiting law holds on the full support, with density6(1−3x). The zero-cost probability is asymptotic to3/n and supplies no limiting atom. Boundedness D/n<=1/3 permits dominated convergence of the tail integrals, giving mean1/9. Differentiating the bivariate series gives the exact first-moment series with numerator2z6 and squared (1−z³) denominator, as stated.

The sampling measure is uniform over ordinary unmarked classes in this fixed quotient family. No labelled, Hamilton-weighted or all-tournament distribution is substituted.

## Independent controls

The standalone Python checker imports no author code. It enumerates sorted integer partitions into three positive parts and applies the literal orbit multiplicity1 or2, rather than enumerating cyclic weight keys. A separate bivariate formal-series multiplication checks the rational function. Univariate coin recurrences independently check every tail coefficient and the first-moment series.

For physical minimum controls it reconstructs the actual graph for each candidate enlarged triple. Self-converseness is tested by enumerating all degree-compatible full permutations via products of permutation classes and then checking every arc, not by assuming that two equal weights are sufficient. It searches increasing added order and compares the first actual successful graph with the formula.

Normal and optimized independent runs agree on:

- 437,844 ordinary classes through order200;
- 6,831 exact bivariate coefficients and6,633 tail coefficients;
- all exact first moments through order200;
- 35 physical starting triples and62 candidate searches, involving42 distinct actual graphs;
- exact mean/order ratios at orders200,1000,10000, recorded only as corroboration of the proved limit.

The copied source additionally reproduces its183,800 classes through150 and76 physical optima with258 actual anti-isomorphism tests in optimized and assertion-enabled debug modes.

## Reproduction and hashes

Run `python3 check.py` and its `-O` form. It has no data dependencies and writes only stdout. Compile the separate copied source checker with C++17 if desired. The all-order bijections and asymptotic proof, not the finite range, establish the result.

- check.py: `085125813a0858155d98045937524a270f0660b8e9ee8bdb67cf446b99fd67bc`.
- independent-normal.json: `3317f70aae176a8a26e867d19df666e58a6e06d940c16b4f8bb50a0391e44d10`.
- independent-optimized.json: `3317f70aae176a8a26e867d19df666e58a6e06d940c16b4f8bb50a0391e44d10`.
