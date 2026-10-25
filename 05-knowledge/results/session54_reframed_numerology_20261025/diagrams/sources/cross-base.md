# Cross-base coefficient support: old clock, new exact intersection sieve

Status: PROVED ELEMENTARY CONSEQUENCES of Session53's inherited coefficient laws, plus independently reproducible finite tests. The general primitive-root and logarithm mechanisms were already proved in the repository. This is a coefficient-specific extraction and extension, NOT a literature originality claim or a Collatz proof.

## Result and exact objects

Let c(n) be the normalized integral square-root coefficients of Cooper's s7 series, as defined in the frozen sources. Session53 establishes: c(n) nonzero modulo3 iff all base9 digits of n lie in {0,1,2}; c(n) nonzero modulo9 iff all digits lie in {0,...,5} and at most one is in {3,4,5}; for positive n, c(n)=2 modulo4 iff n is a power of2, and otherwise c(n)=0 modulo4. Thus the simultaneous exceptional index problem reduces EXACTLY to the ordinary nonnegative exponent sets

E_b = { k>=0 : the complete base9 expansion of 2^k has no digit6,7,8 and at most b digits3,4,5 }, b=0,1.

Define M_r=6*9^(r-1), and S_(r,b) as the residue classes k modulo M_r satisfying the SAME condition on only the lowest r digits of 2^k. For every r>=1 and integer b>=0,

|S_(r,b)| = 2*3^(r-1) * sum_{j=0}^{min(b,r)} binom(r,j).

In particular |S_(r,0)|=2*3^(r-1), |S_(r,1)|=2(r+1)3^(r-1). Their exact periodic densities are respectively 3^(-r) and (r+1)3^(-r). These are NOT densities of full ordinary integer solutions.

## Proof of finite count and recursive extension

The prior clock theorem gives ord_(3^m)(2)=2*3^(m-1), hence exponentiation maps Z/M_r bijectively to (Z/9^r)^*. At the least digit, each of the normal class {0,1,2} and defective class {3,4,5} has two unit choices (1,2 and4,5). At every other digit either class has three choices. For each chosen set of j defective positions there are 2*3^(r-1) units. Summing proves the formula.

Every survivor exponent class has nine lifts from r to r+1. Since exponentiation bijects unit groups and respects reduction, the new base9 digit runs through0,...,8 exactly once. A zero-defect node has three zero-defect children and three additional one-defect children; a node already using its sole defect has only three admissible children. This yields the counts without enumerating all residues. More generally, fixed budget b gives the polynomial-in-r factor above; no coefficient congruence modulo3^(b+1) is claimed for b>=2.

## Inverse-limit geometry, with parity explicitly retained

The exponent completion is Z/2Z times Z_3, not Z_3 alone. Parameterize one parity sheet by k=epsilon+2z, epsilon=0 or1, z in Z_3. The continuous exponential is

Phi_epsilon(z)=2^epsilon*4^z.

It bijects Z_3 with the unit ball 2^epsilon(1+3Z_3). Its inverse is log(x/2^epsilon)/log4; for distinct z,w,

v_3(Phi_epsilon(z)-Phi_epsilon(w))=1+v_3(z-w).

For integer z,w this is LTE; continuity extends it to Z_3. Equivalently the convergent 3-adic log and exponential on1+3Z_3 give the same identity. Hence Phi is a similarity of ratio1/3 and preserves Hausdorff dimension. Put the ordinary 3-adic metric on each sheet and distance1 between sheets.

The digit Cantor set C={sum d_j 9^j : d_j in{0,1,2}} has dimension log3/log9=1/2. Upper bound: at level r it is covered by3^r balls of radius9^(-r). Lower bound: give each digit independent weight1/3. A ball of radius3^(-m) fixes floor(m/2) full digit blocks, so its measure is at most3^(-floor(m/2)), bounded by sqrt3 times its radius^(1/2). The mass-distribution argument gives dimension at least1/2. Restricting the least digit to1 or2 makes two similarity copies of C and leaves dimension unchanged.

Thus the inverse-limit set for b=0 has Haar measure zero and Hausdorff dimension exactly1/2. For b=1, fix the location and value of the optional defect: each such set is an affine translate of a digit Cantor set with one digit fixed. It has dimension1/2 whenever nonempty. Countably many locations and the no-defect set give dimension1/2 by countable stability. The same proof applies to every fixed finite b. Alternatively their cylinder counts have only a polynomial factor. These are geometries of the completed exponent parameter, NOT dimensions or finiteness results about its countable subset of ordinary nonnegative exponents.

## A genuinely ordinary-exponent consequence, without a limit swap

Let I be ANY interval of K consecutive nonnegative integers, K>=1. Choose r minimal with M_r>=K. No exponent residue occurs twice in I, so

|E_0 intersect I| <= 2*3^(r-1) <= sqrt(6K),
|E_1 intersect I| <= 2(r+1)*3^(r-1) <= (r+1)*sqrt(6K).

For r>=2, minimality gives M_r<9K; for r=1 the bound2<=sqrt(6K) is immediate. This proves uniform O(sqrt K) and O(sqrt K log(2K)) bounds and zero upper Banach density. It does not merely exchange a limit with a sieve density. Consequently simultaneous positive mod4/mod3 exceptions below index X number O(sqrt(log X)); for mod4/mod9 the bound is O(sqrt(log X)*log log X), with the usual harmless log(2+.) convention near X=1.

For general fixed b, the analogous elementary bound is O(sqrt K*(log(2K))^b). The lower bound needed for infinitude, or finiteness, is absent.

## What survives the new lens, and what does not

The mod18/mod19/7/63 note's shared clock is structural: ord_27(2)=ord_19(2)=18, and63=2^6-1. Our digit sieve samples the EVEN levels m=2r of exactly this already proved tower. No mod19 phase or Sierpinski-tournament incidence is used in the coefficient proof. In particular the separate63=56+7 orbit-census coincidence remains unexplained and is not promoted to a theorem by sharing integers.

The finite count is NOT special to this digit alphabet. For any D subset{0,...,8}, write d=|D| and u=#{a in D:3 does not divide a}. For a primitive root g modulo9^r, the analogous count is u*d^(r-1). If u>0 the associated nonempty p-adic digit set has dimension log d/log9, transported by the same logarithmic similarity. Thus merely observing a Cantor count cannot explain modular coefficients. The genuine coefficient content lies in why this particular D and this particular one-defect budget occur, supplied by Session53's Frobenius/Lucas identity.

Finite survival never proves global survival: k=6 passes r=1, since64=1 modulo9, but64=(71)_9 fails globally. For EVERY depth r, k=M_r is a positive finite-window survivor because2^M_r=1 modulo9^r. These moving witnesses do not produce a single ordinary exponent surviving all depths. Compactness supplies p-adic exponent limits, not necessarily ordinary integers. Positive Hausdorff dimension and measure zero are compatible; neither proves that E_0 is finite or infinite. E_0 contains0,1,13: 2^13=8192=(12212)_9. The cheap full-digit scan therefore REFUTES the tempting conjecture that0,1 are its only members. Through k=2000, E_0={0,1,13} and E_1={0,1,2,7,13}. The stronger candidate that these finite lists are complete remains unproved and is not promoted by the sieve. A proof would require a new archimedean/full-digit ingredient.

## Prior recovery and source lineage

- sources/prior-transversality-published.md is git-show bytes from published clone HEAD bbe2ae278aeabdf707fdfca3c1bf0d3208564b03, path05-knowledge/results/collatz_procgen_20260922_transversality_foundry.md. Section3 already proves ternary-avoidance counts2^(j-1) and the logarithmic kappa transport. SHA256 c341592c374d56afd25f6fa6aa13f5246a7a4f098927ad60594f68b4ca037f16. Identical to original math/ tracked HEAD1b91742ce55a202b73c80d56561ba56d00b503a6 and its worktree bytes; no modification was made there.
- sources/clock-fractal.md copies path05-knowledge/results/mod18_mod19_seven_sixtythree_fractal_20261006.md. Its SHA256 is f09deca8e7b6b18eb49540a4a27cd3b93d591cccadad2032fc8cad55f0f9993b, identical to git-show bytes from published HEAD bbe2ae278aeabdf707fdfca3c1bf0d3208564b03. Its section1 supplies the clock; sections3 and9 already warn against density/limit overclaims. Verify source manifest for exact bytes.
- sources/mod3-digit-law.md, mod9-digit-law.md, mod4-dyadic-law.md are Session53 frozen REPORT.md files from cooper-prime-dichotomy, cooper-mod9-extension, dyadic-mod4-supplement respectively. They explicitly label inherited proof dependencies. Their modular input ultimately points to commit7818149bf0a1eca778549fdfd1dde34f479421fd, cooper_s7_weight_one_root_frobenius_20261018.md. This packet does not independently reprove that modular input.
- Bounded prior search: results and theorem directories searched for3-adic log, base9 with power/Cantor, cross-base and digit/Cantor. The exact prior found above blocks any novelty claim for the general method. Absence of an exact coefficient-specific match is not a literature search or proof of originality.

## Predictions, tests, stopping boundary

Predictions recorded before holdout computation: every r-level zero-defect node has exactly3 children; one-defect nodes have6 or3 according as the defect remains unused or used; holdout counts at r10 are39366 and433026; all fixed finite budgets have dimension1/2, not an increasing dimension. These are proved structural extensions, tested cheaply by lifting only surviving nodes rather than giant exponent scans.

check.py independently enumerates full unit orbits through r5, checks budgets0,1,2, lifts through r10, tests full ordinary exponents only through2000, and independently computes actual c(n) through4096 by exact integral recurrence. It includes nonprimitive-base, nonunit, parity-sheet, finite/global, and budget-confusion controls. Finite evidence is not promoted to an infinite ordinary-exponent classification. The cheap finite scan is a guard and hypothesis generator, not the main deliverable. No further exponent enumeration is warranted without a sharper full-digit argument.
