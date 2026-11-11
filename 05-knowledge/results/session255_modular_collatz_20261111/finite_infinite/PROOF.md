# A parity bus, prime-power registers, and a closed ternary cycle spectrum

**Session255, 2026-11-11. PROVED IN NOTES:** the scoped elementary theorems below. **FINITE-EXACT:** the accompanying new-claim controls. No canon promotion or literature-priority claim. Positive-integer Collatz, convergence, and ordinary cycle classification remain **OPEN**.

## 0. Inheritance and the construction target

The useful infinite object is not an enormous unlabelled graph. It is one binary input stream shared by finitely many arithmetic registers, together with an exact affine receipt when ordinary integer realizability matters.

Closest inherited mechanism: the depth-dependent ternary digit loss in [the mod6 reverse-tree pieces](../../collatz_mod6_20260917_reverse_tree_pieces.md), Theorem1.2; its three-row observer does not determine arbitrarily deep descendants. Hostile: a finite residue observer forgets an exact ordinary equality, as in [basin and clock observers](../../collatz_basin_observer_20261009.md), section4. Corrected near miss: [measurement independence](../../collatz_measurement_independence_20261005.md) distinguishes atom data, coefficient-cylinder exit, actual descent, and a source-dependent deadline. Least-used sidecar here: the *ordered additive numerator*, not just slope or odd-step count.

Concept board: mod8 multiplication; signed shortcut dynamics; shared-input CRT; ternary prime-power filtration; rational-series compression; integer realization. Anchor = signed modular dynamics. Niche = the exact ternary determinant below. Wildcard = an information-correct product operation transferable to boundary-count problems.

We use both globally signed shortcut maps

    T_sigma(n) = n/2 if n is even; (3n+sigma)/2 if n is odd,
    sigma in {−1,+1}.

Reflection is exact: T_minus(−n)=−T_plus(n). Thus changing the sign of the constant and reflecting the source are linked operations; neither proves a positive-domain convergence assertion. We do not substitute the accelerated odd-only clock.

## 1. Theorem F1: finite branch receipt and exact precision cost

Let w=(b_0,...,b_{k−1}) be any binary word, with s=sum b_i. Define B_0=0 and

    B_{i+1}=3^{b_i}B_i + sigma b_i 2^i.

Then there is exactly one residue c_w modulo2^k whose first k shortcut parities are w. For every integer n in that residue class,

    2^k T_sigma^k(n)=3^s n+B_k,
    c_w=−3^(−s) B_k modulo2^k.

**Proof.** The inverse even branch is y↦2y; the inverse odd branch is y↦(2y−sigma)/3 on the 2-adic integers. Each inverse branch maps onto exactly its designated parity class and shrinks 2-adic distances by1/2. Composing k branches gives one residue class modulo2^k, onto which the unrestricted terminal 2-adic value maps bijectively. Affine composition gives the displayed recursion, which also proves the residue formula. Integer trajectories in that class obey the same first k parities, since each step consumes exactly one available binary digit. No probabilistic assumption is involved. QED.

**Mod8 anchor.** Multiplication by3 has table [0,3,6,1,4,7,2,5], with two odd transpositions 1↔3 and5↔7. But T_plus(1)=2 and T_plus(9)=14 have distinct residues modulo8 although 1≡9 modulo8. The multiplication table determines the numerator, not a quotient after division by2. Computing one output modulo8 requires the input modulo16; k outputs at fixed final precision cost k extra input bits. This is the forward twin of the inherited inverse tree consuming one ternary digit per inverse division by3.

The full infinite parity stream Q_sigma(x) determines one x in Z_2: its nested cylinders have diameters2^(−k). Conversely every x gives such a stream. This proves directly a homeomorphism between Z_2 and the binary one-sided shift space, with Q_sigma T_sigma = shift Q_sigma. This classical elementary mechanism is rederived here, not claimed as a new conjugacy discovery.

## 2. Theorem F2: the finite register with a lossless infinite input tape

For odd q≥3 define

    X_q = inverse limit_k Z/(q·2^k) = Z/q × Z_2.

Here division by2 consumes a binary digit, while it is multiplication by2^(−1) in Z/q. The map induced by signed shortcut dynamics is

    (r,x) ↦ (f_b(r),T_sigma(x)),  b=x modulo2,
    f_0(r)=r/2 moduloq,  f_1(r)=(3r+sigma)/2 moduloq.

Under (r,x)↦(r,Q_sigma(x)) this is precisely a q-state deterministic transducer reading and removing the first bit of an infinite tape. Finite words act by

    r ↦ 2^(−k)(3^s r+B_k) moduloq.

**Proof.** CRT gives the product decomposition compatibly at every level. The parity decision lies in the binary coordinate; the odd coordinate has a unique division by2. F1 gives the conjugacy and composition rule. QED.

**Exact finite interpretation.** Fix r and a length-k input word. Exactly one class modulo q2^k has those two observations, by F1 and CRT. Its natural density among ordinary integers is1/(q2^k). Thus the full-shift tape is not a claim that the parity digits of a specified deterministic integer orbit are independent. It parameterizes all initial cylinders.

**Loss ledger.** Projection to the q-state graph forgets the binary tape. Keeping the tape restores the whole profinite state, not ordinary integrality or positivity. The diagonal embedding n↦(n modq,n inZ_2) is dense, but its positive-integer image is only countable. Density does not give a positive integer with a prescribed infinite itinerary.

## 3. Theorem F3: synchronized CRT is the correct information product

Let q=product q_j with pairwise coprime odd q_j. Write P_b(q) for the row-transition matrix of f_b: entry(r,t)=1 exactly when t=f_b(r). Under the CRT permutation of states,

    P_b(q) = tensor_j P_b(q_j),
    A_q(u) := P_0(q)+uP_1(q)
            = tensor_j P_0(q_j) + u tensor_j P_1(q_j).

This is generally NOT tensor_j(P_0(q_j)+uP_1(q_j)). The latter grants different components independent input bits and introduces spurious trajectories.

**Proof.** Polynomial evaluation and division by the unit2 commute with CRT for each fixed branch. A global branch is one choice shared by all coordinates. Summing only after that choice proves the identity. Independent sums contain mixed-branch tensor terms, which do not belong to the original alphabet. QED.

**Small hostile.** For q=15, sigma=+1, source0, the actual outputs are0 and8. Combining independent choices modulo3 and5 produces0,3,5,8. The outputs3 and5 are fake. Even the one-step information product fails if the common parity bit is discarded.

The minimal modular factors are prime powers, not primes alone: CRT separates distinct primes, while the 3,9,27,... tower retains nilpotent precision. Modulo3 the odd branch is constant; modulo9 its image has three elements. Replacing9 by3 destroys a real register coordinate. For primes p≠3 the odd branch is bijective, while at3 it is not; this is a structural prime/composite and unit/nonunit distinction, not cosmetic table relabelling.

## 4. Theorem F4: rational series for regular modular questions

Let D be a finite deterministic automaton restricting the permitted parity language; retain its state together with r. Attach u to odd edges,1 to even edges, and optional commuting observable weights to any edge. Let M be the resulting finite weighted matrix. For initial row alpha and accepting column beta, the formal generating series is

    sum_{k≥0} z^k alpha M^k beta = alpha(I−zM)^(−1)beta.

It is rational over the field of rational functions in the weights. Coefficients count actual modular paths with multiplicity; two letters landing at the same vertex are still two edges.

**Proof.** Expand the geometric matrix series. Matrix multiplication concatenates paths and multiplies their weights, and summation retains distinct paths. The adjugate formula gives rationality. QED.

This makes an infinite family of *finite observer questions* manageable. It does not make the language of positive integers reaching1 regular. If branch order is needed for arithmetic realization, retain w or the affine receipt (k,s,B); recording only u^s discards B. For example words10 and01 have equal length and odd count but B=1 and2.

## 5. Theorem F5: closed weighted periodic spectrum at every ternary level

For a≥1 put q=3^a and A_a(u)=P_0(q)+uP_1(q). For every k≥1,

    tr(A_a(u)^k) = (1+u)^k − 1 + gcd(2^k−1,3^a).                 (1)

Consequently, as a polynomial identity in Z[z,u],

    det(I−z A_a(u))
      = (1−(1+u)z) product_{j=1}^a (1−z^{2·3^{j−1}}).          (2)

Both identities hold for either sign sigma. Equivalently the weighted dynamical zeta of the compact system X_(3^a), with each odd letter weighted by u, is the reciprocal of (2).

**Proof of (1).** A length-k word w with s odd letters fixes r precisely when

    (2^k−3^s)r = B_k modulo3^a.

If s≥1, its coefficient is a unit modulo3, hence it has exactly one fixed residue at every ternary level, regardless of B or sigma. Summing these words yields (1+u)^k−1. The sole all-even word has B=0 and fixes gcd(2^k−1,3^a) residues. QED.

**Proof of (2).** P_0 is multiplication by2^(−1) on Z/3^a. Its permutation consists of the fixed point0 and, for each j=1,...,a, one cycle of length2·3^(j−1) on the elements of exact additive order3^j. Indeed2 has order2·3^(j−1) modulo3^j: the base order is2 and v_3(4^(3^h)−1)=h+1 by binomial induction, with no smaller exponent because the unit group has that size. Thus

    det(I−zP_0)=(1−z) product_j(1−z^{2·3^{j−1}}).

The formal identity −log det(I−zA)=sum_{k≥1}tr(A^k)z^k/k, obtained by differentiating the determinant identity and comparing the zero constant term, applies over Q[u][[z]]. Equation(1) replaces the fixed eigenvalue1 contribution of P_0 by1+u and leaves the others. Cancelling1−z gives (2); both sides are polynomials, so the formal equality is exact. QED.

**Compact periodic-point interpretation.** Every periodically repeated binary word has exactly one compatible binary coordinate by F1. Its q-state closing residues are precisely those counted by the trace. This gives an actual bijection between compact periodic points of period dividing k and closed length-k labelled paths. No replacement of ordinary periodic integers by these compact points is implicit.

**Corollaries.**

- The level-(a+1) characteristic determinant divided by the level-a determinant is exactly1−z^{2·3^a}. The new precision contributes one cyclotomic block independent of the odd weight and sign.
- For u=1, the adjacency spectrum is2 plus the roots of unity supplied by these blocks, with their algebraic multiplicities. This describes eigenvalues, not a claimed simultaneous diagonalization or Jordan form.
- A word with s≥a odd symbols maps all ternary starting states to one state, because its affine slope3^s2^(−k) vanishes modulo3^a. Thus each odd occurrence pays down one ternary memory digit. This reset property does not imply binary or ordinary-source memory is erased.
- For every nonempty word containing an odd symbol there is a unique compatible 3-adic fixed residue across all a, since the coefficient is a3-adic unit. This is infinite-depth compatibility in the ternary tower, not an integer-cycle theorem.

For example a=1 gives det=(1−(1+u)z)(1−z²); a=2 adds1−z⁶; a=3 adds1−z¹⁸. A family of exponentially larger state graphs has a linear-in-a factored description.

## 5b. Theorem F5b: every ternary enlargement adds only even-branch spectral blocks

Write any odd q as q=3^a m with gcd(m,6)=1, allowing m=1. For m=1 use the single-state matrix A_1=1+u. Let E_n=P_0(n). Then

    det(I−zA_q(u)) / det(I−zA_m(u))
      = det(I−zE_q) / det(I−zE_m)
      = product_{j=1}^a product_{d|m}
          (1−z^{ord_{3^j d}(2)})^{phi(3^j d)/ord_{3^j d}(2)}.     (3)

Thus the ratio is an explicit polynomial independent of the odd weight u and of sigma. The degree of this new factor is q−m. In particular all eigenvalues newly contributed by ternary precision are roots of unity, with their stated multiplicities. This does not assert that the remaining m-register spectrum is simple or that the full operator is diagonalizable.

**Proof.** Fix a word w. If it contains an odd letter, its fixed-point equation has one solution on the3^a register, so synchronized CRT identifies its q-closing states bijectively with its m-closing states. This is true even when there are no m-solutions. For the sole all-even word, the two counts are respectively

    gcd(2^k−1,m) gcd(2^k−1,3^a), and gcd(2^k−1,m).

Consequently the difference tr(A_q^k)−tr(A_m^k) equals tr(E_q^k)−tr(E_m^k) for every k, and is independent of u,sigma. Taking the formal determinant logarithm proves the first equality. Multiplication by2^(−1) on Z/n preserves additive order e for every e|n and splits the phi(e) elements of that order into cycles of length ord_e(2). Hence

    det(I−zE_n)=product_{e|n}(1−z^{ord_e(2)})^{phi(e)/ord_e(2)},

with ord_1(2)=1. Cancelling the divisors e|m leaves exactly (3). Its degree is sum_{j=1}^a sum_{d|m} phi(3^j d)=q−m. QED.

For q=15,m=5, the new factor is (1−z²)(1−z⁴)². Increasing45 over15 adds (1−z⁶)(1−z¹²)². The unrestricted parity bus forces this entire nontrivial spectrum without enumerating either graph.

**A second exact viewpoint: conductor filtration.** On complex functions on Z/q, pullback by f_b acts on the Fourier character chi_l(r)=exp(2pi i l r/q) as

    P_0 chi_l = chi_(l/2),
    P_1 chi_l = exp(2pi i l sigma inv_q(2)/q) chi_(3l/2),

where division by2 in the frequency means its modular inverse, and inv_q(2) is that inverse in the phase. Filter characters by ternary part of their additive conductor. P_0 preserves every grade, while P_1 strictly lowers every positive ternary grade. The grade-zero subspace is exactly functions pulled back from Z/m and carries A_m. Thus the operator is block triangular with odd-branch contribution only below the positive-grade diagonal. This supplies the structural reason for spectral stability: odd steps erase ternary memory, even steps rotate what remains. It is not a claim that all trajectories synchronize on the m-register.

## 5c. Theorem F5c: arbitrary finite drivers preserve the spectral extension

Let D_0,D_1 be any d-by-d matrices over a commutative coefficient ring. They may encode a deterministic or weighted finite automaton, with rejecting transitions omitted. Set

    M_n = P_0(n) tensor D_0 + u P_1(n) tensor D_1.

For q=3^a m as above, there is the exact polynomial identity

    det(I−zM_q) = det(I−zM_m)
      product_{j=1}^a product_{e|m}
        det(I_d−z^ell D_0^ell)^{phi(3^j e)/ell},
    ell = ord_{3^j e}(2).                                    (4)

The entire additional factor depends only on D_0: all odd-driver behavior D_1 remains in the base-m determinant and off-diagonal blocks. Setting d=1,D_0=D_1=1 recovers F5b.

**Proof.** Tensor the Fourier conductor filtration of F5b with the driver space. On the grade-zero quotient the operator is M_m. On every positive ternary grade, the odd pullback strictly lowers conductor and contributes nothing to the diagonal block, regardless of D_1. The diagonal block is E_grade tensor D_0. Characters of additive conductor3^j e form phi(3^j e)/ell permutation cycles of length ell under the even branch. For a cyclic permutation matrix C_ell,

    det(I−z C_ell tensor D_0)=det(I_d−z^ell D_0^ell).

For example this follows by the formal trace identity: the permutation trace vanishes unless ell divides the power, and equals ell when it does. Taking the block-triangular determinant gives(4). The Fourier proof can first be performed over the complex field with independent symbolic entries for D_0,D_1. Both sides are polynomials with integer coefficients in those entries,z,u; their equality is therefore a universal polynomial identity, and specialization proves the assertion over every commutative ring. QED.

**Regular-language control.** For parity words containing no consecutive ones, use two states recording whether the preceding bit was1, with

    D_0 = [[1,0],[1,0]],    D_1 = [[0,1],[0,0]].

The missing odd transition from the second state rejects11. Since D_0 is idempotent of rank1,

    det(I−zM_1)=1−z−uz²,
    det(I−zM_3)=(1−z−uz²)(1−z²).

This holds for either sign and is checked symbolically. More generally the same driver at q=3^a has determinant (1−z−uz²) product_j(1−z^{2·3^{j−1}}).

**What this closes and what it does not.** The extra spectrum under ternary refinement is now known for *every* finite regular constraint, not only unrestricted words. Acceptance generating functions still have the form alpha(I−zM)^(-1)beta. Their numerators, cancellations, and projections depend on the actual odd-driver off-diagonal blocks and initial/final states. A characteristic determinant is not an equality of accepted languages, counts, or ordinary integer realizations. Trace counts close the driver state cyclically; ordinary accepted finite words instead use the specified initial and accepting states. The supplied control computes residue-specific accepted-path series as well as determinants to retain this distinction.

## 6. Theorem F6: exact arithmetic realization sidecar

For a nonempty length-k word, D=2^k−3^s is nonzero and odd. Its unique 2-adic periodic point is the rational number

    x_w = B_k / D.

It has precisely the prescribed parity word, because D is odd and the inverse-cylinder contraction has that unique fixed point. It is an ordinary integer iff D divides B_k, and positive iff additionally B_k/D>0. When these hold, it is a genuine signed shortcut periodic integer; the least period may divide k.

For a general odd modulus q the number of modular fixed states of this word is

    gcd(D,q), if gcd(D,q) divides B_k; otherwise0.

These are necessary finite congruence tests on the integral realization. **All** odd moduli together detect divisibility: if D does not divide B, choose q=|D| (when |D|>1). A fixed finite modulus bank has no such universal completeness guarantee. Even all ternary powers miss every denominator prime other than3, and D is a3-adic unit whenever s>0.

**Decisive hostile.** For sigma=+1 and w=100, B=1,D=5, so x=1/5. It has parity orbit1/5→4/5→2/5→1/5 in Z_2. Every ternary quotient accepts a unique closing path. Modulo5 there is no closing path. The ternary determinant therefore counts this genuine compact periodic orbit, while correctly saying nothing about a positive integer cycle. Positivity is another independent sidecar: w=1 gives x=−1 for sigma=+1 and x=1 for sigma=−1.

## 7. Precise connections to Session254 and forward experiments

Source: [Session254 synthesis](../../session254_four_frontiers_20261111/SYNTHESIS.md), especially its prime-power normal forms, all-row quotients, and geometric conditional multiplication.

**Actual shared-state product map.** For fixed branch word w, the set of modular paths modulo product q_j is bijective under CRT to the Cartesian product of the factor-path sets *with that same word*. Thus the per-word counting functions multiply pointwise. This is the same elementary fibre-product mechanism as Cartesian products of HN path colourings at fixed port colours: independent interiors are combined only after their shared boundary state is fixed. The map here is explicit CRT, the preserved predicate is branch-labelled path/closure, the lost datum on summing first is branch synchronization, and the restoration sidecar is w. No map between Collatz numbers and geometric colourings is claimed. The modulo15 hostile is the cheapest test of the wrong product order.

**Prime-power connection, with limits.** Chan–Pak's mod9/27/81 context response and this ternary transducer both retain finite quotient responses at each precision rather than identifying the prime with its powers. Here reduction3^(a+1)→3^a intertwines every letter, and the exact additional spectral block is computed. This gives a concrete worked model of a compatible tower with a cheap all-depth rule. It is a mechanism transfer, not a claim that the Chan–Pak curvature quotient is a Collatz state or that its positive surgeries follow from our determinant.

**The useful next bets.**

1. Theorem F5c now settles the added ternary spectrum for every regular-language constraint. The next useful target is its accepted-path numerator: determine when initial/final selectors cancel cyclotomic factors, and which odd-driver entries remain observable despite identical added spectra.
2. For q with primes other than3, retain the pair (D modq,B modq), not just s. Search for word families where exact B-divisibility becomes a finite certificate. This attacks an actual lost coordinate.
3. Compile a given finite integer-source family into cylinder receipts and add a proved descent threshold. The present transducer provides exact branch accounting; a uniform positive-source descent theorem would be additional mathematics, not a consequence of spectral radius2.

## Reproduction and scope

Run from the repository root:

    python3 05-knowledge/results/session255_modular_collatz_20261111/finite_infinite/reproduce.py

The script uses exact Python integers/Fractions and SymPy polynomial determinants. It writes `EXACT_RESULTS.json`. It checks all signed binary cylinders through depth9; all closed labelled paths for q=3,9,27,81, lengths1..9 and both signs; symbolic bivariate determinants at q=3,9 and both signs; mixed-prime trace polynomials at q=15,45,21, lengths1..7 and both signs; q=15 versus5 determinant identities at u=0,1,2 and both signs; the two-state no-consecutive-ones driver determinant and residue-specific accepted-path series; the mod8 precision hostile; the synchronized-CRT hostile; and the rational1/5 cycle. These are new-claim controls, not a rerun of inherited large Collatz censuses. The all-a and all-k conclusions are supplied by the proofs, not inferred from those finite universes.
