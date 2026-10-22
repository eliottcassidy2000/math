# Powers of three: exact Collatz cylinders, full finite languages, and a source-budget obstruction

Session29, 2026-10-22. **PROVED (elementary, scoped):** the statements below have complete proofs. **FINITE-EXACT:** accompanying independent arithmetic controls. **OPEN:** universal convergence, eventual entry, and uniform153N annulus coverage. No priority claim. This packet was written outside the repository; no Git action was taken.

## 1. Inheritance and concept board

Closest proved mechanism: `05-knowledge/results/metagraph_trijections_reset_session17_20261020.md`, §7, exact guarded affine-debt composition and finite-state coverage obstruction. The more detailed linked directory is absent from this sparse checkout, so no unseen detail is imported. Closest actual-clock mechanism: `05-knowledge/results/timed_blocks_and_area_clocks_20261020j.md`, §§2–4, ordered source-block insertion, distinct pair ownership, and source-gap registers. Canonical hostile: the arbitrarily long all-valuation-one prefix, not an uncomputed divergent orbit. Corrected near miss: `01-canon/MISTAKES.md`, October20 boundary audit, `(1,1,1,1,3)^k` has valuation clock7k, not8k, and multiplying strict inequalities by a zero cone parameter loses strictness. Least-used relevant sidecar: exponent congruence, retaining the difference between counting exponents and counting integer sources.

Recovered further mechanisms: `collatz_maximal_run_reset_20261019.md` retains an exact exit shell for an expanding word; `collatz_weak_reset_family_20261005.md` prices contraction through actual excess and keeps seed-table hypotheses; `collatz_reset_debt_20261015f.md` refutes the unchanged one-stage7X claim at k23 rather than silently extrapolating lower scales. `THM-4629-twisted-transfer-law-step-count-characters-of-collatz-orbits-hear-the-near-returns-by-congruence-class.md` distinguishes shortcut clock tau, odd count omega, and ordinary clock tau+omega; its infinite spectral/closure conclusions are not dependencies here.

A bounded recovery pass searched356 physically present top-level Collatz/reset Markdown result files, with43 matching files recorded in `recovery_search.json`. This is not a search of every tracked path: the repository is sparse. The earlier Session27 repunit report and catalogue were consulted as historical routes, not reused as Collatz theorems. The supplied text attachment inspires the local-data-versus-global-object question; none of its lattice or website claims is assumed.

Live board:

1. Anchor: powers of3 / exponent residue classes / first reset / missing carry and clock / test n27.
2. Niche: oddparts of3^k+1 / valuation-word cylinders / full finite language / missing infinite-orbit compatibility / binary-lift every short word.
3. Debt: affine maps / upper-bound cone / composition slack / missing legal source and ownership / compare reset with actual descent.
4. Clock: odd-step count versus total divisions by2 / exact cost / sum of valuations / test an intentionally wrong extra bit.
5. Wildcard: local-statistic universality / finite prefix observables / indistinguishability / source-height density / count actual family members in an annulus.

## 2. Conventions and the elementary power-plus-one bridge

For positive odd n set U(n)=(3n+1)/2^a, where a=v2(3n+1). In this algebraic section extend U(1)=1; strict ROOT stopping is addressed below. A word w=(a1,...,ah) has odd-step length h, shortcut clock A=sum ai, and ordinary Collatz clock A+h. These three are not interchangeable.

For k>=1,

- v2(3^k+1)=1 if k is even;
- v2(3^k+1)=2 if k is odd.

This follows from residues modulo8. Thus the oddpart of3^k+1 is exactly U(3^(k−1)). The expression3^k+1 is a one-step numerator on a special input, not the numerator of k repeated3n+1 operations. Without halving, k repetitions give3^k n+(3^k−1)/2; with halving they require the full carry polynomial below.

## 3. Exact valuation-word cylinder lemma

For a positive valuation word w define A0=0,C0=0 and

A_i=A_(i−1)+a_i,   C_i=3C_(i−1)+2^A_(i−1).

Then the candidate composition is

U^h(n)=(3^h n+C_h)/2^A_h.

The sources realizing exactly w are one odd residue class

r_w = (2^A−C_h) (3^h)^(-1) mod2^(A+1).

Proof: requiring the final quotient to be odd gives the displayed congruence. It forces the last predivision numerator to have its prescribed exact valuation. Induct backwards: since the last odd output equals (3x+1)/2^a, its predecessor x is odd whenever it is integral; the final congruence reduced at the preceding powers of2 supplies that integrality and parity. Equivalently, induction on h refines the existing source cylinder by exactly a_h bits. This proves necessity and sufficiency, including the terminal odd bit. All residues are positive odd representatives; adding multiples of2^(A+1) gives every source in the cylinder.

The dependence on A+1 is essential. Using only2^A merges two terminal-parity alternatives and does not certify the exact word.

## 4. Power-of-three finite-language theorem

For L>=3 the powers of3 modulo2^L are exactly the residues congruent to1 or3 modulo8, each occurring once as m runs modulo2^(L−2).

Proof: 3^(2j) is1 modulo8, and multiplication by3 gives the other coset. The elementary valuation identity v2(3^(2j)−1)=3+v2(j) for j>0 proves ord_(2^L)(3)=2^(L−2). The two cosets have precisely that cardinality.

Therefore a word with A>=2 occurs for n=3^m iff r_w mod8 is1 or3. This is equivalent to either:

- its first letter is2; or
- it starts(1,b) with b>=2.

The one-letter word(1) is the remaining shorter special case. Every compatible longer word has one exponent residue modulo2^(A−1), hence density2^(1−A) in the exponent m. In particular every finite tail is available after a first letter2.

This also supplies a direct, logarithmic-number-of-lifts discrete-log algorithm: begin with m=0 or1 modulo2 at level3, then at level j+1 choose m or m+2^(j−2). The valuation identity proves exactly one choice works. The implementation records actual integer exponents, not floating logarithms.

### First-reset classification and a failed payment hypothesis

For m>=1, n=3^m:

- m even: the first valuation is2 and U(n)<n;
- m odd: the first valuation is1 and the second is b=v2(3^(m+2)+5)−1>=2.

For odd m, the first-reset endpoint is(9n+5)/2^(b+1). It is below n iff m≡1 modulo4; it is above n iff m≡3 modulo4. Indeed b=2 exactly in the latter case; b>=3 in the former, with the smallest n=3 already satisfying5<7n. Thus three quarters of exponent residue classes descend by their first reset; the other quarter do not. This is exponent density, not density among integer sources or a total stopping-time statement.

For every b>=2 the second valuation equals b on precisely one odd exponent class modulo2^b. Among odd exponents its frequency is2^(1−b), proved by the power-subgroup theorem, not a stochastic assumption. Subsequent finite valuation tails conditional on any fixed first-reset cylinder have the geometric product frequencies2^(−sum tail).

### No uniform first-descent depth, even for powers of3

For every r>=0, the word(1,2,1,...,1), with r final ones, has source cylinder

n≡−13/9 mod2^(r+4).

It belongs to the power subgroup, so infinitely many powers3^m realize it. The first step expands; the first two steps give(9n+5)/8>n; every later prescribed step expands. Consequently all first r+2 odd-step endpoints exceed the original n. The shortcut clock is r+3 and ordinary clock2r+5.

Thus the first-descent time in odd steps is unbounded among powers of3, despite their uniformly bounded first-reset time. Least exponent residues for r=0,...,8 are3,3,3,3,3,67,67,67,579; these are exact progression representatives, not convergence experiments.

## 5. The user's oddpart family has the full finite Collatz language

Let u_k=(3^k+1)/4 for odd k>=1. This strengthening was proposed by the parent during the session and independently checked here.

For any L>=1, as odd k varies modulo2^L, u_k ranges bijectively over all odd residues modulo2^L. Proof: u_k≡r mod2^L iff3^k≡4r−1 mod2^(L+2). For odd r the right side is3 modulo8, precisely the odd-exponent coset; the power-subgroup order is2^L. The two sets both have2^(L−1) elements. The L=1 case is immediate modulo8.

Combining this with §3 proves:

**Full finite-language theorem.** Every positive finite valuation word w is realized by infinitely many u_k. Its odd exponents k form exactly one residue class modulo2^(A+1). Its density conditional on odd exponents is2^(−A)=product_i2^(−a_i).

This is exactly the finite product law for independent geometric valuation letters, as a natural density on the exponent parameter. It is not a claim of random independent behavior along a particular positive integer orbit.

For an all-ones word of length h, the exponent condition is especially simple:

3^k≡−5 mod2^(h+3), k in one odd class modulo2^(h+1).

Every such u_k increases throughout its first h odd steps, with clocks h and2h. Thus no fixed odd-step first-descent bound holds on this family either. More generally any purported restriction depending only on a fixed finite valuation prefix that holds for every u_k must hold for every positive odd source. A bounded-prefix test cannot distinguish this sparse family from the full finite language.

ROOT guard: a fixed word has only finitely many positive sources hitting1 during its prescribed finite prefix, since each prefix is an injective affine function of the original source. Removing these sources removes finitely many k, leaving infinitely many root-free realizations in each cylinder. Hence the finite-language and hostile-prefix conclusions do not rely on walking around the trivial1-cycle.

**Quantifier boundary.** For every word there is an exponent is not there is one exponent realizing every compatible infinite sequence. In the all-ones case the limiting 2-adic source is−1, not a positive integer. Neither the product frequencies nor arbitrarily deep compatible finite cylinders prove convergence or divergence.

## 6. What this does and does not pay in the153N ledger

The inherited affine debt convention is x_out<=(x_in+delta)/lambda, with lambda>0 and delta>=0. For lambda1,lambda2>0 and delta1,delta2>=0, chronological composition is

(lambda1,delta1) then(lambda2,delta2) -> (lambda1 lambda2, delta1+lambda1 delta2).

The slack lambda−1−delta is the first slack plus lambda1 times the second. This explains closure of the cone delta<=lambda−1, but not legal occurrence or source coverage. It is a weighted cocycle; it is not an unweighted invariant like a Dehn invariant, nor a count of distinct original sources.

There is a sharp concrete source-density obstruction to spending the new exponent frequencies as annulus credit. Consecutive odd-index family values satisfy

u_(k+2)/u_k=(9*3^k+1)/(3^k+1)>8 for k>=3.

For every integer N>=1, the annulus[4N,256N) excludes u1=1. It consequently contains at most two members of this oddpart family: three successive eligible members would have ratio greater than64. Meanwhile its native source budget is153N, and its pair-interaction ledger contains all original old/new and new/new pairs. A density theorem in k has no mass-preserving map to those153N original units. It cannot be substituted for a uniform R' estimate, nor used twice after class refinement.

The result is still useful as a hostile-test generator for any proposed bounded-clock certificate: realize the certificate's missing finite word on a structured explicit family, then retain its exact source, guards, and clock. To become a positive153N result, a separate theorem must map distinct annulus sources to disjoint owned blocks, retain mod15/mod35 registers and arrival times, and pay the unselected interactions. None is supplied by symbolic universality alone.

### An exact mod15/mod35 sidecar, without a false independence claim

For odd k, the pair(u_k mod15,u_k mod35) has period12 in k. At k=1,3,5,7,9,11 its values are respectively(1,1),(7,7),(1,26),(7,22),(1,21),(7,12). This follows by computing3^k modulo420 before dividing by4; the modulus3 coordinate is fixed and the remaining unit order divides12. The six corresponding residues modulo105 are1,7,61,22,91,82.

For a fixed valuation word, k is fixed modulo2^(A+1), hence modulo4. Only three of these six register states are compatible: those k mod12 with the same mod4 value. All three occur infinitely often by the ordinary CRT, since the remaining mod3 choice is free. Thus arbitrary finite valuations coexist with exact, but not independent, mod15/mod35 guards. This is an explicit source-register sidecar. It does not identify the guarded native annulus parent, terminal, ownership, or arrival clock, and does not claim to preserve those missing data.

## 7. Hypothesis/test cycles and strongest survivors

Cycle1: '3^m inputs have quick resets, so the first reset should pay.' Test n27:27→41→31, valuations(1,2), endpoint greater than27. REFUTED. Survivor: exact residue classification, first reset by two steps, and payment iff m even or m≡1 mod4. Revision exposes unbounded later growth via(1,2,1^r).

Cycle2: 'Power inputs form a small symbolic language after the first reset.' Enumerate all positive words of clock at most11; half the nontrivial cylinders occur, precisely those compatible with the two initial guards. REFUTED after the initial guard. Survivor: all finite tails are possible, with exact exponent residues. Stronger parent suggestion then leads to the oddpart family realizing every word.

Cycle3: 'Full finite-language geometric frequencies could supply annulus credit.' Test the map's source-height measure. REFUTED as a mass-preserving inference: at most two family sources lie in any native64-ratio annulus. Survivor: a complete finite-word hostile generator, not a source-budget lower bound.

Implementation near miss: an auxiliary check first expected two solutions after fixing a particular extra upper bit of a second-valuation congruence. It correctly failed. Fixing that bit selects one lift; two solutions would occur only after ignoring it. The corrected check and independent lower-modulus exact-valuation check both pass. No mathematical conclusion depended on the failed check.

## 8. Reproducibility, literature, and boundaries

Run `python3 verify.py > verification.json` and `python3 -O verify.py > verification.optimized.json`; outputs must be byte-identical. Checks are explicit exceptions, not removable assert statements. The run verifies2046 valuation words of summed clock2..11,1023 power-source compatible words, all2046 oddpart realizations, complete power subgroups through modulus2^14,1000 initial power sources,19 growth-prefix power receipts,17 all-ones oddpart receipts, and exact second-valuation congruences through b12. These finite checks audit the implementation, not the all-order proofs.

Classical parity-residue mechanisms are not claimed as new. Primary-source route: C. J. Everett, 'Iteration of the number-theoretic function f(2n)=n, f(2n+1)=3n+2', Advances in Mathematics25 (1977),42–45, https://www.sciencedirect.com/science/article/pii/0001870877900871. The publisher record and Lagarias's author-maintained bibliography route were checked: https://dept.math.lsa.umich.edu/~lagarias/3x%2B1.html. No inaccessible theorem from that paper is needed: the complete cylinder and subgroup proofs used here are above. Search results suggesting generalized3n+3^k dynamics were not imported because that is a different map.

Open next question: what source-sensitive, unbounded-depth descent or ownership invariant survives this full finite-language obstruction? A universal claim based only on forbidding finitely many valuation words cannot exploit the oddpart family. A valid repair needs source size, carry, an unbounded clock, an independently paid resource, or a genuine global event theorem.
