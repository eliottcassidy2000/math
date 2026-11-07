# Session219: Fermat factors, exact fibers and constructive Chan–Pak progress

**COMPLETED: proved arithmetic identities, a sharp computer-assisted all-length grammar theorem, explicit tournaments and independent audits. Chan–Pak completeness and Carmichael’s totient conjecture remain OPEN.** No historical-priority claim.

## 1. The main constructive result: exactly five runs break the Fermat-factor constraint

Put p=641 and q=6700417, so pq=F5=2^32+1 and both primes give2 multiplicative order64. In the inherited physical source grammar, O and Z add vertices while retaining the root; J changes it. A run is a maximal nonempty block of one letter, excluding the final C1-closing J.

**Every O/Z recipe with at most four runs obeys q|H => p|H. Five runs suffice to escape.** The explicit source

    O^31 Z^15 O^9 Z^5 O^9

closes to a72-vertex strong tournament with

    H=1182825546590390812645
     =176530139331685·6700417,     H=275 mod641.

This is an exact minimum in RUN COUNT, not graph order. Each positive run has period64 in its modular action, with residue0 represented by length64, not an empty run. Exhausting the34,087,040 phase recipes with1–4 runs proves the obstruction for ALL lengths. Three different implementations agree. Increasing the final run by64t gives infinitely many five-run escapes, without any source root changes.

A shorter root-changing prefixOOZJOZJJOJJOZ, followed by20+64t O letters and closure, gives

    H_t=4909+6013·2^(20+64t),    order36+64t,
    H_0=6305092397=941·6700417,  H_t=380 mod641.

Full adjacency witnesses for both seeds are provided. An independent weighted Hamilton-path DP on the16 outside vertices proves the COMPLETE36-vertex count and its infinite tail formula. It counts one/two slots for the transitive O block, yielding coefficients16935 and6013. The72-vertex full count uses the inherited exact physical recurrence, with independent adjacency and connectivity checks.

[Five-run theorem](fermat/FIVE_RUN_THRESHOLD.md) · [Phase mechanism and shorter escape](fermat/FERMAT_GRAMMAR.md) · [Independent full-graph count](fermat/FULL_GRAPH_COUNT.md) · [Adversarial audit](run-audit/REVIEW.md).

## 2. Ford–Konyagin inspires exact factorization control

The [joint1999 paper](https://www.ford126.web.illinois.edu/wwwpapers/sigma.pdf) amplifies a known totient-fiber multiplicity k to rk. Its permitted prime suppliers encode exactly r matchings; excluded supports prevent unintended preimages. Ford’s separate1999 Annals paper proves all multiplicities>=2. They are distinct results.

Applying the controlled-decomposition idea to binary Mersenne generators M_n=2^n−1, n>=2, gives a complete normal form:

    the sole relation is M_6=M_2²M_3=63.

Replace every6 by2,2,3. The resulting exponent multiset is unique, by the cited Bang–Zsigmondy primitive-divisor theorem. Thus63^(r−1) has exactly r unordered Mersenne factorizations for EVERY r>=1. At a fixed number of factors, the unordered factorization is unique. The unit M_1=1 must be excluded.

This also proves uniqueness of the unordered exponent pair in the two-run tournament formula H=2+M_(a+1)M_(b+1), a,b>=1. But the only flexible Mersenne relation uses7, a forbidden Hamilton-path count, so it is not itself a legal tournament factorization. The distinction identifies precisely what the arithmetic analogy lacks.

## 3. A sharp all-parameter tournament scaling theorem

Let w(h) be the maximum number of nontrivial strong components in a tournament with H=h. Using the inherited legal blocks3,49,63,343 and only the proved gaps7,21, we obtain the exact optimum:

- w(3^a)=a.
- For b=1,a>=2: w(3^a7)=a−1; a=0,1 are the two gaps.
- For b>=2: w(3^a7^b)=a+floor(b/2).

The upper bound pays one component from the valuation resource a+floor(b/2); the special b=1 case must spend at least two3s to absorb its7. The existing THM-4094 packaging attains the bounds. We credit that construction; the new refinement is sharp optimality throughout the plane.

In particular **w(63)=1 but w(63^k)=floor(5k/2) for every k>=2**. Squaring changes63 from an indecomposable strong count to a product with five strong components. The same number controls the eventual polynomial degree, in vertex order, of the fixed-H tournament-class count. This is not strong-class multiplicity.

[Normal-form proof, smooth-degree classification and explicit source words](ford/PROOF.md) · [Independent construction audit](audit/CONSTRUCTION_REVIEW.md).

## 4. Merge14,26,34,38 without losing the coordinate

[Session218](../session218_smith_repunit_totient_20261107/README.md) supplies the certified bridge and its boundaries:

| Starting value k | Half q=k/2 | phi(k)=phi(q), also Smith base here | First totient repair of target k |
|---|---:|---:|---|
|14|7|6|28=phi(29), multiply by2|
|26|13|12|52=phi(53), multiply by2|
|34|17|16|136=phi(137), multiply by4|
|38|19|18|1216=phi(1217), multiply by32|

The full even-nontotient prefix is14,26,34,38,50,62,68,74,76. The corrected Camion exponent set through64 is4,6,10,12,14,26,34,38,62. These are different predicates:50 separates them. Simultaneous Mersenne/Wagstaff primality at q>5 does prove a genuine common subfamily at k=2q; no converse or infinitude follows.

The scaling operations above act on different coordinates: exponents, phi inputs, missing phi outputs, or a chosen digit base. The new Mersenne theorem acts on factorization descriptions; the new run theorem acts on physical recipes. None identifies their entire sequences.

## 5. A007374 and the Fermat-hole multiplicity law

Write A(m)=#{x:phi(x)=m}. [OEIS A007374](https://oeis.org/A007374) gives the least m with A(m)=n, indexed by n>=2. It is neither a nontotient sequence nor a Camion exponent sequence.

We extend Session218’s membership criterion to full inverse-fiber counts:

    sum A(2^a)x^a=(2−x)/(1−x) product_(F_i prime)(1+x^(2^i)).

F5=641·6700417 therefore causes the exact drop

    A(2^a)=a+2 for0<=a<=31;
    A(2^a)=32 for32<=a<=63.

Doubling2^31 to2^32 decreases multiplicity33→32. No unknown Fermat primality is used. The complete prime-supplier convolution counts every preimage of2^a q, including the exceptional q² supplier when q is a Fermat prime. Second differences recover supplier births, which mere zero/nonzero membership forgets.

The two F5 factors cannot replace the missing pure dyadic generator: phi(641)=2^7·5 while phi(6700417)=2^7·3·17449. A multigraded Euler product preserves these different odd sectors. This is why equal order64 does not imply equal totient behavior.

[Exact fiber formulas, general Fermat-hole theorem and deconvolution](totient/FERMAT_FIBER_TRANSFERS.md) · [Independent arithmetic audit](audit/DYADIC_REVIEW.md).

## 6. What now helps a Chan–Pak proof

The combined contribution is concrete: a sharp complexity obstruction was found and repaired by an explicit richer recipe; arithmetic factorization fibers were classified exactly; and existing smooth-target packaging was proved optimal. Ford–Konyagin’s useful lesson is to control intended AND unintended sources, retaining their support or endpoint data.

The next global obstacle remains direct construction in the prime and7-times-prime lanes identified by THM-4094. Multiplicative packaging cannot manufacture a prime target. Modular coverage alone cannot prove integer coverage. Our five-run theorem concerns a particular divisibility test, not completeness of five-run words or of the full grammar. Exact totient multiplicity is not a map to tournament-class multiplicity.

The concurrent bounded-root-change preregistration at6d0e88f0e is preserved and supplies no theorem dependency. It counts final-D word descriptions, whereas this packet separately tracks C1 profiles, run count, actual graphs, multiplicative decompositions and strong-component degree.

## Reproduce and validation limits

From this directory run

    python3 verify_all.py

Requires Python3 and a C++17 compiler;25 processes replay all substantive certificates. Included independent checks cover177,100 bounded Mersenne products,439 attained smooth-resource targets,245 producer plus130 independent inverse-fiber comparisons, complete small-target totient sieves, all64 pure-dyadic exponents,34,087,040 run-phase recipes by independent methods,728 small literal-graph path counts,18 larger prefix checks and the full36-vertex weighted count. The first ten A007374 minima are recovered by complete small-target computation; no new OEIS term is claimed. Individual reports state exact universes and distinguish all-parameter proofs from finite regressions.

SHA256SUMS freezes the packet. Publication uses small verified checkpoints and exact remote readback; the shared dirty checkout and concurrent work are preserved. There is no global spectrum proof, totient-conjecture proof, or proof that only five Fermat primes exist.
