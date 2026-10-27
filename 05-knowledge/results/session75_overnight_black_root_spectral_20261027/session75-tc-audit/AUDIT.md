# THM-4812 {7,21} presentation-fibre audit

## Verdict and pin

**The all-order exclusion `tc(T)=H(T)/|Aut(T)| not in {7,21}` passes an independent extension audit, conditional on the inherited n<=11 theorem THM-4790, OCF THM-002 and the classical group/bound theorems identified below.** No counterexample was found. Three compressed/incorrect proof sentences in the frozen source need the explicit repairs below. The source's original large order-12 permutation-type scan was **not** rerun; a separately justified, much smaller orbit-window computation replaces that gate.

Frozen source commit: `87d61f67673a6b10cbafff1bf9e1377f3c41e56c`.
Read-only source checkout: `/workspace/scratch/07d5532f2fce/session75-repo`.
Own packet: `/tmp/session75-tc-audit`; independent window packet: `/tmp/session75-tc-window-audit`.
No shared-checkout modifications, installs, commits, pushes or remote messages were performed by this audit.

`019882e00ccf2738bc39fc4f998daabcfb5ed7fc` was inspected only as a separate diff. Its new part(g), generalized567-case search and exact small-fibre class counts are **not included** in this acceptance. Parts(b)-(d) retain the old wording there.

## Scope distinction

- `H(T)` is the number of Hamiltonian paths. Its all-order exclusion of7 and21 is the inherited THM-1370, not the new result audited here.
- `tc(T)=H(T)/|Aut(T)|` counts path presentations modulo automorphisms. THM-4790 previously established its exclusion only through order11 (and records a prior independent audit).
- This audit independently closes the new all-order extension for `tc`, with a replacement order12 gate and reproduction of the83-case high-order window.
- This does **not** prove that every other odd value occurs for either invariant.
- Part(e)'s almost-all-prime touchability is already existing S73 work, not a new uniqueness proof. No claim that2 and5 are the only prime untouchables is accepted here.

## Analytic audit

### Orbit product

For any vertex partition, OCF's positive terms supported within the blocks give `H(T)>=product H(T[B])`. For the **actual** `Aut(T)`-orbits, restriction is a faithful embedding into the product of the full induced-tournament automorphism groups. Thus division gives `tc(T)>=product tc(T[O])`; each induced orbit tournament is vertex-transitive. It is not legitimate to substitute an arbitrary partition in this last vertex-transitivity assertion. The code's later prescribed subgroup orbits may be smaller than actual Aut-orbits, but that only over-enumerates: the actual group is among the enumerated groups, and the final full-Aut count decides candidates.

An alternative justification of the H product bound is repeated deterministic merging of vertex-disjoint directed paths while preserving each path's order; the final word determines the input paths. OCF remains needed for the triangle-packing bound.

### Vertex-transitive small orders

Automorphism groups of tournaments have odd order, so transitive degrees are odd. At prime degree a Sylow prime element is a full cycle, giving a circulant.

For degree9 choose a Sylow3 subgroup containing a Sylow3 subgroup of a point stabilizer: its point orbit has size9. Embed it in the standard Sylow subgroup `C3 wr C3` of `S9`. Its intersection with the base `F3^3` is nonzero, and invariant under cyclic coordinate shift. Every nonzero invariant submodule contains the diagonal line (the shift is a single unipotent Jordan block). A lift moving the three blocks either has order9 and is a regular9-cycle, or has order3 and, with the diagonal3-cycle, generates a regular `C3^2`. Thus every degree9 vertex-transitive tournament is Cayley on `C9` or `C3^2`. This removes any dependence on completeness of the source's15-regular-tournament input.

Independent finite corroboration: all50 subgroups of the standard81-element Sylow group were enumerated;13 are transitive,9 contain regularC9 and the remaining4 contain regularC3^2.

At degree15, Feit–Thompson and the prime-power degree theorem for solvable primitive groups force a block system of sizes3 or5. With five3-blocks the quotient isC5, the nonzero base kernel is a shift-invariant submodule of `F3^5`; with three5-blocks it is a submodule of `F5^3`. The respective cyclotomic factors are irreducible (orders of3 mod5 and5 mod3 are4 and2), so each kernel is diagonal, augmentation or full. Transporting coordinates with a complement of prime order makes its action pure block shift. The diagonal/full cases contain regularC15. The augmentation cases contain exactly the stated groups `A3 semidirect C5` or `A5 semidirect C3`. Their projection on every two distinct blocks is surjective, so between-block arcs are uniform, yielding exactly the8 orientations of `R5[C3]` or the8 orientations of `C3[R5]`.

Independent pull-Held–Karp plus a simple exact permutation backtracker checked every526 Cayley connection sets:

| group | sets | minimum tc |
|---|---:|---:|
| C3 |2|1|
| C5 |4|3|
| C7 |8|9|
| C9 |16|39|
| C3^2 |16|39|
| C11 |32|1729|
| C13 |64|94949|
| C15 |128|162165|
| C17 |256|792878221|

The small values appear only for C3, R5 and P7 (plus K1). For degree15, source Python normal/-O replay and an independently constructed16-substitution input both give `197030475/1215=162165` and `193215375/375=515241`. `n!<2^64` for n<=17, so these C++ H counts cannot overflow. Brute permutation H checks cover all Cayley inputs through order7.

### Greedy triangle bound: repair required

Put `C=m(m^2-1)/24` and `E=D+1=(3m^2-19)/8`. The product for k disjoint triangles is a valid lower bound only while every factor `C-iE` used is positive; it must stop at the first nonpositive factor and never resume after signs turn positive again. The independently implemented rational sum with that truncation exceeds `23*3^((m-1)/2)` at all192 odd m from19 through401. Its smallest ratio is about2.6168, at19.

The source's claim that its floor-term logarithmic margin increases at every m is false: the reported margin is32.77469 at401,31.97459 at403 and25.52175 at419. The decrease from401 to403 is exact because `403^60 < 3*401^60`.

A complete replacement tail is elementary. Let `k=floor(m/20)>=20`. Then

`C-kE >= (11m^3+37m)/480 >= 11m^3/480`,

so the single positive term satisfies

`2^k i_k >= (11m^3/240)^k/k! >= ((11/12)m^2)^k >= ((1100/3)k^2)^k`.

Since `m<=20k+19`, the required right side is at most `23*3^(10k+9)`. The exact base inequality

`440000^20 > 23*3^9*177147^20`

proves the result at k20, and `(1100k^2/3^11)^k` increases for all integer k>=20 (its base is already>1). This proves the entire tail without floating-point acceptance or false pointwise monotonicity.

### Finite order window: repairs required

For allowed exceptional counts `(r,p)` define

`B_v(n)=max 5^r*21^p*3^floor((n-5r-7p)/3)`

subject to `3^r*9^p<=v` and nonnegative remaining order. For v7 the pairs are(0,0),(1,0); for v21 they additionally include(2,0),(0,1). Once n>=10, all such pairs fit, hence `B_v(n+3)=3B_v(n)` while Busch gives `f(n+3)=5f(n)`. Three residue bases13/14/15 close the v7 tail;23/24/25 close the v21 tail, with21 separately excluded. A finite check to199 alone would not prove the infinite tail.

For the remaining v21 orders13..20,22 there are exactly24 configurations. Their maximum index cap is6. The source's statement "odd and<=6, hence1 or3" omits5. Exact enumeration repairs this: the only configuration with cap>=5 is n13 with two C3 blocks and one P7 block, product group order189;5 does not divide it. All configurations whose product group order is divisible by5 have cap<5.

An index3 subgroup of an odd-order group is normal: the coset action has odd transitive image inS3, hence imageC3 and kernel equal to the subgroup. All such subgroups are kernels of nonzero linear functionals on the3-primary abelianization, modulo nonzero scalar. The source generators give precisely these kernels. Transitive projections are explicitly checked.

## Finite window reproduction and independent order12 replacement

The separate window audit reconstructed the complete Cartesian product groups and every eligible subgroup/orbital input independently; normal and `python -O` jobs are byte-identical. It replaced nauty (unavailable) with an exact orbit-stabilizer automorphism counter, validated against brute permutations on all33,867 labelled tournaments through order6 and structured controls. H arithmetic was widened to128 bits and each final H checked against uint64; module block-run counts were cross-checked with Held–Karp as below.

| gate | configurations | subgroup jobs | strong generated tournaments | full-Aut candidate calls | module/HK cross-checks | hits |
|---|---:|---:|---:|---:|---:|---:|
| source13..20,22 |24|83|185642|48908|93646|0|
| independent replacement n12 |6|45|46968|10104|45384|0|

Every one of the83 normalized output lines matches the source's stored result, including73 module and10 non-module cases. The order12 gate has33 module and12 non-module cases. Repeated generated labelled tournaments across subgroup cases are expected; these are **not** counts of distinct isomorphism classes or distinct tournaments over the full union.

At order12, `f=375`; for tc21 the only configurations `(K1,C3,R5,P7; |P|, index cap)` are:

- `(3,3,0,0;27,1)`
- `(0,4,0,0;81,4)`
- `(5,0,0,1;21,1)`
- `(2,1,0,1;63,3)`
- `(1,2,1,0;45,2)`
- `(2,0,2,0;25,1)`.

Thus only indices1 and3 occur. For tc7 only `C3^4` with index1 remains, already contained in this union. The same complete orbital-enumeration proof applies, so the alternate gate excludes both7 and21 at order12. It does not claim to reproduce the source's18-million-scale order12 scan or its3,491,072 Aut calls.

The source order12 type lemma was also reviewed: its group arguments are sound with elementary Sylow support details supplied, and its40 hyperplane fact was independently replayed. This is supporting analysis, not a substitute claim that the old finite scan ran.

## Primary literature and inherited dependencies

- [Dixon, The maximum order of the group of a tournament, Theorem1](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/682531E83B0A3AA042CFC5852474E914/S0008439500067461a.pdf/div-class-title-the-maximum-order-of-the-group-of-a-tournament-div.pdf): directly checked the `3^((n-1)/2)` automorphism bound and its odd-order/solvable-group scope. Its proof also states the elementary-abelian regular normal subgroup and prime-power degree structure for solvable primitive groups.
- [Busch, A Note on the Number of Hamiltonian Paths in Strong Tournaments, Corollary1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/): directly checked all three residue-class minima used as f(n).
- THM-4790's n<=11 exclusion and prior independent audit, and THM-002's all-order OCF, were inspected but their historical finite computations/external OCF proof were not rerun here.

Part(f)'s strong-fibre finiteness argument was checked algebraically: at most `log_3 v` orbit pieces have tc>1; their total order is at most `M(v)log_3 v`; Dixon yields `|Aut(T)|<=3^(n/3)v^(M(v)/6)`, and Busch gives the stated order bound. The repaired greedy tail proves `M(v)=O(log v/log log v)`. This is finiteness of strong isomorphism classes only; unrestricted fibres are infinite by stacking tc1 components.

## Reproduction and remaining boundaries

Run `./reproduce.sh` in this packet for portable normal/-O exact analytic/group checks, independent Cayley census and degree15 controls. The window packet has its own scripts, manifests and exact-aut controls. No external install is required; Python3 and a C++17 compiler suffice.

Remaining boundaries: source old order12 scan not rerun (replaced); n<=11 and OCF inherited; generalized part(g) not audited; no uniqueness claim for prime untouchables; no full odd-spectrum completeness claim. These are explicit scope boundaries, not unresolved gates in the repaired new exclusion proof.
