# Session94 independent source, port, clock, and enumeration audit

2026-10-29. Source pin: `88c969273f2b6bc4ea2a8c5d29bb36c558eb59d2`.

**VERIFIED / PROOF-REVIEWED**, with exact scope below. This audit imports no producer Python modules. It does not claim to re-prove every inherited theorem or audit an unrestricted class atlas. All code is standard-library Python; all checks raise explicit exceptions and survive `-O`.

## Reproduction

```
python3 check_independent.py
python3 -O check_independent.py
python3 check_sources.py
python3 -O check_sources.py
python3 check_final_classification.py
python3 -O check_final_classification.py
```

Normal and optimized outputs agree byte for byte: **2,137 + 125 + 122 = 2,384 always-active checks**. The count is not an estimate of proof strength. All three scripts are self-contained. They do not require stored receipts, producer scripts, or an external tournament bank.

## 1. Sharper constructive marking loops

The proposed bound is sound within the precise square-root re-marking operation menu:

\[
c_k=(3^k-2k-1)/4,\qquad \text{steps}\le 3c_k.
\]

A centralizing tree automorphism has a portrait coefficient at each internal node, with `a(-w)=-a(w)` and zero coefficients on all-zero prefixes. Choose one representative of each remaining pair by requiring its first nonzero digit to be1. These give exactly `c_k` free trits. Descending-depth execution of paired local rotations realizes the supplied portrait. The code checks all243 elements at depth3, all3 elements at depth2, and200 fixed-seed samples at depth4.

Each canonical paired node rotation is a literal based triangle. For an arbitrary pair, a centralizer conjugator with nonzero portrait entries strictly above the node sends its canonical prefix to the target and has trivial suffix sections; thus conjugacy produces the desired local generator exactly. Every such conjugacy through depth4 is checked as a full permutation, not just by group order.

Conjugating the three marking endpoints is not execution of the conjugator as an additional path. Hence three legal re-marking edges suffice per nonzero coefficient. Full endpoint descriptions, constructing their permutations, and applying them to named ports remain charged separately. This is an upper bound, not an optimality theorem, and not a bound on black/native tournament moves.

The portrait bijection also proves, under uniform `C_k`, that this prescribed protocol costs `3 Binomial(c_k,2/3)`, with mean `2c_k`. This is neither minimal-distance distribution nor an arithmetic-source law. Exact small-level cost multiplicities are independently checked.

## 2. Fixed points are insufficient for weighted subset attachments

Independent literal automorphism enumeration gives

```
A_2 fixed points: 0:62, 3:12, 6:6, 9:1.
C_3 fixed points: 3:124, 9:86, 15:24, 21:8, 27:1.
```

Uniform group sampling has `E X_2=1`; uniform sampling from the four values in its support has mean `9/2`. Replacing source multiplicity by support already breaks the first moment.

All512 subsets of the nine-vertex action are independently partitioned into orbits. The size inventory is

```
1,5,14,30,46,46,30,14,5,1;
```

there are192 binary-coloring orbits. The fixed-point-only triple formula gives `86/3`, not30. Internally permuted selected vertices require cycle data.

A stronger hostile is also independently enumerated: regular `C9` and `C3×C3` actions have the same entire fixed-point distribution, hence all ordered-port moments, but their three-subset inventories are10 and12. Even every moment order does not universally reconstruct the cycle index.

The producer's cyclic-wreath cycle-index recurrence and centralizer substitution `p_i→p_i²` were proof-reviewed. This duplicates each cycle; it does **not** double cycle lengths or independently resample paired sides. The polynomial theorem for fixed subset size r, degree r and leading coefficient `(d−1)^r`, follows by the stated lower-coefficient induction. Signed palettes are valid as presented because all cycle lengths are odd; arbitrary even-cycle transfer would require the full weight substitution.

## 3. CRT menus and a sparse full-carrier certificate

The odd-period ring normal-form proof was reviewed: the common divisor of `F_(n−1)−1,F_n` has no odd prime factor by the norm obstruction, and has exactly one factor2 when `3|n`. The resulting ring is `Z/L_n` or `F4×Z/(L_n/4)`. Prime-power factors, not just reduced residue fields, are necessary.

An independently implemented cyclic-word generator verifies bijectivity and full `{+1,reflection}` partition refinement at n9 and15, for observed moduli19 and11,31,341. The `{rotation,reflection}` reflected-double minimality proof is valid and operation-relative. Word reflection is not implicitly treated as ring conjugation.

For n25, the independent checker verifies the six legal-word witness pairs excluding every proper nontrivial divisor block of cyclic order167761. Combined with the elementary translation-coset theorem and the inherited odd reader bijection, this proves the finite primitivity/full-future-carrier statement without enumerating all167761 words in this audit. It does **not** imply all-period primitivity. The prime-power control checks order55 for37 modulo121 and order5 after reduction modulo11, together with the lost nilpotent11.

A wording caution communicated to the producer: equal fibre **cardinalities** along quotient orbits follows from bijections; arbitrary source weights need their own invariance or transport rule. No weighted source is silently made uniform.

## 4. Arithmetic valuation/tree factor and clocks

The exact inverse-branch cylinder argument was reviewed for `3n+1`, `3n−1`, and unrestricted Antihydra odd return. The punctured odd2-adic domain is essential: singular sources and their backward iterates must be removed. The older literal whole-space homeomorphism statement is not inherited unchanged. This puncture was already repaired in the Session64 dynamics/audit reports; current propagation of that correction is not a newly discovered error.

Under conditional odd Haar, valuation digits have law `P(a=j)=2^-j`. Parity gives critical0/3 offspring, but its clock is not independent. With `a=1+b+2G`, the residual G is independent of offspring; the whole a is not.

Independent slot-queue enumeration through five internal vertices verifies ternary tree counts `1,1,3,12,55,273`. For `N=3m+1` and total valuation `D=4m+1+2h`, the claimed exact mass is

\[
\frac{1}{3m+1}\binom{3m+1}{m}\binom{3m+h}{h}\,2^{-(4m+1+2h)}.
\]

A tree has m internal vertices and2m+1 leaves. Weak compositions of h among N vertices give the second binomial. Consequently `E[D|N]=2N−1/3`. This proof preserves breadth-first order, residual marks and the unused tail when claiming source recovery.

The finite-generation mixed-clock formulas were separately derived. If `C_k=Cov(Z_k,D_k)` and `V_k=Var(D_k)`, root conditioning gives

```
C_(k+1)=C_k+2/3+4k,
V_(k+1)=V_k+2+8k²+8k/3.
```

Thus the producer's closed formulas follow. The positive source1 has valuation2 forever, despite almost-sure Haar tree extinction. The negative cycle `−5↔−7` has phase words `(1,2)…` and `(2,1)…`, respectively dying immediately or never emptying the tree queue. These are not conclusions about positive-integer Collatz convergence. Standard, Terras and Antihydra clocks remain distinct.

## 5. A000568-adjacent pointed composition

The all-order marked SCC recurrence was proof-reviewed using the unique total order of SCCs of the union of at least one tournament layer. The fixed-singleton Burnside channel formula, differentiation of `A(1−S)=1`, and occupied-component recurrence are valid with occupied terms using at least one port. Union-strong need not mean each layer is strong.

Independent labelled Burnside enumeration—not the producer's class-bank method—through order5 gives:

| n | zero through three ordered-port counts | corresponding strong counts | two ports same/different SCC |
|---|---|---|---|
|4|4,12,44,172|1,4,16,64|28/16|
|5|12,48,224,1088|6,24,116,572|170/54|

At order2 the same/different counts are2/2. Summing occupied-component counts recovers the total, but separate marginal port counts cannot determine occupancy.

The prime-child-slot scalar jet method was also proof-reviewed. For `K=P∘F`, differentiation in `p1` changes only `F1`; all higher Adams components have zero `p1` derivative. Ordinary Bell-polynomial inversion is triangular with unit diagonal `A1^j`. Multiplication by `F^j` fills the j ordered distinct named child slots. This is not vertex-pointing: the strong order4 tournament has three root child slots but four vertices. The method recovers aggregate slot inventories, not individual prime tournament classes or a complete symmetry action.

## Scope and remaining obligations

No fatal discrepancy was found in the reviewed current claims. The independent scripts cover the exact finite ranges stated above; larger producer censuses remain producer finite checks, not independent tests by attribution. Higher-degree slot computations and arbitrary-depth proofs are supported by proof review, not an invented exhaustive replay. The inherited canonical tournament decomposition and all-odd legal-word bijection remain explicit dependencies.

## Additional extensions reviewed before closing this audit

**Saturated templates.** The active-depth compression argument is sound: retain first-nonzero tag depths and jointly folded suffix LCA depths. For r selected vertices, each occupied tag pair contributes at most one initial depth plus m−1 branching depths, so at most r total distinct active depths. Deleting inactive levels preserves side tags and all first differences; zero insertion is inverse on **orbits**, not literally on source subsets, using the inherited complete folded invariant. Consequently the binomial-basis template counts are nonnegative and Newton evaluation avoids iteration through huge depth. This does not remove the template's remaining joint LCA/orientation data.

**Divisor jets.** For odd d, `partial_pd Psi_i(F)=Psi_i(partial_p(d/i)F)` when i divides d, and zero otherwise. No extra factor i appears. This validates the producer's divisor-triangular first-derivative recurrence. The `1/d` in the direct Burnside derivative arises from canceling cycle multiplicity against the partition denominator. These are rational cycle weights, not integer class counts. In particular the interpretation of `H9[x3]` uses cycle type `(9,3)`, not `(9,1,1,1)`; the latter would make the nine-vertex orbit a module. The detailed 11/12 prime-class census interpretation remains attributed to its source, not independently re-enumerated here.

**Nine-candidate inverse reader.** Independent proof review sharpened the proposed44-candidate golden carry box to9. Put `delta=F/alpha`, with real coordinate x in `[0,1]` and conjugate y in `[-phi,1]`. Its coefficient offsets satisfy

```
delta_b=(x-y)/sqrt(5) in (-1,2),
delta_a=(phi*y+phi^(-1)*x)/sqrt(5) in (-2,1].
```

For the integral carry `aa+bb*phi` centered at `(A+B*phi)/L`, integrality gives only

```
aa=floor(A/L)+{-1,0,1},
bb=floor(B/L)+{0,1,2}.
```

The producer incorporated this9-candidate bound. The audit independently checks the resulting period55 sparse certificate using embedded legal words and Horner evaluation; it imports neither decoder nor producer verifier. Exact prime factorization and all22 proper nontrivial divisor blocks are checked. Thus, conditional on the explicitly inherited all-odd reader bijection, the translation/reflection action on **312,119,004,989** states is primitive. Verification uses44 word occurrences, not a state census. This is a finite explicit action theorem, not an all-period conjecture or a prime-degree group-classification claim.

The recurring distinction is precise: **group-uniform moments, invariant object weights, arithmetic source measures, and faithful actions are different data**. The new recurrences connect them only through specified maps and menus. No major open conjecture is settled by these recurrences alone.


## Final bounded frozen-source audit: classification and loop protocol

The frozen CRT continuation and parent loop packet received one final pass, with no further extensions.

- Reopened the primary [Jones paper](https://arxiv.org/pdf/1209.5169), Theorem1.2 case k=0, and its explicit CFSG dependency. Its exceptions are exactly prime affine degree, projective repunit degree, and degrees11/23 as used here.
- `check_final_classification.py` independently excludes projective degrees by **dimension-wise exact binary search**, rather than the producer's base-wise exhaustive loop. For each dimension d>=3 with `2^d−1<=N`, it locates the first integer base whose repunit reaches N and verifies strict overshoot. Higher dimensions are impossible. Dimension2 is excluded because N−1 is even and not a power of2. Exact trial division verifies every factor and full products for both degrees.
- A separate two-state path DP counts reflection-fixed words, obtaining377 and514229. Its forced-final-zero boundary is the reflected cycle's edge between the last paired vertices. Both reflection transposition counts are even; the full translation cycles are even too.
- Thus the stated identifications `G25=A167761` and `G55=A312119004989` follow from independently checked sparse primitivity, the inherited reader theorem, Jones/CFSG, and these exclusions/signs. They do not follow from reflection counts alone. No gigantic group was enumerated.
- Parent `loop_protocol/check_loop_protocol.py` was inspected and replayed in normal/optimized modes: **127,370 producer checks**, identical outputs. These checks include actual anti-involution endpoints for every paired node through k4, reversed triangles, descending portrait composition, and the exact uniform-source binomial cost. These producer checks are additional to—not included in—the2,384 independent audit checks.
- The frozen nine-carry decoder and Jones checker were also replayed normally and with `-O`, with identical outputs. The period55 certificate was already independently verified from embedded words without importing the decoder.

One wording correction was requested: the final two-transitivity obstruction should concern **loopless relations/on distinct ordered pairs**; diagonal relations remain invariant. This does not affect the impossibility of an invariant tournament on the same alternating action. No remaining fatal mathematical discrepancy was found within the reviewed scope.


## Frozen editorial addendum: positive-odd auxiliary extinction

The final dynamics corollary and root README were reviewed. The density-one proof is sound: for each finite division-clock cutoff, stopped words are a finite prefix-free family and their disjoint positive-odd residue cylinders have natural density equal to their conditional Haar mass. These finite-union densities approach1. The full auxiliary-extinction set contains every finite union, so its lower natural density is1; no unsupported exchange of infinite union and natural density occurs.

The exact hostile word `(2,1,1,1)` has packet `(P,B,A,r)=(81,103,5,41)`. Direct arithmetic independently confirms `41→31→47→71→107`. On every positive source in the cylinder `41 mod64`, its endpoint exceeds its start by `(49x+103)/32`. Hence this density result is explicitly about auxiliary queue extinction, not native descent or Collatz convergence. Source1 remains an exceptional nonextinct tree with a cyclic native orbit.

The README correctly separates aggregate recurrences, faithful carriers, source measures, clocks and classification dependencies. One minor editorial correction was requested: saturated-template coefficients are **nonnegative**, not necessarily positive. The live Antihydra sentence repair is reported by the publisher at commit `77683b7ced760e4fff7a3e2a431a417caf46a925`, propagating the already-proved Session64 punctured-domain correction; this audit does not claim a newly discovered theorem or independent remote mutation. No shared ledger edit was performed by this auditor.

This addendum closes the requested bounded review; no further mathematical extension was undertaken.
