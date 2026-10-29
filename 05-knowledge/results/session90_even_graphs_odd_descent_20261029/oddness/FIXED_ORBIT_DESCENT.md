# Oddness as fixed-orbit descent, with a deliberately forgetful path observer

Session90, 2026-10-29. **CITED / PROVED classical mechanism; PROVED elementary corollaries; FINITE-EXACT new comparisons.** No external-priority claim, new moonshine theorem, or major-open-problem closure.

## 1. Inheritance and the useful change of definition

Rather than treating “odd” as a property of a displayed number, ask whether the **group being divided out has order coprime to the symmetry being imposed**. The exact object is a compatible pair of actions, not a numerical resemblance.

Closest inherited mechanisms:

- `THM-409-self-converse-perspective-flip-involution.md`: converse acts canonically on ordinary vertex orbits, although an arbitrary anti-automorphism need not be an involution.
- `THM-643-blue-black-line-parity-structure.md`, `THM-644-gridsym-fiber-law-and-anti-redei.md`, and the later proved `THM-647-anti-redei-theorem.md`: reversal-twisted Hamiltonian paths, ordinary path-orbit parity, and the anti-Rédei theorem. THM-644's old conjectural status is superseded by THM-647, not a fresh open problem.
- Session89, `encoding/ODD_SELF_CONVERSE_LIFT.md`: paired consecutive tournaments correspond to odd self-converse tournament classes; retaining an anti-involution changes the automorphism group to its centralizer.
- Session32, `paley/PALEY_MONSTER_BRIDGE.md` and `monster/MONSTER_BABY_PRIME_PARTITION.md`: a genuine Monster-local quotient action and a genuine level-structure/Frobenius obstruction.

Canonical hostiles are the even-group examples below, the C3 joint-alignment example, and Paley7 with one mirror-pair arc flipped. Corrected near misses: a fixed-orbit theorem is not a full groupoid equivalence, and a count observer is not an attachment carrier. The least-used sidecar recovered here is the compact pair-mask path DP in `05-knowledge/results/hollow_frame_blue_extrema_20261020/replay_subset_dp.cpp`.

## 2. The fixed-orbit theorem

Let a finite p-group A act by automorphisms on a finite group G with p not dividing |G|. Let X be a finite G-set equipped with a compatible A-action:

`a(g x) = a(g) a(x)`.

Write `C=C_G(A)`, `X^A={x:a(x)=x for all a}`, and let A act on the ordinary orbit set X/G. Then the natural map

`X^A / C  ->  (X/G)^A`

is a bijection.

**Proof.** An A-stable G-orbit has cardinality dividing |G|, hence prime to p. Every nontrivial A-orbit in it has cardinality divisible by p, so it contains an A-fixed point. This proves surjectivity. If fixed x and y lie in the same G-orbit, their transporter

`Trans_G(x,y)={g:g x=y}`

is a coset of G_x, has cardinality prime to p, and is A-stable by compatibility. The same argument gives an A-fixed transporter, an element of C. This proves injectivity. No solvability theorem or classification of finite simple groups is required for this p-group case.

This is a special case of classical coprime-action/Glauberman fixed-point theory, **not a new theorem of this session**. Bibliographic primary record: George Glauberman, *Fixed points in groups with operator groups*, Math. Z.84 (1964),120–125, [DOI](https://doi.org/10.1007/BF01117119). A directly inspectable later primary presentation is Michael Burkhart, *Fixed point conditions for non-coprime actions*, [arXiv:2308.12286](https://arxiv.org/abs/2308.12286), introduction. Our elementary proof above is sufficient for every use here; no stronger non-coprime assertion is imported.

### Exact information boundary

For fixed x, the fixed object's stabilizer is `C ∩ G_x = C_{G_x}(A)`. It generally differs from G_x. Consequently the theorem preserves the **set of classes**, not their full original automorphism groups, orbit sizes, or groupoid masses. It is natural under maps compatible with both actions. Apply it to a space of joint tuples or relational structures to retain joint information; applying it separately to each marginal does not recover the joint orbit.

The two coprimality failures are independent:

1. G=C2 acts regularly on a two-point X; A=C2 acts by the same swap, commuting with G. X/G has one fixed class, but X^A is empty. Surjectivity fails.
2. G=C4 acts on X=G/C2. Let A=C2 invert G and act on X by inversion, which is the identity on this two-point quotient. Both points are fixed. But C_G(A)={0,2} acts trivially on X, giving two fixed-point orbits over one fixed G-class. Injectivity fails.

## 3. Constructively repair a transporter, without searching the group

Specialize to odd |G| and an involution alpha acting on G and X. If fixed x,y satisfy g x=y, put

`h = g^{-1} alpha(g)`, `a = sqrt(h)`, `c = g a`.

Then h belongs to G_x and alpha(h)=h^{-1}. In an odd-order finite group, h has a unique square root in its cyclic subgroup, and uniqueness gives alpha(a)=a^{-1}. Thus

`alpha(c) = alpha(g) alpha(a) = g h a^{-1} = g a = c`,

and a fixes x, so c x=y. The repaired transporter centralizes alpha and respects all data already built into x and y.

For a given permutation g and a given involution implementing alpha, this correction is linear in the permutation degree in the ordinary array-operation model: each odd cycle of h is rotated by `(length+1)/2`. It requires no enumeration of G. It does **not** find the initial transporter g or a fixed representative of an arbitrary orbit. Those search tasks can still be expensive.

This corollary is the same square-root repair principle as the companion carrier lane's conjugator of two anti-involutions, applied to an arbitrary fixed relational object. It is an explicit operational use of classical descent, not a polynomial-time graph/tournament isomorphism claim.

### Shared-port hostile

Let G=C3 act by translation on X=Y=C3, and let alpha negate each coordinate. Each marginal quotient X/G and Y/G is a singleton fixed class. But diagonal G-orbits on X×Y are indexed by `d=y−x in C3`; alpha sends d to −d. Only d=0 is fixed. Thus:

- marginal fixed class counts: 1 and 1;
- joint alignment classes: 3;
- fixed joint alignment classes: 1.

Oddness supplies a fixed representative **after the joint class is specified and shown stable**. It does not conjure the lost relative alignment from the two marginal classes. For any regular G-torsor pair with chosen alpha-fixed origins, the same calculation identifies fixed relative alignments with C_G(alpha), rather than all of G. This is the relevant boundary for LCA ports, clone mergers, and later named guards.

## 4. Odd self-converse tournaments and a fast, smaller observer

Let W have order 2n+1, let tau be a supplied anti-involution, and put G=Aut(W), C=C_G(tau). A tournament automorphism group has odd order: an involution would exchange two vertices and reverse their arc. The compatible involution on Hamiltonian paths is

`P -> reverse(tau(P))`.

The theorem gives a bijection from tau-fixed paths modulo C to converse-fixed ordinary path orbits. Since ordinary automorphisms act freely on Hamiltonian paths, if h_tau counts the fixed paths, then the blue presentation count is exactly

`b(W) = h_tau(W) / |C|`.

This recovers and makes explicit the orbit bookkeeping in THM-643/644. THM-647 already proves h_tau odd by pairing all nonfixed Hamiltonian paths and using Rédei's theorem. Do not call this new anti-Rédei parity.

### The pair-mask recurrence

Tau fixes one vertex r and exchanges n pairs. Every fixed Hamiltonian path is uniquely

`v_1,...,v_n,r,tau(v_n),...,tau(v_1)`,

where the first half picks one vertex from each pair. It suffices that the first half follows arcs and v_n points to r; all remaining arcs follow by anti-invariance.

For a set S of used pairs and a signed vertex v in one of those pairs, let D(S,v) count such partial first-half paths ending at v. Initialize `D({pair(v)},v)=1` and update

`D(S union {pair(w)},w) += D(S,v)`

whenever pair(w) is unused and v points to w. Sum D(all pairs,v) over v pointing to r. For n=0 the value is 1.

There are at most `2n*2^(n−1)` endpoint states and O(n² 2^n) elementary transitions, using O(n 2^n) integer cells. Bit costs depend on the integer sizes. A straightforward full Hamilton-path subset DP on 2n+1 vertices uses O(n² 4^n) transitions. These are different observables: the folded DP computes h_tau, **not H(W)**. No measured runtime ratio or complete-enumeration speedup is asserted.

The old `hollow_frame_blue_extrema_20261020/replay_subset_dp.cpp` already contains this pair-mask recurrence. A tempting 3^n signed-subset state is unnecessarily large: once an old pair has been used, its chosen sign cannot affect any future continuation except through the current last vertex. The contribution here is recovering the smaller existing mechanism, porting it to the consecutive-pair lift, and exposing the centralizer-dependent boundary.

### An exact eight-object hostile family

Take Paley7 on F7, x→y iff y−x is in {1,2,4}, with tau(x)=−x. Independently reverse any subset of the three internal tau-pair arcs {1,6},{2,5},{3,4}. The folded recurrence never reads these arcs, so **all eight objects have h_tau=9**. Direct permutations independently give:

| Number of reversed pair arcs | Number of masks | H | |Aut| | |C| | b=h_tau/|C| |
|---|---:|---:|---:|---:|---:|---:|
|0|1|189|21|3|3|
|1|3|159|1|1|9|
|2|3|129|1|1|9|
|3|1|99|9|3|3|

Thus a coordinate can be free to discard for one count while remaining essential for full structure and even for its symmetry-normalized version. This eight-object comparison is FINITE-EXACT; the all-order independence of internal pair arcs follows directly from the first-half path description.

### Why two square-matrix counts agree globally but not fibrewise

For uniformly random normalized n×n matrices M, the mean of `h_tau(lift(M))` is exactly `n!`. There are `2^n n!` candidate first halves. Each requires n independent fair-bit constraints: n−1 edges on distinct unordered pair-pairs, plus the last vertex's root incidence. Each cross edge is either one matrix bit or the XOR of the two bits for that unordered pair-pair, hence fair; the involved coordinate blocks are distinct. Thus each word occurs in exactly a `2^(-n)` fraction of the matrices. This proves

`sum_M h_tau(lift(M)) = n! 2^(n²)`.

A pair class with centralizer C has `n!/|C|` normalized labelled matrices, while its blue presentation fibre has `h_tau/|C|` elements. Therefore the sum recovers the inherited blue mass `sum_classes b=2^(n²)`. It does **not** identify those two fibres class by class. For Paley7, n=3, the matrix fibre has2 elements but the blue fibre has3; for the transitive seven-vertex class those sizes are6 and1. This provides a small exact obstruction to a class-preserving “same square bits” trijection. A general set bijection based only on total cardinality is not ruled out.

## 5. The Monster comparison: genuine common mechanism, different hypotheses

Session32 recovered a genuine local action: for a Monster 7A element g, its cyclic normalizer modulo the Held kernel is `N_M(<g>)/He ≅ 7:6 ≅ AGL(1,7)`. Its square-multiplier subgroup 7:3 is Aut(Paley7); the nonsquare coset reverses orientation. This is an affine/coset action, not conjugation on the seven elements of <g>. The quotient discards the Held kernel and does not identify a tournament with the Monster. The explicit subdirect-product description is in Breuer–Magaard–Wilson, [*Some steps in the verification of the ordinary character table of the Monster group*](https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/doc/ctblm.pdf), §4.13; that source also records the 2A centralizer 2.B in §3. These imported facts consume only the cited local action and centralizer, not any VOA or Hamilton-path claim.

The Baby Monster analogy is marking an involution and retaining its centralizer. But the Monster has even order and two involution classes, 2A and 2B. For example, in M×C2 the elements (1,t), (a,t), (b,t), with a in 2A and b in 2B, are involutions in the same outer coset but not conjugate under M. There is no odd-kernel uniqueness theorem there. The literal mark's conjugacy sector is indispensable. In the tournament signed group, by contrast, the ordinary kernel has odd order and all anti-involutions are conjugate. That is exactly why the Session89 class lift can forget its mark while changing its stabilizer.

A second inherited bridge concerns **Frobenius-fixed objects after retaining extra level data**. Aricheta, [*Supersingular Elliptic Curves and Moonshine*](https://arxiv.org/html/1809.07421), Theorems1.1 and1.3, relates supersingular rationality with level to modular-curve genus and Fricke elements in Monster centralizers. Session32 independently computed that primes29,41,59,71 have all supersingular j-vertices rational, but at level two each has a single conjugate pair of nonrational coarse isogeny edges. This is an explicit example of a fixed base object with nonfixed extra structure, analogous to the C3 joint-alignment hostile in §3. It is not an application of our odd-group theorem: no suitable odd quotient action is supplied.

## 6. What happens to {7,21}, and what remains open

Keep four different statistics separate:

1. H(T), the number of ordered Hamiltonian paths;
2. H(T)/|Aut(T)|, ordinary presentations;
3. h_tau(T), paths fixed by a supplied anti-involution;
4. b(T)=h_tau(T)/|C_Aut(T)(tau)|, converse-fixed presentations.

The all-order omissions7,21 for H are inherited `THM-1370-h-spectrum-omits-7-21-all-n.md`, with its explicit 2016 priority correction. The proof uses strong-component multiplicativity and strong Hamiltonian-path floors. It is not a consequence of Monster primes or fixed-orbit descent. Paley7 has respectively `(189,9,9,3)` for these four statistics, so its vertex count7 and group order21 cannot be substituted for any of them. Moreover the literal lift census below attains h_tau=7 repeatedly. This alone blocks transporting the H-exclusion into the folded observer.

The useful unifying rule is operational: **coprime fixed-orbit descent can remove an auxiliary symmetry mark from class existence, but its centralizer and joint alignment must remain wherever later operations observe them.** The folded path observer illustrates a different compression, justified by a Markov recurrence rather than by class reconstruction. These are compatible mechanisms, not one universal “free bit” principle.

Open directions: characterize the optimal observer state under specified allowed updates; exploit centralizer-orbit aggregation of the pair-mask DP without losing transition multiplicities; compare the resulting weighted state construction with the paired class generator. Nothing here solves graph or tournament isomorphism, spectrum completeness, LRC(14), moonshine, or a major open conjecture.

## 7. Reproducible evidence

Run `python3 check_fixed_orbit_descent.py` and `python3 -O check_fixed_orbit_descent.py`; outputs must agree. Standard library only, no repo import, no network, explicit active checks. The final run executes30,760 checks.

- 391 compatible actions of odd cyclic groups of orders≤31 on all their quotient torsors, with 7,991 explicit transporter repairs.
- Both even-group failures above.
- Paley7 discovered by all7! permutations: Aut21, centralizer3; all128 vertex subsets, 12 ordinary subset orbits, 8 fixed orbit classes,16 fixed subsets;84 repaired transporters.
- Paley7 paths:189 total,9 ordinary orbits,9 tau-fixed paths,3 centralizer orbits;63 paths in the three blue ordinary orbits.
- All531 labelled consecutive-pair matrix lifts at n=0,1,2,3. Folded counts agree with literal Hamilton-word enumeration. All4,165 internal-pair-flip controls preserve the folded count.
- Exact C3 joint control and all eight Paley mirror-pair flips.
- Lift n=3 folded histogram:1:48,3:64,5:96,7:192,9:112. These are labelled matrix inputs, not unlabelled tournament class counts.

Source pin for all repository references: `3eb49fa63d7f` (Session90 startup). Published paths named above are under `01-canon/theorems/` unless otherwise stated; Session32 means `05-knowledge/results/session32_square_monster_20261023/`, Session89 means `05-knowledge/results/session89_consecutive_tournament_pairs_20261028/`.
