# Excess-two sources, consecutive Mersenne factors, and a physical edge repair

Status: **PROVED elementary constructions; INHERITED historical classifications; FINITE-EXACT experiments; OPEN global coverage.** No external priority claim. All paths in this note are repository-relative.

## 1. Provenance ledger and normalization

The requested sequence is **not the quasiperfect-number sequence**. Write sigma(n) for all divisors, s(n)=sigma(n)−n and Delta(n)=sigma(n)−2n.

- `07-reflections/perfect-aliquot-fixed-point-carrier-s645.md` and `05-knowledge/hypotheses/HYP-2221-perfect-aliquot-fixed-point-carrier.md`: the early fixed-point/defect-shell carrier; Delta=0 perfect, Delta=1 quasiperfect, Delta=−1 almost perfect. This is idea provenance, not a global existence theorem.
- `05-knowledge/results/ogg_triangular_chowla_20260927.md`, §5, and experiment `04-computation/experiments/ogg_triangular_chowla_20260927.py`: remove the unit divisor, obtaining Chowla's map c(n)=sigma(n)−n−1. Its fixed points are quasiperfect numbers, none known. The elementary parity argument gives square or twice-square; the stronger odd-square condition is inherited there as Cattaneo, not re-proved or externally checked here.
- `05-knowledge/results/aliquot_divisor_overlap_20261028/THEORY.md`: the factor-pair question s(n)=1+a+b, ab=n, is equivalent to tau(n)=4 or Delta(n)=2. The two branches are disjoint. Its bounded scan through 1,000,001 gives 20,104,464,650,1952,130304,522752. Its original tau≤10 classification and incoming tau12 extension are distinct contributions.
- `05-knowledge/results/session83_aliquot_divisor_grammar_20261028/FIXED_EXCESS_AND_FIBERS.md` and `LOW_TERM_ATLAS.md`: guarded prime-extension/low-arity history.
- `05-knowledge/results/dyadic_module_spectrum_20261026/ALIQUOT_OVERLAP.md`: the eleven-proper-divisor classification gives 650 and1952, recovering the fifth term without treating it as an extrapolated sequence.
- `05-knowledge/results/session141_metric_aliquot_square_paths_20261102/aliquot_twenty/PROOFS.md`, §5: complete all-size theorem with tau(n)≤17, exactly the five requested integers. Its `sourcepins/index.json` preserves earlier sources at commit c190e11b7795ed4e2d2f3aeb667e86d08238cdbf. The next example130304 has tau18. This is inherited, not a new classification.
- `05-knowledge/results/session258_endpoint_seed_20261111/SYNTHESIS.md`: endpoint cut law; a+2 newborn count occurs only as1→3 or3→5; strong-core seed classification; exact multiplier q versus odd slope c versus output N.
- `05-knowledge/results/session258_endpoint_seed_20261111/arithmetic/RESEARCH.md`, §§2–3: physical Carmichael twins, two-run shifted Mersenne products, and why source multiplicity does not fill an empty fiber.
- `05-knowledge/results/session259_arithmetic_graph_growth_20261111/SYNTHESIS.md`: sum/product transport and strict divisibility topology are arithmetic structures, not a tournament realization theorem.

For20 specifically, s(20)=22 and c(20)=21. Thus20 is mapped to its successor by the Chowla map, not fixed by it. For every listed term, c(n)=n+1. In particular c(20)=21 is a forbidden tournament H-count. This is a decisive hostile to any claim that Chowla outputs of excess-two numbers automatically lie in the Hamilton-path spectrum.

## 2. Board and new bridge

Anchor: exact Hamilton-path coverage. Niche: factor-pair compression and excess shells. Wildcard: removing one directed edge from a cyclic triangle blowup.

Live coordinates: (i) divisor source and prime guards; (ii) the neighboring target n+1; (iii) physical cyclic blowup; (iv) marked arc and path classes; (v) endpoint readiness. Closest inherited mechanism: two-run Mersenne-product graph families. Hostile:20→21. Corrected near miss: a near count does not supply its missing neighbor. Least-used sidecar: the actual edge used by each path, rather than just total H.

### Proposition A: the consecutive Mersenne bridge

For every integer a≥2 put

    p_a=2^(a+1)−3,  n_a=2^a p_a,  M_a=2^a−1.

Then

    n_a+1=M_a M_(a+1).

When p_a is prime, n_a is an excess-two number. Conversely the inherited two-prime classification says every excess-two number with at most two distinct prime factors has exactly this form with p_a prime. The prime guard is essential for this arithmetic assertion; it is unnecessary for all graph constructions below.

The identity is expansion. Under the prime guard,

    sigma(n_a)=(2^(a+1)−1)(p_a+1)=2n_a+2.

The terms20,104,464,1952 correspond to a=2,3,4,5. The interleaved650=2·5²·13 is not in this two-prime family. The formal a6 member8000 fails the prime guard (p6=125); its graph below still exists.

### Proposition B: a literal graph at n_a+3

Let T(a,b) be the tournament on transitive blocks A of size a, B of size b, and a singleton z, with A→B→z→A. Then, for all a,b≥1,

    H(T(a,b))=(2^a−1)(2^b−1)+2.

Proof: if z is first or last there is one path in each case. If z is interior, every path has the unique form

    A_left, B_left, z, A_right, B_right,

with each subsequence increasing in its transitive block, B_left nonempty and A_right nonempty. Choose the nonempty A_right and B_left independently. This gives exactly (2^a−1)(2^b−1) interior paths. No labels or endpoint conditions were discarded.

Thus T(a,a+1) is a strong tournament of order2a+2 with H=n_a+3, for every a≥2. For excess-two source n_a the Chowla target n_a+1 lies precisely two below this explicit strong family. This is an exact map, not a numerical analogy:

- source: dyadic exponent a, with optional prime guard;
- target: explicit T(a,a+1);
- preserved numerical relation: H=c(n_a)+2 under the guard;
- lost information: divisor lattice and primality; neither affects the graph;
- required sidecar for coverage: an actual operation realizing H−2, or a separate graph at c(n_a).

It does not show c(n_a) is impossible or possible. In particular105 is not excluded merely because the suggested construction gives107.

## 3. A new marked-edge response identity

### Proposition C: boundary crossing flip

Assume a,b≥2. In T(a,b), let u be the last A vertex and v the first B vertex. Reverse only u→v. The resulting tournament T*(a,b) satisfies

    H(T*)=H(T)+2−2^(a−1)−2^(b−1)
          =4AB−3A−3B+5,
    A=2^(a−1), B=2^(b−1).

Proof: only Hamilton paths using the changed adjacency can be gained or lost. In the original graph, the paths using u immediately followed by v number

    L=2+(A−1)B+A(B−1)=2+2AB−A−B.

The2 counts z at either end. When z is internal, uv lies in the left block or right block. The subset choices in Proposition B give the two displayed products.

For paths gained after reversal, contract the consecutive pair vu into a marked symbol w. All A\{u} vertices dominate w, w dominates all B\{v} vertices, and w cannot be immediately adjacent to z: neither z→v nor u→z exists. There are two end-z paths because both residual blocks are nonempty. For interior z, w lies either left or right of z; each choice requires a nonempty residual A-right and residual B-left. Thus

    G=2+2(A−1)(B−1)=2AB−2A−2B+4.

Therefore H(T*)−H(T)=G−L=2−A−B. QED.

The graph remains strong: a residual A vertex and a residual B vertex still connect all three blocks; u reaches B\{v}→z→A and v reaches z→A→u. The assumptions a,b≥2 are load-bearing.

### Corollaries useful for constructive coverage

1. The flip decreases H by2 exactly at a=b=2, giving11→9. For all a≥2 on the consecutive family b=a+1 its decrement is3·2^(a−1)−2≥4, never2. The most natural marked crossing therefore cannot perform the missing-unit repair.
2. For every a,b≥2, write p_a=2^(a+1)−3 and p_b=2^(b+1)−3. Then

       4H(T*(a,b))−11=p_a p_b.

   This is a physical shifted-product chamber with an explicit constructor and O(a+b) vertices. There is **no primality assumption**.
3. Taking a=b yields actual strong counts

       H=(p_a²+11)/4,

   including9,45,213,933 for a=2,3,4,5. The odd-square relation now follows from an actual path bijection. It is not a use of the hypothetical existence of quasiperfect odd squares.
4. At b=a+1 the repaired counts are n_a+5−3·2^(a−1), namely19,97,445,1909. These provide a second explicit strong family near the excess-two source, on the opposite side of its Chowla target.

The existing two-run chamber is4H−8=4M_aM_b. The new marked-edge chamber is4H−11=p_ap_b. They use different affine factor coordinates. This gives a cheap target-dependent arithmetic constructor: test whether4N−11 factors into two numbers each of the form2^r−3, r≥3; if so, recover the two block sizes and emit T*. This is a sufficient, sparse coverage rule, not an iff for all tournaments.

## 4. Exact experiments and remaining directions

`check.py` is a standard-library-only reproduction engine. It checks the explicit graphs by subset Hamilton-path DP, the smallest cases independently by permutation enumeration, the complete one-edge neighborhoods of T(a,a+1) for2≤a≤5, and the arithmetic identities without using sigma's prime formulas at composite p.

FINITE-EXACT: in those four complete one-edge neighborhoods, no reversal yields n_a+1. The smallest positive decrease is respectively4,10,22,46, attained by the marked boundary arc. This originally motivated the following **scoped conjecture, now PROVED by this session’s [all-edge theorem](../LOCAL_RIGIDITY.md)**: for every a≥2, every single-edge reversal of T(a,a+1) which decreases H decreases it by at least3·2^(a−1)−2. The experiment concerns this explicit family only, not all realizers of its H-count. Proposition C proves attainability; the later all-edge theorem proves minimality over every other arc.

Next high-value experiment: derive the full (i,j)-marked crossing response, retain rank within A/B, and classify which two-edge switches can hit n_a+1 for a≥3. A positive answer would construct an infinite aliquot-adjacent family without proving Chan–Pak. A negative answer would identify a precise local rigidity and motivate root changes. Another direction: combine the4N−11 chamber with ordinal products of already physical counts; product closure is lawful, whereas assuming an unattained Mersenne factor is physical is not.

The historical quasiperfect problem, universal odd-slope seed existence, all excess-two classification, and Chan–Pak remain OPEN. No verification-first detour or unsupported broad equivalence is being substituted for these constructive additions.

## 5. Positive closure: the entire neighboring family is decided

The parent recovered an important additional mechanism: **THM-4748**, `01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md`, statement4. For every k≥4 its Fibonacci fan F_k has H=2^k−1. Explicitly, the chain1,…,k+1 is transitive, and an apex0 loses exactly to floor(k/2), k+1−floor(k/2), and k+1. This is an inherited all-k physical constructor, not an inference from numerical coverage.

### Theorem D: complete classification for the dyadic formal family

For a≥2, the integer

    N_a=2^a(2^(a+1)−3)+1

is a tournament Hamilton-path count **if and only if a≥3**. No primality hypothesis is needed. In particular, among the two-prime excess-two sources,20 is the only source whose Chowla successor is absent from the spectrum.

Proof. For a≥4 take the ordinal sum F_a→F_(a+1). Every Hamilton path must finish the first block before entering the second, and therefore their counts multiply. Proposition A gives H=N_a. This is an explicit graph of order2a+5.

For a3, the seven-vertex tournament with outgoing bitmasks

    [116,45,80,101,10,84,18]

has exactly105 Hamilton paths. Both subset DP and exhaustive permutation enumeration in `check.py` certify the finite base case. Labels are0,…,6 and bit j in row i means i→j. For a2, N2=21, excluded by inherited THM-1370, `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`. QED.

This closes the positive branch that the one-edge probe left open. Source replacement is essential to the present construction: at a3 the new graph has7 rather than8 vertices; for a≥4 the constructed ordinal sum has2a+5 rather than2a+2 vertices. Failure of an edge move in the nearby strong graph is therefore fully compatible with actual realization. The constructed large graphs are not asserted strong or endpoint-ready.

### Proposition E: every even excess-two source has a Mersenne factor

Let n=2^a m, a≥1 and m odd, satisfy sigma(n)=2n+2. Then

    n+1=(2^(a+1)−1)·sigma(m)/2,

and sigma(m)/2 is an odd positive integer. Proof: multiplicativity gives (2^(a+1)−1)sigma(m)=2(n+1); n+1 is odd, so sigma(m) has exact2-adic valuation1. Divide by2.

Consequently, if a≠2 and sigma(m)/2 has a physical tournament realizer, the Chowla successor n+1 has one by ordinal sum: the Mersenne factor is3 for a1 and is supplied by F_(a+1) for a≥3. For a2 the factor is7, so this argument genuinely fails; it is not legitimate to substitute a graph with H7. A missing or forbidden half-sigma factor also blocks this particular construction, not the final product in general. The two-prime case has sigma(m)/2=2^a−1 and is precisely Theorem D's consecutive-Mersenne specialization.

For650 the decomposition is651=3·217. For the other four requested values it is21=7·3,105=15·7,465=31·15,1953=63·31. These factor sides explain why the seven-vertex105 base case cannot be replaced by its naive forbidden-factor decomposition.

The stronger question “every even excess-two source except20 has realizable Chowla successor” is **OPEN in this packet**, not a consequence of Theorem D. The next arithmetic task is to compile physical realizers of the half-divisor-sum coordinate, retaining its factorization and avoiding unjustified use of7 or21 factors.

### The exceptional650 is also closed constructively

Take the transitive chain1,…,9 and an apex0 that loses exactly to2,5,6,9. Call this graph A217. The inherited one-apex formula in THM-4748 gives

    H(A217)=1+(2+16+32+256)−(4+8+64+2+8+4)=217.

Join a directed triangle before A217. Its13-vertex ordinal sum has H=3·217=651, the Chowla successor of650. Both counts are independently checked by subset DP in `check.py`, which saves literal row masks in `receipt.json`.

Therefore the entire requested five-term source list has a complete successor answer:20 maps to forbidden21;104,464,650,1952 map to physically constructed105,465,651,1953. The all-even-excess-two extension above remains open. This final construction closes the exceptional three-prime entry rather than discarding it from the bridge.
