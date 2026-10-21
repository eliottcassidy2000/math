# Session24 arithmetic recovery: 341, the bracket primes, and signed tournament structure

Date: 2026-10-21. Read-only inheritance from `math/`; no Git operations or repository edits. This packet is a research note plus independent exact checks, not a canon promotion. No claim of literature novelty is made.

## Executive result

1. The original `{2,3,11}` is an exact **odd-square-bracket escape set**. It is not originally a multiplicative-order set. Its later appearances in integer Collatz cycle periods, Mills-type steps and the eleven-cell golden model have different mechanisms and must not be identified.
2. The original `341=11*31` is a member of the Collatz inverse trunk `N_j=(4^j-1)/3`, and its supposed failure to produce primitive primes was already refuted. Both 11 and 31 are primitive divisors at exponent five.
3. A much newer, genuinely structural 341 appears in `golden_signed_relation_geometry_20261019.md`: `Z[phi]/(phi^15-1, phi-3^4) = Z/341`. Golden time is multiplication by **81**, not 4.
4. Joining these actual maps produces a complete small tournament family at composite order341. The trunk affine map `R(n)=4n+1` is conjugate to multiplication by4 under `L(n)=3n+1 mod341`. On the common carrier, `<4,81>=QR_11 x QR_31`, of order75. Exactly **16 labelled connection sets / four tournament isomorphism classes** preserve both clocks. Eight sets give lexicographic Paley products; eight give **prime (module-free) tournaments**, and a graph-invariant certificate proves the full automorphism group imprimitive in all cases. The two prime forms have full group exactly `(C11 semidirect C5) x (C31 semidirect C15)` and a three-vertex regular-subset certificate.
5. The independent Fano recovery maps the squareclasses on primes `{2,3,11}` to seven vertices. The inherited Hamming-coset condition then picks16 of128 line orientations, all Paley7 with189 Hamiltonian paths; all112 other line orientations have171. This is a real map but does **not** distinguish `{2,3,11}` from any other triple of primes.

The new useful common principle is: **retain every support/conductor stratum and its sign action, not just the units or the shared order.** At341, multiplication by2 and by3 preserve many skew sections on units but no complete tournament, because a different CRT zero stratum carries a flipped orbit in each case.

## 1. Inheritance and correction ledger

Paths below are relative to `math/`. Every status is local to the stated claim.

### A. Original bracket meaning (PROVED elementary, recovered)

`05-knowledge/results/procgen_numerology_20260923_odd_square_brackets.md`, §§1–2:

- `B_n=((2n-1)^2,(2n+1)^2]`.
- A prime p and its double share `B_n` only if `2(2n-1)^2<(2n+1)^2`, hence n<=2.
- In `(1,9]` the primes are2,3; in `(9,25]` the only qualifying prime is11. Thus the set is exactly `{2,3,11}`.

`05-knowledge/results/procgen_brackets_20260924_pairings_transitions.md`, §1.4 sharpens it:

- The same prime escape set holds for every multiplier `r in (25/13,25/11]`.
- The complete pairs `(k,p)` with k>=2 integer and kp in p's own bracket are `(2,2),(3,2),(4,2),(2,3),(3,3),(2,11)`.
- `rho(p)=(next odd square)/p` is the faithful scalar threshold. Ordering it is only a transitive ranking; it does not itself define a useful cyclic tournament.

**Correction:** `01-canon/MISTAKES.md`, entry “2026-09-24 Collatz tree/brackets prompts”, and the same pairings note §1.1 distinguish the real-interval cutoff from the integer cutoff. The former floor formula can overshoot by one; the latter is `max{m:r<=(2m+1)^2/((2m-1)^2+1)}`.

### B. Original341 (PROVED elementary, recovered)

`05-knowledge/results/collatz_mod6_20260917_wild_typing.md`, §3 (especially lines64–74); and `05-knowledge/results/collatz_mod6_20260917_scaffolding_audit.md`, §6:

- `N_j=(4^j-1)/3=R^(j-1)(1)`, `R(n)=4n+1`, and `3N_j+1=4^j`.
- `N_5=341=11*31`, with `ord_11(4)=ord_31(4)=5`. Both factors are new primitive divisors at j=5. “341 stalls primitive prime generation” is **REFUTED**, not an open conjecture.
- For j>=3, `ord_(N_j)(2)=2j`; base2 pseudoprimality is equivalent to `6j | 4^j-4`, while base4 pseudoprimality is equivalent to `3j | 4^j-4`.
- **Corrected near miss:** these criteria do not always coincide. At j=4, N_4=85 is base4-pseudoprime but not base2-pseudoprime; the missing coordinate is `v_2(N_j-1)=2`.

This packet needs no general Zsigmondy or Cipolla citation: the facts used at341 are proved or checked directly. The historical note credits those classical results rather than claiming them as new.

### C. Newer341 (PROVED elementary in the source note; not independently audited canon)

`05-knowledge/results/golden_signed_relation_geometry_20261019.md`, §3, especially lines69–100:

`O=Z[phi]`, `phi^2=phi+1`, and

`O/(phi^k-1,phi-3^r) = Z/gcd(3^(2r)-3^r-1,3^(rk)-1)`.

At `(k,r)=(15,4)`, the gcd is341, and phi maps to81. This is the largest cyclic quotient for that specified generator identification. It is a map of rings/time labels, not an integer Collatz convergence result.

The same note §4 supplies an important hostile to extrapolation: modulo121, tripling has order5 while the golden lift37 has order55. The order11 agreement does not survive prime-power lifting; `Z/121` also differs from the121-element split ring `F_11 x F_11`.

### D. Which eleven is actually Paley (scoped canon, recovered)

`01-canon/theorems/THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md`, statements1,3,4:

- For the standard (not Terras) Collatz map C on its period-five rational points, the golden parity reader kappa satisfies `kappa(Cx)=4kappa(x)` in F_11.
- `<3>=<4>=QR_11`; the55-element affine group is regular on55 unordered pairs and has the Paley11 orientation as an arc orbit.
- `4*3=1 mod11`, so tripling is **backward**, not forward, golden time.

`01-canon/MISTAKES.md`, **MISTAKE-584**, overrides historical prose: all-odd-period bijectivity was not proved; the recorded computation is finite-exact. The named paper's P is multiplication by5, Q by4. The period10 projection has a factor2 twist. None of those corrections is erased here.

`01-canon/theorems/THM-4591-integer-cycles-of-the-terras-map-to-period-301993-are-the-five-pythagorean-cycles.md`, Remarks, separately records periods1,2,3,11 in its finite-exact signed-integer census. It explicitly warns that the cycle-period11 is not an explanation of Golay11 or the eleven-cell torus. `05-knowledge/results/seam_mills_20260925.md`, §5, separately proves the only prime solution of `q=p^3+3` is p=2,q=11, by parity. These do not explain the bracket theorem.

### E. Signed symmetry and Fano (recovered mechanisms)

- `01-canon/theorems/THM-4725-no-uniform-dilation-rule-for-prime-circulants-signed-set-burnside-and-the-p-mod-4-dichotomy-with-the-paley-rule.md`, §2, Theorem B: a section fixed by an action exists exactly when no orbit flips the commuting antipodal involution. Its stated canon status is **PROVED (pre-audit)**. Only this elementary signed-set observation is used; no dependence on the file's unrelated Chebotarev claim.
- `05-knowledge/results/arithmetic_braids2_20260917_fano_code.md`, §§1–2: rank-four line-incidence map, Hamming coset, syndrome correction, and finite-exact171/189 split. This is older work, not a new result of Session24.
- `01-canon/theorems/THM-1745-leaf-graded-arborescence-filtration-721-shadow.md`, §1: free automorphism action on Hamiltonian paths and on two-leaf spanning out-arborescences. The path divisibility below is also immediate directly.

## 2. Transfer I: two arithmetic clocks become a signed tournament family

### 2.1 Exact source, map and preserved predicate

Source A is the trunk affine dynamics modulo341. Since3 is invertible,

`L: Z/341 -> Z/341`, `L(n)=3n+1`, satisfies `L(Rn)=4L(n)`.

Source B is the golden ring quotient from §1C, with phi-time sent to81. Thus the **common state carrier** is `A=Z/341`, equipped with commuting multipliers4 and81. This is not merely a coincidence of their periods.

Target is a Cayley tournament T_S on vertices A, with

`x -> y iff y-x in S`, `S disjoint (-S)`, `S union (-S)=A\{0}`.

The orientation gauge is the choice of one side in every antipodal difference pair. There are no ties for distinct vertices. The preserved predicate is that both clock maps and all translations are tournament automorphisms. Pulled through L, the first clock is literally R.

**Lost information:** integer size, positive sign, denominator, actual Collatz branch bit, and first-hit stopping time. None is claimed recoverable from the tournament. Retain those separately for any dynamical use. The quotient also does not inherit field arithmetic, prime-degree theorems, or prime-circulant type classifications.

### 2.2 General conductor-stratified count (PROVED elementary)

For odd m and a unit a, write X_d for differences of additive order d, for each divisor d>1 of m. It has phi(d) elements, and multiplication by a has orbit length `ord_d(a)` on X_d.

A fixed skew section exists iff no multiplication orbit meets its negative. This is equivalent to `ord_p(a)` odd for every prime p|m, and equivalently `ord_m(a)` odd. Indeed an even order modulo an odd prime creates a half-order value -1 on the additive-order-p stratum; conversely odd global order cannot carry x to -x. For completeness, if ord_p(a)=d is odd, then a^d=1 modp and the binomial theorem gives a^(d*p^(e-1))=1 modp^e. Thus all prime-power orders, and their global least common multiple, remain odd.

If it exists, its number is

`2^( sum_(d|m,d>1) phi(d)/(2 ord_d(a)) )`.

Proof: the distinct a-orbits pair freely under negation, and one orbit of every pair is chosen. This is exactly the fixed-section clause of the recovered signed-set principle, now retaining all additive conductors.

At m=341, the orders at d=11,31,341 are:

| clock | d=11 | d=31 | d=341 | full fixed sections |
|---|---:|---:|---:|---:|
| 2 | 10 | 5 | 10 | 0 |
| 3 | 5 | 30 | 30 | 0 |
| 4 | 5 | 5 | 5 | 2^34 |
| 81 | 5 | 15 | 15 | 2^12 |

**Hostile to units-only reasoning:** on the300 units, clock2 has `2^15` skew sections and clock3 has `2^5`; nevertheless neither extends to a full invariant tournament. Explicit obstructions:

- `2^5*31 = -31 mod341`, on the additive-order11 stratum;
- `3^15*11 = -11 mod341`, on the additive-order31 stratum.

The two involutions are32 and309 modulo341. Neither is the global -1=340. Looking only for -1 in the unit multiplier subgroup therefore misses the obstruction.

### 2.3 Joining the two lawful clocks (PROVED elementary)

Under CRT,4=(4,4) and81=(4,19). The ratio81/4 has coordinates(1,28), with28 of order15 in F_31. Therefore

`K=<4,81>=QR_11 x QR_31`, `|K|=5*15=75`.

On differences, there are four unit orbits of size75, two nonzero-first/zero-second orbits of size5, and two zero-first/nonzero-second orbits of size15. Negation pairs these eight orbits into four pairs. Hence there are exactly16 joint-clock connection sets.

The full exhibited group is

`G=A semidirect K`, `|G|=341*75=25,575`.

It preserves both CRT coset partitions. Its action is transitive but imprimitive. This subgroup observation ALONE says nothing about the full automorphism group. The separate invariant proof in §2.6 controls the full group.

### 2.4 Four normal forms, and the prime/module distinction

Let (p,q) be either (11,31) or(31,11), and let chi_p,chi_q be the two Paley characters. On CRT coordinates(a,b), for delta=+1 or-1 define the sign of a nonzero difference by

- `s(0,b)=chi_q(b)`;
- `s(a,0)=delta chi_p(a)`;
- `s(a,b)=chi_p(a)` if a and b are both nonzero.

The section is s=+1. Independent signs in front of chi_p and chi_q provide four labelled gauges per normal form. These are precisely the16 sets: on the four nonzero CRT quadrants an antipodal Boolean sign must be one of `+chi_p,-chi_p,+chi_q,-chi_q`; the two axes supply the remaining choices.

For delta=+1, this is `P_p[P_q]`, the lexicographic Paley product. Its fibres of sizeq are modules. For delta=-1, it is an **axis-reversed product**: reverse the P_p tournament in each fixed-b slice, while leaving the mixed-coordinate arcs unchanged.

**PROVED for all distinct odd primes, not just341 (structural closure, verified here).** The axis-reversed product is prime for any two skew Cayley tournament factors on Z/p and Z/q. The proof does not use the Paley characters:

1. Every Cayley tournament on Z/(pq) is strong: its connection set contains1 or-1, giving a Hamiltonian cycle.
2. In a strong tournament, the maximal proper modules are disjoint. Intersecting modules have a module as their union. If two maximal proper modules overlap without inclusion, their union must therefore be all vertices. Write their nonempty regions as `X=A\B`, `I=A intersect B`, `Y=B\A`. Module homogeneity forces `X -> I -> Y` and `X -> Y`, or all reversed, contradicting strong connectivity.
3. Every vertex lies in a maximal proper module, so these modules form a canonical invariant partition. In a vertex-transitive decomposable tournament their common size is strictly between1 andpq.
4. Under regular translations the block M through0 is a subgroup: for m in M, `M+m` meets M at m, so it equals M. Thus M is additively closed and, being finite, is a subgroup. At orderpq the only nontrivial possibilities are the two CRT axes.
5. Neither axis is a module of the axis-reversed product. A vertex(0,b), b!=0, distinguishes vertices(a,0) by the positive and negative first-coordinate differences, so `Z/p x {0}` is not a module. A vertex(a,0), a!=0, distinguishes(0,0) from every(0,b), b!=0, because the former difference uses the reversed sign and the latter the unreversed sign. Thus `{0} x Z/q` is not a module either. Contradiction.

This applies to arbitrary skew Cayley factors at distinct prime orders, not only p,q=3 mod4 Paley factors. The two-axis subgroup list uses distinct primality. The companion `COPRIME-EXTENSION.md` extends the proof to arbitrary coprime odd factor orders by checking every rectangular subgroup. It also gives a C3 x C3 diagonal-module hostile, so coprimality cannot simply be dropped.

**Subgroup-constancy test.** More generally, the same canonical-module argument proves that a cyclic Cayley tournament is prime iff no nontrivial proper additive subgroup is a module. A subgroup H is a module precisely when each nonzero coset of H is wholly in S or wholly outside S. At orderpq there are only two such subgroups, so two linear scans of the connection vector decide primeness. The checker compares these axis-coset tests against every pair-closure result for all16 cases.

**Independent FINITE-EXACT control at341:** all16 gauges were tested by minimal-module closure for every anchored pair{0,d}, d=1..340. Translation covers every pair of vertices. The two reversed normal forms always close to341; the two aligned forms have proper closures of sizes31 and11. See §4. The exact tests predated the structural proof and remain as independent evidence.

All four normal forms are nonisomorphic. For each oriented edge x->y record the number of common out-neighbours. The following profiles count one representative per difference, so multiply every frequency by341 for the full isomorphism invariant:

| dominant p | delta | prime? | common-out-neighbour value:frequency |
|---|---:|---|---|
| 11 | +1 | no | 77:155, 162:15 |
| 11 | -1 | yes | 77:80, 79:75, 152:15 |
| 31 | +1 | no | 82:165, 167:5 |
| 31 | -1 | yes | 82:90, 84:75, 137:5 |

The multiplier group sends each four-gauge set transitively to itself; the four distinct profiles separate them under arbitrary isomorphism. Thus there are exactly **four isomorphism classes**, not just four multiplier types.

**Why this is useful to primitive-action proofs:** prime/indecomposable tournament does not imply that an exhibited transitive automorphism subgroup acts primitively. The two prime examples above carry a concrete imprimitive group with two visible block systems. The full-group claim requires the extra proof below. Conversely an invariant block system need not consist of tournament modules. Any transfer from a block-system argument to a tournament decomposition needs that extra predicate.

### 2.5 What reaches Hamiltonian counts (PROVED, sharply limited)

G acts freely on ordered Hamiltonian paths, because fixing such a path fixes every vertex. Consequently, for every member of the family,

`25,575 | H(T)`.

The inherited two-leaf free-action proof gives the same divisibility for the total number of two-leaf rooted spanning out-arborescences. These are honest structural invariants transported from the actual action. No value of H at341 has been computed, no equality of H between the four forms is asserted, and no extremality or pure-blue conclusion follows.

### 2.6 Full automorphism group and a signed-regular subset certificate

This upgrade is NOT inferred from the exhibited subgroup. For each normal form, the maximum common-out-neighbour count in the table occurs exactly when the two vertices have the same dominant CRT coordinate. The checker verifies this equivalence for every nonzero difference. Thus a genuine graph invariant recognizes the disjoint union of p cliques K_q. Every automorphism preserves that partition, proving that the **full** automorphism group is imprimitive.

For the reversed form, between any two blocks the majority direction is precisely the P_p orientation: q(q-1) arcs have that direction, and the q opposite arcs are exactly the matching of equal b-coordinates. Any automorphism therefore induces sigma in Aut(P_p) on blocks. Write its within-block permutations as tau_a on the b-coordinates. Preserving each exceptional equality matching forces `tau_a(b)=tau_a'(b)` for all a,a',b, hence one common tau. Internal block orientation forces tau in Aut(P_q). Conversely each such pair is an automorphism. Therefore

`Aut(T_reversed)=Aut(P_p) x Aut(P_q)`.

A separate complete backtracking search in `paley_automorphisms.py` checks the two factor groups without importing any prime-degree classification theorem. Fixing0 gives exactly the square multipliers:5 maps at11 and15 at31. Translation then gives full orders55 and465. The searches visit52 and572 nodes and retain every candidate consistent with assigned directed pairs. Hence, at341,

`Aut(T_reversed)=G`, `|Aut|=55*465=25,575`.

The three-corner subset `U={(0,0),(1,0),(0,1)}` has trivial stabilizer. Its two occupied rows have sizes2 and1, fixing0 and1 in the first factor. Its two occupied columns likewise fix0 and1 in the second. An affine map fixing two distinct points is the identity. CRT identifies U with `{0,155,187}` in Z/341; direct enumeration of all25,575 affine automorphisms confirms its stabilizer is only(1,0).

Here the **signed action is the ordinary permutation action on subsets** by automorphisms and anti-automorphisms. There is NO set complementation. Negation is an explicit anti-automorphism, and every anti-automorphism differs from it by an automorphism, so `Aut^±(T)=G union (-1)G`, of order51,150. This full signed group preserves both CRT partitions. If an anti-automorphism stabilized U, its coordinate maps would fix both distinguished rows and both distinguished columns individually. But a tournament anti-automorphism cannot fix two distinct vertices, because it would have to reverse their edge while fixing its endpoints. Thus no anti-automorphism stabilizes U. Direct enumeration of all51,150 signed affine permutations independently confirms that only the identity stabilizes `{0,155,187}`. U is genuinely signed-regular, and its signed orbit has size51,150. Its size3 is minimal: every2-set {u,v} is stabilized by the anti-reflection x->u+v-x, and every singleton has an ordinary stabilizer containing a conjugate of K, of size75.

## 3. Transfer II: the bracket primes label a Fano/Hamming repair mechanism

This is a recovered mechanism specialized to the requested prime triple, not a special theorem about its escape property.

### Exact object and map

Use the squareclasses supported on E={2,3,11}. Assign basis labels2->001,3->010,11->100 and map products to parity vectors of exponents. The nonzero labels are

`2,3,6,11,22,33,66`.

Multiplication followed by square removal is XOR. A Fano line is `{x,y,x+y}`; in integer labels the seven lines are

`{3,11,33}`, `{2,11,22}`, `{6,11,66}`, `{2,3,6}`, `{3,22,66}`, `{2,33,66}`, `{6,22,33}`.

Every pair belongs to exactly one line. Choose a cyclic orientation on every line; these seven choices orient all21 pairs and produce128 tournaments. All are regular since each vertex lies on three cyclic lines.

### Gauge, syndrome and preserved target

Index points and line normals by nonzero vectors in F_2^3, with `L_a={x:a dot x=0}`. Changing basis signs by a bit vector g changes line signs by

`(Bg)(a)=sum_(x in L_a) g(x)=S+a dot v`.

The image has dimension4 and is the `[7,4,3]` Hamming code; the kernel has dimension3. In the recovered Cayley–Dickson convention the valid gauge orientations are the affine coset `96+im(B)`, not the zero-containing linear code itself. Their relative syndrome is `sum_a (m(a)+96(a)) a`.

A nonzero syndrome a is repaired uniquely by reversing all three arcs of line L_a. It is not an arbitrary single-edge flip. Independent checks here reproduce the inherited results:

-16 zero-syndrome orientations, all isomorphic to Paley7, each H=189;
-112 nonzero-syndrome orientations, each H=171;
-the repair increases H by18 in this bounded family;
-two line errors can decode to a different valid gauge, so this is not arbitrary-error recovery.

### Loss ledger and stopping boundary

The map preserves squareclass multiplication, labelled Fano incidence, and the precise sign-gauge predicate. It forgets exponent heights and magnitudes, prime order, and the odd-square-bracket escape predicate. Three arbitrary independent primes produce exactly the same mechanism. The prime labels do not select an orientation, the incidence alone does not compute H, and H is not generally invariant under changing gauge. Here the equality189 is an independently verified property of this particular coset.

The useful transfer to a tournament proof is the **syndrome sidecar**: before reading symmetry or path counts from a seven-line quotient, retain its three-bit compatibility defect. The price of the quotient is explicit, and the cheapest hostile is a two-line error. A literal connection from bracket escape to this defect has **not** been found.

## 4. Exact verification and reproduction

Files:

- `check.py`: standalone standard-library checker, no repo imports;
- `receipt.json`: exact output;
- `receipt-optimized.json`: identical output from `python3 -O`;
- `paley_automorphisms.py` and `paley_automorphisms.json`: independent exhaustive factor-automorphism check;
- `COPRIME-EXTENSION.md` and `coprime_extension.py`: proved coprime-odd-order extension, exhaustive45/15-vertex controls, and a9-vertex noncoprime hostile;
- `probe_joint.py`: earlier exploratory version, not required for the final receipt.

Run:

`python3 session24-arithmetic/check.py > session24-arithmetic/receipt.json`

`python3 -O session24-arithmetic/check.py > session24-arithmetic/receipt-optimized.json`

`cmp session24-arithmetic/receipt.json session24-arithmetic/receipt-optimized.json`

745 explicit checks; failure checks use exceptions, not assertions. The module-closure and Fano computations also compare independent algorithms.

### Why module closure proves the finite claim

Given a pair M={u,v}, if a vertex z outside M has different orientations toward two members of M, every module containing M must also contain z. Repeatedly adjoin all such distinguishers until stable. The stable set is itself a module and is contained in every module containing the original pair; it is therefore their minimal closure.

A tournament is prime iff every distinct pair closes to the whole vertex set. Translation invariance permits u=0. The code checks all340 choices d and all16 sections, not just sampled pairs. A bitset queue implementation is independently checked by naive set/adjacency expansion on representatives of all eight K-orbits for every section. The orbit reduction is only used for the independent cross-check, not the primary complete test.

The Fano path census is checked by subset dynamic programming and all7! ordered permutations for every one of128 orientations. The Paley identification of all16 gauge orientations is separately checked over every permutation; the base gauge has exactly21 such identifications, agreeing with the known affine action without assuming its full automorphism group in advance.

## 5. What is new, what is recovered, and what remains open

**Recovered:** the original bracket definitions;341 trunk and corrected pseudoprime criteria; newer golden341 quotient; signed-set fixed-section mechanism; Fano/Hamming map and171/189 census; existing free-action divisibility.

**New synthesis in this packet:** the common carrier with literal trunk conjugacy and golden quotient; complete joint-clock341 family; units-only sign obstructions on opposite CRT axes; four-class classification; the all-distinct-prime primeness proof and independent341 controls; full-group CRT block certificate; three-corner signed regular subset; explicit primitive-action hostile; the resulting25,575 path divisibility. This is a checked research result, not a literature-priority claim or canon promotion.

**Honest frontier:** the module-partition argument closes primeness for all distinct-prime skew Cayley factors, and §2.6 pins down the two prime341 full automorphism groups. Hamiltonian counts at341 remain uncomputed. The companion note proves the broader coprime-odd-factor primeness theorem; no universal noncoprime-factor theorem or universal full-automorphism formula is claimed. The Fano construction does not explain why the bracket primes are exactly{2,3,11}. None of these finite carriers proves Collatz, LRC14, or a global tournament classification.
