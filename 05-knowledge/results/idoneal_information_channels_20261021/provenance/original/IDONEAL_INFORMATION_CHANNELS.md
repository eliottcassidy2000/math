# Idoneal class information, compatible registers, and actual path reflections

2026-10-20/21. **PROVED elementary / CITED classical genus theory / FINITE-EXACT form controls.** No new idoneal-completeness proof, Collatz convergence theorem, or six-colour forcing theorem is claimed. The positive maps retain their specified group operations; the obstructions retain actual source labels and typed clocks.

## 1. Completeness is a separate input, not a shortcut through this proof

The classical status is carefully distinguished from this repository's adopted input. Brandhorst–Sonel–Veniani, *Idoneal genera and K3 surfaces covering an Enriques surface*, Remark1.5, records65 known idoneal numbers, completeness under GRH, and at most two possible additional numbers without that assumption. The older one-extra assertion applies to squarefree exceptions and must not silently replace the bound for all orders.

The repository separately adopts openai/math #003. The October5,2026 *Quasi-Riemann Hypothesis* preprint explicitly states a zero-free half-plane theorem and derives completeness of the65-number list from the absence of Landau–Siegel zeros, class-group computations, and conductor reduction. We directly inspected the public TeX source (Git blob2000faddbbebac5de0ecfe0b962534ea61a852d5). This note does not independently re-audit that analytic proof. It records completeness as an inherited consequence of the adopted input, not as something proved by the constructions below.

None of our class-group, reflection, register, or Collatz examples needs completeness of the list. They concern explicit small quadratic orders and are unconditionally certified by reduction of positive binary forms.

Primary sources:
- [Kani, Idoneal Numbers and Some Generalizations, Theorem6 and genus discussion](https://mast.queensu.ca/~kani/papers/2011K1.pdf).
- [Brandhorst–Sonel–Veniani, Remark1.5](https://arxiv.org/html/2003.08914v7).
- [Elsenhans–Klüners–Nicolae, class groups of small exponent](https://arxiv.org/abs/1803.02056), whose stated computation/completeness theorem retains its zero-free hypothesis.
- [OpenAI, Quasi-Riemann Hypothesis, public source](https://github.com/openai/math/blob/main/preprints/The-Quasi-Riemann-Hypothesis-October-5-2026/build/paper2.tex).
- Repository inheritance: THM-4566-hilbert-class-field-planes-are-the-idoneal-planes.md, golden_tail_idoneal_20261017c.md §6, escape_renewal_idoneal_20261017d.md §§6–7, digit_idoneal_difference_20261018b.md. Those already distinguish orders from fields and genus codes from discarded phase.

## 2. An exact information-loss theorem

Let A be the proper ideal class group of the quadratic order of discriminant D=-4d. Classical genus theory identifies its genus group with A/A², where A²={a²:a in A}. The idoneal condition is exactly A²={1}, equivalently every class has order at most2. This is an order statement, not merely a statement about the maximal order of Q(sqrt(-d)).

For any finite abelian group A, each fibre of pi:A->A/A² has cardinality |A²|. For a uniformly selected class X, elementary Shannon entropy therefore gives

 H(X | pi(X))=log_2 |A²|.

Proof: conditional on a coset, the uniform class is uniform on its |A²| elements; average this constant conditional entropy. In particular the genus channel loses no class information iff d is idoneal. No asymptotic or distributional claim about integer sources is involved.

This is a precise meaning of “information flow”: the hidden subgroup is explicitly known, and its disappearance is equivalent to the target property. A residue register or an actual-source flow requires its own kernel analysis; cardinality resemblance does not transfer this theorem.

## 3. Source63 exposes the conductor and the hidden bit

For d=63, D=-252=6²(-7). Thus the quadratic order has conductor6 in the field Q(sqrt(-7)). Its complete list of primitive reduced positive forms is

 (1,0,63), (7,0,9), (8,-6,9), (8,6,9).

Reduction gives one representative per proper class. The last two forms are distinct inverses; the first two are self-inverse. Therefore the class group has order4 and an element not of order2: A is cyclic of order4. Its square subgroup has size2, and the genus channel loses exactly one bit. In contrast the maximal order of discriminant-7 has class number1. Forgetting conductor6 erases all four order classes and would give the wrong idoneal verdict.

An explicit genus character is chi=(-3/n) for an integer n represented primitively and coprime to252. The first two forms have chi=+1; the last two have chi=-1. This can also be seen directly modulo3: the first two forms represent nonzero squares modulo3, while the last two represent twice such squares. Hence the two source classes in each genus cannot be reconstructed from chi alone.

No class-group composition is inferred from a numerical count alone: order4 together with the displayed non-self-inverse class rules out C2xC2 and forces C4. The standard inverse-class rule is (a,b,c)->(a,-b,c).

## 4. Source273 gives a literal lossless graph code

For d=273, D=-1092 is fundamental. Its complete reduced primitive forms and one positive represented integer n coprime to1092 are:

| form | n | ((-4/n),(-3/n),(-7/n)) |
|---|---:|---|
| (1,0,273) |1| (+,+,+) |
| (2,2,137) |137| (+,-,+) |
| (3,0,91) |103| (-,+,-) |
| (6,6,47) |47| (-,-,-) |
| (7,0,39) |67| (-,+,+) |
| (13,0,21) |73| (+,+,-) |
| (14,14,23) |23| (-,-,+) |
| (17,8,17) |17| (+,-,-) |

All eight forms are self-inverse, using the standard reduced-form boundary identifications when |b|=a or a=c. Thus A is C2³ and273 is idoneal. The prime-discriminant decomposition is (-4)(-3)(-7)(13)=-1092. Classical genus theory makes the displayed three characters homomorphisms; the explicit eight distinct rows prove that together they give a group isomorphism A->F2³ after replacing each minus sign by1.

This yields a concrete, operation-preserving graph map, not an arbitrary labelling of eight objects. On vertices0,1,2,3 encode the three bits(b1,b2,b3) by graph edge-orbits

 {01,23}:b1,  {02,13}:b2,  {03,12}:b3.

These three edge-orbits are the three perfect matchings of K4. Their eight unions are exactly all regular simple graphs on four labelled vertices: one empty graph, three matchings, three complementary four-cycles, and the complete graph. Thus the target is a recognizable whole graph family, not an arbitrary eight-element selection. Symmetric difference is closed on this family. Call the graph G_a. Then

 G_(ab)=G_a symmetric-difference G_b,
 rho(G_a)=G_a=G_(a^-1),  rho=(03)(12).

Let eta be the actual class represented by(6,6,47), whose character code is(1,1,1). Then

 G_(eta a)=complement(G_a).

The positive graph operation matching complement is translation by this specified nontrivial class, NOT class inversion. Inversion is the identity on an idoneal class group and therefore cannot correspond to fixed-point-free graph complementation. The full source and three chosen genus characters make every map reversible.

Apply the ordered-graph staircase map of the reflection-resolution packet to get eight explicitly marked five-vertex tournaments. Graph complement becomes reversal of the nonpath arcs; rho becomes converse followed by reversal of the fixed Hamiltonian order. All eight selected paths are blue. Graph regularity in this ordered chart does not mean tournament regularity: the empty graph gives the transitive five-tournament, which is not regular. This does not assert that their entire underlying isomorphism classes are points of symmetry, that their other paths are blue, or that arithmetic classes are distinct unmarked tournament classes.

The same construction works abstractly for any elementary2 class group of rank r with chosen character basis: embed its bits in reflection-invariant graph coordinates. If complement is also to be retained within the image, choose a nonzero class eta and send it to the all-one vector. A sufficiently large reflection-fixed space permits this. Rank0 has no such nonzero class and requires an enlarged carrier. None of these choices is canonical merely from the count of classes.

## 5. Compatible mod15/mod35 registers are not a genus channel

The actual odd sources63 and273 have exactly the same compatible registers:

 63 mod15=273 mod15=3,
 63 mod35=273 mod35=28.

Both are primitive B-chain sources in the inherited sense: neither is5 modulo8. Yet63 is not idoneal and273 is idoneal, by the complete form calculations above. Thus even the joint mod105 register does not determine the idoneal predicate on this native source family. It also does not determine a primitive transport step.

For the unstopped primitive map P, use a=2 on s=1 mod8 and a=1 otherwise, then repeatedly remove B(m)=4m+1 ancestors from values5 mod8, adding2 to the braid clock per removal. Direct calculation gives

 P(63)=95 with braid clock1,
 P(273)=51 with braid clock4, since(3*273+1)/4=205=4*51+1.

The successor residues modulo105 are95 and51, so no deterministic map on the joint register alone realizes P. These are single unstopped primitive transitions, NOT first-arrival clocks to a specified frontier. If a source is already in that frontier its stopped first-arrival clock is zero instead. The fully odd map U has different data: U(273)=205 with valuation2. These clocks must not be exchanged.

No claim is made that the two sources have equal charges, the same terminal, or interchangeable legal children. Their failure to share those data is the point: a CRT register is a guard, not the source itself.

## 6. Idoneality is not a Collatz invariant or a two-state factor

For the actual odd map U(n)=(3n+1)/2^v2(3n+1),

 U(3)=5, U(7)=11, U(17)=13.

The small reduced-form check certifies3,5,7,13 idoneal and11,17 non-idoneal. Thus idoneality can be lost and gained. More sharply,3 and7 share the same idoneal flag but their successors do not. Consequently no deterministic function of the idoneal flag alone intertwines U. No ROOT issue is involved:3 and7 are both nonROOT sources.

This rules out the proposed simplest information channel, not every possible arithmetic connection. The positive class-group/graph code in section4 preserves precisely stated algebraic operations; it supplies neither a decreasing Collatz potential nor source-authenticated annulus cancellation. The65-number completeness input does not repair either missing implication.

## 7. Verification and dependencies

Run python3 verify_idoneal_channels.py and python3 -O verify_idoneal_channels.py, comparing receipts. The script independently enumerates primitive reduced forms by the exact bound a<=sqrt(|D|/3), validates every displayed form and representative, checks all three-character codes under many coprime representations, checks the graph/tournament operation table, and replays the literal primitive and odd steps. It also recovers the65 idoneal candidates up to1848 by the self-inverse-form test. This last bounded computation is not a global completeness proof.

The only non-elementary structural inputs are classical reduction, inverse-class identification, and genus-character theory. Their roles are explicit. The analytic adopted #003 input is unnecessary for every new map and hostile in this note.
