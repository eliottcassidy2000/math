# Session44 synthesis: actual guards, shared ports, and the coordinates that survive

Status: proof and exact-certificate synthesis. Publication verification is recorded separately; this document does not turn a selected observer packet into a complete native wiring theorem.

## 1. The principal result is a distinction between guard strengths

The original electrical phase is not an arbitrary clock residue. On the actual converse-merged tournament metagraph it is

    phi_x(T,P) = x_[T] - x_[K_P(T)] mod8,
    Lx=0 mod8, x_root=0.

The complement partner K_P keeps the actual Hamiltonian path arcs and reverses the other arcs. The original node guard uses owned sums of these phases with actual multiplicities.

Three requirements must remain separate:

- preserving one chosen character in each port;
- preserving every harmonic character in each port;
- preserving the aggregate node currents for every character.

The second implies the first. The third also retains WHERE each current lands, so there is no general implication from universal scalar phase equality to preservation of the full node-current vector. For the fixed-source order-eight domain below, target injectivity supplies that missing implication. Our explicit packets retain these scopes. Equality of three-bit outputs is not equality of target ownership or of native realizability.

A same-modulus hostile is the weighted four-cycle with conductances01=03=1 and12=23=2. Its grounded mod8 harmonic vectors are0000,0040,0424,0464. Ports0→1 and0→3 have equal phase for every such vector, but the character0424 gives current vectors(4,4,0,0) and(4,0,0,4). A scalar sum at the original source is not the whole network current.

## 2. An exact all-size electrical quotient, with an unusually strong order-eight realization

For any integer modulus m and zero-sum charge d,

    d.x=0 for every grounded harmonic x
       iff d belongs to im(L mod m).

The correct additive charge quotient is K/mK, where K is the graph's critical group. An owner-marked Smith transformation, not just the invariant factors, gives its coordinates. Extra boundary constraints change the annihilator and must be included explicitly.

At all moduli simultaneously, two vertex charges agree precisely when the vertices are joined by a path of literal unit bridges. An elementary integral-flow proof is included. This is a standard bridge-contraction phenomenon specialized and proved here, not a priority claim.

The complete native order-eight atlas is especially rigid. Checked harmonic characters modulo8 separate all3528 merged nodes except

    {0,3527}, {53,3492}, {109,3523}.

Each exceptional pair has an integer unit-leaf image witness, so no modulus can separate it. These are exactly the three literal bridges. Crucially, no pair can occur as two different targets at the same native source: a leaf can only be reached from its parent, and that parent has no self-loop.

Therefore, on ALL2^21 normalized native presentations, fixing the source and every original mod8 phase is equivalent to fixing the merged target owner. Independent reconstruction checks all438216 exported pair rows and all3528 source neighborhoods.

For a candidate perturbation, retaining every original harmonic vector means kernel inclusion, not necessarily equality of old and new harmonic banks. This does not recover a literal Hamiltonian path or complete wiring. The universal charge module has3369 binary bits, whereas its vertex evaluations have only3525 equivalence classes. A twelve-bit label for those vertex classes is not an additive replacement for the charge module.

See [electrical quotient](electrical-duality/UNIVERSAL_ELECTRICAL_GUARDS.md) and [independent audit](audit-electrical/AUDIT.md).

## 3. Actual coordinate dependence becomes smaller when observations are adaptive

The original four switch sources have409,407,523,527 actual Hamiltonian paths. All four require six fixed indexed owner clocks to determine the original chosen phase. But the exact optimal worst-case adaptive query counts are:

| Source | Chosen phase | Full target / universal phase bank | Literal path |
|---:|---:|---:|---:|
|1238|5|5|5|
|2242|4|4|4|
|2246|4|4|4|
|2315|4|5|5|

The policies are explicit, with independent lower-bound recursion, not only successful decision trees. Uniform-path expected cost is a separate objective; a mean-optimal policy can use five queries even when the minimax optimum is four. The source2315 row shows that compressing the output really can help adaptively even when its fixed-coordinate footprint stays unchanged.

With the actual first and last owners already supplied, all210 nonempty endpoint fibres in these four sources need at most three additional adaptive queries for phase, target or path. This is an exact statement about those fibres, not an all-order query bound. Queries can branch on up to eight owner labels; they are not automatically ternary trees.

See [adaptive observers](actual-guards/ADAPTIVE_OBSERVERS.md), [endpoint supplement](actual-guards/ENDPOINT_OBSERVERS.md), and [audit](audit-adaptive/AUDIT.md).

## 4. Genuine coupled moves survive the strongest per-slot guards

A selected two-slot packet at source1238 has two compatible states. Each slot keeps its actual endpoints and target class; both slots keep fixed shared arcs(6,7) and(5,4). Two free shared owner variables exchange0 and1 at clocks1 and4 in both slots together. A one-slot change violates the shared-interior equality, but the simultaneous two-slot move is valid.

The four literal paths and complement partners are distinct. Direct target isomorphisms prove preservation of EVERY target-dependent electrical character in each slot, at every modulus. This is not merely a common rotation of a source automorphism orbit; that source has trivial automorphism group.

There is also a genuine necessary three-slot assignment cycle at source2242 and target3442. Actual indexed-clock pins produce eligible sets{p0,p1},{p1,p2},{p0,p2} within the complete three-path target fibre. Exactly two injective assignments exist. They differ in all three slots, so pairwise assignment moves cannot connect them. All fixed shared arcs and per-slot universal guards survive. This is a ternary cycle, not evidence of a recursive tree.

The general fixed-pin mechanism is alternating matching: an unused literal port allows vacancy propagation, while fully occupied alternating cycles require simultaneous reassignment. For dynamically shared owner values, an exact one-hot integer compiler retains every occurrence, integer congruence lift and port-capacity row. Conformally primitive kernel moves connect its bounded fibre, but no small move bank or efficient all-size algorithm is promised.

Global caveat: newly selected ports must be unoccupied in the unchanged remainder, or an explicit global injective reassignment is necessary. The full canonical atlas already contains its paths. These certificates do not alter its complement involution or invent new metagraph edges.

The positive move is not confined to order eight. Its complement-target isomorphism fixes both core endpoints. Prepending and/or appending arbitrary transitive source/sink chains preserves that isomorphism, fixes the padded owners, and preserves the shared arcs at their shifted clocks. This gives an actual sharp two-slot guarded move at EVERY order n>=8, still outside source automorphism orbits. Nonempty padding makes this a decomposable family, not a theorem about regular or primitive tournaments. The proof does not assume that numeric phase values5 and4 extend to higher-order networks, or that whole padded target fibres retain their original sizes.

See [coupled packets](coupled-moves/THEOREM.md), [independent native audit](audit-coupled/AUDIT.md), and [all-size stable padding](stable-padding/THEOREM.md).

## 5. Shared kernel lifts have a precise gluing obstruction

For authenticated local lift cosets r_i+K_i and the homomorphism D recording all actual shared-port mismatches, compatible lifts exist exactly when

    delta-D(r) belongs to D(product K_i).

When nonempty, the global lift set is a torsor for the kernel of that restricted map. Its size need not equal the product of the local lift counts: the exact count is the restricted kernel size.

The actual C9-cell control has81 independent local kernel choices but only27 compatible joint lifts at a fixed quotient coordinate. Equal residues modulo3 can still hide incompatible owner values0 and3 modulo9. The missing data are literal alignment coordinates, not a missing numerical count.

This is an all-size algebraic criterion CONDITIONAL on authenticated finite-abelian charts and the actual mismatch map. It does not declare arbitrary native attachment relations to be cosets.

## 6. Connected quotients can conceal disconnected native lifts

For whole accepted ternary torsor fibres and authenticated invertible transports, choose a spanning-tree gauge. Components of the supplied lifted graph are precisely the orbits of the based affine cycle transports and internal kernel moves. In the translation case, rank d in an r-dimensional fibre gives3^(r-d) components. Fixed sections and connected components answer different questions.

An actual padded two-C9-cell tournament of order20 supplies27 Hamiltonian observers with fixed endpoint owners. Common cell rotation projects to a connected ternary quotient, but has three native nine-state components. The cycle wrap translates the retained fibre by(1,1), leaving a missing kernel direction. The actual one-cell rotation by3 supplies that direction and joins the components. Every move preserves the complement-target class by an explicit isomorphism, so no invented order20 electrical potential is needed.

With specified probabilities the quotient has gap8/9, while each native C9 component has gap16/27 times(1-cos(2pi/9)); the global unrepaired gap is zero. Another exact connected PSD family keeps quotient gap8/9 while its native gap tends to zero. Connectivity alone does not bound the transition weights or mixing clock.

These are supplied-move statements. Partial relations, unequal fibres, missing transports and additional native edges must not be silently completed into this cover interface.

See [holonomy compiler](lift-holonomy/README.md) and [independent audits](audit-holonomy/AUDIT.md).

## 7. What really lies behind 7/15,8/15 and the ternary fractions

Incoming work already explains the order-four pair through the same raw conductance4/5 and different retained masses(3,4) versus(2,6). This session proves exactly when a source compression discards only zero modes, recovering the safe order-four compression rather than assuming it from the shared fraction.

At order five, the full and compressed nontransitive coupling gaps have algebraic degree seven. The complement-paired compressed gap has degree five. Independent characteristic polynomials, finite-field irreducibility and exact root isolation certify these facts. Thus the rational Berggren gap-coordinate tree cannot encode general native gaps.

A different, genuinely native ternary process is the cyclic-cell continuation transfer I+Q. Conditional first-step probabilities have Jacobsthal denominators

    J_r=(2^r-(-1)^r)/3.

For remaining length r>=3, apart from1/2, their unordered pair is((J_r-1)/(2J_r),(J_r+1)/(2J_r)). Neither7/15 nor8/15 occurs. They are individually Farey-adjacent to1/2, not to each other. The consecutive Jacobsthal denominator ratios J_(r+1)/J_r obey a generalized continued fraction whose transfer determinant is minus two rather than a unit.

Their exact Berggren address exists, but its depth grows exponentially in the actual remaining path length; a ternary heap address grows doubly exponentially. An address is not a native execution clock. The polynomial transfer I+zQ retains the integer displacement that a phase residue forgets.

See [spectra and continuation fractions](fractions-ternary/README.md) and [audit](audit-fractions/README.md).

## 8. The meaningful mod-nine alternative is a different representation, not an error

At modulus8, ternary Fourier inversion works in a split root extension because three is a unit. At modulus9, the same averaging loses information. A tensor Newton basis with an integral determinant-one evaluation matrix gives lossless interpolation over ANY coefficient ring, with no modular division.

For one cyclic shift T, D=T-I satisfies D^3=-3TD. Over Z/(3^k), its exact regular-representation nilpotency index is2k+1. For integers k,r>=1, on the faithful regular table representation of r commuting ternary coordinates, the entire augmentation ideal has exact index

    2(k+r)-1.

Arbitrary modules inherit the upper bound, but need not attain it. In this faithful regular representation, the last nonzero ideal has precisely three elements, generated by3^(k-1) times the sum of all translations. It records total table mass modulo3. This is a coupled coefficient phenomenon, sharper than adding separate nilpotency bounds.

The older binary deck law and the new ternary law are genuinely related through finite-difference representations. Neither converts the original electrical phase to a different modulus. Newton coefficient indices are not Fourier frequencies; replacing one support-rank rule by the other fails on an explicit rank-one table.

See [Newton compiler](ternary-newton/THEOREM.md), [multivariable theorem](ternary-newton/MULTIVARIABLE_FILTRATION.md), and [independent review](audit-newton/MULTIVARIABLE_REVIEW.md).

## 9. A missing harmonic coordinate already forbids two native rectangle menus

The two incoming black rectangles use sources1238 and2246 and target pairs(881,2762) and(2650,2843). Each passes the chosen primitive character modulo8, and even modulo16 or64 respectively. That did not mean it passed the whole mod8 bank.

For four distinct endpoints and differences u,v, a rectangle has perturbation c(uv^t+vu^t). It preserves every original mod8 harmonic vector iff cu and cv both lie in im L. Here all three relevant charge differences have EXACT order four: authenticated mod4 characters scaled by two separate2d, while explicit full vectors verify Lz=4d mod8. No claimed completeness of a proposed character basis is needed.

A paired black move has c=2a, so the exact electrical condition is that amplitude a be even. But the full actual source1238 path fibres to2762 and2843 each contain only one path. Distinct-port capacity limits each menu to|a|<=1. Therefore neither menu admits a nonzero universally mod8-safe amplitude. Increasing aggregate black-edge weight does not manufacture a second literal path.

There is also an all-size necessary coordinate rule. For a fixed source difference u, disjointly supported zero-sum net target vector v, and coefficient2^b, a genuine harmonic source gap2^s times an odd unit modulo2^k forces every v_j to be divisible by2^max(0,k-b-s). In this case all net target amplitudes must be even. With unit net capacities, only incidence-neutral trades remain on that source pair. Other source pairs, larger capacities and cancellations yielding zero net incidence are not excluded.

Even target amplitudes alone are not sufficient: the source-row obligations remain. The supplement gives the exact general criterion and a counterexample to dropping those rows. Its preservation claim is inclusion of the original harmonic bank in the new one, not equality of both banks.

See [rectangle bank and capacities](rectangle-bank/THEOREM.md) and [independent direct-edge audit](audit-rectangles/AUDIT_RECTANGLE_BANK.md).

## 10. The remaining all-size obligation

The most useful unification is not that every object is secretly the same ternary tree. It is that quotient compatibility needs its missing sidecar:

- a target charge needs an owner-marked coordinate;
- a guard observation needs actual input dependence;
- a local attachment needs its common kernel lift;
- a connected base needs its cycle action on those lifts;
- an assignment needs distinct available literal ports;
- a spectral coordinate needs stationary mass and the actual clock.

An all-size native wiring theorem still needs a complete authenticated move inventory across different observer orbits, actual guarded eligibility and global injective port allocation, plus quantitative control if a mixing statement is intended. The new criteria and native witnesses close concrete subproblems; they do not replace these remaining obligations.
