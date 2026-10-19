# An integer Pell half-step joining Gold equality parameters to the eleven-point Collatz carrier

**2026-10-19. PROVED explicit affine and finite conjugacies.** This extends the [symbolic Gold family](gold_two_chain_fork_20261019.md) and [Gold-budget/Golay bridge](gold_budget_golden_golay_20261019.md). It does not transfer a general conjecture across these carriers.

## 1. Equality parameters have an integral square root of their recurrence

The curved Gold boundary is

    Q(m,n)=n²-mn-(m-1)(m-2)=0.

The structural theorem classifies all solutions with m,n≥2 and m≤n as the orbit of (2,2) under

    R(m,n)=(2m+3n-3, 3m+5n-6).

Define the integral affine map

    S(m,n)=(3-m-n, 3-m-2n).                              (1)

Direct expansion proves Q(S(m,n))=Q(m,n) and S²=R. Thus there is an exact one-step map on the full integer equality curve whose square is the positive-dimension recurrence.

The positive admissibility chamber is not preserved by S. For m,n≥2, its image has m'≤-1 and n'≤-3; S² returns to the admissible boundary. In particular the S orbit begins

    (2,2), (-1,-3), (7,10), (-14,-24), (41,65), ...

The negative pairs are formal integer solutions of Q=0, not dimensions of actual posets. Two steps preserve the structural family; one step does not.

## 2. Norm coordinates make the half-step transparent

Put

    X=5m-6,     Y=2n-m.

Then

    X²-5Y²=-20Q(m,n)-4.

On the equality curve, X²-5Y²=-4. In Z[φ], φ=(1+√5)/2, the maps have the simple form

    S: X+Y√5 → -φ²(X+Y√5),
    R: X+Y√5 →  φ⁴(X+Y√5).                            (2)

Indeed S sends (X,Y) to (-(3X+5Y)/2, -(X+3Y)/2). These halves are integral on the parameter lattice, as the original integer formulas (1) show. The corresponding affine fixed point is (6/5,3/5).

The Gold-budget obligation g(b)=(1-2b)/(1-b) has fixed points α=(3-√5)/2 and β=(3+√5)/2. Its cross-ratio z=(b-α)/(b-β) obeys z(gb)=-βz=-φ²z. Thus (2) realizes the same one-step multiplier directly on the Pell norm coordinate. This is an equality of specified algebraic actions, not an assertion that successive conditional comparison probabilities equal g iterates.

## 3. Modulo eleven, the Pell curve is exactly the nonzero golden carrier

Work over F_11 and choose √5=4. Define

    u=X+4Y=m+8n-6,
    v=X-4Y=9m-8n-6.

On Q=0,

    uv=-4=7 mod11.                                    (3)

Hence u is nonzero. Conversely, every u∈F_11* determines v=7/u and then uniquely

    X=(u+v)/2, Y=(u-v)/8,
    m=(X+6)/5, n=(Y+m)/2.                             (4)

All denominators here are units modulo eleven. Thus Q=0 over F_11 has exactly ten points, in explicit bijection with F_11*. Equations (1)–(2) become

    u(S(m,n))=2u(m,n),
    u(R(m,n))=4u(m,n).                                (5)

The S action is one ten-cycle; the R action is two five-cycles. The admissible positive integer equality orbit starts at u(2,2)=1 and reduces to precisely

    1,4,5,9,3,

the quadratic-residue subgroup H. Its other five modular states do not lift to any integer equality pair with m,n≥2 and m≤n, by the complete positive-boundary classification. They do occur on the alternating negative formal branch of the S orbit. This is an explicit local/global admissibility distinction, not an assumption that every modular state is realizable by positive dimensions.

For the projective budget coordinate, put

    b=(9u-5)/(u-1),

with b=infinity when u=1. This maps the ten Pell states to P¹(F_11) minus the two fixed budget points 5 and9. By the already proved cross-ratio identity, it intertwines S with g and R with g². The positive equality sequence follows the five-cycle

    infinity → 3 → 10 → 4 → 0 → infinity

under g².

## 4. Exact Collatz dictionary and the signed syndrome step

Let κ be the inherited standard-Collatz labeling in [THM-4592](../../01-canon/theorems/THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md), restricted to its ten nonzero period-five rational points. Its labels are

    1/13,2/13,4/13,8/13,16/13 ↔ 3,9,5,4,1,
    -5,-14,-7,-20,-10       ↔ 8,10,7,6,2.

The standard map C uses 3x+1 on odd rational 2-adic integers and x/2 on even ones. It multiplies κ by4. Consequently the map

    Ψ(m,n)=κ⁻¹(m+8n-6 mod11)

satisfies Ψ(R(m,n))=C(Ψ(m,n)) on the full modular equality curve. In particular the actual positive-dimensional Pell boundary maps cyclically to

    16/13 → 8/13 → 4/13 → 2/13 → 1/13 → 16/13.

Let D be the permutation of this ten-point carrier defined by κ(Dq)=2κ(q). Then D²=C, and Ψ(S(m,n))=D(Ψ(m,n)). Its orbit interleaves the two Collatz five-cycles:

    16/13, -10, 8/13, -5, 4/13, -14, 2/13, -7, 1/13, -20.

D is a square-root permutation on this finite carrier, not a new integer Collatz branch formula. In particular D is not numerical negation.

Extending multiplication by2 to the eleven coordinate labels by fixing0 gives exactly the signed Golay coordinate action from the earlier packet: 2=-3² mod11. On reciprocal syndrome modules it is

    (z_+,z_-) → (z_-^9,z_+^9),

and on the canonical single-syndrome decoder it is G_dec(z)=J_dec(z^9), with G_dec²(z)=z^81. Thus the same half-step has three explicit carriers: formal Pell equality parameters, projective budget residues, and a signed permutation representation on Golay syndromes.

The ten-point Pell curve carries the multiplicative cyclic subgroup. Translations u→u+1 do not preserve F_11*: they can hit0. Therefore the full affine group of order55 or110 is not being asserted to act internally on this nonzero-norm conic. Adding the fixed0 coordinate enlarges the carrier and is essential when using the full affine representation.

## 5. What this connection can and cannot do

The occurrence of eleven is now tied to the exact equality-family recurrence, not merely to a common count of cells. Positive Gold dimensions select one of the two modular Collatz cycles, while an integral square root of the recurrence alternates between positive and negative formal parameter chambers. The field-unit action explains the shared period structure.

The Gold inequalities themselves are established by the real conditional-shuffle and monotone-fibre proofs. Their modular images cannot decide those inequalities. The Collatz carrier here is a fixed finite rational periodic set; it is not all positive integers. The Golay module is a representation of label permutations, not a map from extension counts to a code correcting a failed proof.

A useful next target is to identify which actual comparison-strategy transformations, if any, realize these unit actions while preserving positive dimensions and exact count inequalities. The two-step equality-family recurrence is one such structural parameter map; the one-step half-map explicitly exhibits the positivity obstruction.

Reproduce the targeted identities with `python3 05-knowledge/experiments/gold_pell_collatz_halfstep_20261019.py`. The companion JSON records the exact ten-state dictionary and initial integer orbit. No prior broad census is rerun.
