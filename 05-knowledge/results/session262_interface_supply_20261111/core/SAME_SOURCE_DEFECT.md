# Two exact defect branches at the same labelled source graph and base cost

**PROVED**, assuming the rooted response kernel proved in [PROOF.md](PROOF.md).
This is a constructive all-n family with its initial case, not density or a
claim of arbitrary increments. Chan–Pak remains OPEN.

## 1. The universal budget becomes a fixed-source comparison

For the bare apex z in our Q[A,B] construction, A_z={0}. The cut graph G_z
has one edge from0 to each out-neighbour of0 in the old source. The quotient
row Q₀=198 points to the three singleton vertices1,2,6 and every vertex
of the module at port7. Consequently, with b=|B|,

    G_z=K_(1,b+3) plus isolated vertices,
    d(G_z)=1+2^(b+3).                                  (1)

A star dominating set either contains its centre and an arbitrary leaf
subset, or consists of all leaves without the centre. Isolated vertices
are forced. This proves(1) without counting any old path family.

The inherited inequality says that the source's exact born increment
J=U₀−E₀ must pay 2^(b+3), leaving quantitative slack

    σ=U₀−E₀−2^(b+3).                                  (2)

For positive tail length t, root0 also beats the t tail vertices, so

    d(G_z)=1+2^(b+3+t),
    σ=2^t(U₀−2^(b+3)),          t≥1.                  (3)

The original base count remains K in both cases. The zero-tail endpoint
correction is visible in the exact defect itself, not merely in notation.

## 2. Equal-order modules freeze both G_z and K

For each n≥3 choose

    A_n=TT_n,
    B_n=C₃ ⇒ TT_(n−3).

The second module is a genuine triangle at n=3; an empty ordinal suffix
has its usual neutral role and does not invoke an absent Q port. Both Q
ports remain nonempty. Exchanging A_n,B_n preserves K by the inherited
symmetric kernel. Because their orders are equal, the vertex intervals
and all out-neighbours of0 are also unchanged. Thus the TWO source graphs
G_z are literally the same labelled graph, not just graphs with equal d.

Let x=2^n,z=3^n. Ordered r-slot covers allowing empty slots, denoted Q_r,
satisfy Q_r(X⇒Y)=Q_r(X)Q_r(Y), by independently concatenating each slot.
Their elementary triangular relations give

    A_n: H=1, C₂=x−2, C₃=z−3x+3;
    B_n: H=3, C₂=3x/2−6, C₃=11z/9−9x/2+9.

Here Q₂(C₃)=12 and Q₃(C₃)=33. The rational-looking coefficients are
integers at every n≥3. Inserting these physical profiles gives the common
base count

    K_n=111x²/2−9x+27.                                (4)

The two positive-tail capacities are

    U_AB=15x²+682xz/9+45x/2−44z/9−24,
    U_BA=15x²+93xz+63x−12z−24.                        (5)

Their endpoint corrections are E_AB=21x and E_BA=21x/2. Hence

    U_BA−U_AB=(155x−64)z/9+81x/2 >0,
    J_BA−J_AB=(155x−64)z/9+51x >0.                   (6)

This is an exact, positive and explicit choice of response at every n,
with the same source domination graph and the same old Hamilton count.
No change to the source inequality's scalar input (K,d) detects the choice.
The rooted response coordinates do detect it and determine every output.

## 3. Actual initial case and uniform operations

At n=3 the source order is12, and

    K=3507, d=65,
    (U_AB,E_AB,J_AB)=(17352,168,17184),
    (U_BA,E_BA,J_BA)=(21204,84,21120).

The bare outputs therefore have counts20691 and24627, with slack17120
and21056 respectively. At positive tail length1 the outputs are38211 and
45915, and all later tails are K+2^t U with their specified capacity.

Increasing n adds one transitive vertex to A and one unselected sink to
the suffix of B. Swapping the complete modules is an involution. The
class and its equal-order guard survive these source operations for every
n≥3. Adding further tail vertices reuses the chosen source. A single
source graph is not asserted to regenerate after its output is fed into
another unrelated quotient.

The source quotient is strong even though B is decomposable. Its lift and
the apex/tail outputs remain strong by the argument in PROOF.md. Thus this
fixed-source contrast is not obtained merely by putting triangle shells
outside the whole source.

## 4. Constructive significance and exact limit

The universal inequality makes the compulsory invoice explicit: K is the
old count and d−1 is the mandatory part of the new count. This family keeps
BOTH invoices fixed while changing the exact excess through a real
count-preserving source exchange. Thus the constructive task naturally
lives over a *fibre* of (old count, domination graph), with rooted response
as its controllable coordinate. The inequality neither forbids this
movement nor selects a desired point in the fibre.

What is now supplied uniformly: two explicit, reusable response branches
in each n-fibre, physical transitions between fibres, a target formula,
and the initial n=3 case. What is still missing: a sufficiently rich set
of branches to hit arbitrary target counts, or a return that turns this
growing branch gap into a fixed binary digit. The rapidly growing gap in
(6) is retained as a limitation rather than described as interval coverage.

## Reproduce

    python3 -B core/same_source_defect.py

Four literal controls at n=3,4 and t=0,1 compare both orientations with
(4)–(6), check the labelled cut graph and base count, and record all outputs
in `same_source_receipt.json`. The all-n statement is the displayed profile
calculation, not extrapolation from these controls.
