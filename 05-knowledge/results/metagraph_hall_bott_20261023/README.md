# Marked harmonic charges, linking units, and the useful Hall/Bott boundary

2026-10-23. **PROVED:** the elementary all-size charge/phase extraction and
the stated algebraic obstructions below. **FINITE-EXACT:** the order-eight
charge, its 105/256 self-link, and the fixed-phase rejection of the declared
four-terminal port switches. **CITED:** graph monodromy pairing and the
specified physics/Clifford facts. **OPEN:** a complete bounded-size port
signature or all-order tournament-wiring classification.

The practical result is a cheap new guard on actual ports. A previously
authenticated harmonic lift assigns a residue to each merged node. Every
line inherits the difference of its endpoint residues. The candidate wiring
must conserve these weighted differences at every node. For the surviving
order-eight four-terminal switches, just three bits of the native assignment
reject every feasible nonzero switch without recomputing a Smith form.

## Inheritance and concept board

The closest proved mechanism is the [higher dyadic filtration](../metagraph_order8_electrical_20261023/README.md):
its last primitive lift survives to modulus 256. The canonical hostile is
the [same-port-capacity wiring counterfeit](../metagraph_order8_frontier_20261023.md),
which preserves local margins while changing global electrical data. The
corrected near miss is treating an eight-dimensional carrier as Bott
periodicity, or a group order as a complete linking invariant. The least-used
sidecar is the *odd unit* of a linking form, visible in the repository's
older knot/Goeritz calculations and factorial-block unit register.

The board is: actual port incidence; harmonic phase; finite linking form;
characteristic refinement; marked stabilization; graded matrix words.
The anchor is a realizability guard, the niche is a phase finer than Smith
multiplicities, and the wildcard is the QCD/root-lattice analogy with its
precise loss of word ownership.

Recovered mechanisms, with their current scope:

| Earlier lane | Exact surviving mechanism | Transfer here |
|---|---|---|
| [Knot decider / Goeritz pairing](../UNKNOT1-DECIDER-STATUS-AND-DESIGN-klein-S691.md), corrected by MISTAKE-318 | A cyclic first homology group need not represent the required linking norm; a sufficient local move is not a complete decision procedure | Keep the pairing unit after retaining Smith factors; successful port rejection is not completeness |
| [Marked Bott packet](../collatz_bott_marked_periodicity_20261007d.md) | Stabilized matrices retain operations only with a marked corner; a stable class loses carries | Keep actual node/port charges when passing to a lattice invariant |
| [Fano/Clifford carrier](../creation_fano_20260925.md) | Primitive marked vectors preserve divisibility; a norm or eight labels do not | Primitive harmonic potentials certify exact torsion order |
| [THM-3191, factorial-block exterior Clifford law](../../../01-canon/theorems/THM-3191-factorial-block-exterior-clifford-law-and-global-carry-smith-profile.md) | Valuation thickness and prime-free unit are distinct coordinates | The dyadic rank tower gives thickness; the odd linking unit is an additional register |
| [THM-4762, cat-map modular data](../../../01-canon/theorems/THM-4762-the-arnold-cat-map-is-the-commutator-of-t-and-s-pisano-bounds-in-modular-data.md), corrected by MISTAKE-665 | The modular commutator depends on discriminant-form data, not only the group or lattice signature | A possible downstream invariant only after the quadratic/characteristic sidecar is supplied |
| [THM-868, score-lattice membership](../../../01-canon/theorems/THM-868-e8-bridge-score-lattice.md) | At order eight the score-deviation half-vector lies in the stated E8 coset | This is a score embedding, not a derivation of the new 256 torsion or of Bott recurrence |

These are independent uses of markings, not identifications of their native
categories. The current [Bott packet](../collatz_bott_marked_periodicity_20261007d.md)
explicitly separates the historical informal identification of several
eights. MISTAKE-484 distinguishes the theta series of Z8 and E8. We retain
the score-membership statement and do not inherit the old numerical analogy.

## 1. The all-size inverse-free charge rule

Let A be any nonsingular symmetric integral matrix, and put

    K(A) = Z^r / A Z^r,
    b_A([c],[d]) = c^T A^(-1) d modulo Z.

This is a well-defined perfect symmetric pairing. Replacing c by c+Az
changes it by z^T d, an integer. If c pairs integrally with every d, then
A^(-1)c is integral, so [c]=0; nondegeneracy implies perfection for finite
groups. For a root-reduced graph Laplacian this is the standard graph
monodromy pairing; see [Shokrieh, Theorem 3.4 and Proposition 3.7](https://arxiv.org/html/0907.4764).
We claim a new application/certificate here, not novelty of that pairing.

For every positive m there is an explicit isomorphism

    ker(A mod m) -> K(A)[m],
             x -> [A x / m].                              (1)

The input is x modulo m. Changing its representative by mz changes the
output by Az. The map is injective since Ax/m=Az implies x=mz. If m[c]=0,
then mc=Ax for some integer x, proving surjectivity. If x is primitive
modulo a prime p and m=p^j, its image has exact order p^j: k times the
image is zero precisely when kx/m is integral.

Write y=Ax/m. The pairing can be read without an inverse:

    b_A([c],[y]) = c^T x/m modulo Z,
    b_A([y],[y]) = x^T A x/m^2 modulo Z.                    (2)

For a graph, keep the root at potential zero. A directed port from i to j
has charge e_i-e_j, so its phase is (x_i-x_j)/m modulo Z. Reversing the
port negates the phase. This is an exact map from an actual endpoint-marked
line to a character value, not an arbitrary coloring of the lines.

If u=x^T A x/m is a unit modulo m, the cyclic subgroup generated by [y]
is nondegenerate and splits orthogonally. Indeed for any [c], solve
ku = m b_A([c],[y]) modulo m. Then [c]-k[y] is orthogonal to [y]. The
intersection is zero because u is a unit. This retains the *unit* of the
cyclic linking form as well as its order.

The same construction works for several lifted potentials X: the charge
columns AX/m have Gram matrix X^T A X/m^2 modulo Z. If the integer matrix
X^T A X/m is nonsingular modulo p, their p-primary block is nondegenerate.
Thus the natural extension of the dyadic rank tower is a matrix of unit
residues with its charge markings. Merely storing ranks loses these units.

## 2. The actual order-eight charge has self-link 105/256

The frozen native [electrical receipt](../metagraph_order8_electrical_20261023/receipt.json)
contains a primitive potential x modulo 256 on all 3,528 merged nodes,
rooted at node 3,527. The new consumer reads that potential and the actual
full edge export independently of the original lifting engine. It checks

    A x = 0 modulo 256,
    x^T A x = sum_edges w_ij (x_i-x_j)^2 = 6,819,178,752.

Consequently y=Ax/256 is an exact order-256 charge and

    u = x^T A x/256 = 26,637,417,
    b_A(y,y) = 105/256 modulo Z.                           (3)

This is an exact new occurrence of 105 with an identified carrier. It is
not an intrinsic preference for that numerator: multiplying y by 27 gives

    105 * 27^2 = 1 modulo 256,

so the same cyclic summand has a generator with self-link 1/256. The four
normalizing multipliers modulo 256 are 27, 101, 155, 229. The invariant
unit squareclass is 1 modulo 8, whereas the specific 105 depends on the
chosen lift/generator. The rank tower had retained only the denominator.

This squareclass is well-defined for the unique highest cyclic factor of
the present group. After splitting off the nondegenerate order-256 factor,
the complementary 2-primary group has exponent four. Replacing its generator
by ag+h, with a odd and h in that complement, changes the numerator modulo
256 to a^2 u plus a multiple of 64. Since odd squares are one modulo eight,
the unit modulo eight is unchanged. This argument uses the large exponent
gap; it is not asserted for arbitrary repeated highest factors.

## 3. A fixed character rejects the remaining four-terminal switches

The concurrent electrical packet studies the balanced update

    A_t = A + 2t (u v^T + v u^T),
    u=e_1238-e_2242, v=e_2246-e_2315.

Its source/port capacities permit t=-1,0,1,2. We use the inherited exact
potential, whose entries at those four vertices are

    (97,18,65,212).

Thus u^T x=79 and v^T x=-147, and the Kirchhoff defect is

    (A_t-A)x = t*(-294,294,158,-158)

on those four vertices, zero elsewhere. Modulo 256 the t=1 defect is

    (218,38,158,98),

and every allowed nonzero t still gives a nonzero defect. The same fixed
marked character therefore rejects all three alternatives. The concurrent
electrical lane sharpened this to an exact three-bit guard: both 79 and
-147 are odd, so the defect vanishes modulo eight exactly when 4 divides t.
Among t=-1,0,1,2 only t=0 passes. Modulus eight is the smallest binary
modulus that rejects all three: the t=2 matrix perturbation is divisible
by four, so no harmonic-potential test modulo two or four can distinguish
it from the native matrix. This minimality is for this specified scalar
switch family and class of modular harmonic tests, not all port ambiguities.

The four retained residues modulo eight are (1,2,1,4), and the t=1 defect
is (2,6,6,2). This tests four differences after the native certificate has
been supplied, rather than repeating the global dyadic calculation for
each candidate. The 256-level lift supplies the character; its three-bit
reduction is enough for this adversary.

**Scope:** the native character is a known necessary invariant of this
fixed target network. It is not derived from the coarse capacities alone.
A different network may legitimately possess a different character. This
is a target-specific extra guard, not a standalone characterization of
tournament wiring. To make the rule uniform in order, one must derive or
transport suitable marked potentials from the tournament/port structure.

The all-size rule being tested is simply

    for each node i: sum_j w_ij (x_i-x_j) = 0 modulo m.       (4)

Equivalently, charge phases satisfy conservation at every node. Retaining
which port contributes which difference is essential: an unowned histogram
cannot impose (4). The saved phase histogram uses the export's increasing
endpoint gauge; it is not claimed invariant under arbitrary relabeling.

## 4. A minimal hostile: the complete rank tower still misses a unit

Consider two genuine three-vertex integer-conductance graphs with reduced
Laplacians

    A_1 = [[9,-1],[-1,1]],    A_3 = [[3,-1],[-1,3]].

The first has eight parallel root-to-vertex edges and one leaf edge. The
second is a triangle of conductances 2,2,1. Both have Smith factors (1,8),
and their complete dyadic dimension sequences are (1,1,1,0). Their nonzero
one-dimensional terminal forms over F2 agree as well. Yet the first marked
generator has self-link 1/8, the second 3/8. These cyclic linking forms are
not isometric: every odd square is one modulo eight. This is a structural
hostile to replacing linking data by the entire dyadic rank sequence.

It is not a same-margin tournament counterfeit; the concurrent actual-port
counterfeit is the separate strong control. Conversely even the full marked
finite pairing does not reconstruct a graph: every unit-weight tree has
trivial critical group. The two three-vertex rooted trees with reduced
matrices [[2,-1],[-1,1]] and I_2 have different rational inverse responses
but identical trivial quotient pairings. The parent lane supplied an even
sharper same-margin hostile: the labelled paths 0-1-2-3 and 0-2-1-3 have
the same labelled degree vector (1,2,2,1). With root 3 their reduced
Laplacians are unimodular, hence their embedded image lattices are both
Z^3 and all marked quotient characters are trivial. They nevertheless
have different labelled wiring. A full marked finite pairing plus node
margins is therefore still not a universal wiring classifier.

## 5. What the fractional Hall analogy genuinely transfers

For Abelian Hall effective theories an integral symmetric K matrix and an
electromagnetic charge vector t encode inverse-matrix response: the filling
factor is t^T K^(-1)t, mutual quasiparticle statistics involve K^(-1), and
genus-g degeneracy is |det K|^g. Tong's treatment also retains additional
geometric data such as the shift. See [Tong, chapter 5, pp.168-170](https://www.damtp.cam.ac.uk/user/tong/qhe/five.pdf).

The precise algebraic transfer is K=A and marked quasiparticle charge
vectors equal to port/divisor charges. It preserves the finite inverse
pairing and its characters. The metagraph provides no physical Hall fluid,
gap, electromagnetic t, or metric/geometric shift. No experimental phase
or Hall conductance follows from the tournament graph alone.

There is an extra parity boundary. If A is even, x^T A^(-1)x/2 modulo one
descends to its discriminant group. With odd diagonal it need not: A=[1]
already gives 0 and 1/2 for two representatives of the same zero class.
A characteristic covector c, with c_i=A_ii modulo two, repairs this:

    Q_c([x]) = (x^T A^(-1)x - x^T A^(-1)c)/2 modulo Z.     (5)

Replacing x by x+Az changes (5) by

    z^T x + (z^T A z - z^T c)/2,

an integer. Its polarization is b_A. This is a quadratic refinement with
the stated characteristic choice; it is not necessarily homogeneous.
Homogeneous modular-data uses require the corresponding additional choice
and hypotheses. This blocks an automatic application of THM-4762 to an odd
metagraph lattice. The finite bilinear/marked-phase guard (1)-(4) needs no
such extra choice and is the useful immediate transfer.

## 6. QCD triples suggest an incidence algebra, then reveal its loss

QCD has gauge group SU(3): quark colors form its three-dimensional
fundamental representation and gluons its eight-dimensional adjoint
representation. These are representation dimensions, not three graph
colors or eight tournament vertices. See [Particle Data Group, QCD review,
section 9.1](https://pdg.lbl.gov/2015/reviews/rpp2015-rev-qcd.pdf).

There are two exact algebraic maps worth extracting.

First, send a directed edge i->j to the root e_i-e_j. Summing over all arcs
gives the score-deviation vector 2s_i-(n-1). On three vertices the sum is
zero exactly for a directed cycle. The six possible roots form A2, the
root system of SU(3). But every regular tournament of larger odd order
also sums to zero. This map preserves net flow and loses cycle placement,
so the full endpoint incidence must be retained. The *undirected* triangle
has reduced Laplacian A2's Cartan matrix [[2,-1],[-1,2]], but reversing or
transitively orienting its edges does not change that Laplacian.

Second, retain matrix units E_ij for arcs, working in the complexified
matrix Lie algebra sl_n(C), not the compact real algebra su(n). For
distinct i,j,k,

    [E_ij,E_jk]=E_ik.

This preserves composable path endpoints. However the ungraded generated
Lie algebra is maximally lossy: every strongly connected tournament on n
vertices generates sl_n. A simple directed path gives E_ij by iterated
commutators; strong connectivity supplies both E_ij and E_ji, whose bracket
gives E_ii-E_jj. All off-diagonal units and traceless diagonals follow.

Therefore the useful operation is to retain bounded *graded words*, their
distinct-vertex conditions, and actual ports, instead of replacing them
by the generated algebra. This matches the present rooted five-vertex
incidence program: contractions on shared vertices are the missing datum
after separate neighborhood summaries agree. It does not supply a new
fixed five-vertex completeness theorem.

## 7. Bott, E8, and the correct notion of stabilization

As ordinary algebras, real Clifford stabilization gives Cl_(0,k+8)
isomorphic to Cl_(0,k) tensor M_16(R); the graded formulation additionally
retains the grading on the matrix factor. Complex periodicity is two.
These algebra/module and stable-topology statements are not an eight-step
graph return. See [Atiyah-Bott-Shapiro, sections 4-6](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/abs.pdf).

The reusable structural idea is *which data survive stabilization*. If U
is an integral unimodular symmetric matrix, A direct-sum U has the same
finite bilinear pairing as A: U^(-1) is integral, so its quotient is zero.
For even lattices an E8 summand preserves the discriminant quadratic form
and adds eight to signature. This is a precise reason a discriminant-only
observer loses a rank-eight summand. It is not the same theorem as Bott
periodicity, nor an explanation of order-eight tournament wiring.

An elementary graph hostile is even stronger: attaching any number of
unit leaves at the selected root adjoins an identity block and leaves the
finite pairing on the old ports unchanged. A boundary electrical response
also ignores dangling trees after their internal potentials equilibrate.
Thus the stabilization sidecar is the complete node/port embedding and
integer response, not merely a rank or a period. The marked-corner device
from the old Bott packet is the right design pattern.

## Connection contract and next test

| Source -> target | Preserved predicate | Lost data | Required sidecar / cheapest test |
|---|---|---|---|
| Native lifted potential -> critical charge -> port phases | Exact torsion order, linking values, nodewise modular conservation | Conductance magnitudes and full wiring | Actual endpoints and weights; apply (4) to a candidate switch |
| Rooted metagraph Laplacian -> Hall-style finite pairing | Algebraic inverse-pairing/charge formulas | Physical realization and geometric response | Charge vector, characteristic choice, physical model if making physical claims |
| Tournament arcs -> A_(n-1) roots | Net flow/score deviations | Which cycles own the flow | Rooted motif/port incidence; compare regular tournaments |
| Arc matrix units -> generated sl_n | Existence of composable paths | Grading, ownership, almost all strong-tournament distinctions | Keep graded distinct-vertex words and their boundary labels |
| Integral lattice -> stabilized discriminant pairing | Finite linking data | Unimodular summands and their vertices | Mark the original coordinate corner and all new ports |

The next precise target is to select a *small bank* of marked harmonic
potentials which annihilates the remaining feasible port-switch lattice.
This is a finite exact linear/congruence problem at each supplied order:
each candidate switch is tested by its defect on each retained potential.
The present bank of one modulus-eight potential closes the declared
four-terminal family. Uniformly deriving such banks from rooted motif
records, without already knowing full wiring, remains open.

## Reproduction

    python -B 05-knowledge/results/metagraph_hall_bott_20261023/check.py
    python -B -O 05-knowledge/results/metagraph_hall_bott_20261023/check.py

The standard-library consumer performs 438,333 generating checks (438,334
on frozen comparison), hashes all inherited inputs, computes the energy
both as an edge sum and as a sparse dot product, checks every Kirchhoff
congruence, checks primitive order and the four normalizers, and compares
independent 2-by-2 inverse/Smith/quotient controls. It does not call or import
the original dyadic engine. The receipt's LF SHA256 is
`ad1d2bbbb2856daa3414dad34aca3763608a2d34be3484436c0e17f4d978e84f`.
The parent lane supplies the independent full-port export audit. The new
receipt does not recertify all tournament ownership from scratch.
