# Linking forms beyond the dyadic tower

23 October 2026. **PROVED:** the lattice/type distinctions, exponent-two
linking boundary, characteristic-refinement switch law, and odd-scaling
separation below. **FINITE-EXACT:** all 64 order-five fixed-path presentations,
all 42 connected disjoint paired-black two-edge switches, and the terminal
order-eight potential witness. **NOT CLAIMED:** a new tournament realization,
an all-order short reconstruction signature, or a physical interpretation.

**Incoming reconciliation:** before freezing, read-only recovery of
`bc3c88ec033` found the parallel `metagraph_hall_bott_20261023` and
`metagraph_linking_lifts_20261023` packets. They already establish the
inverse-free charge map, native order-eight summand `⟨1/256⟩`, its
highest-unit squareclass, the characteristic-covector convention, and
small equal-rank-tower/different-linking controls. Those are prior
repository mechanisms, not new claims of this packet. Their native
generator has numerator 105; ours has 153, with
`9²·153 ≡105 (mod 256)`, so the cyclic forms agree. This does not assert
that the actual saved node potentials are related by multiplication by
nine; their node conventions and lower-order components can differ.
The new finite increment here is the exhaustive order-five quadratic
switch census and its four positive refinements. The odd-scaling control
also makes the comparison explicitly of the *entire labelled tower*,
not merely its dimension sequence or unmarked terminal forms. The n8
terminal computation is an independent replay in the Session38 bank.

## Inheritance and the useful increment

Read-only inheritance used repository HEAD
`e41a2c82b3` (the incidence-wiring bank reconciliation), including the preceding
one-bit order-eight repair and marked odd-cell work. The closest mechanism is
the Session38 intrinsic dyadic tower in
`session38_incidence_ports_electrical_20261023/electrical/PROOF.md`.
The canonical hostile is its order-five switch with tree counts 41096 and
36680. The corrected near miss is identifying the critical group with torsion
in the graph's own homology, or treating mod-two tower forms as the full
integral discriminant form. The least-used sidecar is a characteristic
covector, needed to define a quadratic refinement when the lattice is odd.

Anchor / niche / wildcard: reconstruction constraints beyond the dyadic
tower / integral cut and cycle lattices / exact finite quadratic phases.
The concept board consists of labelled ports, loop-blind Laplacians,
cut discriminants, dyadic filtrations, odd unit classes, and characteristic
refinements. Every connection below specifies what is retained and lost.

Two precise results emerge:

1. The original order-five hostile is invisible to its **entire 2-primary
   linking pairing**, not just the Smith group or the tower. A rooted
   characteristic quadratic refinement detects four other feasible graph
   switches that the same tower misses, but still misses the original one.
2. The full native order-eight graph has a cyclic discriminant summand
   `⟨1/256⟩`. Odd scaling by three leaves its entire labelled dyadic tower
   unchanged but changes this linking summand to `⟨3/256⟩`, nonisometric even
   after allowing all other 2-primary summands. Scaling is a colored-network
   control, not a same-fibre actual-port counterfeit.

## 1. Graph topology, cut lattices, and critical groups have different types

Expand integral edge multiplicities as parallel literal edges. Let B be an
oriented vertex-by-edge incidence matrix; a loop has zero column. For a
connected graph with n vertices and m literal edges,

    H_1(G;Z) = ker B = Z^(m-n+1),     H_0(G;Z) = Z.

In particular graph homology has no torsion here. Write Z(G)=ker B for the
integral cycle lattice and C(G)=im B^T for the integral cut lattice, both
inside the standard unimodular edge lattice Z^E. Their real spans are
orthogonal complements. Since B maps Z^E onto the zero-sum vertex lattice,
with kernel Z(G), it induces

    Z^E / (Z(G)+C(G))  ≅  Z^V_0 / im(BB^T)  = K(G).

The finite critical group is therefore a lattice gluing quotient, not graph
homology torsion. Delete one vertex row from B; its rows form a basis of
C(G), with nonsingular Gram matrix A equal to the grounded Laplacian.
The cut discriminant C(G)^*/C(G) is coker A, equipped with

    b_A([x],[y]) = x^T A^(-1)y modulo Z.                 (1)

It is well-defined because adding Az changes the pairing by the integer
z^T y. It is nondegenerate: pairing integrally with every integer y forces
A^(-1)x to be integral. A basis change transports the discriminant form,
so (1), up to isometry, does not depend on the grounded cut basis.

A loop supplies a primitive orthogonal coordinate of the cycle lattice.
It disappears from the gluing quotient and from every invariant computed
only from A, including all linking and quadratic data with an A-defined
sidecar. **Boundary:** if the native total degree d_v is also supplied, loop
multiplicity is recovered as (d_v-A_vv)/2. With colors, use separate color
degrees. Even A plus these degrees cannot identify individual path-labelled
port matchings inside an unchanged aggregated edge cell.

No 2-cells, 3-manifold, physical flux space, or quantum dynamics are supplied
by this graph construction. Such an interpretation would require an actual
additional map and its choices. The finite phase calculations below are
ordinary exact invariants of quadratic modules.

## 2. Why linking cannot repair the original order-five hostile

For any nonsingular integral symmetric A, put V=ker(A modulo 2). The map

    phi_A : V → (coker A)[2],     u ↦ [A u_hat / 2]      (2)

is an isomorphism. It is independent of a change u_hat→u_hat+2z. It is
injective because Au_hat/2=Az implies u_hat=2z. If 2x=Az then z is harmonic
modulo two and x=Az/2, proving surjectivity. Under (2),

    b_A(phi_A(u),phi_A(v)) = u_hat^T A v_hat/4
                          = q_1(u,v)/2 modulo Z.       (3)

Thus q_1 always determines the restricted linking pairing on the
2-torsion subgroup. If the whole 2-primary group has exponent two, it
determines the entire 2-primary linking form. The exponent hypothesis is
needed only for this last sentence, not for (2).

At order five the native q_1 is nonsingular on dimension three. Both the
native graph and the inherited counterfeit

    remove paired-black 1--4, 3--7; add 1--7, 3--4

have exactly the same labelled q_1 and critical 2-group (Z/2)^3. Hence their
2-primary linking forms agree under the canonical harmonic identification.
Their tree counts remain 41096 and 36680, so the full critical groups are
already separated by their odd orders. This packet makes no stronger claim
about matching full critical groups.

## 3. A declared characteristic refinement detects four other switches

A characteristic covector is an integer c with c_i≡A_ii modulo two. Define

    q_c([x]) = (x^T A^(-1)x - x^T A^(-1)c)/2 modulo Z. (4)

Changing x to x+Az changes (4) by
x^T z+(z^T A z-z^T c)/2, an integer. Its polarization is b_A. For an even
lattice, c=0 is the natural homogeneous choice. For an odd lattice a choice
is required; it is not valid merely to divide b(x,x) by two without one.
Here root 9 and the grounded degree covector c=diag A are explicitly
retained. Choosing diag A again after an arbitrary basis/root change need
not transport the same refinement. This is a sidecar-dependent invariant.

For a paired-black switch removing ab,cd and adding ac,bd, all four vertices
distinct, write alpha(u)=u_a-u_d and beta(u)=u_b-u_c. The literal Laplacian
changes by 2(xy^T+yx^T), where x=e_a-e_d and y=e_b-e_c. Its diagonal and
hence the chosen c are unchanged. Substitution into (4) at phi_A(u) gives

    Delta q_c(u) = alpha(u) beta(u)/2 modulo Z,         (5)
    Delta q_1 = alpha⊗beta + beta⊗alpha over F2.        (6)

When (6) vanishes, either one functional is zero, or alpha=beta. In the
latter nonzero case (5) is the extra linear character alpha/2. This is the
precise information missed by the bilinear half-form, not a new group
factor or a deeper nonzero tower layer.

The exhaustive census has 42 connected switches: 25 change q_1 and 17 do
not. Of those 17, four change the quadratic value distribution and thirteen
leave the entire labelled quadratic refinement unchanged. All seventeen
have tree valuation exactly three, independently checked by Bareiss
determinants, and hence complete 2-primary group (Z/2)^3.

A positive example removes paired-black 4--7 and 5--6, then adds 4--5 and
6--7. Its tree count is 34360. The counts of quadratic values
0,1/4,1/2,3/4 are respectively

    native:     1,1,3,3;      counterfeit: 3,1,1,3.

Their exact Gauss sums sum_x exp(2πi q_c(x)) are -2-2i and 2-2i. Different
value counts prove nonisometry of the quadratic refinements without any
classification theorem or floating-point phase computation. These are
same-degree colored graph switches, not claims of another tournament
realization, and no endpoint-signature transport is certified here. The
original 41096/36680 hostile remains identical even under this refinement.

## 4. Independent replay of the native order-eight cyclic linking summand

Inherit the independently audited Session38 decomposition

    K(G)[2^∞] = (Z/2)^3358 ⊕ (Z/4)^4 ⊕ Z/256.

The frozen terminal potential v has root coordinate zero, an odd
coordinate, and Av divisible by 256. The new checker independently replays
every coordinate of Av from the complete compressed literal edge ledger.
It obtains

    v^T A v = 14960072960,
    v^T A v /256 = 58437785 ≡153 modulo 256.

Consequently g=[Av/256] has exact order 256 and b(g,g)=153/256 modulo Z.
The odd numerator makes its cyclic subgroup nondegenerate, hence an
orthogonal direct summand. Since 13²·153≡1 modulo 256, replacing g by 13g
normalizes its form to ⟨1/256⟩. Its orthogonal complement has exponent at
most four, by the inherited group decomposition.

For any other element of order 256, write h=a g_0+r using this orthogonal
splitting, with a odd, ord(r)≤4, and b(g_0,g_0)=1/256. Then
256 b(h,h)≡a²≡1 modulo eight: the r contribution is a multiple of 64.
The same conclusion holds without first choosing an orthogonal expression,
because cross terms also vanish modulo eight. Thus the highest-order unit
class modulo eight is an invariant of the **entire unmarked 2-primary
linking form**, not just a specially selected generator.

## 5. Odd scaling preserves the whole tower but changes the linking form

For an odd positive t, set A'=tA. For every k,

    {z: A'z divisible by 2^k} = {z: Az divisible by 2^k},
    q'_k(u,v) = t q_k(u,v) = q_k(u,v) over F2.

Hence the complete labelled tower (V_k,q_k), not only its ranks, is
identical. This statement concerns the mathematical tower, not the full
integer currents retained as computational certificates.

Use the potential-based identification [x]↦[t x] between 2-primary parts.
It is injective since tx∈tAZ implies x∈AZ, and is surjective on 2-primary
parts because the added group index is odd. Under this map b scales by t.
For t=3 the n8 generator g'=[3Av/256] has self-linking 203/256; after the
same multiplier 13 it has self-linking 3/256. The highest-order unit class
is now three modulo eight, so the two 2-primary linking forms are not
isometric. An alternative localized-identity convention scales by 1/3 and
gives 51/256 instead; these descriptions differ by an odd-square change of
generator and give the same unit class three. The receipt uses g', not the
localized-identity convention.

This is a literal odd-conductance scaling of the actual n8 graph. It keeps
colors, support, connected components and black parity but changes fibre
degrees and capacities; it is **not** an admissible actual-port wiring
counterfeit under fixed native margins. It proves the narrower structural
claim that the complete dyadic tower does not determine the 2-primary
linking form. The elementary scalar controls A=[4] and A'=[12] give the
same separation on Z/4 with forms ⟨1/4⟩ and ⟨3/4⟩.

## 6. Reproduction and audit boundary

The checker imports no earlier research producer. It independently rebuilds
all order-five presentations and canonical labels, then all 42 switches.
Fraction-free determinants, exact rational inverses, all 2-torsion
pairings, characteristic representative-independence, quadratic
polarization and both switch formulas are checked directly. The n8 part
replays the saved terminal integer lift from every frozen edge record; it
does not recompute the earlier complete 3528-node tower or a huge inverse.

    python3 check.py --wiring PATH_TO_SESSION38_WIRING
    python3 -O check.py --wiring PATH_TO_SESSION38_WIRING

The default compares the frozen receipt without rewriting it. Explicit
`--write` regenerates it. Required external input is `edges.tsv.gz` from
the Session38 full-wiring packet, checked by SHA256; the terminal witness
is included here. The source of that witness is the k=8 entry of
`electrical/n8_complete_tower_receipt.json` in the Session38 electrical
packet. Both normal and optimized runs must agree, with no removable
assertions. `probe.py` is exploratory and is not part of the checkpoint.

**Open boundaries:** no same-native-port-margin n8 counterfeit with the
entire tower unchanged and different linking has been constructed here;
no intrinsic root-free characteristic refinement is asserted; no
Laplacian-based invariant authenticates literal port identities within
unchanged aggregated cells. All-size reconstruction still needs its
authenticated incidence sidecar.
