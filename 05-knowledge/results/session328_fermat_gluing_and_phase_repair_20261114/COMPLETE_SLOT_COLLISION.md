# A Fermat triangle can erase the whole slot difference

The ordinary gluing N=H+2^j L has a higher-slot transfer that tells which source choices really have the same cost. Here there is a literal pair whose **entire optional-slot polynomial** coalesces after gluing the smallest Fermat block, even though their original Hamiltonian counts and rooted responses differed.

## Exact optional-slot gluing

For a rooted A at v and integer z≥1, reserve one designated slot for v. Delete v and split the reserved slot into the pre-v piece P and post-v piece Q, leaving z−1 ordinary optional slots. Let L_z(A) count these covers if only the outgoing join is enforced (Q empty or v→first(Q)). Let E_z(A) count covers with Q empty and the incoming join enforced (P empty or last(P)→v). Write q_A=Q_(A−v)(z). Also write J_A=Q_A(z)/z for the fully lawful covers with this designated slot; define the dual incoming relaxation R_A and first-at-v count S_A similarly.

Take A and a right block B sharing v, with all cross arcs A−v→B−v. Outside the reserved slot the two pieces concatenate in that order. Inside it the four cases, according to whether B contributes a nonempty pre-v piece and whether A contributes a nonempty post-v piece, give

    Q_out(z)/z = (L_A−q_A)(R_B−q_B) + (J_A−E_A)q_B
                 +q_A(J_B−S_B) + E_A S_B.

For example with both pieces present, their joins to v are the only restrictions; with neither present, A must end and B must start at v. The four selectors are disjoint and recover the original split covers.

Now take B=R_j, the order-(j+2) near-transitive block with extreme edge reversed, rooted at its first vertex, j≥1. Its old chain is transitive of length j+1, with v beating all but the last. The only legal nonempty B pre-piece ends at that last vertex. Thus R_B−q_B=J_B−S_B=(z+1)^j: put the last vertex in the pre-piece and assign the other j to any of the z+1 split slots. When pre is empty the only disallowed post-piece is the singleton last vertex, leaving (z−1)^j choices for the other slots, so S_B=z^(j+1)−(z−1)^j. Consequently

    Q_out(z)= z^(j+1) Q_A(z) + z (z+1)^j L_z(A) − z (z−1)^j E_z(A).    (1)

At z=1 this is the Hamiltonian seam H+2^j L. At z=2 it is the exact operational cost 2^(j+1)Q_A(2)+2·3^j L_2−2E_2. For equal-Q inputs, two distinct gluing lengths can both retain equality iff L_2 and E_2 both match. A single observed equal clock does not suffice. These formulas are a cover proof, complementary to the signed odd-cycle transfer in [Fermat clock memory](FERMAT_CLOCK_MEMORY.md).

## A literal cancellation and its all-length difference identity

Use the following order-eight strong tournaments, both rooted at vertex 6 (zero-based outgoing masks):

    A=(6,12,80,149,163,15,59,103),
    B=(2,52,73,131,45,77,155,55).

They share (Q_A(2),L_2,E_2)=(5536,6163,433), hence gluing either with the same R_j has exactly the same slot budget at every j. The scalar profiles (H,L_1) are (349,860) and (355,857). Different sources therefore cost the same before and after the operation. With the triangle R_1 they produce the *same* H=2069, Q2=58256 and Q3=735021. The outputs are nonisomorphic even up to converse: sorted outdegrees are (1,1,4,4,4,6,6,6,6,7) versus (1,1,3,5,5,5,6,6,6,7), and complementing either cannot give the other (either complement has maximum 8, while each original maximum is 7).

In fact let D(z)=Q_A(z)−Q_B(z). Exact slot DP at z=0,…,n on the source and first two glued graphs gives these polynomial identities:

    D(z)=2 z²(z²−4),
    ΔQ_(j=1)(z)=0,
    ΔQ_(j=2)(z)=−z D(z).

Degrees are at most the actual graph orders, so agreement at those independent integer values is exact. Substituting j=1,2 into (1) and eliminating the two rooted differences forces ΔL_z=−D(z)/2 and ΔE_z=D(z)/2. Thus for **every** common length j≥1 the full response difference is

    ΔQ_j(z)= z D(z) [z^j−((z+1)^j+(z−1)^j)/2].

At j=1 the block erases the entire difference. At j≥2 the Hamiltonian difference is 3·2^j−6, whereas the two-slot difference stays zero and the three-slot difference has the opposite sign. In particular j=11 gives H difference 6138=31·198, so the two H values are complementary first-phase partners mod396 while paying the same literal Q. This last j is near-transitive, not a genuine Fermat-number size; the identical-polynomial cancellation at j=1 is a genuine Fermat triangle.

These two sources have different Q_(A−v)(2), namely1184 and1116. Keeping the old root for a further common seam does not retain the same cost; matching only the observable output at the first step is not a state-closed induction. There is nevertheless a further exact move: reroot each ten-vertex output at vertex8, the first new vertex. Their rooted Q2 arms then agree, and appending another R_1 gives strong equal-Q2 sources with H=6853 and6859, Q2=403552. This is a verified two-step route with a purposeful root change; no claim of repeating it indefinitely is made here.
