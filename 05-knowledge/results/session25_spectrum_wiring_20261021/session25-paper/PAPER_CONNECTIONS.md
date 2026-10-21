# Session 25: exact acceleration, phase retention, and signed symmetry

## Source and status

The attached paper is Marco Dos Santos, Vivek Veeriah, Alex Havrilla, Andrew Bolt, Goran Zuzic, Satinder Singh, and Tom Zahavy, *At the Limits of Computability: Steps Towards Autonomous AI Discovery on the Busy Beaver Frontier*, arXiv:2610.22603v2, revised 19 October 2026. Source: https://arxiv.org/abs/2610.22603v2 . The PDF's linked proof repository is https://github.com/google-deepmind/cryptids . Both destinations were opened and verified. Local attachment: `2610.22603v2.pdf` (Library identity retained). The independent review rendered Appendix A.9 on printed/PDF page 18 and visually checked both rule families against the extraction.

The paper reports eleven new BB(6) machine accelerations, with formal proofs, plus further formalizations of previously studied machines. It does not resolve these eleven machines' halting, classical Collatz, or a tournament-spectrum conjecture. Definition 2.1, page 4, expressly makes its Cryptid classification provisional and broader than conventional community usage. Section 6, page 10, explicitly calls possible classical Collatz implications speculative. The paper is not an S2/S3 symmetry paper.

These notes prove elementary structural statements and check the printed A.9 reparameterization. They do NOT independently compile the paper's Lean files or certify its raw-machine implementation. The distinction matters: the correctness of the raw Turing-machine acceleration is reported from the paper, whereas the coordinate checks below are independently reproduced.

## The exact transferable mechanism

Section 3.4, page 8, fixes the raw machine and start configuration. Acceptance requires that each macro rule is proved against that machine and that the start and every nonhalting output match a rule input. Rules can be repaired, but the target machine cannot be changed. A claimed compression therefore needs local rule correctness AND global closure, not just successful samples.

Here is a precise sufficient mathematical interface. Let F:X->X be deterministic raw evolution with an absorbing halt state h_X. Let M include an absorbing macro halt state h_M, let e:M->X be an encoding, let G:M->M be a total domain-closed macro map, and let tau:M->positive integers. Terminal symbols may be adjoined to complete partial transition systems; require G(h_M)=h_M, e(h_M)=h_X, and choose tau(h_M)=1. Require

    F^{tau(m)}(e(m)) = e(G(m)),

with halt reflection (e(m) is halt iff m is halt), a verified raw prefix from the designated initial state to e(m0), and domain closure. Induction then gives

    e(G^k(m0)) = F^{sum_{j<k} tau(G^j(m0))}(e(m0)).

Because positive integer times have unbounded cumulative sum and raw halt is absorbing, the raw and macro trajectories halt together. This is a stuttering simulation, not necessarily a one-step conjugacy, and need not be injective. If raw durations are not available, reachability acceleration does not provide runtime or time-weighted density formulas. The printed A.9 rules supply initialization time but not a full physical-time polynomial.

For a proposed quotient q:M->Q with Q=q(M) (the actual image), a deterministic induced map exists if and only if q(G(x))=q(G(y)) whenever q(x)=q(y). Proof: necessity follows from commutation qG=Gbar q; sufficiency defines Gbar(q(x))=q(G(x)), independent of representative. Terminal status must also be constant on fibers. For preservation of the specified exact raw-duration label of every macro step, tau must factor through q as well. This is a requirement for the clock-labelled quotient, not a claim that endpoints uniquely determine raw durations. This is the criterion to test before any state label is discarded.

For a nondeterministic metagraph the corresponding condition is preservation of the outgoing transition relation (or multiset if multiplicities matter), including edge colors, endpoints, and ports. Scalar path counts, automorphism orders, and even degree data are not that transition relation.

## Appendix A.9: an exact binary-phase encoding, not a phase quotient

Machine: `1RB1RC_1LC1LC_1RE0LD_1LB1RF_---1RA_0RD1RD`.
The appendix defines two tape configuration families C0(a,b,c) and C1(a,b,c), whose first block lengths are 2a+2 and 2a+3 respectively. Their merger is the bijection

    phi(i,a,b,c)=(2a+i,b,c), i in {0,1}.

Its inverse reads i=A mod 2 and a=floor(A/2). No bit is lost. With all parameters nonnegative, the unified system is:

    C(a,0,0) -> Halt
    C(2k,0,1) -> C(5,2,k+1)
    C(2k+1,0,1) -> C(2,0,k+3)
    C(2k,0,c+2) -> C(2,k+3,c)
    C(2k+1,0,c+2) -> C(5,k+3,c)
    C(a,b+1,0) -> C(a+3,0,b)
    C(a,b+1,c+1) -> C(a+6,b,c).

Start is C(5,2,0), reached after 100 raw TM steps. Direct substitution of 2a+i into all old rule branches proves phi G_old=G_new phi. The supplied checker independently verifies all 59,582 states with i in {0,1}, a,b,c in 0..30.

The tempting map (i,a,b,c)->(a,b,c) is NOT a dynamical quotient. For input (a,b,c)=(0,0,1), phase0 gives C1(2,2,1), while phase1 gives C0(1,0,3); forgetting the output phase still gives different triples. Thus neither the phase toggle nor the naive phase-forgetting projection preserves these transitions. The number two does not create an S2 symmetry.

## Recovering the actual previous S2 x S3 construction

Source: `math/05-knowledge/results/arithmetic_seams_20260921_dynamics.md`, sections 2-3; also its synthesis companion. Three distinct rational projective points have a unique projective cyclic interpolant. In a suitable lifted basis:

    B = [[0,-1],[1,1]], J = [[1,1],[0,-1]],
    B^3=-I, B^6=I, J^2=I, JBJ=B^-1.

In these displayed coordinates take v=(1,0). Its six B-orbit vectors lie above the three marked projective points, and J fixes v and sends B^k v to B^(-k)v. The exhibited group acting on this standard orbit is

    <B,J> = D6 of order12 = <B^3> x <B^2,J> ~= S2 x S3.

For a general nonzero rational or real v=(x,y), the fixed displayed J need not preserve its B-orbit: v=(1,2) has Jv=(3,-2), outside that orbit. Instead P_v=[v,Bv] has determinant x^2+xy+y^2>0 and commutes with B. The reflection J_v=P_v J P_v^-1 fixes v and satisfies J_v B J_v=B^-1, giving the corresponding order12 action on that orbit. Thus the reflection must be transported with the chosen orbit/basis; it cannot silently be held fixed. The earlier repository proof makes precisely this basis-dependent choice.

Projection forgets the central sign and has kernel {I,-I}; it is genuinely equivariant for the projective S3 action. But only the cyclic order3 subgroup preserves the directed three-cycle. J reverses arrows/time; it does not commute with the cycle map. In the lifted order12 group exactly the six rotations commute with B, independently checked by the script.

This yields a real structural bridge to tournaments: Aut(C3)=C3, whereas its automorphism/anti-automorphism group is S3. Ordinary tournament automorphisms cannot include any nontrivial involution: an involution swaps some pair u,v and would have to preserve and reverse the arc between them. Hence finite tournament automorphism groups have odd order (Cauchy's theorem). Neither S2 nor S3 is an ordinary tournament automorphism group. They can occur as signed symmetries, actions on decorations, or quotients of different objects, with the action stated explicitly.

The concrete six-vector example for the rational cycle of x^2-29/16 is also in the earlier note; agreement of a degree1 projective map with the degree2 polynomial on three points is not a global conjugacy. No map from this six-point construction to actual classical Collatz dynamics has been established here.

## Consequences for complete metagraph wiring and spectrum work

The session24 synthesis reports a block-contiguous path lift T->T[C3] commuting with converse reversal and preserving signed node status. It also records the explicit failure to commute with fixed-path complementation kappa:

    H(kappa(C3[C3]))=2983, H(kappa(C3)[C3])=27.

Thus that node lift alone cannot transport complete colored adjacency. This is the same logical boundary as the A.9 phase counterexample: the representation may preserve some invariants while failing the transition equation needed for the intended conclusion.

A useful certificate for the session25 wiring lane should therefore carry the actual marked path/port state, target node, converse action, and edge multiplicity. Test each proposed compression against the quotient criterion above. A proof of that criterion plus closed rules would be a genuine transfer of the paper's method; a matching occurrence of 2,3,6,7,21 is not a transfer theorem.

For Hamiltonian-count spectrum the paper offers a discovery-and-certification method, not an existing theorem transferring BB transitions to H(T). No structure-preserving map was found from any printed cryptid recurrence to the tournament-count spectrum. Likewise there is no claimed bridge proving a new Hadwiger-Nelson lower bound. These are explicit scope limits, not evidence of impossibility.

## Verification and provenance

`verify_connections.py` and `verification.json` contain reproducible finite controls: 59,582 A.9 states, the binary-phase counterexample, exact group order12, and centralizer order6. The all-parameter reparameterization proof is by the seven guarded branches, not inferred from the census. The signed group proof is algebraic, not inferred from enumeration. No repository files were modified by this lane.


## Independent-review lineage

This is a reviewed copy of the producer note, not an edit of its original. Corrections make the orbit-dependent reflection, quotient image, absorbing halt conventions, clock labels, and A.9 source page explicit. The original and reviewed hashes are recorded in the separate review manifest.
