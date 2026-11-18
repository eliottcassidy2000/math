# Index-p projective actions and the rigidity required by a +2 selector

**CITED classification, elementary proved selector obstruction, and literal tournament controls.** The exceptional labels5,7,11 do not by themselves supply Chan–Pak sources.

## Where the three actions genuinely arise

For primes p>=5, PSL(2,p) has a subgroup of index p exactly at p=5,7,11, with subgroups A4,S4,A5 respectively. This follows from Dickson's subgroup classification as stated in McCullough–Wanderley's research paper, Theorem2.2: https://math.ou.edu/~dmccullough/research/pdffiles/traces.pdf . An index-p subgroup has order(p²-1)/2, coprime to p. Its order is too large for a dihedral subgroup; its p-free part is too large for an affine subgroup. There is no proper subfield for prime p. The exceptional alternatives have orders12,24,60, giving precisely5,7,11. The replay explicitly constructs these three subgroups and their index-p coset actions.

These groups have even order and their faithful degree-p actions cannot preserve any tournament orientation. A nontrivial involution swaps at least one vertex pair; the directed arc between that pair cannot be preserved. By Cauchy's theorem, every tournament automorphism group therefore has odd order. The finite Platonic actions can organize undirected frames or a family of orientations, but are not automorphism groups of a single tournament.

This agrees with the repository's distinction between finite modular quotients and actual objects carrying their action: see THM-2768, `01-canon/theorems/THM-2768-modular-c2-c3-quotients-to-a4-s4-and-bass-serre-cycle-ranks.md`. No identity between that Bass–Serre construction and tournament counts is claimed.

## A rigorous filter on the target-selected source

Let G be ANY automorphism subgroup fixing a tournament's marked root. Every retained-root binary attachment extends G by fixing each new vertex. G acts freely on Hamiltonian paths: an automorphism fixing an ordered spanning path fixes every vertex. Thus |G| divides the ordinary count after EVERY word.

If two such outputs differ by2, |G| divides2. Since |G| is odd, |G|=1. Therefore:

> A physical source admitting a two-word +2 repair must be rooted-rigid.

For a sole-exit source, the entire rooted profile h,e,t,d is divisible by |G|. The weighted-start objects are ordered spanning path-cover pairs, so the same free-action proof applies. Its near-floor defect h+2eta-3d is consequently divisible by |G|. A nontrivial rooted symmetry is incompatible with defect±2.

This is a necessary structural filter, not a claim that rooted rigidity guarantees the target profile. It sharpens the inverse problem: the smaller source chosen from the target must satisfy the exact marked-count equations AND break every root-fixing symmetry.

## Literal odd-action controls at5,7,11

Use the cyclic5 tournament with differences{1,2}, and the quadratic-residue tournaments at7 and11. Their natural affine odd actions are respectively C5, C7:C3 and C11:C5. The selected deletion vertex0 has stabilizers1,3,5 in these explicit subgroups. Append a new sole-exit vertex beating only0, then root there. The resulting data are:

|deletion order|deletion H|endpoint eta|new source h|defect h+2eta-3H|
|---|---:|---:|---:|---:|
|5|15|3|33|-6|
|7|189|27|423|-90|
|11|95095|8645|239305|-28690|

The p=7 and p=11 defects cannot become±2 through any common rooted continuation because the surviving stabilizers have orders3 and5. At p=5 there is no such stabilizer obstruction, but the actual defect still misses±2. Absence of the obstruction is not sufficiency.

A useful counting identity supplies the table without guessing: if the deletion tournament Q is vertex-transitive of order p, ordinary count D and optional ordered two-cover count P2, then its weighted start at each vertex is(P2-D)/p and its endpoint count is D/p. The new sole-exit source has

    h=D+(P2-2D)/p, eta=D/p,
    defect=P2/p-2D.

Thus its ±2 readiness is the exact equation P2=2pD±2p. This is an additional physical profile equation, not an implication of the exceptional projective action.

The broader lesson is operational: preserve symmetry while it compresses enumeration, then break the root stabilizer before asking for a primitive +2 correction. This explains an obstruction and a required design change; it does not explain the special nonrealizability of7 and21 by the numerical coincidence of labels alone.

Replay: `python3 finite_edit_recontact/psl_rigidity_replay.py`. It constructs the matrix subgroups and coset actions, verifies affine tournament automorphisms, and independently counts the three sole-exit sources and their two-word repairs.
