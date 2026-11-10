# Session240: the missing 8- and 12-valued orientation coordinate

Status: **PROVED** general matching and induced-action statements below; **FINITE-EXACT** enumerations at 3, 7, 11, 19 and all 32 circulant connection sets at 11. No new global Hamiltonian-path extremality or Chan–Pak coverage theorem is claimed.

## Inheritance and question

Inherited source: `05-knowledge/hypotheses/HYP-9357-paley-hamiltonian-path-maximality-and-the-kantor-galois-exceptional-primes-3-7-11.md`, with correction `01-canon/MISTAKES.md`, MISTAKE-651. The latter repairs the Kantor citation and the passage from difference-set equivalence to design isomorphism; it does not itself assert that tournament and design groups coincide. The source correctly distinguishes the index-p design actions from the tournament. Prior exact design computation is `04-computation/experiments/tournament_pseudoprimes_20261021/audit/a20_designs.py` (nauty). Group identifications PSL(2,7), PSL(2,11) and stabilizers S4, A5 are inherited, not new identifications from group order alone.

Concept board: (1) Hamiltonian paths and endpoint counts; (2) Paley tournaments; (3) neighbourhood incidence designs; (4) skew point-block matchings; (5) homogeneous spaces and induced actions. Anchor: identify exactly what information separates 21/55 from 168/660. Niche: test the matching fibre against p=19. Wildcard: recover the familiar p+1 action from the exceptional p-point action.

## 1. Exact restoration sidecar (PROVED)

Let a finite incidence structure have n points V and n distinct blocks B. A **compatible skew matching** is a bijection f: V -> B such that

- x is not incident with f(x);
- for x != y, exactly one of y in f(x) and x in f(y) holds.

Then x -> y iff y in f(x) is a tournament T_f, whose out-neighbourhoods are exactly B. Conversely, every tournament on V with precisely this set of out-neighbourhoods gives a unique such f. This is a bijection, not an analogy.

Write G for the side-preserving design automorphism group. Because blocks are distinct, each point action alpha uniquely determines the block action beta. G acts on compatible matchings by

    (alpha,beta) . f = beta o f o alpha^{-1}.

The stabilizer G_f is exactly Aut(T_f): preserving the matching is precisely beta(f(x))=f(alpha(x)), which is precisely preservation of tournament arcs. Thus forgetting the point-block correspondence loses a G-set of possible tournament orientations. In general it need not be a single orbit. It is not a torsor: stabilizers need not be trivial.

For the canonical Paley matching f(i)=B_i, its stabilizer is the diagonal condition alpha=beta. This turns a vague instruction to 'remember block labels' into an exact computable coordinate.

## 2. New finite-exact fibre calculation

Take V=F_p and B_i=i+QR_p, so i -> j iff j-i is a nonzero square. Dependency-free exhaustive incidence backtracking and an independent matching search give:

| structure | design group | tournament group | all compatible matchings | one G-orbit? |
|---|---:|---:|---:|---|
| Paley 3 | 6 | 3 | 2 | yes |
| Paley 7 | 168 | 21 | **8** | yes |
| Paley 11 | 660 | 55 | **12** | yes |
| Paley 19 | 171 | 171 | **1** | yes |
| rotational 7 | 14 | 7 | 2 | yes |
| rotational 11 | 22 | 11 | 2 | yes |

The substantive extra statement is **all** compatible matchings have been found, and they constitute one orbit. Merely computing 168/21 or 660/55 would not establish this: other orbits could have existed.

For 7 and 11 the matching action is faithful and 2-transitive: the stabilizer of a matching has orbit sizes (1,7), respectively (1,11), on all matchings. Hence the incidence design recovers a concrete 8- or 12-point faithful action of the inherited PSL group from its exceptional 7- or 11-point action. The canonical matching distinguishes one of those 8 or 12 points and leaves exactly the tournament group. We do not claim an explicit coordinate isomorphism with the projective line here.

The original point stabilizer has size 24 or 60 and orbits (1,6) or (1,10). After matching restoration its size drops to 3 or 5 and its orbits split into (1,3,3) or (1,5,5). Thus the restored sidecar restores the two arc orientations that the 2-transitive design action merges.

**Stronger hostile test at order 11.** All 2^5=32 labelled skew connection sets were enumerated, not merely the Paley and interval examples:

- 20 have design group 11, tournament group 11, and 1 matching;
- 10 have design group 22, tournament group 11, and 2 matchings;
- 2 have design group 660, tournament group 55, and 12 matchings (QR and its complement).

Every case has one matching orbit. These are labelled connection-set counts, not isomorphism-class counts. In this complete universe, 12-valued restoration is isolated to the Paley pair. At 19 the Paley design itself uniquely recovers the tournament, so neither the 8/12 phenomenon nor p+1 matchings extrapolates.

## 3. Hamiltonian-path bridge and its exact limit

For a single matching f, let H_f be its set of directed Hamiltonian vertex sequences. A design element outside G_f need not send H_f into itself. Indeed 168 does not divide 189, and 660 does not divide 95095. Their remainders are 21 and 55.

The correct G-set is the disjoint union

    X = {(f,h): f a compatible matching, h a Hamiltonian path of T_f}.

If all matchings form one orbit, then X is the induced G-set G x_{G_f} H_f. The G-action is free: fixing a pair fixes the ordered path, hence all n points, hence the design automorphism. Moreover X/G is canonically in bijection with H_f/G_f. To prove this, transport the matching of a pair to f; any two such transports differ by G_f. This also proves injectivity and surjectivity of the quotient correspondence.

At 7 and 11 the exact counts are:

| p | H(T_f) | matchings | size of X | free G-orbits = free tournament orbits |
|---|---:|---:|---:|---:|
| 7 | 189 | 8 | 1512 | 9 |
| 11 | 95095 | 12 | 1141140 | 1729 |

Thus the larger design group is legitimately available only after enlarging the object by precisely the missing matching coordinate. The factor cancels: this gives **no improved divisibility denominator** for H(T_f), and no new forbidden/attainable raw Hamiltonian-path value in the Chan–Pak problem. It does provide the exact action and obstruction for attempts to import exceptional symmetry into that problem.

## 4. Endpoint witness: the lost coordinate changes a real observable

Independent subset path DP, with initial vertex 0, gives for Paley 7 exactly 1 path to each square endpoint and 8 to each nonsquare endpoint. For Paley 11 the counts are 628 and 1101. These are not equal, despite 2-transitivity of the bare design. They are a concrete obstruction to treating the design action as an endpoint-preserving symmetry of the fixed tournament.

The affine square-multiplier tournament group acts regularly on each of the two ordered distinct-endpoint orbits. Consequently the free path orbit totals split as 1+8=9 and 628+1101=1729. Path totals are 21(1+8)=189 and 55(628+1101)=95095. These endpoint refinements are checked here; H totals and tournament orders were inherited and are independently reproduced.

## 5. Reproduction, completeness and limitations

From repository root:

    python3 05-knowledge/results/session240_seven_eleven_symmetry_20261110/paley/verify.py
    python3 05-knowledge/results/session240_seven_eleven_symmetry_20261110/paley/audit.py

Both commands were also run with `python3 -O`: normal and optimized outputs are byte-identical. Every validation uses an explicit `check` helper raising `RuntimeError`, not an optimization-sensitive assertion. A deliberately false check was verified to raise under `-O`.

Outputs: `output.jsonl` (JSON records followed by the pass marker) and `audit.out`.

Incidence search assigns every point permutation one coordinate at a time. Blocks are partitioned by their incidence bit-pattern on already assigned source points and candidate target points. A branch is rejected only when corresponding cell sizes differ, a necessary condition for any completion. At full depth distinct blocks have distinct signatures, producing the unique block action. Every completed action is independently checked against all block images. Thus the search neither assumes affine form nor uses inherited group orders as a filter.

Matching search assigns every block bijection, pruning only the required diagonal and already decided skew conditions. Hence it is exhaustive. The orbit obtained by transporting the canonical matching is explicitly compared with the entire matching list. `audit.py` independently scans all 7! point permutations and all 7! block matchings at 7, and exhausts all 32 order-11 circulants. The 19 computation exhausts matchings and group permutations via backtracking but does not compute Hamiltonian-path counts. No prime-power or all-prime conclusion follows.

Novelty here means new within this Session240 investigation: the explicit complete matching fibres, their faithful 2-transitive actions, the endpoint obstruction, and the all-32 control census. We make no literature-priority claim.
