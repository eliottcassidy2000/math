# Sharp strong eleven-symmetry gap, extending the inherited seven-symmetry theorem

Session240, 2026-11-10. **PROVED from inherited orbit-profile and parity facts plus the FINITE-EXACT census below.** This packet specializes an existing general prime-action theorem; it does not claim a new general theorem or literature priority. The verifier is `probe.py`, with identical normal and optimized outputs `probe.json` and `probe.optimized.json`.

## 1. Objects, inheritance, and scope

For a finite tournament T, H(T) counts ordered directed Hamiltonian paths. Let H(empty)=1. All symmetry below means an actual automorphism of the actual tournament, not a formal quotient or a divisor of H alone.

Inherited sources, identified by ID and path:

- THM-5099, `01-canon/theorems/THM-5099-value-refined-scc-recurrence-and-sharp-seven-symmetry-gap.md`, and `05-knowledge/results/session239_sequence_engines_20261110/proof_bridges/PROOF.md`, section 4: the general prime-action gap reduction and sharp seven-symmetry gap 2163. Its independent audit is `05-knowledge/results/session239_sequence_engines_20261110/audit/strong_gap_audit.md`.
- THM-5094, `01-canon/theorems/THM-5094-cyclic-orbit-profile-cost-and-sharp-sevenfold-path-floor.md`: cyclic orbit-profile product bound and its ordinal-sum equality characterization. The p=11 core minimum 92411 already appears in `05-knowledge/results/session238_planar_chanpak_20261110/cyclic_cost/PROOF.md` and its exact census. It also occurs much earlier in THM-132, `01-canon/theorems/THM-132-z11-ocf-alpha1-not-monotone.md`. This minimum is inherited and independently replayed, not newly discovered.
- THM-3121, `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`: the path-cover substitution kernel. Its ordered two-cover specialization is reproved below to fix conventions.
- LEM-003, `01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md`: automorphisms act freely on ordered Hamiltonian paths. Redei parity, as used in THM-5099, says H(T) is odd.

Closest mechanism: separate one moving prime orbit from two or more, then retain the two-cover sidecar. Canonical hostile: transitive singleton padding preserves H and prime symmetry at unbounded order if strongness is dropped. Corrected near miss: a subgroup denominator cannot be replaced by the full automorphism-group order. Least-used sidecar: ordered two-path covers, which vary at p=11 even though they were constant at p=7.

Research board: cyclic orbit profiles; strong quotient triangles; induced path monotonicity; ordered two-cover costs; finite multiplier orbits. The connection is from an order-p action on T to a cyclic core R and its triangle inflation. It preserves an induced lower-bound witness but forgets all other vertices and most ambient arcs. The missing numerical coordinate is F2(R), not just H(R). The cheapest decisive test enumerates all 32 eleven-circulants, with no assumed classification of arbitrary eleven-vertex tournaments.

## 2. Two-cover identity and general prime gap (inherited, reconstructed)

Let R be a tournament and let C3[R,1,1] denote its triangle inflation with arcs R -> x -> y -> R. Define the **ordered**, nonempty two-cover count

    F2(R) = sum over nonempty proper S subset V(R) of H(R[S]) H(R[V(R)\S]).

Then

    H(C3[R,1,1]) = 3 H(R) + F2(R).                       (1)

Indeed, a spanning directed path using R in one consecutive run has one of the three patterns R,x,y; x,y,R; y,R,x. These contribute 3H(R). Otherwise the unique possible run pattern is R1,x,y,R2. Both R-runs are nonempty, and choosing their ordered vertex partition and directed paths gives exactly F2(R). No other run pattern is possible because x and y each occur once and their adjacency to R is prescribed.

For an odd prime p define

    mu_p = min H(R),
    gamma_p = min (3 H(R) + F2(R)),

where R ranges over all tournaments fixed by a specified p-cycle, equivalently all 2^((p-1)/2) difference-set orientations on Z/p.

**Inherited prime-gap theorem (THM-5099).** If T is strong and admits an automorphism sigma of order p, then either |T|=p, in which case T is a p-circulant, or

    H(T) >= min(gamma_p, mu_p^2 + 2p).                  (2)

Here is a full reconstruction of the reduction. The sigma-orbits have size 1 or p; at least one has size p. Let s be the number of moving orbits.

If s>=2, THM-5094 gives H(T)>=mu_p^s. Equality would force an ordinal sum of the orbit blocks. There are at least two blocks, so equality contradicts strongness. Both H(T) and mu_p^s are odd multiples of p: divisibility follows from the free path action on T and on every p-circulant, and oddness is Redei parity. Therefore

    H(T) >= mu_p^s + 2p >= mu_p^2 + 2p.               (3)

The additive 2p is not a consequence of the product bound alone; it also uses its equality characterization, parity and divisibility.

If s=1 and |T|>p, let R be the moving orbit. All other vertices are fixed by sigma, so each sees R uniformly. Thus R is a module. Contract R to r to obtain a strong quotient Q: projecting directed paths gives directed walks, and removing loops and repetitions preserves reachability. Q has more than one vertex and cannot have exactly two, hence has at least three.

Every vertex r of a nontrivial strong tournament lies in a directed triangle. Its out- and in-neighborhoods A,B are nonempty. If there were no r->a->b->r, every cross edge would point B->A, preventing A from reaching r or B. Choose such a triangle. Its lift gives an induced C3[R,1,1] in T.

For completeness, H is monotone under induced tournament inclusion. For each old Hamiltonian path insert a new vertex v into an available slot: prepend if v beats the first vertex; otherwise append if the last beats v; otherwise the path starts at a vertex beating v and ends at one beaten by v, so some consecutive pair a,b has a->v->b. Fix a deterministic choice of slot. Deleting v recovers the old path, proving injectivity. Iterate to obtain

    H(T) >= H(C3[R,1,1]) >= gamma_p.

This proves (2). This insertion argument proves non-strict monotonicity only; no equality classification is silently extracted from it.

## 3. Complete finite inputs and the new sharp elevenfold consequence

The cyclic gauge is vertices 0,...,p-1, with x->x+d for d in D, where D chooses exactly one of i,p-i for 1<=i<=(p-1)/2. Mask bit i-1 replaces i by p-i. The universe has 2,4,8,32 labelled orientations for p=3,5,7,11 respectively, with no further filters.

All 32 eleven-circulants have the following joint values. Multiplicity counts orientations for the specified cycle, not unlabelled isomorphism classes.

| H(R) | F2(R) | H(C3[R,1,1]) | multiplicity |
| ---: | ---: | ---: | ---: |
| 92411 | 1571526 | 1848759 | 10 |
| 93027 | 1574958 | 1854039 | 10 |
| 93467 | 1577598 | 1857999 | 10 |
| 95095 | 1586310 | 1871595 | 2 |

Thus mu_11=92411, gamma_11=1848759 and mu_11^2+22=8539792943 > gamma_11.

**New specialization.** Every finite strong tournament T with an automorphism of order 11 satisfies exactly one of these size alternatives:

1. |T|=11 and H(T) belongs to {92411,93027,93467,95095};
2. |T|>11 and H(T)>=1848759.

The second bound is sharp. Take R with D={2,3,4,5,10}, mask 1. Then C3[R,1,1] has 13 vertices and H=1848759. It is strong: the cyclic quotient supplies paths between blocks and between any two vertices within R through x,y. Rotation of R, fixing x,y, is an automorphism of order 11. This establishes sharpness without any bound on the order of the ambient tournament.

For comparison, replayed inherited data give (mu_p,gamma_p)=(3,15),(15,135),(175,2163) at p=3,5,7. At p=7, all eight F2 values equal 1638; six cores have H=175 and two have H=189. At p=11, F2 is not constant. The observed matching minimizers of H and F2 at p=11 is a finite fact, not a theorem for all primes.

**Equality scope.** Within the 32 triangle inflations, equality holds precisely for masks

    1, 8, 9, 12, 14, 17, 19, 22, 23, 30.

These ten cores constitute a single multiplier orbit under x->a*x, a nonzero modulo 11. The verifier constructs that orbit explicitly. Therefore their triangle inflations are all isomorphic. No converse assertion that arbitrary circulant isomorphisms must be multipliers is needed.

For an arbitrary ambient equality case H(T)=1848759, the chosen order-11 automorphism must have exactly one moving orbit by (3). Its core must be one of the ten minima, and any triangle-through-core inflation selected in the reduction has exactly the same H as T. These are necessary conditions. This packet does **not** prove that equality forces |T|=13 or classify all possible equality cases. A strict induced-inclusion theorem under appropriate hypotheses would be an additional ingredient, not a consequence of the insertion injection alone.

## 4. Consequences, hostiles, and limits

By Cauchy's theorem, if a strong T has H(T)<1848759 and |T|!=11, then 11 does not divide |Aut(T)|. Equivalently, any fixed-permutation strong histogram at a profile containing an 11-cycle vanishes below this threshold except at the single orbit profile (11). For a general permutation with order divisible by 11, take a power of order 11 before applying the result.

Strongness is essential: placing any number of singleton SCCs before or after a minimum eleven-circulant preserves H=92411 by ordinal-sum multiplicativity, and preserves its order-11 rotation. The verifier checks one such padding at both p=7 and p=11. A strong tournament with one moving p-orbit cannot have exactly p+1 vertices, since its contracted two-vertex quotient would not be strong.

The bound concerns raw H and an actual order-11 subgroup. It gives H/11>=168069 in the non-core case, not H/|Aut(T)|>=168069 for unrestricted automorphism groups. It is a symmetry obstruction, not a spectrum existence or completeness theorem. Chan-Pak completeness and LRC(14) remain OPEN. The general reduction, the sevenfold gap, the two-cover kernel and mu_11 are inherited. The newly supplied content is the exact elevenfold two-cover/inflation table and its sharp all-order strong-gap consequence. No comprehensive literature novelty search has been conducted.

## 5. Reproduction and checks

From the repository root:

    python3 05-knowledge/results/session240_seven_eleven_symmetry_20261110/gaps/probe.py > /tmp/s240-gap.json
    python3 -O 05-knowledge/results/session240_seven_eleven_symmetry_20261110/gaps/probe.py > /tmp/s240-gap-O.json
    cmp /tmp/s240-gap.json /tmp/s240-gap-O.json
    cmp /tmp/s240-gap.json 05-knowledge/results/session240_seven_eleven_symmetry_20261110/gaps/probe.json

Every check raises RuntimeError on failure and survives optimization. Exact integer subset DP enumerates each Hamiltonian path by its visited subset and final vertex. It returns all induced-subset path counts, whose complementary products give F2. A separate DP on each full inflation verifies (1) for all 46 orientations. The code also checks tournament axioms, strongness, rotation invariance, oddness and p-divisibility.

As an algorithmically independent path, uncached depth-first enumeration of directed simple paths checks one core representative of every explicitly constructed multiplier orbit, and directly counts all 1848759 directed Hamiltonian paths of the sharp 13-vertex inflation. It does not import another result or use the two-cover identity or subset-count DP. Explicit multiplier relabellings cover every core, and all rows are nevertheless independently processed by the subset DP. Positive controls are H(transitive three-tournament)=1 and H(C3)=3. Hostiles test non-strong padding. Complete output includes every D, mask, H, F2, inflation value, multiplier orbit and all minimum masks.
