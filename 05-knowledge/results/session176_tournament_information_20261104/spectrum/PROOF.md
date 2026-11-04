# A persistent unit edge gate and a complete nonuniform flow circuit

Session 176. **Proved explicit constructions and exact counting identities; unrestricted spectrum completeness remains open.** No historical novelty or new globally attained count is claimed. The starting material is Session 175's endpoint interfaces, orphan formula and finite-clone obstruction. Ordinary H=7 and H=21 remain globally impossible. The infinite prime and 7p lanes are not closed here.

## 1. An actual unit gate at arbitrarily large orders

Outgoing-neighbor masks define two tournaments on vertices 0,1,2,3,4:

    L = (16,17,3,7,12),    U = (16,25,3,5,12).

Their sole differing edge is 3→1 in L versus 1→3 in U. Replace vertex 0 by a transitive module X of r≥1 vertices, with its external edges inherited from vertex 0. Call the resulting tournaments L_r and U_r. The changed edge is outside X, so they still differ by **one literal edge reversal**, not by a simultaneous switch of r edges.

**Theorem 1.** Both tournaments are strong and

    H(L_r) = 6·2^r−1,    H(U_r) = 6·2^r+1.

Thus the strong |ΔH|≤2 response graph contains a two-state interval at every order r+4, with unbounded center 6·2^r. Examples are {11,13}, {23,25}, {47,49}, {95,97}, {191,193}, {383,385}, {767,769}, {1535,1537}, {3071,3073}.

Strongness is immediate from the shared Hamilton cycle

    0→4→3→2→1→0.

Inflating a vertex of a strong quotient to any nonempty tournament preserves strongness: routes between modules give access to every vertex, and routes leaving and returning to a module connect its vertices.

Here is a complete finite symbolic proof of the count, valid for all r. For any fixed seed Q with marked vertex v, order the other m vertices by a permutation p and choose a nonempty subset of the m+1 gaps to contain runs of X. A gap after a must satisfy a→v; a gap before b must satisfy v→b. Every unoccupied internal gap must satisfy its direct edge a→b. Let a_k count legal permutation/gap patterns with k occupied gaps. Then

    H(Q[X_r]) = Σ_(k=1)^(m+1) a_k · k! S(r,k),                 (1)

where S is the Stirling number of the second kind. Indeed distribute the r labelled, linearly ordered module vertices onto the k ordered nonempty runs; within each run their order is forced. These are exactly the surjections onto k positions, numbering k!S(r,k). This bijection includes every possible departure from and return to the module.

For our two seeds the entire nonzero pattern table is below. Entries are (a_1,a_2); all a_k for k≥3 vanish. All omitted permutations contribute zero.

| Other-vertex order | L | U |
|---|---:|---:|
| 1324 | (0,0) | (1,0) |
| 1342 | (0,0) | (1,1) |
| 1432 | (2,1) | (2,1) |
| 2134 | (0,0) | (1,0) |
| 2143 | (2,1) | (2,1) |
| 2431 | (1,1) | (0,0) |
| 3142 | (2,1) | (0,0) |
| 3214 | (1,0) | (1,0) |
| 3421 | (1,1) | (1,1) |
| 4213 | (0,0) | (2,1) |
| 4321 | (2,1) | (2,1) |
| **Total** | **(11,6)** | **(13,6)** |

There cannot be three module runs: every nonterminal run must be followed by vertex 4, the only outside vertex that v beats. Two such nonterminal runs would repeat vertex 4. Finally 1!S(r,1)=1 and 2!S(r,2)=2^r−2, which prove the theorem.

## 2. Return-to-source flux, and why the unit gate does not recurse freely

The same table works if X is **any** nonempty tournament A. Let h=H(A), and let

    q_2(A) = Σ_(∅≠S⊊V(A)) H(A[S]) H(A[V(A)\S]).

This counts two ordered nonempty path runs covering A. The two-run skeleton coefficients remain six, so

    H(L[A]) = 11h+6q_2(A),
    H(U[A]) = 13h+6q_2(A).                                (2)

Equivalently, defining F_A(2)=2h+q_2(A), these are 6F_A(2)−h and 6F_A(2)+h. The one-pass stream changes by 2h while the entire two-run return stream balances exactly. This is actual tournament realization of a two-register information-flow circuit, not arbitrary selection of a numerical path profile.

**Unit-gate rigidity.** The reversal has ΔH=2H(A). Hence it remains a unit gate exactly when A is transitive. A transitive tournament has H=1. Every nontransitive tournament contains a directed triangle, whose three paths extend injectively under successive vertex insertion, giving H≥3. Thus recursively substituting a nontransitive counting gadget amplifies the step to at least six. The apparent independent binary control is not available.

Even within transitive inflation the identities

    H(L_(r+1))=2H(L_r)+1,    H(U_(r+1))=2H(U_r)−1

must not be advertised as unrestricted binary instructions. The sign is tied to the structural gate state; changing it resets the choice to one of the same two numbers at the new scale. The centers double, so the intervals from Theorem 1 do not overlap across orders.

## 3. Relation to the hard 7p lane

For r=3k with k≥1, the upper endpoint is

    H(U_(3k)) = 6·8^k+1 = 7·((6·8^k+1)/7).

Modulo 49, 8^k≡1+7k, so the displayed quotient is 1+6k modulo 7. Thus v_7(H(U_(3k)))=1 unless k≡1 modulo 7. The k=3 example is 3073=7·439, with prime cofactor 439 (trial division through 19 suffices). The r=3 endpoint 49 instead has valuation two.

This supplies infinitely many strong outputs with odd 7-adic valuation one, but does **not** prove infinitely many prime cofactors, cover every 7p, or cover every prime. It is a finite-template inflation family and therefore remains subject to Session 175's polylogarithmic prime-coverage bound. Its contribution is a persistent unit **edge** response, not a complete spectrum algorithm.

## 4. A growing nonuniform interface circuit that counts all orphans

A separate construction removes the fixed-template restriction. Choose a transitive backbone b_1→...→b_n, k labelled feedback vertices Z with any tournament among them, and orient every backbone-feedback edge independently. The feedback number k may grow. This specifies actual tournaments, including nonuniform interfaces, rather than merely admissible-looking registers.

Every Hamilton path determines a unique permutation (z_1,...,z_k) of Z and k+1 increasing backbone runs placed before, between, and after these vertices. Some runs can be empty. Consequently the following finite-state circuit counts H exactly.

Fix a feedback permutation. Process b_1,...,b_n. Each of the k+1 runs has state:

- 0: empty;
- 1: nonempty, with current last vertex beating the run's right feedback boundary, or with no right boundary;
- 2: nonempty, with current last vertex not beating that boundary.

Start with all runs empty and weight one. For the next backbone vertex b, choose exactly one run j. If its state is empty and it has left boundary z_j, require z_j→b. Append b and replace that run's state by 1 or 2 according to its right boundary. All other states are unchanged. An occupied run can always receive a later backbone vertex, even if its old state was 2, because the backbone is transitive and its previous last vertex ceases to touch the right boundary.

At the end discard states containing 2. For each empty internal run between z_j,z_(j+1), require z_j→z_(j+1). Empty boundary runs impose no condition. Sum surviving weights over all k! feedback permutations.

**Theorem 2.** This circuit returns the exact Hamilton-path count, using at most k!·3^(k+1) states and O((n+1)·k!·(k+1)·3^(k+1)) integer additions.

**Proof.** Given a Hamilton word, all vertices in one backbone run increase, and assigning each backbone vertex to that run reconstructs the unique processing history. Starting a run checks its left edge; its final state checks its right edge; an empty internal run checks the feedback-feedback edge. These are all edges in the resulting word. Conversely every accepted assignment constructs a valid Hamilton word and has a unique feedback permutation. The two maps are inverse. The state bound is direct. ∎

Arithmetic cost excludes integer bit complexity. The extra 1 covers final-state summation, including an empty backbone; final feedback-edge legality tests and adjacency preparation fit within the same stated operation bound under constant-time adjacency access. k! separate copies are a deliberately simple circuit, not a claimed minimal automaton. With k=n_total−1 and a one-vertex backbone this covers all tournaments, so expressiveness alone is not a compression or completeness theorem.

Crucially no step assumes that deletion of a feedback vertex leaves a Hamilton path. Consecutive feedback vertices and all possible splits of backbone runs are included. Paths that are orphans under deletion are counted automatically. This fixes precisely the loss identified by Session 175's insertion examples.

## 5. What would close the remaining gap

The persistent gate proves that single-edge ±2 control can coexist with unbounded strong counts. The circuit proves exact realizability and complete counting for genuinely growing nonuniform interfaces. Neither supplies a path through response components with growing overlapping H-intervals. To complete the spectrum by this route one still needs a structurally closed collection of circuit states whose realized response components overlap and contain every required prime and 7p target. An abstract register update is insufficient, and substituting arbitrary modules into the unit gate demonstrably fails.

## 6. Reproduction and scope of checks

Run `python3 check.py` and `python3 -O check.py`; the JSON outputs agree exactly. No software installation or network data is needed.

- All 2,730 two-feedback signatures with backbone sizes 0 through 5 are checked against subset endpoint DP. A further 64 seeded three-feedback signatures are checked, for 2,794 streaming comparisons total.
- 682 of the two-feedback signatures, all with backbone size at most 4, are independently counted by literal vertex permutations.
- Formula (1) is checked 1,024 times: every labelled four-vertex seed, every marked vertex, and each transitive module size 1 through 4.
- Formula (2) is checked for all 75 labelled nonempty tournaments of orders at most four as substituted modules, including strongness of both outputs.
- The two gate profiles are reconstructed by exhaustive legal skeleton patterns. Strongness, sole edge difference, and the closed count are checked independently for each r=1,...,9.

The unrestricted statements rely on the bijective proofs, not extrapolation from these finite checks. `explore.py` includes the deterministic profile search that first found the gate. No claim is made that its finite search is an unrestricted spectrum census.
