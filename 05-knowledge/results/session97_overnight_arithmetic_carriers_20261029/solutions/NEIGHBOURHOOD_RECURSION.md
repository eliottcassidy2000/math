# Tournament solutions: well-founded recursion still needs attachment maps

**Status:** CITED classical construction; PROVED elementary dependency-depth and transport deductions below; FINITE-EXACT hostile/control checks. No literature-priority claim. This is a decision-carrier result, not an unrestricted enumeration formula.

## Source recovered

Scott Moser and Daniel Allcock, *Binary Relations from Tournament Solutions, and Back Again*, version April 25, 2013, [author-hosted primary PDF](https://web.ma.utexas.edu/users/allcock/research/tournaments.pdf), Definitions 2.2–2.3 and 3.1, Propositions 1–2. The author's page lists it as submitted; no publication status beyond that is inferred.

Their input is an isomorphism-equivariant neighbourhood rule assigning a proper vertex subset N_y to each vertex y. For a tournament solution S, form the challenge relation x→y exactly when x belongs to S(T[N_y]). Its top-set is the union of source strongly connected components. Repeating this operator gives a unique fixed solution, independent of the initial solution. Proper-subset descent supplies the induction. Choosing the in-neighbourhood gives the tournament equilibrium set. General neighbourhood rules need not satisfy Condorcet criteria; the paper states additional sufficient hypotheses.

## A sharper computational parameter (elementary deduction)

Fix the rule N on every induced tournament encountered. Set h(empty)=h(singleton)=0 and, for larger T,

    h(T) = 1 + max({0} union {h(T[N_y]) : y in V(T)}).

Let F_N be the challenge/top-set operator. For every initial nonempty, isomorphism-equivariant tournament solution S, F_N^h(T)(S)(T) is independent of S and is the unique recursively defined value. In fact every iterate k>=h(T) has this value. Here exponent h(T) means the number of global operator iterations before evaluating at T, not h iterations independently restarted at each child.

**Proof.** Empty and singleton answers are forced. Each child has depth at most h(T)-1. At iteration k>=h(T), induction fixes every named child winner set at iteration k-1, hence fixes every challenge edge, and hence fixes its top-set. The same induction defines and proves uniqueness of a fixed solution. Since each child is proper, h(T)<=|T|-1 for nonempty T. This is a dependency-DAG improvement over a uniform size bound, not a claimed sharp runtime lower bound.

An implementation can evaluate only the reachable induced-subset DAG, memoizing repeated states. If that DAG has D labelled states of size at most n, assembling adjacency and computing SCCs costs O(D n^2), excluding the neighbourhood-rule cost. This is a conditional algorithmic bound: D may be exponential. It is not an A000568 speedup or a polynomial-time TEQ claim.

## Five-vertex hostile: every child type agrees, final answer sizes differ

Take the rigid tournament on 0,...,4 with edges

    0→1, 2→0, 0→3, 4→0, 1→2,
    1→3, 1→4, 2→3, 2→4, 3→4.

Choose singleton neighbourhoods, written as the selected vertex for y=0,...,4:

    A = (2,0,1,0,3),
    B = (4,0,1,2,3).

Every selected vertex beats its owner. Every child tournament is the same singleton, every local solution returns its only vertex, and all five child-size and child-isomorphism-type rows agree. Nevertheless A assembles a source directed triangle 0→1→2→0 with outward tail 0→3→4, while B assembles the source directed five-cycle 0→1→2→3→4→0. Thus the final selected sets are {0,1,2} and {0,1,2,3,4}. Even the output cardinality is not recoverable from the child-type bank.

This is a valid pair of global equivariant rules, not merely arbitrary named data: the host has trivial automorphism group, so transport A or B to its entire isomorphism class unambiguously, and use ordinary in-neighbourhoods on every other class. Both global rules use proper subsets. The example compares two rules whose erased local type data agree; it does not claim TEQ itself has two answers on this host.

## The exact retained carrier and relation to prior work

For each child canonical representative C_y, retain an actual isomorphism

    i_y: C_y → T[N_y].

Cache the winner subset W(C_y), then assemble challenge edges from

    i_y(W(C_y)) × {y}.

A different canonicalizing isomorphism changes i_y by a child automorphism; equivariance makes W(C_y) invariant under that automorphism. Consequently this transported set is well-defined. Thus caching is safe **with** embeddings, and demonstrably unsafe if only child class and winner count remain. The exact sufficient parent-level decision carrier is the assembled challenge relation; it need not reconstruct all tournament edges. No claim of universal minimality is made.

Source→target: embedded induced tournaments → challenge graph → selected vertices. Preserved predicate: fixed-solution membership. Destroyed information: unselected edges may disappear. Required sidecar: each child's named winner-image, or an embedding determining it. Cheapest hostile: the five singleton children above.

Targeted repository inheritance at source pin `5d6eb91dd76eee86b91e04d856daa32c32b3d61d`:

- `05-knowledge/results/joint_port_germs_20261028/PORTS.md`, J1: a tournament register encodes an ordered tuple of ports; the actual attachment names matter. The present construction gives a different, classical recursive decision operation for which those embeddings cannot simply be dropped.
- `05-knowledge/results/balanced_choice_network_quotients_20261021/PROOF.md`, §§1–2,6: a particular choice-network quotient preserves selected observables but is explicitly not the physical tournament metagraph quotient. Here likewise the solution quotient is operation-specific; a social-choice solution is not a metagraph equivalence class.
- Filename search for retentiveness/tournament-solution/choice-function/revealed-preference routes found no direct prior solution-theory packet. This is a bounded search report, not proof the repository never treated this literature.

The useful unifying pattern is **well-founded recursion plus equivariant transport**: recursion eliminates dependence on arbitrary seed data, but does not eliminate dependence on attachment maps. Convergence and informational sufficiency are separate obligations. Nothing here identifies Leech lattices with tournament solutions, nor equates a source SCC with a black metagraph core.

## Reproduction

Run `python3 check.py` and `python3 -O check.py`; both reproduce `receipt.json`. Standard-library checks exhaust all 1+2+8+64+1024=1099 labelled tournaments of orders 1–5 for the ordinary in-neighbourhood rule, two initial solutions, and the depth bound. The separate hostile check verifies the host's rigidity, legal in-neighbour selections, and final sets of sizes three and five. Explicit exceptions remain active under `-O`. The all-order conclusions rest on the proofs, not this finite census.
