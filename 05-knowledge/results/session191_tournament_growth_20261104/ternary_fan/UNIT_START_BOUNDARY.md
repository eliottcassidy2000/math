# The unit-start interface is intrinsically a one-cycle interface

**PROVED, with an exhaustive five-name certificate.** Every strong tournament of order at least three with a unit-start vertex is uniquely Hamiltonian. Moreover, no tournament with exactly thirteen Hamiltonian paths, at any order, has a unit-start or unit-end vertex. The latter is an actual attainable Hamilton-path count, not one of the global gaps 7 and 21.

## 1. The cycle budget of a marked start

Let hc(T) count directed Hamiltonian cycles modulo cyclic rotation, and s_v(T), e_v(T) count directed Hamiltonian paths starting and ending at v. For every tournament,

    hc(T) <= min(s_v(T),e_v(T))  for every vertex v.       (1)

Indeed, cut a Hamiltonian cycle immediately before v to obtain a path starting at v, or immediately after v to obtain one ending at v. Each map is injective: the removed arc is the unique closing arc of that path. For a strong tournament of order at least three, Camion's theorem supplies at least one Hamiltonian cycle. Thus s_v=1 or e_v=1 forces hc=1. More generally, a start budget b forces hc<=b.

Consequently every strong output of the two-port fan that retains its unit-start vertex lies in the uniquely Hamiltonian stratum. This does not make the construction invalid; it identifies the precise structural class controlled by its interface.

## 2. All strong realizations of H=13

Up to isomorphism there is exactly one strong tournament with thirteen Hamiltonian paths. It has five vertices and outgoing rows

    D = (2,12,9,17,7).

Its two Hamiltonian cycles, each written from vertex zero, are exactly

    0 -> 1 -> 2 -> 3 -> 4 -> 0,
    0 -> 1 -> 3 -> 4 -> 2 -> 0.

Its start and end vectors are

    s=(2,2,2,3,4),       e=(4,3,2,2,2).

The verifier enumerates all 1,099 labelled tournaments of orders one through five. Direct Hamilton-word enumeration and a separate predecessor-subset recurrence agree on every H and end-count vector. The displayed core start vector is counted directly; padding start vectors are checked separately by applying the recurrence to the converse tournament. Precisely 120 labelled five-name tournaments have H=13; exhaustive canonicalization over all 5! permutations identifies exactly the displayed class. Direct cycle enumeration fixes vertex zero, so it counts cycles once modulo rotation, not once for each starting vertex. The explicit thirteen paths and both cycles are generated in `thirteen_receipt.json`.

For larger strong orders, Busch's bound H>=5^((n-1)/3) excludes H=13 as soon as n>=6, since

    5^5 = 3125 > 13^3 = 2197.

Thus the five-name census and this strict cutoff prove the complete all-order strong classification. This use of Busch is cited, not a new proof of his lower bound.

## 3. Nonstrong padding cannot restore a unit endpoint

The strong components of a tournament are totally ordered: C_1->...->C_t. Every Hamiltonian path traverses them in that order, and restriction/concatenation gives

    H(T)=product_i H(C_i).

A nontrivial strong component has at least three Hamiltonian paths, by cutting one Hamiltonian cycle at each vertex. Since 13 is prime, an H=13 tournament therefore has exactly one nontrivial strong component, with H=13, and every other strong component is a singleton. By Section2, every such tournament has the exact form

    TT_a -> D -> TT_b,       a,b>=0.                  (2)

This is an all-order structural classification, not merely a set of sampled paddings.

If a=0, the nonzero start counts are exactly the five coordinates of s above; the b suffix vertices have start count zero. If a>0, the first source is the only possible starting vertex, and its start count is 13. Dually, if b=0 the nonzero end counts are the five coordinates of e above; if b>0 the last sink is the only possible endpoint, with count 13. Hence no tournament of count13 has any endpoint coordinate equal to one.

The verifier separately checks the stated start and end formulas on all sixteen paddings 0<=a,b<=3. Those checks support the implementation; the component argument proves (2) and the endpoint conclusion for unbounded a,b.

## 4. Exact implication for a global construction

A proposed construction whose *every realized tournament* is required to carry a unit-start or unit-end mark cannot realize the full odd spectrum outside {7,21}, because it cannot realize13. Keeping strong unit-start marks throughout a recursive grammar confines its strong states of order at least three to hc=1.

This does not disprove a grammar with multiple state types, a separate non-unit seed13, transitions that temporarily or permanently leave the unit-start class, or a construction with larger cycle/start budgets. It also does not assert that unit-start absence characterizes either global gap. The practical consequence is that a global proof needs an interface broad enough to cross cycle strata, rather than treating a preserved unit port as a universal invariant.

## 5. Reproduction and source scope

    python3 verify_thirteen.py
    python3 -O verify_thirteen.py

The script is self-contained, uses only the Python standard library, and takes no external data input. Both modes retain the same receipt and byte-identical stdout. Checks use explicit exceptions, not removable assertions.

The cycle-cut injection and component argument are elementary; no novelty is claimed for them or against existing small tournament classifications. Camion's Hamiltonicity theorem and Busch's lower bound are classical inputs. The latter is A. H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC13(2006), N3, https://doi.org/10.37236/1141, used with the same source scope as the Session190 and second-generation packets. The contribution here is the exact interface boundary, its all-order thirteen-path certificate, and the clarified consequence for this session's construction strategy.
