# Response-weighted modular inventories

**Status:** PROVED elementary aggregation and response identities, conditional only on the inherited canonical modular decomposition theorem; FINITE-EXACT executable family census. No claim of a Chan–Pak proof, complete prime catalogue, or global polynomial-time tournament enumeration.

## Inheritance and contribution

The canonical decomposition is already present in `../tournament_composition_classification_20261026/CLASSIFICATION.md` (paths in this paragraph relative to the results directory). The ordinary species equation, prime cycle-index census, and strong-component split are inherited from `a000568_intransitivity_families_20261026/FAMILIES.md`, not new here. The repaired child-colour stabilizer and invariance of the full ordered path-cover profile under cyclic-triple fibres are inherited from `tournament_composition_classification_20261026/fibre/RESULT.md` (F1 and F4). The existing walk-resolvent line includes `01-canon/theorems/THM-3248-q4-paired-owner-stirling-compiler.md`; our innovation is **aggregating whole isomorphism-class inventories by their compositional responses without confusing response equality with actual child isomorphism**.

Anchor: count ordinary classes and Hamiltonian-path distributions together. Niche: carry reciprocal automorphism mass through the same response quotient. Wildcard: a finite weighted transducer computes each response operation, while a separate Burnside inventory accounts for hidden state multiplicities. The canonical hostile is the 81/27 automorphism boundary, not an arbitrary failed numerical example.

Live concepts: canonical strong modules; ordered path covers; Burnside fixed-child counts; labelled species mass; interacting seams. A quotient is useful only while the allowed operation respects its response equivalence.

## 1. Response congruence and separator compiler

Let q_r(T) count ordered r-tuples of nonempty directed paths that partition V(T). Set q(T)=(q_1,...,q_n); q_1=H and q_n=n!. These are NOT the unordered path-cover numbers p_r=q_r/r!.

For a fixed tournament Q on k ports, with adjacency matrix A and D=diag(x_1,...,x_k), define formal series

    R_Q(x) = 1^t D (I-A D)^(-1) 1,
    W_Q(x,z) = z R_Q(x)/(1-z R_Q(x)).                 (1)

For each child T_i define the coefficient contraction L_i(x_i^r)=q_r(T_i) for 1<=r<=|T_i|, and L_i(x_i^r)=0 otherwise, including r=0. This is a linear functional, NOT a ring homomorphism. Then

    sum_l q_l(Q[T_1,...,T_k]) z^l
      = (L_1 ... L_k) W_Q(x,z).                     (2)

**Proof.** A directed path in the substituted tournament splits uniquely into maximal runs within its modules. Its port word follows directed Q arcs, hence is counted by R. An ordered list of l such words is counted by R^l, with z marking the number of global paths. Across distinct global paths, even identical successive port labels remain separate runs. For a fixed port-word list containing r_i occurrences of port i, its runs are filled bijectively by an ordered r_i-path cover of T_i, ordered by their appearance in the global path list. Multiply q_{r_i}, sum, and obtain (2). Conversely split any global cover into its maximal module runs. This proves both directions, with no factorial missing. QED.

Equivalently W=z 1^t D (I-(A+zJ)D)^(-1)1 by the rank-one resolvent identity. Every transition contributes z+A_ij; a diagonal transition MUST be cut, and a legal off-diagonal transition may be cut or joined. This provides the code's three-port dynamic-programming kernel. Formal power-series inverses are legitimate because D has zero constant term. Degree bounds permit exact truncation before contraction.

**Corollary (contextual congruence).** Equality of full q profiles is preserved by every uniform tournament-substitution context. This is stronger than equality of H. No sequence of those contexts can expose hidden state bits encoded by equal-q children. Nonuniform seams are outside this statement.

## 2. Exact ordinary-class pushforward

Let a child family be closed under isomorphism and let m(a) be the number of its isomorphism classes with response a. Let G=Aut(Q). Q must be a prime canonical quotient, so its substitution children are the canonical maximal proper strong modules. (A linear quotient is handled by the ordered strong-component construction below.) Write F_Q(a_1,...,a_k) for (2).

For a target response b, the number of ordinary substitution classes is

    M_Q(b) = (1/|G|) sum_{g in G}
      sum_{a:Cyc(g)->Responses}
        [ F_Q(lift_g a)=b ] prod_{c in Cyc(g)} m(a(c)).       (3)

Here a labels permutation cycles of g; lift repeats that response at each port in the cycle. Sizes can be marked or restricted before summation.

**Proof.** G acts on tuples of ACTUAL child isomorphism types. Canonical decomposition identifies its orbits with ordinary classes of substitutions. The output response is constant on each orbit by relabelling invariance. A tuple fixed by g must assign the same actual child class to every port of each cycle. Once the response a of that class is prescribed, there are m(a) choices per cycle, not m(a) raised to the cycle length. Apply Burnside on the invariant fibre of output response b. QED.

This is not a claim that response buckets classify graphs. It is an exact pushforward of an existing classifier, with ordinary-class multiplicities as an essential sidecar. Using just a cycle index without the cycle's PORT POSITIONS can also lose F_Q information for general Q: permutations of the same cycle type need not induce the same response interaction. Retain actual group action, or a proven finer symmetry reduction.

### Reciprocal-automorphism mass sidecar

Let u(a)=sum_{[T]:q(T)=a}1/|Aut(T)|. Then

    U_Q(b) = (1/|G|) sum_{a_1,...,a_k}
           [F_Q(a_1,...,a_k)=b] prod_i u(a_i).              (4)

For a representative tuple t the decorated substitution has automorphism order |Stab_G(t)| prod_i |Aut(T_i)|. Its orbit has |G|/|Stab_G(t)| ordered tuples, which proves (4) after summation. Consequently n! U(b) counts labelled tournaments in response bucket b. This does not recover the distribution of individual automorphism orders. The 81/27 example shows why that distinction matters.

### Linear roots

Every tournament has a unique ordered list of strong components, all arcs directed forward. Thus use the unambiguous first-component grammar

    A = S + S A,
    S = vertex + sum_{allowed prime Q} Q[A,...,A]/Aut(Q).   (5)

Both ordinary multiplicity and reciprocal-Aut mass multiply at S A and add across alternatives. Nontrivial series children are not unrestricted A at every linear node: that would count different parenthesizations repeatedly.

The response operation for an ordinal sum has the explicit form

    q_k(A -> B) = sum_{i,j,ell: i+j-ell=k}
        q_i(A) q_j(B) k! / ((i-ell)! (j-ell)! ell!).        (6)

Choose i unordered A paths and j unordered B paths, join ell pairs, and order the resulting k paths. This gives (6), proves its symmetry in the two responses, and explains why the graph's ordered component sequence can be forgotten by q but not by class multiplicity.

## 3. A closed infinite-family implementation

Let C be the class of tournaments whose canonical prime quotients are all C3. Include arbitrary ordered strong-component sums and recurse indefinitely. This is an exact infinite family, not a census of all tournaments. Every n has members; many primes of order >=5 are absent.

For ordinary scalar counts A(x) of nonempty objects and S(x) of strong objects,

    A = S/(1-S),       S = x + (A(x)^3 + 2 A(x^3))/3.       (7)

For reciprocal-Aut mass series the corresponding equations are

    A_u = S_u/(1-S_u), S_u = x + A_u^3/3.                 (8)

These inherited species identities are used only as independent scalar checks on the new RESPONSE-refined implementation.

The response F_C3(a,b,c) is symmetric under all S3 permutations. Rotation follows from C3 automorphisms. A transposition reverses the quotient orientation; reversing all global paths, and using q(T*)=q(T), preserves the result. This does NOT enlarge the ordinary graph action C3 to S3.

Therefore enumerate an unordered triple of response types and accumulate:

- a,b,c all distinct: classes 2 m(a)m(b)m(c); mass 2u(a)u(b)u(c).
- a,a,b with a!=b: classes m(a)^2 m(b); mass u(a)^2 u(b).
- a,a,a: classes (m(a)^3+2m(a))/3; mass u(a)^3/3.

Each row contributes to the same computed response F_C3. These formulas contain the essential difference between an observable-state symmetry and true object symmetry. If m(a)=2, the last row is FOUR ordinary classes, not one response orbit and not 8/3. Taking a=q(C3->1)=q(1->C3) recovers exactly the inherited 81/27 boundary.

The implementation dynamically inventories all q buckets at each size, separately for A and S, with both m and u. Canonical-form computation is absent from the production loop. It computes the full H histogram by summing bucket multiplicities with equal first coordinate.

## 4. Efficiency and boundaries

At a k-port prime, replacing child class counts C_n by response counts R_n reduces the tuple engine from class tuples to response tuples. This is an exact structural reduction, not a claim that R_n is polynomially bounded. Coefficient dimensions and integer bit lengths still grow. The automaton kernel has finitely many states at any fixed size bound; the unbounded exact profile is NOT a finite-state object.

For the C3 kernel through total degree N, there are O(N^3) count triples, three last-port states and O(N) separator-degree coefficients: O(N^4) integer cells/operations with the direct recurrence. A single child-triple contraction uses at most O(n1*n2*n3*(n1+n2+n3)) scalar operations. Across the inventory, cost is dominated by the number and sizes of distinct response triples, not the total number of ordinary classes represented. The tables report measured compression, not a comparison with the world's fastest graph enumerator.

**Loss ledger.** Source: canonical decorated module trees. Target: (q,m,u) inventories. Preserved: ordinary class counts, full ordered path-cover distributions, H distributions, labelled counts, and all uniform substitution operations. Lost: recovery of individual graphs, individual automorphism orders, port-rooted endpoint information, adjacency-sensitive cross-module seams. Sidecar for the last loss: rooted response data, developed in the interacting-automata lane of this session. Neither this quotient nor H alone proves Chan–Pak, yields all prime tournaments, or makes a large generated subfamily exhaustive.

## Reproduction

From repository root:

    python3 -B 05-knowledge/results/session268_composable_tournaments_20261111/modular/census.py --max-n 16 --out /tmp/modular-receipt.json

`receipt.json` stores the measured run, all H distributions through the run bound, and full q buckets through n=12. Timing is environment dependent. Runtime checks use explicit exceptions, not assertions, so Python -O keeps them active. `independent_controls.py` additionally compares the actual labelled universe at n<=5 with the inventory buckets, and records a scalar-H loss witness. No historical theorem is wholesale re-audited.

## 5. Exact contextual equivalence (joint session corollary)

The spectrum lane's physical tomography theorem (`../spectrum/PROOF.md` §2) supplies the missing converse. For T of order n, let

    P_d(T)=H(C3[TT_d,TT_d,T]), 1<=d<=n-1.

It proves P_d is linear in q_1,...,q_{d+1}, with coefficient (d!)^2 on q_{d+1}. Thus H(T) together with these n-1 probes determines q(T). Combined with §1:

**For two tournaments of the same order, equality of full q profiles is equivalent to equality of Hamiltonian-path counts in every uniform substitution context.**

The forward implication is our congruence. For the reverse, equality in the displayed contexts gives triangular recovery of q. All probes have at most 3n-2 vertices, so arbitrary context quantification has an explicit bounded separating set. If the first difference is q_r, the d=r-1 context separates by ((r-1)!)^2 times that difference. Credit for physical tomography belongs to the spectrum lane, not to a rerun here.

The independent control finds the first scalar-only loss inside this family at n=6:

    q_A=(9,150,756,1704,1800,720),
    q_B=(9,126,684,1656,1800,720).

Concrete representatives are A=C3[TT_3,1,1]->1 and B=C3->C3. Both H=9, but C3[child,1,1] gives H=177 and 153. Equal scalar H therefore cannot replace the full response in inventory composition, even before introducing nonuniform seams.
