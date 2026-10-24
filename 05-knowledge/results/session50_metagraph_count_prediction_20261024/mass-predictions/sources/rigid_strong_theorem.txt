# Finite pure-blue atoms and the complete rigid-strong alphabet

Status: PROVED deductions from the CITED published theorems and the strong-cycle reduction; FINITE-EXACT circulant exclusion through order31. See AUDIT.md for independent checks. The complete pure-blue strong alphabet is not resolved by this packet; a separate classification may strengthen these intermediate bounds.

## 1. Objects, inherited mechanism, and exact count conventions

H(T) counts directed Hamiltonian vertex words. HC(T) counts directed Hamiltonian cycles up to cyclic rotation, not rooted cycles. Aut(T) acts freely on H(T), so the fixed-path tiling count is tc(T)=H(T)/|Aut(T)|. A strong tournament means strongly connected; the singleton is included as a trivial strong atom where stated.

The strong-cycle note [strong-cycle reduction](../merged_metagraph_purity_and_strong_core_20261020/PURE_BLUE_STRONG_CYCLE_REDUCTION.md) proves that every nontrivial strong pure-blue T has odd order and is a regular circulant; moreover the successor of EVERY Hamiltonian cycle is in Aut(T). That result uses purity of all Hamiltonian-path reversals, not just one chosen cycle. Our closest inherited canon is THM-2454, [THM-2454 center law](../../../01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md), whose center law is proved but whose alphabet was finite-exact only through13. THM-2477, [THM-2477 earlier rigidity result](../../../01-canon/theorems/THM-2477-rigid-strong-classification-through-sixty-two.md), previously left a rigidity gap above62. Section3 closes that gap without using its nonprimitive group argument.

Hostile control: regular cyclic C5 has H=15, Aut=5, tc=3 and is pure-blue, so nonrigidity is not nonpurity. Nontrivial strong pure-blue and strong rigid are different predicates.

## 2. Universal finiteness, with no invented asymptotic cutoff

CITED: Bill Cuckler (single author), *Hamiltonian Cycles in Regular Tournaments*, Combinatorics, Probability and Computing16(2) (2007),239–249, DOI https://doi.org/10.1017/S0963548306008170. The published abstract at https://www.cambridge.org/core/journals/combinatorics-probability-and-computing/article/abs/hamiltonian-cycles-in-regular-tournaments/5CB605A35014A30F60F740D74199D3E7 states a uniform lower bound HC(T) >= n!/(2+o(1))^n for every regular tournament. The indexed author manuscript Theorem1 states: for each epsilon>0 there is N(epsilon) such that all regular T of order n>N(epsilon) satisfy HC(T)>=n!/(2+epsilon)^n. Its N is not explicitly supplied here.

CITED: John D. Dixon, *The Maximum Order of the Group of a Tournament*, Canadian Mathematical Bulletin10(4) (1967),503–505, DOI https://doi.org/10.4153/CMB-1967-048-9. Theorem1 in the primary PDF proves |Aut(T)| <= 3^((n-1)/2) for every tournament. Primary PDF: https://www.cambridge.org/core/services/aop-cambridge-core/content/view/682531E83B0A3AA042CFC5852474E914/S0008439500067461a.pdf/div-class-title-the-maximum-order-of-the-group-of-a-tournament-div.pdf . Its proof invokes Feit–Thompson through solvability of odd-order groups.

For strong pure-blue T, the audited injection gives

    n!/(2+o(1))^n <= HC(T) <= |Aut(T)|-1 <= 3^((n-1)/2)-1.

Fix epsilon=1. For n>N(1), n!/3^n grows faster than 3^((n-1)/2), giving a contradiction for all sufficiently large n. Thus the orders of strong pure-blue tournaments are bounded. Since there are finitely many tournaments at every finite order, there are finitely many strong pure-blue isomorphism classes. This argument is qualitative; it does not turn the unnamed N(1) into a numerical cutoff.

### 2a. Inverse-pair sharpening

For any tournament T of order n>=3, let F be the set of automorphisms that are a single n-cycle. Inversion pairs F into distinct pairs {g,g^{-1}}. The arc between v and g(v) has the same orientation for every v because g is a transitive automorphism. Exactly one of g and g^{-1} therefore gives a directed Hamiltonian successor. If T is pure-blue, every Hamiltonian successor is an automorphism, so

    2 HC(T) = |F| <= |Aut(T)|-1.

This equality does not claim that every automorphism is an n-cycle. The injection plus pairing is stronger than the bound needed for qualitative finiteness.

### 2b. An explicit alternative from Thomassen, not from Cuckler's asymptotic

CITED primary publisher abstract, retrieved through its indexed page: Carsten Thomassen, *Hamilton Circuits in Regular Tournaments*, in *Cycles in Graphs*, North-Holland Mathematics Studies115 / Annals of Discrete Mathematics27 (1985),159–162, DOI https://doi.org/10.1016/S0304-0208(08)73005-6, primary page https://www.sciencedirect.com/science/article/pii/S0304020808730056/pdf?md5=2edc6467c1fcfb38ffb347da2a134065&pid=1-s2.0-S0304020808730056-main.pdf . The exact statement in the publisher abstract is: every strong tournament of minimum outdegree at least3k+3 has at least4^k k! distinct Hamiltonian circuits. The full proof was not retrieved here; this is a cited published theorem, not an independently reproved theorem. The exponent printed as uppercase K in the indexed abstract is corroborated as lowercase k in the standard survey statement; there is only one bound parameter.

For a pure-blue strong tournament of odd order n>=13, put k=floor((n-7)/6). Then k>=1 and (n-1)/2>=3k+3, so

    2*4^k*k! <= 2 HC(T) <= 3^((n-1)/2)-1.

If n belongs to {6k+7,6k+9,6k+11}, then (n-1)/2<=3k+5. At k=21,

    2*4^21*21! > 3^68.

The ratio of the left/right expression when k is increased by1 is 4(k+1)/27>1 for k>=21. Hence every odd n>=133 is impossible. With the audited odd-order reduction, every strong pure-blue atom has order at most131. The same exact bound separately excludes n=121,127,129. The slightly weaker original injection gives order at most137 and separately excludes127,133,135. These are elementary exact-integer comparisons, not asymptotic numerology.

The effective bound is useful only after the Thomassen statement and its count convention are accepted as a published input. Qualitative finiteness in Section2 is independent of this extra source.

## 3. Rigid strong tournaments are completely classified

**Theorem.** The only nontrivial strong tournament with H(T)=|Aut(T)| is C3. Including the singleton, the rigid strong alphabet is exactly {1,C3}.

Proof. Since Aut(T) acts freely on Hamiltonian paths, equality says that the Hamiltonian paths form one orbit. Camion's theorem supplies a Hamiltonian cycle. The path obtained by deleting one cycle edge closes back into a cycle. Whether a Hamiltonian path's last vertex beats its first is automorphism invariant, so EVERY Hamiltonian path closes. In particular T contains no Hamiltonian bypass: a Hamiltonian path whose first vertex beats its last.

We use the following published tournament result, not the more general digraph result or its exceptional digraph classification.

CITED: Abdelhamid Benhocine and A. Pawel Wojda, *On the existence of specified cycles in a tournament*, Journal of Graph Theory7(4) (1983),469–473, DOI https://doi.org/10.1002/jgt.3190070413. Primary publisher abstract: https://onlinelibrary.wiley.com/doi/abs/10.1002/jgt.3190070413 . For every tournament of order n>=7 and every p in {2,...,n}, it guarantees D(n,p), defined by reversing p-1 consecutive edges of a directed n-cycle. For p=2, this is exactly a Hamiltonian bypass. Thus n>=7 is impossible.

For the remaining small orders we avoid any exception list. Let C=(v_0,...,v_{n-1}) be a Hamiltonian cycle. Its consecutive cuts are Hamiltonian paths and lie in one Aut orbit. The map sending one cut word coordinatewise to the next is exactly v_i -> v_{i+1}, so the full n-cycle successor is an automorphism. A tournament automorphism has odd order (an even-order element powers to an involution that reverses the edge of a transposed pair). Consequently n is odd, excluding4 and6; there is no strong order2 tournament. This successor identifies T with a circulant, so it is regular.

At n=5, the four skew connection sets on Z/5Z are {1,2}, {1,3}, {2,4}, {3,4}. Multiplication by a nonzero residue identifies all four with {1,2}. In that presentation the word (0,2,3,4,1) is a Hamiltonian path, with steps2,1,1,2, and its first vertex0 beats its last1. This is a bypass, again impossible. Finally C3 has exactly three Hamiltonian paths and three automorphisms, so it is rigid. QED.

This theorem closes the classification left open in THM-2477, by an existence obstruction rather than a comparison of lower path counts and upper group orders. It does not classify pure-blue strong atoms: C5 remains the essential nonrigid pure-blue witness.

### 3a. Independent checks of the potentially confusing order-five names

`probe_cycles.py` exhausts all120 vertex permutations. The regular cyclic tournament has scores(2,2,2,2,2), H=15, Aut=5, HC=2, and five bypass paths. The exceptional tournament T(5) in Darbinyan's later bypass theorem has arcs

    0->1,1->2,2->3,3->0,0->4,2->4,4->1,4->3,0->2,1->3.

It has scores(3,2,2,1,2), H=15, Aut=3, HC=3, and no bypass. These are not the same tournament. None of the proof above depends on Darbinyan's exceptional theorem.

For provenance only, Darbinyan's arXiv1405.0002 Theorem12 appeared at CSIT2015, DOI10.1109/CSITechnol.2015.7358250. It is also restated as Theorem5.4 in Darbinyan–Karapetyan, Discrete Applied Mathematics216 (2017),537–549, primary PDF https://noad.sci.am/Content/136097/1-s2.0-S0166218X16300610-main.pdf . The broader theorem is not needed here.

## 4. All-order counting consequence, without guessing the remaining atoms

Let b_s count isomorphism classes of pure-blue strong tournaments of order s, including b_1=1, and set B(z)=sum b_s z^s. By Section2, B is a finite polynomial; by the audited parity theorem it has only odd powers. The classification in Section3 gives the rigid-strong polynomial R(z)=z+z^3 exactly.

By THM-2454's proved center law, every pure-blue tournament is uniquely an ordered stack consisting of a left word of rigid strong atoms, its converse-reversed right word, and an optional central pure-blue strong atom. Thus, including the empty object,

    P(z) = (1+B(z))/(1-R(z^2)) = (1+B(z))/(1-z^2-z^6).

This is a global rational ordinary generating function even before B is explicitly classified. It is a statement about isomorphism classes, not labelled counts. The denominator is now known exactly; only the finite numerator remains unknown. In particular all even-order pure-blue counts are the Narayana-type values c(m) at order2m, where

    c(0)=1, c(m)=0 for m<0, c(m)=c(m-1)+c(m-3) for m>=1.

Equivalently c(0)=c(1)=c(2)=1. The odd count at2m+1 is sum over odd s of b_s*c((2m+1-s)/2). Since deg B<=131, the full counts obey p_n=p_{n-2}+p_{n-6} for every n>=132. The even formula holds at every even order, not merely eventually; the new eventual odd recurrence begins at n=133.

The same rigid-atom classification independently gives the OGF 1/(1-z-z^3) for ALL rigid tournament classes (not necessarily self-converse): strong components form an arbitrary ordered word in {1,C3}, and H and Aut multiply over the uniquely ordered strong decomposition. Pure-blue even classes are the palindromic subfamily, not all rigid classes.

## 5. Finite-exact and symbolic obstructions

The exact experiment examines every labelled circulant presentation on Z/nZ for odd n=3,5,...,15, namely all2^((n-1)/2) skew connection sets. There is no isomorphism filter and no sampling. At each n>=7, every presentation has a non-automorphic Hamiltonian-cycle successor among the first two cycles found by lexicographic DFS anchored at0. A bad positional reversal of a cut path is then supplied and independently checked directly from connection differences. An additional guarded exhaustive sweep covers every skew connection set at each odd n=17,19,...,31; each also has a certificate among its first two cycles. Its complete order-wise universes contain256,512,1024,2048,4096,8192,16384,32768 presentations. Thus the complete strong pure-blue alphabet through31 is {1,C3,C5}, using the all-order circulant reduction rather than a full tiling-cube search. The n31 run completed with a2,000,000-node per-presentation cap and45-second whole-run cap; all presentations completed, maximum nodes1,937,950. An earlier200,000-node attempt stopped without certifying the whole order, and that incomplete result is retained separately rather than silently relabelled. All final certificates are checked by a separate verifier. Positive controls n=3,5 exhaust all cycles and have no bad successors.

For the entire carousel family S={1,...,m} on Z/(2m+1)Z, m>=3, there is a symbolic obstruction:

    C=(0,1,...,m-2, m,m+1,...,2m-2, 2m,m-1,2m-1).

This is a Hamiltonian cycle: its increments belong to {1,2,m}. The successor sends the arc0->m-1 to1->2m-1, whose difference2m-2 exceeds m, so the latter is not an arc. Hence every carousel of order at least7 is impure. The restriction m>=3 is necessary; C5 is pure.

The independent subset-DP controls reproduce carousel HC at n=3,5,7,9,11,13,15 as1,2,17,222,5109,166562,7604433. The lexicographic C3[C3,C3,C3] control has HC=207 and Aut=81; it violates even the weaker cycle-count inequality. These examples show that vertex transitivity, circulancy, regularity, and large automorphism groups are not sufficient for purity.

## 6. Connection contract and residual question

Source: all Hamiltonian paths and their reversal symmetries. Map: take all cuts of each cycle, compose successive reversals, and retain the successor permutation. Preserved predicate: purity implies every successor is an automorphism. Lost information: arbitrary nonclosing Hamiltonian paths. Restoration sidecar: explicit bad path reversal or a bypass for the rigid one-orbit case. Decisive tests: C5 distinguishes pure from rigid; a carousel witness distinguishes circulant from pure; the nonregular Darbinyan exception distinguishes two order-five names.

The residual of this packet is classification of possible odd pure-blue strong atoms beyond31; this is a boundary of the present methods, not a claim about the status of separate classification work. With the CITED bound in Section2b, the only possible new orders are the47 odd sizes33,35,...,119,123,125,131. In particular B(z)=z+z^3+z^5 plus unknown terms supported on those sizes. The theorem does NOT assert that {1,C3,C5} is the complete pure-blue alphabet. The short carousel splice does not yet extend to arbitrary skew connection sets.

Reproduce: `python3 probe_cycles.py`. Output and small per-presentation witnesses are in `probe_cycles.out` and `cycle_certificates.json`; extended witnesses are `circulant_n17_certificates.json.gz` through `circulant_n31_certificates.json.gz`. The guarded extension commands are listed in README.md. No claim relies on a floating-point computation.
