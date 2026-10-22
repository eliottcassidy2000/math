# Independent hostile audit: fixed feedback strata

Verdict: PASS for the all-order arguments within their stated scope. One recommended wording clarification about unconditional endpoint probabilities is recorded below. Review date: 2026-10-22. No repository or Git writes.

## Bound sources

- `session26-fixed-apex/proof.md`: SHA256 `5f202c71ac19cbea002fce39425529d04e8dd4d1a2a5c04196355b4b681ba458`.
- `session26-fixed-apex/check.py`: SHA256 `01d366ef41807407ede0823a442610a3bc2f05ea0cedc4b9390f0d6c2af1be1d`.
- `session26-fixed-apex/results.json`: SHA256 `65e8ee9f63e33e40e8395c9dd7010cd8a7e7550e8aa90a6086c27bd21623f991`.
- Attachment-mean dependency: `session26-attachments/proof.md`, SHA256 `8d4f75f0bf851b9544e046f6dc8e2b11dda7fedc71e1980d9dc716eede24587d`. This review directly checked the needed maximal-core-run argument and leading coefficient; the dependency's other conclusions are not needed here.

## All-order proof audit

1. A competing deletion set differs in exactly r core/base vertices. Its retained core vertices have mutually independent cross-edge restrictions to the m-r remaining base vertices. Each has exactly m-r+1 allowed insertion profiles among 2^(m-r) possibilities. Required core-core consistency can only lower the probability. Counting exchanged subsets yields precisely the displayed union bound. For fixed k its r=1 term is k m^2 2^(1-m); every r>=2 term is O_k(m^(2r)2^(-rm)), hence lower order. This is not a bound uniform in growing k.

2. The property of having a unique deletion k-set and distinct attachment profiles is invariant under isomorphism. A generic automorphism must preserve that set, then fix the transitive base pointwise, then fix its distinct-profile core. Thus the generic ordinary class multiplicity is exactly k!, and every exceptional class consumes at least one exceptional presentation. Collision probability is at most binom(k,2)2^(-m). These facts justify class asymptotics without assuming all classes rigid or all deletion sets unique. Every tau<k class is exceptional when m>=1: extend a smaller feedback set to distinct k-sets.

3. For a fixed ordered pair of core vertices, the m two-step witnesses use disjoint edge pairs, giving failure probability (3/4)^m. A union bound puts all core vertices in one SCC except with probability at most k(k-1)(3/4)^m. On that good event the two endpoint conditions are necessary and sufficient for strong connectivity. Their probability is (1-2^(-k))^2 UNCONDITIONALLY when m>=2. Consequently the difference between this probability and the strong probability is bounded by the core-SCC failure probability. This is enough; no conditional independence is needed or assumed. Recommended clarification: explicitly say “unconditionally” in the source and include this absolute-difference inequality. The k=1 case has no core-SCC error. The endpoint statement deliberately excludes m=1.

4. A generic anti-isomorphism preserves the unique core set and reverses the unique base order. Its square is an automorphism and so is identity. Its action on the km cross coordinates is therefore an involution, each orbit requiring opposite bits; fixed coordinates are impossible. At most 2^(km/2) cross assignments exist for each of O_k(1) core possibilities. Exceptional self-converse classes are bounded by the exceptional-presentation count. Dividing both bounds by L/k! tends to zero. Strong class counts have the same positive fixed-k leading scale, so the conclusion survives restriction to strong classes. Converse pairing then gives exactly (N+SC)/2.

5. The action on choices of deleting k-set and its ordered core labels is free: a stabilizer fixes the core pointwise and the ordered transitive complement pointwise. Hence k! f(T)/|Aut(T)| is the exact presentation count, also for nonrigid classes and nonunique sets. A Hamiltonian word determines a core permutation and an assignment of every base vertex to its k+1 intervening slots. Base order in each slot is forced, proving H<=k!(k+1)^m. The excess weighted mass is sum_(f>1)(f-1)H/|Aut|, bounded above using the number of nonunique presentations and this H bound. No exceptional-class multiplicity is discarded.

6. The attachment expectation is unconditioned. Maximal core runs with r=k give the top exponential coefficient k!/4^k; all remaining channels have base at most k. Averaging over the finitely many core tournaments preserves the uniform O_k(k^m) error. Multiplying by L/k! and then subtracting the multiplicity excess gives exactly the displayed V coefficient and relative error. It does not prove a strong-conditioned Hamiltonian-mass formula and does not replace each class's H with an average.

7. Blue mass bounded by SC times max H has relative error O_k(2^(-km/2)+m^2 2^(-m)). In the inherited tiling convention H/|Aut| is the number of Hamilton-word orbits, since the action on words is free. Reflection's fixed tilings can only come from self-converse classes and are no more numerous than their H/|Aut| counts, so the stated rough bound is sufficient. Nonfixed reflection orbits have size two. Once the selected stratum has no internal edge, these give the boundary count (V-B)/2. This last interpretation uses the established reflection/incidence convention, not an independently changed graph definition.

8. The all-edge lower inequality also has a direct check: partition a Hamilton word into floor(n/3) disjoint consecutive triples. The two consecutive arcs stay forward in both tournaments; the third flips. Exactly one of the two tournaments has a cyclic triple on each block. Their two feedback sets together must hit all blocks, giving tau+tau'>=floor(n/3). Thus two tau<=k endpoints are impossible for n>=6k+3, including after isomorphism/converse identification. This is an all-order statement, not inferred from the small checks.

9. For the upper inequality, the exceptional consecutive Hamilton edges induced on a maximum transitive subtournament form a linear forest. A vertex cover of each component with ell vertices uses floor(ell/2) vertices, so their total is at most floor((n-tau)/2). Off this cover, all remaining edges reverse and the retained subtournament is transitive. The resulting deletion count is tau+floor((n-tau)/2)=floor((n+tau)/2). This works for every path and tournament, including empty forest and singleton components; sharpness is not asserted. Combining with tau<=k gives the stated neighbor interval.

## Independent executable checks

`independent.py` uses adjacency bitsets, source-removal transitivity DP, Hamiltonian subset DP, and full permutation canonicalization. It does not import the source checker. Per-class exact multiplicities, rather than merely their aggregate mass, are verified for every presentation with 2<=n<=5 and 1<=k<n. Each presentation checks H's upper bound; each stratum checks uniqueness and profile-collision probability bounds and the weighted-mass correction inequality. This includes the nonrigid C3 and nonunique deletion witnesses.

The second enumeration fixes the directed Hamiltonian path 0->1->...->n-1 and enumerates every choice of nonconsecutive arcs through n=7. This represents every tournament/path pair up to simultaneous relabelling. Both transfer inequalities hold for all 33,868 such representatives (1,1,2,8,64,1024,32768 across n=1,...,7). Explicit RuntimeError guards remain active under Python optimization. Full exact output is in `results.json`.

Finite controls support, but do not replace, the all-order arguments above. The results concern fixed feedback strata only. No unrestricted A000568 asymptotics, growing-k uniformity, global tournament classification, strong-conditioned attachment law, or literature-priority claim is inferred.

## Final editorial delta and binding

The final proof was read in full and its SHA256 independently recomputed as `0ee3eed0b3d2ef55926ad210c79ba67efa1014491e5e0c7841fa555cee87fa62`. PASS is rebound to this final source. Its revised title, explicit unconditional endpoint probability and absolute-error inequality, and direct references to the proved attachment and transfer dependencies preserve all mathematics. The endpoint wording recommendation is resolved. The attachment dependency has the exact reviewed hash listed above; the mass conclusion is a theorem using that proved dependency, not a conjecture conditional on an unproved mean. The reflection/incidence convention remains the specified convention. Publishing the lower-transfer dependency at the linked sibling `session26-blue/RESULT.md` is a layout requirement, not a mathematical hypothesis.

Source-checker replay matches the supplied results exactly. Independent normal and optimized Python executions have byte-identical output.
