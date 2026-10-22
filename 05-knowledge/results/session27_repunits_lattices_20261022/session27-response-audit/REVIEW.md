# Independent Session 27 response audit

Verdict: PASS for the stated identities, finite equivalence theorem, finite census through six, and score-root hostile. No universal reconstruction conclusion is supported or claimed. No repository or Git writes were made.

## Proof review

1. The one-port polynomial correctly partitions paths into new-vertex-first, last, and internal cases. Internal deletion creates ordered spanning two-path covers, including absent-arc bypasses. The k=2, mask (0,1) hostile refutes insertion-only counting. The Boolean coefficient representation is unique.
2. Converse reversal transposes C and exchanges S/E, yielding exactly g_u=2(E_u−S_u)+sum_v(C_uv−C_vu). Its sum is zero. Vertex transitivity makes S/E constant and makes the row/column sums of C equal. This proves converse blindness without arc transitivity. Hamming-weight-only dependence requires the stronger arc-transitive hypothesis used in the proof. The formal divergence-kernel interpretation is correctly distinguished from realizable tournament data.
3. The two-port block cancellation is valid: after relabelling the added vertices, only their mutual arc changes. Counting adjacent blocks gives D_uv=P_u−P_v−J_uv. Zero-mask calibration gives P_u=F(e_u,0)−F(e_u), since the initial zero port is a universal source. Then B=C+C^T, J=P_u−P_v−D, C=(B+J)/2, S=(P−L)/2, E=L+S−rowsum(C). The battery count is exactly 1+2k+3 binom(k,2); diagonal C is zero and needs no measurement. This recovers a sidecar, not the tournament or the complete two-port response.
4. Differentiating a Hamiltonian-path monomial at x_uv cuts it into the required two paths. This works even on absent arcs. No monomial uses both directions of an edge, so reversing an actual u→v changes H by C_vu−C_uv. The transitive-three hostile correctly shows that J is not an orientation oracle.
5. The inherited slot construction has common nonzero initial row and common transition matrices; only the final vector uses core adjacency. Thus the difference has dimension D=2k!3^k, not 2D. Reachable spaces stabilize by length D−1, including the case of full dimension first reached at D−1; invariant completion then proves equality for all words. Independent-row expansion and rational membership elimination in equivalence.py are sound, and actual row/word pairs are retained separately from elimination vectors. Repeated transition destinations retain multiplicity. Empty core is separately handled.
6. Census groups by the actual one-port polynomial before refining with exact two-port values. Discarding singleton groups is valid because they are already distinguished. Recompiled census through six reproduces census.out exactly. The k=7 exploratory output is explicitly incomplete and is not evidence of completed separation. No general injectivity, minimal realization, or novel general automaton algorithm is claimed.

## Independent tests and reruns

- Source SHA256SUMS: all 14 entries verified before and after checks.
- independent.py uses full permutation enumeration for extended Hamiltonian counts, and an independently organized subset/left-path/right-path construction for ordered covers. It does not call the supplied Hamiltonian DP or supplied cover counter.
- All 75 labelled cores of orders 1–4: 1,098 one-port evaluations, complete two-port tomography batteries, and 410 single-edge-flip checks passed.
- 100 deterministic controls compare the supplied common realization directly to independent permutation enumeration (k=0–4, attachment lengths bounded respectively by 7,5,4,3,3).
- Normal and optimized independent outputs agree; explicit exceptions keep checks active under optimization.
- Supplied audit.py rerun: 33,866 polynomial-vs-DP checks and 4,869 automaton-vs-DP checks; positive/negative equivalence controls pass.
- Supplied two_port.py rerun: all 97 audit cores plus cyclic controls of orders 3,5,7 pass. The order-seven cyclic core has H=175 and is not the Paley-seven core.
- Supplied structure.py rerun: 8 cyclic-module-flip pairs at order four; 80 such pairs plus 12 regular-converse pairs at order five; 1,099 converse-criterion checks.
- Root two_port_current.py rerun: exactly 10,650 unordered-pair checks, no failures. This count is 2*1 + 8*3 + 64*6 + 1024*10. These are pair checks, not 10,650 distinct cores.
- The stronger deleted-pair hypothesis remains refuted already at k=2: D_01=0 while ±2H(empty)=±2. Do not promote it.

## New exact-score/E8-root hostile

The construction itself was independently verified, not only its JSON values. Start with the tournament on Z/7 having residue differences {1,2,4}; reverse the directed triangle 0→1→3→0. Each core remains regular with outdegree 3. Append vertex 7 with u→7 exactly for u in {0,1,2}.

Full permutation enumeration gives:

- Original core: H=189; appended tournament: H=657.
- Triangle-flipped core: H=171; appended tournament: H=643.
- Both labelled score vectors are exactly (4,4,4,3,3,3,3,4).
- Both centered vectors are exactly (1/2,1/2,1/2,−1/2,−1/2,−1/2,−1/2,1/2).

Each eight-vertex count enumerates all 40,320 permutations, independently of root's subset DP. The centered vector has integer-even coordinate sum zero, all coordinates half-integral, and squared norm 2. It therefore lies in the standard E8 lattice and is a root. This is the same labelled lattice point with two different Hamiltonian counts, stronger than equal norm or a theta-shell coincidence. Consequently H does not factor through this labelled score/deviation map. It says nothing against two-port reconstruction.

## Provenance

SOURCE_SHA256SUMS records the exact reviewed inputs, including proof.md, its supplied manifest, all supplied Python/C++ sources, root hostile/current sources and JSON, the root dependency probe_responses.py, and inherited slot proof/counter. SHA256SUMS records this audit's deliverables. Reproduction outputs and independent.json are in this directory. No external novelty search was undertaken; the audit verifies mathematics and the supplied computational scope.
