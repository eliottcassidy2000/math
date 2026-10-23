# Fano line flips: exact continuation state, fourth-order blindness, and cover memory

Status: PROVED elementary coding/observer consequences; FINITE-EXACT base Hamiltonian identity independently certified two ways. No literature-wide novelty claim. The Hamming-code structure and 171/189 split are inherited from `05-knowledge/results/arithmetic_braids2_20260917_fano_code.md` and audited `collatz_mod6_20260922_paley_fano_octonion_design.md`, not new discoveries here. This note extracts their operation algebra and sensitivity consequences.

## Inheritance and target

Closest mechanism: the repo's Fano incidence code and THM-2176's operation-sensitive continuation quotient. Hostile: equal current Hamiltonian counts need not have equal responses to a specified line reversal. Corrected near miss: a graph quotient can preserve every local move and the desired scalar without lifting closed walks to closed walks. Least-used sidecar: the three-bit syndrome and its sixteen-element deck group.

Source: the 128 cyclic orientations of the seven lines of one fixed labelled Fano plane. Operation: reverse **all three arcs** of one specified line. Target: the Hamiltonian-path count after arbitrary words of these operations. These are not single-arc reversals, and not the fixed-Hamiltonian-path chord complement defining the user's merged metagraph. The construction supplies a different, explicitly defined move graph.

## 1. Coordinate convention and independently checked base identity

Use vertices Z/7 and lines L_i={i,i+1,i+3}, i=0,...,6. The baseline orientation is the Paley tournament: x→y iff y−x is 1,2,4 modulo7. Each line is a directed triangle. A bit m_i reverses that line. Every pair belongs to exactly one line, so this defines a tournament T_m for every m in F2^7. All vertices have score3: each lies on three cyclic lines and receives one outgoing arc from each.

Identify the vertex labels with nonzero F2^3 vectors in the order (1,2,4,3,6,7,5). The corresponding line normals are

(a_0,...,a_6)=(4,1,3,7,6,5,2).

Define s(m)=XOR of the a_i for which m_i=1. The inherited, independently rechecked identity is

H(T_m)=171+18·1_{s(m)=0}.                                      (1)

Exactly16 masks have syndrome0, and112 have nonzero syndrome. The checker verifies (1) first by subset dynamic programming on all128 tournaments, then by a separate expansion over all5040 vertex words. This finite identity is not promoted to arbitrary regular tournaments or arbitrary Fano-like geometries.

## 2. The minimal state depends on the allowed intervention

**Theorem A (labelled continuation).** Two masks have the same H after every word of specified line reversals iff their syndromes agree. The minimal deterministic continuation state therefore has eight classes, although current H has only two values.

Proof. Reversal i sends s to s+a_i, so equal syndromes have equal future values by (1). If s≠t and one is0, the empty word distinguishes them. If both are nonzero, reverse the unique line with normal s: its new syndrome is0, while the other's is t+s≠0. This one-step context distinguishes their counts. Every syndrome occurs, so eight states are necessary and sufficient. QED.

This is an exact context-dependent minimality result. It does not imply that an arbitrary numerical encoding needs eight real coordinates: eight deterministic states require only three bits, and linear realization is a different question.

**Theorem B (averaged continuation).** If every step chooses a line independently and uniformly, the two H-level sets themselves form an exact Markov quotient, with transition matrix

[[0,1],[1/7,6/7]]

in the order (189,171). The probability of reaching a Paley state on the next step from any non-Paley state is1/7. Its first hitting time is geometric with mean7; mean return time from a Paley state is8. This is first return to the *Paley set*, not to the original labelled tournament.

Proof. The seven directions a_i exhaust the nonzero syndromes. From0 all seven moves leave; from s≠0 exactly one reaches0. The strong conditional statement holds after every non-hitting history, giving the geometric law. QED.

Thus scalar compression fails for controlled moves and succeeds for this randomized operation. No invariant is 'complete' without specifying its consumers.

## 3. The complete-graph quotient and its lost holonomy

The line-flip graph is Q7. The syndrome map is a regular graph covering Q7→K8: from any mask the seven neighboring syndromes are each other syndrome exactly once. The fiber is the [7,4,3] Hamming code C=ker(s), of size16, acting by translations.

Every triangle in K8 lifts to an **open** three-edge path in Q7. Its edge differences are a,b,a+b, and its final displacement is the nonzero weight-three codeword supported on those coordinates. The7 possible displacements are the Fano line words; each occurs for8 of K8's56 triangles. They generate all C. Thus the graph is locally faithfully covered while a closed quotient walk can change the labelled orientation inside its fiber. Keeping a deck-coordinate restores the lift.

Proof of the generating assertion without enumeration: the weight-three words express e_a+e_b=e_{a+b} modulo their span. Therefore any vector is equivalent to the coordinate expression of its syndrome in a fixed F2-basis. A zero-syndrome vector lies in that span. The other assertions follow directly from distinct nonzero columns. QED.

This is a concrete graph-cover analogy to affine cochain obstructions in THM-3398. It is not the same carrier: one has finite binary deck translations, the other integer gaps, strict widths, and a CRT lift. What transfers is the necessity of cycle closure, not an LRC conclusion.

## 4. The entire Hamiltonian variation is fourth-order

Character orthogonality applied to (1) gives the exact Walsh expansion

H(T_m)=693/4+(9/4) Σ_{b≠0} (-1)^{Σ_i (b·a_i)m_i}.             (2)

Every nonzero b has b·a_i=1 on exactly four of the seven normals. Hence (2) has a constant term and precisely seven degree-four characters, and no degree1,2,3,5,6,7 terms.

Consequences:

1. Under uniform independent line bits, conditioning on any assignment of at most three specified bits leaves E[H]=693/4. This is blindness of those *conditional averages*, not equality of all pointwise first derivatives.
2. Some four-bit assignment changes that conditional mean. Take the support of one nonzero character. No other distinct weight-four character is contained in it, so the conditional mean becomes693/4 ±9/4.
3. Let x_i be the independent mean of (-1)^{m_i}. The expected Hamiltonian count is the multilinear polynomial

   E_x H=693/4+(9/4) Σ_{b≠0} ∏_{i:b·a_i=1}x_i.

   At x=0 its gradient, Hessian and all third derivatives vanish, but its fourth derivative tensor does not. It is a saddle: on the four-coordinate support of one term, choosing product positive or negative gives opposite deviations arbitrarily close to0. This demonstrates a precise limitation of Hessian-only diagnostics, without identifying this polynomial with the unrelated rank88 Gaussian Hessian.
4. If A is Q7 adjacency, A(H−693/4)=−(H−693/4). Equivalently Σ_i H(T_{m+e_i})=1386−H(T_m). For uniform discrete flips E[H_t|m]=693/4+(−1/7)^t(H(T_m)−693/4). For rate1 flips per coordinate, the continuous-time centered expectation decays as e^(−8t).

These are exact spectral consequences, not fitted asymptotics.

## 5. All-r extension of the mechanism, not of the tournament identity

For r≥2 put n=2^r−1 and use every nonzero vector of F2^r as a column. The syndrome map F2^n→F2^r gives a graph cover Q_n→K_(2^r), with fiber size2^(n−r). Its code indicator has character expansion

1_{s=0}=2^(−r) [1+Σ_{b≠0} χ_{(b·a)_(a≠0)}].

Each nonconstant term has degree2^(r−1). Consequently the indicator is invisible to conditional averages fixing fewer than2^(r−1) individual bits, although r parity measurements recover its full continuation state under specified coordinate flips. Its centered adjacency eigenvalue is−1. These follow by the same proof for every r, not finite extrapolation.

There is **no** asserted tournament construction at r≥4 whose Hamiltonian count equals this code indicator. The generalization concerns the coding/observation mechanism. Values7,8,16,128 arise here from dimensions, kernels and quotient sizes, not a new common numerical invariant connecting all repository topics.

## 6. Reproduction and stopping boundary

Run `python3 check_fano_observer.py > fano_observer.json` and also under `python3 -O`. The standard-library implementation has848 explicit checks, tests all128 orientations, all5040 independent word expansions,379 partial assignments, all56 quotient triangles and all21 pairs of distinct nonzero syndromes. Its r=2..5 instances supplement the all-r proof; they are not its proof. Normal and optimized outputs agree.

Cheapest further target: compute analogous continuation rank for other explicitly defined regular-tournament triangle-reversal carriers. Do not transfer the eight-state classification to arbitrary triangle flips, all regular tournaments, or merged metagraph wiring without proving operation compatibility.
