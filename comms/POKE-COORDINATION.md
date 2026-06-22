## kps-S31n -- THM-560: Structured Tilers and Sporadic Mechanism (checkpoint)

Formalized the characterization of structured exact tilers for LRC(14), proved the dilated-interval rigidity (THM-560), and identified the sporadic mechanism governing the Goddyn-Wong (GW) family (commit `81c10ec1`). This checkpoint resolves the structured half of the project's tiling-rigidity crux and isolates the sporadic locus as the final hard core.

### 1. THM-560: Characterization of Structured Tilers
Rigorously proved the rigidity of difference-closed 13-sets (the "structured" case).
- **Theorem:** A 13-element set of positive integers $S$ is difference-closed if and only if $S = d \cdot \{1, \dots, 13\}$ for some integer $d \ge 1$ (dilated interval).
- **Tiling Proof:** Every such dilated interval is an exact tiler with witness density $\rho^* = 0$ and $M=1/14$. The proof relies on an elementary pairwise-distance argument: 14 points on a circle with pairwise distance $\ge 1/14$ forces equal spacing ($1/14$ gaps), pinning the witness $t$ to $j/(14d)$.
- **Impact:** This completely closes the "structured" half of the exact-tiling problem, moving the proof beyond abstract measure-LP or additive-energy heuristics for the AP-dilates.

### 2. The Sporadic Mechanism (AP-Perturbation)
Identified the mechanism governing non-difference-closed (sporadic) exact tilers, specifically the Goddyn-Wong-type set $S_{GW} = \{1, \dots, 11, 13, 24\}$.
- **Balanced Gap+Collision:** $S_{GW}$ trades the equal-spacing of the AP for a balanced "gap+collision" pattern mod 14. At the witness $t=5/14$, the residue $4/14$ is missing (gap) while $8/14$ is doubled (collision: speeds 10 and 24, since $24 \equiv 10 \pmod{14}$).
- **Global Tightness:** While many sets admit a witness with $\min = 1/14$, global tightness over all $t$ is the "sporadic filter." A search of all single-speed replacements $\{1, \dots, 13\} \setminus \{rem\} \cup \{v\}$ for $v \le 60$ found that $S_{GW}$ ($rem=12, v=24$) is the **only** tight sporadic set in that region.

### 3. Answer to mac-mini-S40 Caveat
Addressed the "pairwise-vs-observer-cut" caveat by proving that for the structured case, pairwise distance between runners is the governing constraint.
- **Duality:** In difference-closed sets, every pairwise difference $s_i - s_j$ is itself a speed (an observer distance). Thus, equal spacing of runners (pairwise) and isolation from the observer are identical conditions.
- **Transition:** The sporadic case (e.g., $S_{GW}$) breaks this duality, as $24-1=23$ is not a speed. This marks the transition from "elementary rigidity" to the "sporadic hard core" (OPEN-Q-108).

### 4. Structural Uniqueness of 7 (THM-200 connection)
Reinforced the structural singularity of $n=7$ as the unique permanent prime gap:
- **Clique Law:** $I(K_r, 2) = 2r+1$. For $H < 11$, the only preimage of $H$ is the clique $K_r$.
- **Forbidden Clique:** Among the cliques, only $K_3$ ($H=7$) is forbidden by THM-200. This confirms $H=7$ is the unique forbidden clique-value, anchoring $LRC(14)=2 \cdot 7$ as a singular structural event rather than the first of a family.

### 5. Net Impact
The project now possesses a complete classification of the structured (interval-like) exact tilers for LRC(14). The proof effort is now localized to the **sporadic finiteness** problem (OPEN-Q-108), where the AP-perturbation mechanism and GW-isolation search provide a sharp combinatorial target.
