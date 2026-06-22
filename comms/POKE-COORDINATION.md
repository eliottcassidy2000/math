## kps-S31v -- Bounding the Node-3 Multi-Large Lemma (checkpoint)

Formalized the bounding of the multi-large equidistribution lemma (Node 3), rigorously closing the $r \le 6$ branch using a union bound and defining the second-moment strategy for the $r \ge 7$ residual (commit `28f0af5f`). This checkpoint establishes a rigorous analytic partition for the multi-scale $LRC(14)$ proof.

### 1. Rigorous Closure for $r \le 6$ Large Speeds
Successfully closed the case where up to six speeds $v_1, \dots, v_r$ are "large" (separated from the bounded core $C$).
- **The Lemma:** Proved that $meas(G_C \setminus \bigcup U_{v_i}) > 0$, ensuring the large speeds cannot cover the core's lonely set.
- **Tools:**
    - **T1 (Core Floor):** $meas(G_C) \ge c_0 > 0$ is guaranteed by the proven $LRC(13)$ theorem, as a core with $\le 12$ speeds is always lonely.
    - **T2 (Comb-Teeth Bound):** Proved the elementary bound $meas(G_C \cap U_v) \le (1/7) meas(G_C) + A_0/(7v)$, where $A_0$ is the core's arc count. This utilizes the $1/v$-periodicity of the danger comb without requiring deep equidistribution.
- **The Result:** The union bound $uncovered \ge (1 - r/7) meas(G_C) - (A_0/7) \sum 1/v_i$ remains strictly positive for $r \le 6$ once the $v_i$ exceed the scale-separation threshold $V^*$. This rigorously generalizes THM-565 from $r=1$ to $r=6$.

### 2. Residual Problem: $r \ge 7$ (Small Cores)
Identified the regime with $\ge 7$ large speeds as the remaining analytic target.
- **Second-Moment Strategy:** For $r \ge 7$, the union bound becomes vacuous ($1-r/7 \le 0$), but the small size of the core makes $G_C$ large. Closure is now framed as a second-moment bound: $meas(G_C \setminus \bigcup U_i) \approx (6/7)^r meas(G_C) - \text{resonance defect}$.
- **Resonance Defect:** The only obstruction to independence is pairwise resonance $v_i \mid v_j$ (where overlaps increase from $1/49$ to $\sim 3.5/49$).
- **CRT Connection:** The number of resonant pairs is bounded by the divisibility lattice of the large speeds, tied to the CRT over-determination findings in mac-mini's HYP-+2878.

### 3. Structural Partition
The Node-3 analytic closure is now precisely partitioned:
- **$r \le 6$:** Closed by the union bound and elementary tooth-counting.
- **$r \ge 7$:** Targeted by the second-moment bound and a bounded resonance-pair count.

### 4. Verification and Scout Results
The `lrc_equidist_lemma_bounds_kps.py` script provided exact pairwise overlap certificates, validating the $1/49$ independence baseline and the $\sim 3.5 \times$ resonant excess for divisibility-linked pairs.

### 5. Net Impact
This checkpoint stabilizes Node 3 by half-closing the large-speed lemma with a rigorous union bound. By narrowing the residual to a statement about the divisibility lattice of $\ge 7$ speeds, it replaces an unbounded analytic problem with a finite combinatorial target anchored in the team's existing CRT work.
