> **HISTORICAL SYNTHESIS / GLOBAL GIRTH CLAIM SUPERSEDED (2026-08-25).**
> THM-264 is finite-exact only through order six.  Order seven refutes
> `alpha_1>=3 => girth(Omega)=3`; global no-intermediate girth and global
> complement acyclicity is refuted at order nine.  See MISTAKE-507.  Later claims
> in this file that explain all tournaments through that dichotomy are
> speculative, not current canon.

# The Binary Skeleton

**Session:** kind-pasteur-2026-03-21-S18b
**Arising from:** the finite THM-264 census, the six-patterns reflection, and a careful reading of prior work

---

## The Discovery

The original session observed `girth(Omega(T)) in {3,infinity}` through the
exhaustive range `n <= 6` and conjectured it globally.  That extrapolation is
**OPEN**, not a theorem; its proposed `alpha_1` threshold already fails at
order seven (THM-343).

This is not an isolated fact. It is one bone in a skeleton that runs through the entire project.

---

## The Skeleton: 26 Binary Phenomena

A careful scan of the repository reveals that binary dichotomies — quantities that take exactly two values, or properties that hold everywhere or nowhere, with no intermediate state — appear at every level of the theory. Here are the major ones, organized by what kind of binary they are.

### Type A: Always/Never (Universal Binary)

These are properties that hold for ALL tournaments or NONE:

1. **H(T) is always odd** (Redei) — H mod 2 = 1, always. The deepest binary.
2. **beta_2(T) = 0 always** — Path homology in degree 2 vanishes universally. Not "usually 0" — ALWAYS 0, for every tournament ever tested (~47,000).
3. **Every tournament is weakly connected** — Consequence of Redei. No disconnected tournaments exist.
4. **tr(M) = H at odd n, tr(M) = 0 at even n** — The transfer matrix trace is either the full path count or exactly zero, determined by a single bit (parity of n).

### Type B: Mutual Exclusion (Seesaw Binary)

These are pairs of properties where exactly one holds:

5. **beta_1 * beta_3 = 0** — 1-holes and 3-holes never coexist. The im(d_2) mediator allocates dimensions exclusively to one or the other.
6. **Adjacent-odd seesaw: beta_{2k-1} * beta_{2k+1} = 0** — Generalizes the beta_1/beta_3 exclusion to all adjacent odd Betti numbers. Verified through n=8.
7. **rank(d_2) in {n-1, n} only** — Two values, nothing between. The value n-1 occurs iff beta_1 = 1.

### Type C: Sharp Threshold (Phase Transition Binary)

These are properties that switch from universal to possible at a precise n:

8. **beta_3 onset at n=6** — Zero for all tournaments at n <= 5. First nonzero at n=6. Requires 2 disjoint 3-cycles = 6 vertices.
9. **beta_4 onset at n=7** — Zero below. First at n=7 (Paley T_7: beta_4 = 6).
10. **beta_5 onset at n=8** — Zero below. First at n=8 (~0.2% prevalence).
11. **alpha_1 = 3 impossible at n <= 6, achievable at n >= 7** — The common-vertex constraint breaks at n=7.
12. **Per-path identity: works at n <= 5, fails at n >= 6** — The mu-triviality boundary.
13. **I(Omega, x) real roots: always at n <= 8, can fail at n = 9** — The claw-free boundary.
14. **Omega is interval graph: always at n <= 5, fails at n = 6** — Graph class boundary.
15. **Kneser girth transition: K(5,2) has girth 5, K(6,2) has girth 3** — The Petersen breaks at n=6.

### Type D: Forbidden Values (Gap Binary)

These are specific values that can NEVER be achieved:

16. **H = 7 is impossible for all n** — The atomic forbidden value. Proved via cycle-forcing.
17. **H = 21 is impossible for all n** — The compound forbidden value (3 x 7). Proved via poisoning graph.
18. **H = 7 and H = 21 are the ONLY permanent gaps** — All other odd values are eventually achieved. Binary: permanent vs temporary gap.

### Type E: Structural Dichotomy (Classification Binary)

These are properties that split tournaments into exactly two classes:

19. **Finite Omega girth census** — THM-264 proves values `{3,infinity}` only for `3 <= n <= 6`; the global statement is **OPEN**.
20. **Self-complementary vs non-SC** — SC tournaments have involution anti-automorphisms; non-SC don't. beta_1 > 0 requires SC at n=6.
21. **Hereditary maximizers: odd n only** — At odd n, regular maximizers are hereditary. At even n, NO maximizers are hereditary.
22. **Pfaffian classification at n=6** — beta_1 > 0 implies |Pf| in {1,3}; beta_3 > 0 implies |Pf| in {7,9}. Perfect separation.
23. **Blueself/blackself tiling: no mixing** — No isomorphism class contains both types.

### Type F: Modular Binary

24. **H mod 4 = 1 + 2*alpha_1 mod 4** — Binary switch: alpha_1 even -> H = 1 mod 4; alpha_1 odd -> H = 3 mod 4.
25. **det(B) = Pf^2 at even n, det(B) = 0 at odd n** — Even/odd parity creates the fundamental spectral dichotomy.
26. **Finite complement census** — complement acyclicity holds through `n=6` but is **REFUTED** globally by THM-264's order-nine `K_3` hostile.

---

## What the Skeleton Means

These 26 binary phenomena are not 26 unrelated facts. They are manifestations of a single underlying principle.

### The Principle: Tournament Structure Is Maximally Decisive

A tournament is a complete oriented graph. Every pair of vertices has a definite winner. There are no ties, no abstentions, no partial outcomes. This maximal decisiveness at the pairwise level propagates upward through every algebraic structure built on tournaments:

- **Homology is maximally decisive:** beta_2 = 0 always (the second boundary map is always exact). The seesaw beta_1 * beta_3 = 0 means: the homological "energy" is allocated to one degree or the other, never shared. This is decisiveness in homological dimension.

- **Girth supplied a finite binary-looking signal:** only girth three or
  infinity occurs through order six.  Whether intermediate conflict girth
  occurs later is open, so completeness alone does not yet explain it.

- **H is maximally decisive modulo 2:** Always odd. The parity is locked. This is Redei's theorem — the most fundamental binary property of tournaments.

- **Score determination is maximally decisive:** The root cycle profile determines H exactly at n=5, then fails by exactly 4 (one quantum) at n=6. The failure is not gradual — it's a quantum jump at a sharp threshold. Decisiveness to the last quantum.

### Why No Intermediate States?

The original text attributed every missing intermediate state to the slogan
**tournament completeness forces saturation**.  That is not a proof.  Three
cycles need not form a conflict triangle: THM-343's order-seven examples have
three canonical odd cycles and `Omega=K_1 disjoint_union K_2`.  Any global
girth argument must control the full cycle-intersection hypergraph, including
distinct directed cycles with the same vertex support.

This is exactly like percolation on the complete graph. In the Erdos-Renyi model G(n,p), there is a sharp threshold for the appearance of a giant component. Below threshold: many small components. Above threshold: one giant component. The transition is sharp — no "medium-sized components" persist.

Tournament conflict is percolation on a complete substrate. The substrate's completeness forces the percolation to be maximally sharp: the "giant component" (triangle in Omega) appears the instant there are enough vertices.

### The Zero-One Law for Tournament Properties

The Erdos-Renyi zero-one law states: for any first-order property of random graphs, the probability that G(n, 1/2) satisfies the property tends to 0 or 1 as n -> infinity. There is no property that holds with limiting probability 1/2 or any other intermediate value.

Tournament theory has an analogous zero-one structure, but it's deterministic rather than probabilistic:

- beta_2 = 0 with probability 1 (universal, deterministic)
- beta_1 * beta_3 = 0 with probability 1 (universal, verified)
- Omega girth in `{3,infinity}` throughout the exact `n<=6` census; the
  all-order version remains open

These are not statements about random tournaments. They are statements about ALL tournaments. The binary skeleton is not a statistical phenomenon — it is algebraic.

---

## Where the finite `{3,infinity}` signal suggested analogies

The following are historical prompts.  They do not promote the finite Omega
census to an all-order dichotomy.

### 1. The Kneser Girth Transition (Found by opus-S119)

K(n,2) — the Kneser graph on 2-subsets of [n] — has girth:
- K(3,2): infinity (no edges)
- K(4,2): infinity (matching, no cycles)
- **K(5,2): 5** (Petersen)
- K(6,2): 3 (first triangle from 3 disjoint pairs)
- K(n,2) for n >= 6: 3

Thus `K(n,2)` has girth infinity for `n=3,4`, five at `n=5`, and
three for `n>=6`.  This makes the Petersen graph a useful hostile for any
putative conflict-graph representation, but THM-264 does **not** exclude the
Petersen graph from occurring as some `Omega(T)` at larger order.

### 2. The Transfer Matrix Trace (tr(M) in {H, 0})

The transfer matrix M satisfies tr(M) = H at odd n, tr(M) = 0 at even n. This is another {value, zero} dichotomy. The parity of n acts as a binary switch that either enables or completely annihilates the trace. No partial trace exists.

The trace identity is exact where separately proved; its resemblance to the
finite girth census is only an analogy.

### 3. The Seesaw (beta_1 * beta_3 = 0)

This is a {left, right} dichotomy rather than {finite, infinite}, but the structure is the same: two mutually exclusive states with no interpolation. A tournament's homological "energy" is entirely in degree 1 or entirely in degree 3, never split.

The proposed cycle-level dichotomy remains open beyond the finite census.

### 4. The H = 7 Impossibility

H = 7 requires (alpha_1 = 3, alpha_2 = 0). But alpha_1 = 3 with alpha_2 = 0 means 3 pairwise-conflicting cycles — i.e., a triangle in Omega, i.e., girth(Omega) = 3. AND it requires no other cycles. But girth = 3 with alpha_1 = 3 forces alpha_1 >= 4 (the triangle implies a shared vertex which generates a 4th cycle). So H = 7 is impossible BECAUSE the girth-3 phase is too entangled: once you enter it, you can't stay at alpha_1 = 3.

**H = 7 is impossible because the binary phase transition overshoots.** You want exactly 3 cycles in the girth-3 phase, but the girth-3 phase FORCES >= 4 cycles. The threshold is too eager.

### 5. The Partition Function Bridge (opus-S127)

Opus discovered that tournament H and stabilizer code distance d are both evaluations of the same partition function type at different points:
- H = I(Omega, 2) — evaluation at x = 2 (counting regime)
- d = lim_{x->0} of a related function (minimum weight regime)

The partition-function/code-distance paragraph is conjectural.  Girth alone
does not monotonically determine an independence polynomial, and no typed map
from these tournament conflict graphs to the cited stabilizer constructions
was supplied.  Retain it only as a prompt for a controlled pair of graph
families with all side data fixed.

### 6. Ecological Intransitivity (External Analogy)

Kerr et al. (Nature 2002) found that bacterial competition (E. coli with colicin-producing, colicin-sensitive, and colicin-resistant strains) exhibits a sharp phase transition:
- **Local dispersal:** Rock-paper-scissors cycles persist indefinitely. Biodiversity maintained. (Analogous to girth = 3: cycles everywhere.)
- **Global dispersal:** One strain dominates. Hierarchy established. (Analogous to girth = infinity: no competition cycles.)
- **No intermediate regime:** There is no stable state with "some cycles but not everywhere." The spatial structure either sustains cycles or doesn't.

This is the {3, infinity} pattern in biology. The tournament conflict graph of the ecosystem either has triangles (local rock-paper-scissors) or is acyclic (global hierarchy). The spatial range of interactions acts as the parameter that switches between phases, just as alpha_1 acts as the switch for Omega.

### 7. The Ramsey Connection

Ramsey theory's most basic result: R(3,3) = 6. In any 2-coloring of the edges of K_6, there is a monochromatic triangle. Equivalently: among 6 people, either 3 are mutual friends or 3 are mutual strangers. No intermediate state.

The former Ramsey paraphrase was false: an `Omega` with at least three
vertices need not contain a triangle, as THM-343's `K_1 disjoint_union K_2`
hostile shows.  THM-264 currently supplies only the finite `n<=6` census.

The Ramsey number R(3,3) = 6 appears in our setting as well: n = 6 is where the Kneser girth transitions from 5 to 3, where beta_3 first appears, where the per-path identity first fails, and where alpha_2 first becomes possible. Six is the tournament Ramsey threshold.

---

## The Deeper Question

Why is tournament structure so relentlessly binary?

I think the answer is that tournaments are **complete.** Every pair has a definite orientation. There are no "missing" edges, no ambiguity, no partial information. This completeness is the source of all the binary phenomena:

- Completeness constrains cycle intersections, but it does not by itself make
  the conflict graph dense; the exact threshold claim is refuted.
- Completeness makes every path extendable (Redei: H >= 1 always), locking the parity.
- Completeness makes the path homology exact in degree 2 (beta_2 = 0: every apparent 2-hole is actually a boundary).
- Completeness makes the score sequence highly constraining (the root cycle profile almost determines H).

In a general digraph, all of these binary properties fail. Beta_2 can be nonzero. Girth can be 4, 5, 6. H can be even (or zero). The binary skeleton is **specific to tournaments.** It is the algebraic signature of complete pairwise comparison.

This suggests a meta-theorem: **every "soft" structural property of general digraphs becomes "hard" (binary) when restricted to tournaments.** The completeness constraint rigidifies the entire structure, collapsing continuous spectra into discrete alternatives.

If a repaired meta-theorem exists, it must state its carrier and quantifiers
separately for each invariant.  The finite girth signal cannot currently be
grouped with the proved parity statements as one theorem.

---

## Concrete Predictions

If the "completeness implies binary" principle is correct, it predicts:

1. **The chromatic number of Omega(T) should have a narrow range** — probably {1, 2, 3, ...} but with a sharp threshold separating low from high, analogous to girth.

2. **The clique number of Omega(T) should equal the independence number minus some universal constant** — because completeness constrains the two to be closely related.

3. **Any graph invariant of Omega(T) that takes a continuous range for general graphs should take a discrete (possibly binary) range for tournament conflict graphs.**

4. **Removing completeness (studying general digraphs) should break ALL the binary properties simultaneously.** HYP-327 already confirms this: removing even ONE edge from a tournament can make beta_2 > 0.

5. **Weighted tournaments (where arc "strength" varies continuously) should exhibit a phase transition from binary to continuous as the weight variance increases.** At zero variance (standard tournaments), all properties are binary. As weights become non-uniform, the binary structure should soften.

These predictions are testable. The fourth is already partially verified (HYP-327). The others are new.

---

*The durable lesson is methodological: a striking finite binary pattern is a
reason to search for the missing coordinate and a hostile boundary case, not
to infer an all-order law.  Here the missing coordinates were canonical
directed-cycle multiplicity and the full cycle-intersection hypergraph; the
Petersen-as-conflict question and intermediate-girth question remain open.*
