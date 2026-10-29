# Audit and scope

Same-agent adversarial and alternate-path audit,2026-10-28. No subagent review
is claimed. Proofs remain in this research packet; no scarce theorem ID was
reserved and no unreviewed stub was promoted into canon.

## Exact reproduction

`check_bridges.py` runs with the Python standard library both normally and
under `python3 -O`; explicit runtime checks are not disabled by optimization.
Both runs give2509 checks and identical JSON. The script rewrites its own
`receipt.json` deterministically.

- Every one of55 geometric point pairs is checked in Z[zeta5]. An independent
  cyclic-autocorrelation check uses divisibility by Phi5, avoiding the first
  routine's polynomial reducer.
- All15 native arcs, all55 Paley arcs and their translation closure are
  explicit. The mirror fixes0,1,7;25 Paley arcs keep and30 reverse orientation.
- Plain permutation breadth-first closure, without a group library, gives660
  elements and55 diagonal point/block actions. Separately, the explicit12-frame
  action is faithful and identifies the generators with projective matrices.
- The12-frame universe is complete for the stated group orbit. Agreement
  scores are15 twice and6 ten times. The1->2 probe has squared length4-phi,
  deliberately failing the unit-contact condition.
- The period7 A29 result does not rely on the exploratory SymPy order. The
  saved12-factor word has cycle type(3)(25)(1); its25th power is(0 20 23).
  Prime degree, generator parity, and the cited Jordan theorem finish the
  argument. No inference to every later Lucas prime is made.
- Threshold controls include330 positive rational cyclic-gap sources,
  sampled homotopy parameters, nonuniform thresholds, equalities and zero
  coordinates. They are controls, not the proof of the all-parameter theorem.
- Word counts and their signed six-period law are checked for3<=n<=17.
  The15-object symmetry hostile uses reflection fixed counts1 versus3.
- Late integration through b07ec51b20 adds exact checks of both cyclotomic
  plane drawings and their golden Gram eigen-equations. A native star Gram
  is constructed and rationally row-reduced for1<=m<=5, supporting the direct
  all-size proof that fixed homotopy type does not bound return rank.

## Independent conceptual check of the gap topology

The explicit homotopies admit a different check by a convex-cover nerve.
For a maximal face F of K let

    U_F={g in Delta:g_i<a_i for every i outside F}.

These relatively open convex sets cover B_K(a). An intersection over faces
F_1,...,F_r is nonempty iff their intersection F is nonempty: when F is empty
the strict inequalities contradict equal total mass; when F contains j,
the vertex e_j works. In the strictly positive gap simplex use a point
sufficiently close to e_j instead. Every nonempty intersection is convex.
The same nerve describes the cover of K by its maximal simplices, whose
intersections are simplices on F. The standard finite good-cover nerve
argument therefore reproduces the homotopy type obtained by the explicit
maps. The constructive proof in GAP_TOPOLOGY is the primary proof and also
provides equivariance without choices of cover refinements.

## Hostiles and information audit

1. **No unmarked LRC substitution.** The all-runner space has cycle-complex
   topology; one specified runner gives a contractible space. Integer
   velocities, collision transport and the marked target have not been encoded
   in the homotopy equivalence. LRC(14) stays open.
2. **No arithmetic-isometry substitution.** The same finite labels can carry
   golden rotation, Paley translations or plane distances, but these operations
   preserve different relations. The geometric phi reads as8, while the
   golden reader uses4. Incoming Session91 proves the corresponding ambient
   golden matrix is not an isometry of any positive-definite flat metric.
3. **No scalar oddness substitution.** Odd period, odd cardinality, trivial
   Royle sign, odd automorphism order and converse fixed-class counts are
   separately typed. Lambda5 itself is sign-even with even automorphism order.
4. **No positive-summand conclusion.** Its15-edge signed copy count vanishes
   in every self-converse eleven-vertex tournament, despite a positive native
   embedding in Paley11. The pairing count is algebraic, not an empirical
   cancellation estimate for a Collatz source row.
5. **No general compact closure.** The exact A29 hostile blocks extrapolation
   of the eleven-point biplane mechanism to the next Lucas prime.
6. **No new packing optimum.** The plane Lucas realization is not the physical
   eleven-square contact graph or the metric torus. Prior packing sources and
   their current validation scope are routed through Session91.

## Inheritance and priority checks

- Read THM-4592 with the later Session74 all-period proof before deriving a
  weaker finite reader statement. The current three-word even exception is
  respected. Incoming Session91's toggle and norm refinements are linked.
- THM-370 supplies the static gap criterion; the new step is a continuous
  threshold-space equivalence and the explicit marked-target comparison.
- Brown2004 already displays the exact pentagonal biplane mirror recovered
  here. PSL(2,11), its order660 and the biplane are credited as classical.
- The generic plane cube realization was found in the stronger inherited
  cubelike theorem and is not claimed as a new result.
- The Monster4/5/6 filtration is inherited with MISTAKE-668's cover/index
  correction. MISTAKE-081's false self-converse halving identity is not used.
- No historical reflection, finite census or RESERVED statement is a proved
  dependency. No demonstrated live mathematical error was found requiring a
  new MISTAKE entry; false candidate extensions were stopped and recorded here.

The figure uses floating coordinates only for display. All geometric truth
claims use the exact checker. Its rendering script needs matplotlib; the
mathematical checker does not.
