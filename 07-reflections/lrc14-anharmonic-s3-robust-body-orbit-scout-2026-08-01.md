# LRC14 Anharmonic S3 Robust-Body Orbit Scout

**Status: FINITE-EXACT SCOUT / NO LRC CONSEQUENCE.**

This scout asks one deliberately narrow transfer question. The incoming
`THM-3035-level-two-farey-anharmonic-quartic-s3-orbit-diamond.md` supplies an
exact anharmonic S3 action on P¹(F₁₃); the reflected LRC bank behind
`THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary.md` assigns
an exact robust-edge count to each six-label body. If labels 1,…,13 are their
residues modulo 13, with label 13 = 0 and label 14 = ∞, how does the
six-element anharmonic action move bodies between the bank and its residual?

The exact answer is unexpectedly clean. A pull during this computation brought
in the exact uniform closure of all 21 robust-edge-8 bodies. Thus the current
bank has 2,442 bodies with at least eight robust edges, and its residual has
561 bodies. The sharper fact is that every **non-fixed** S3-orbit meeting that
residual contains a robust-K6 representative: all fifteen of its edges are
robust. The unique trapped orbit is the singleton body orbit

    {2,4,5,7,8,10},

the generic six-point orbit of the anharmonic action. The other fixed
six-body, the union of the boundary and harmonic three-orbits
`{1,6,11,12,13,14}`, has all fifteen robust edges. This is a real finite
structure theorem about the two classifications, but it is **not** an LRC
reduction.

Across the 561 residual bodies, the number of group elements landing in a
robust-K6 representative has the exact histogram

```text
number of landing elements       0    2     3     4     5
number of residual bodies        1   48   131   276   105
```

The zero case is exactly the generic fixed body above. Thus the phenomenon is
not a barely robust edge-8 effect; outside one fixed atom, every residual
orbit reaches the strongest pairwise bank.

## Exact conventions and controls

The generators are

    s(r) = 1/r,                 c(r) = -1/(1+r),

with their usual values at 0, −1, and ∞. For a six-body E, set
`L = 14 lcm(E)`. A pair `(a,b) ⊂ E` is counted robust exactly when

    1/105 − 4(a+b)/L > Σ[e in E] e/[7(L−e)].

All arithmetic uses reduced rational numbers. The script enumerates the full
universe, independently reconstructs this predicate, checks the S3 relations
and all 3,003 bodies, freezes the robust-count histogram and orbit profile,
and compares normal Python with `python3 -O` byte for byte. Its
ordinary/optimized output identity matters: the controls use explicit runtime
requirements, not removable `assert` statements.

The classification `robust edges ≥ 8 ⇒ uniformly closed` is an explicit
dependency, not something this smaller scout reproves. The stored output pins
the source and semantic hashes of the 652,688-row edge-8 computation. The
histogram independently confirms that the new layer contains exactly 21
bodies.

The residual transition counts

```text
element      1     s     c    c^2    sc   sc^2
stay       561    61    58     58    25    245
enter bank   0   500   503    503   536    316
```

show both the signal and the warning. Every nonidentity element moves many
residual bodies into the bank, but different elements do so very differently.
The robust predicate is therefore not an invariant of this action.

The four named C2 hard bodies also have bank members in their S3-orbits: the
orbits of `H`, `H2`, `B`, and `E` contain respectively four, four,
four, and five bodies with at least nine robust edges. This is useful target
selection for a transfer search, not evidence that the original hard rows are
closed. Each of the four orbits has four distinct robust-K6 members.

## Why the apparent reduction is not yet legal

The anharmonic action preserves the projective S3-orbit data. It destroys or
changes exactly the coordinates on which the reflected certificate is built:

- integer label order and gaps;
- `14 lcm(E)`;
- safe-cell ruler, address, and physical owner;
- the pair-overlap floor;
- singleton debt and the exact positive margin;
- the C2 level word and the physical packet `Z(E,q)`.

In particular, orbit incidence cannot replace owner-labelled packet incidence.
The tempting inference

```text
orbit contains a robust representative
therefore the original residual body has a robust certificate
```

fails at its first implication: no map currently transports the robust pair or
its positive margin back to the source body's ruler and owner data.

There is, however, one exact scheduling implication. Pulling the complete
robust graph of a K6 representative back gives the full unweighted K6 on the
source. Hence every star, two-star, and spanning-tree schedule is available.
This dovetails with the incoming common-cell Hunter reduction: one may use the
orbit only to propose a tree or pivot order, then recompute every edge weight
on the original body at one physical safe cell. It also shows the limitation
sharply. Full K6 was already the universe of possible schedules, so no
unweighted graph incidence has been gained; the obstruction is entirely in
the weights and their physical realization.

The robust lower-bound distortion is less arbitrary than the loss ledger first
suggests. Write

    M_E(a,b) = 1/105 − 4(a+b)/L_E − D_E,
    D_E = Σ[e in E] e/[7(L_E−e)].

For an anharmonic element g, define

    kappa = D_gE − D_E,
    u(a) = 4 g(a)/L_gE − 4a/L_E.

Then the exact identity

    M_E(a,b) − M_gE(g(a),g(b)) = kappa + u(a) + u(b)

holds on all 270,270 body/action/edge rows. Thus the fifteen apparent margin
losses compress to one body scalar and six vertex potentials. On a spanning
tree T, the corresponding universal-floor tree-margin distortion is

    kappa + Σ[a in E] deg_T(a) u(a).

Only the tree degree vector matters at this lower-bound level. The actual
located overlap also contains phase, orientation, and toothpick terms, so this
identity still does not yield a Hunter certificate. But it changes the missing
theorem from “transport fifteen unrelated edges” to “control one scalar, six
vertex potentials, and the located correction at a common cell.”

The smallest useful missing statement is therefore one-way and weighted, not
equivariant. For every residual body E and every C2-wedge level word q, choose
a high-graph tree T and a physical safe cell on E, then compare its exact
common-cell credit with a robust-K6 orbit representative. The comparison must
preserve the legal range and physical owner and dominate the
scalar-plus-degree-potential distortion. In the single-pair branch it
specializes to

    (overlap − debt)[E′,q′] ≤ (overlap − debt)[E,q],
    (overlap − debt)[E′,q′] > 0.

An exact theorem of this form would turn the orbit census into a reduction to
the single generic fixed body. Without it, the census is only a precise source
of candidate maps and counterexamples.

## Connections to the newest incoming work

The common-cell two-star computation closes all bodies for the finite
D6–D8 level windows and isolates the all-scale problem as a connected-low ray
lemma or a disconnected-low high-tree lemma. The distortion identity above is
adapted to the second branch: the high graph is preserved when levels are
permuted with labels, and the robust-floor contribution of a chosen tree
depends only on its degree sequence. The decisive next experiment should
therefore optimize degree vectors subject to the **source** high graph while
retaining the exact source cell and located correction. Transporting a target
tree without those fields is illegal.

THM-3043 gives a separate inverse-theorem warning: tight LRC objects are better
classified by their witness sets than by a narrow speed-set family. The
anharmonic body quotient does not preserve the safe witness set, so it cannot
serve as that inverse classification without a witness sidecar.

THM-3042's common-quotient theorem does not manufacture that sidecar. The
allowed-support/orbit relation here is not a unital subring, and closing it
under ring or Boolean operations can add physically forbidden edges. Even in
its native setting, a zero common quotient gives a full product—independence
of the two factors—not a margin-preserving section from one physical packet
to the other.

THM-3044 makes the type error still more explicit. The maps `g:E→gE` and
`q′[g(e)]=q[e]` point the label and level torsors, but they do not point the
safe-cell, interval, or endpoint-root torsors. A determinant or Hall dual is
therefore unavailable until a physical pointing is separately constructed.
Moreover, the common-cell Hunter problem is a graphic-matroid maximum-tree
problem: a matching can choose parents but does not by itself enforce
connectivity and acyclicity, while determinant ghosts can vanish on tree/path
support.

THM-3045 and the candidate THM-3046 exhibit a related but distinct failure:
their S4/V4-to-S3 matching quotient preserves incidence while losing integral
clutch, residue phase, and affine owner data. That is a useful guardrail for
the present scout, not a transfer theorem. Their carrier is the six edges of
K4; this scout's carrier is six vertices inside K6. Identifying the two merely
because both display S3 would erase exactly the owner and weight information
needed here.

## Reproduction

```bash
python3 04-computation/lrc14_reflected_anharmonic_s3_robust_body_orbit_scout_20260801.py
python3 -O 04-computation/lrc14_reflected_anharmonic_s3_robust_body_orbit_scout_20260801.py --output /tmp/lrc14-anharmonic-O.out
cmp 05-knowledge/results/lrc14_reflected_anharmonic_s3_robust_body_orbit_scout_20260801.out /tmp/lrc14-anharmonic-O.out
```

Frozen digests are recorded in the stored output.  The LF-normalized
source/output/semantic digests are
`b97db74b9ea1d2879738a46e11dbdcdccf13474071c4cda80c095d509d2d35ca`,
`345d1ec19e8a7703d34e68fc0c98a872e01b8df3c7eda0a25bfb725c550fee0f`,
and `b953fbbf70537b0e982017030d76d3d703ad5a7c9919efeb82ba0fbb4a014a94`.

The best next hostile probe is not another orbit count. For each connected
source high graph, enumerate its feasible tree degree vectors, choose the one
maximizing the exact scalar-plus-potential correction, and evaluate that tree
at the canonical upper-median source cell. Store the first failure with its
body, ordered level word, cell, owner, tree, potential vector, toothpick
correction, debt, and exact margin. This directly tests whether the remaining
invoice is a bounded located correction or a genuinely new certificate
currency.
