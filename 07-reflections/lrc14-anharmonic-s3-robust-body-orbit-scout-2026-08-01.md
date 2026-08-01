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
561 bodies. Every **non-fixed** S3-orbit that meets that residual also meets
the current bank. The unique trapped orbit is the singleton body orbit

    {2,4,5,7,8,10},

the generic six-point orbit of the anharmonic action. The other fixed
six-body, the union of the boundary and harmonic three-orbits
`{1,6,11,12,13,14}`, has all fifteen robust edges. This is a real finite
structure theorem about the two classifications, but it is **not** an LRC
reduction.

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
closed.

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

The smallest useful missing statement is one-way, not equivariant. For every
residual body E and every C2-wedge level word q, find a closed orbit
representative E′ and pull its exact certificate back into the physical packet
of (E,q). The pullback must preserve the legal range and physical owner, map
every certificate atom to a valid source atom, and preserve strict positivity
of every required margin. In the single-pair branch this specializes to

    (overlap − debt)[E′,q′] ≤ (overlap − debt)[E,q],
    (overlap − debt)[E′,q′] > 0.

An exact theorem of this form would turn the orbit census into a reduction to
the single generic fixed body. Without it, the census is only a precise source
of candidate maps and counterexamples.

## Reproduction

```bash
python3 04-computation/lrc14_reflected_anharmonic_s3_robust_body_orbit_scout_20260801.py
python3 -O 04-computation/lrc14_reflected_anharmonic_s3_robust_body_orbit_scout_20260801.py --output /tmp/lrc14-anharmonic-O.out
cmp 05-knowledge/results/lrc14_reflected_anharmonic_s3_robust_body_orbit_scout_20260801.out /tmp/lrc14-anharmonic-O.out
```

Frozen digests are recorded in the stored output. The semantic digest is
`a9b174728c0749b12a5194d428f4ffc47c0556025ae05690c48b76991f0fd112`.

The best next hostile probe is not another orbit count. It is to choose the
highest-margin robust representative of each of `H`, `H2`, `B`, and
`E`, then attempt an explicit safe-cell/owner pullback and record the first
coordinate that fails. That test directly measures whether the missing theorem
needs only an owner sidecar or a genuinely different certificate currency.
