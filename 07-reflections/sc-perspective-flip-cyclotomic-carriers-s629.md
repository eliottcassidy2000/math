---
source: codex-2026-06-04-S629
status: theorem + finite audit + speculative synthesis
tags: [tournaments, self-converse, perspectives, unit-distance, cyclotomic, HYP-2205]
---

# SC perspective flips and cyclotomic carriers

The Claude snippet was useful, but only after separating three things that were
trying to occupy the same word "21":

```text
H=21     : forbidden Hamiltonian-path count / OCF evaluation.
n=21     : Moser unit-distance carrier size.
21 edges : triangular/Harborth lattice scalar echo at n=11.
```

Those are not the same object.  They are three shadows of carriers that share
the primitive-cube-root angle only after their side channels are collapsed.

## The SC Correction

The cleanest theorem from this session is THM-409:

```text
For a self-converse tournament T, Anti(T) is a coset over Aut(T),
and edge reversal induces a canonical involution on Aut(T)-vertex orbits.
```

This answers the user's question about "how nodes in an SC class get swapped
with each other, back and forth" with a correction:

```text
vertices do not have to swap back and forth;
rooted perspectives do.
```

S629's exact atlas catches the mistake.  Through `n=5`, every
anti-automorphism found has the expected fixed point plus transpositions
pattern at odd `n` or transpositions at even `n`.  But at `n=6`, the atlas finds
anti-automorphisms with vertex cycle type `(6,)`.  So an anti-map can rotate
six labelled vertices before returning.

The reason this does not break complement symmetry is that `sigma^2` is an
automorphism.  Once vertices are quotiented by `Aut(T)`, that square vanishes,
and the rooted-perspective map is a genuine involution.

The mathematical object is therefore:

```text
Aut^+- (T) = Aut(T) union Anti(T),
```

a `Z/2`-graded permutation group or, equivalently, an anti-automorphism coset
torsor acting on rooted perspectives.

## The Old Perspective Curiosity Repaired Again

S586 already corrected the old count:

```text
P(3)=4=U(4),
P(4)=12=U(5),
P(5)=48<U(6)=56.
```

The mistake was thinking the root alone was extension payload.  The extension
operation was root-blind.

S629 gives the self-converse version of the same warning.  Complement merging
is observer-blind unless the rooted perspective flip is retained.  The fixed
SC class is not a point.  It is a class equipped with a conjugation action on
its root perspectives.

That is exactly the pattern LRC keeps teaching: if the predicate is
observer-coupled, the quotient must keep the observer/root/fiber payload.

## H=7 and Unit Distance

The honest relation is:

```text
Phi_3(2)=7,
Phi_3(4)=21,
3*Phi_3(2)=21.
```

The unit-distance side sees the `60` degree/Eisenstein angle as edge-count
geometry.  The tournament side sees `7` and `21` as forbidden H evaluations.
The same integers appear, but in different semirings:

```text
unit distance : additive edge carriers with spine/bulk/frontier labels
tournament H  : Hamiltonian-path / OCF evaluation with strong-component labels
```

So `u(5)=7` is not a contradiction to H=7 impossibility.  It is:

```text
7 = 4 unit-spine edges + 3 tile/bulk edges.
```

Likewise the `21` lattice echo at `n=11` is:

```text
21 = 10 spine edges + 11 bulk edges,
```

while the exact unit-distance value is already `u(11)=23`.

The useful reading is not equinumerosity.  It is a warning: if a quotient makes
these objects equal by scalar alone, it has forgotten the side channel that
proves or blocks the predicate.

## What Ends After n=7

There are two meanings of "the Eisenstein integer pattern ends after n=7."

If it means the first closed centered hex shell, then yes: `1+6=7` is the last
moment where the construction is just "take the center and all first-shell
Eisenstein neighbors."  After that, each new point is a frontier choice, and
the relevant data are boundary ears, gains, traceability, and direction
support.

If it means the triangular/Eisenstein lattice stops being useful after `n=7`,
then no.  The triangular/Harborth machine continues to produce strong rows and
the S599z lattice traceability theorem follows it through `n<=28`.  But S626,
S627, and S628 show why this is no longer enough: non-lattice Moser carriers,
point-order gauge artifacts, and spine/bulk side channels are now the proof
surface.

So the right phrase is:

```text
the first-shell Eisenstein picture ends at n=7;
the retained-carrier Eisenstein/Moser recursion begins there.
```

## The n=21 Carrier

THM-408 makes `n=21` concrete:

```text
P_2^-:
  vertices = 21
  unit edges = 57
  unit spine = 20 edges
  bulk = 37 edges
```

This is not `H=21`.  It is a 21-vertex carrier that refuses the scalar-collapse
reading.  The recursion is:

```text
add one full Moser slab:
  +8 vertices
  +27 unit edges
  +8 spine steps
  +19 bulk edges
```

That `27=3^3` is a real cube-root-of-unity smell, but the exact proof object is
the slab word and its unit-shell bridge `(-1,1,0,0)`.

The next serious `n=21` question is:

```text
classify 57-edge 21-cores by endpoint-compatible ears.
```

Each ear should carry:

- gain packet;
- unit-spine compatibility;
- direction-support mask;
- canonical orbit class;
- deletion resilience;
- totally-unfaithful obstruction label.

This is the unit-distance analogue of keeping observer-source, endpoint owner,
and pressure labels in LRC.

## Roots of Unity Without Overclaiming

The primitive cube root of unity is genuinely present in:

- the triangular lattice angle;
- the Eisenstein integer norm;
- `Phi_3(x)=x^2+x+1`;
- the LRC `C=2n-1` triadic burden rows such as `C=21` and `C=27`;
- the complement/converse conjugation metaphor for SC tournaments.

But the `Cl_2(pi/3)` / `1.014` bridge is not yet a theorem here.  I would keep
it as a search light, not a load-bearing beam:

```text
proved: Sawin-style unit-distance arithmetic carriers amplify beyond n^(1.014)
proved: THM-409 gives the SC conjugation/perspective object
proved: THM-408 gives the Moser slab unit spine
open : a theorem equating the exponents or volumes
```

This is still valuable.  It tells us what to instrument:

```text
count norm/conjugation-fixed side channels,
not only scalar H values or scalar edge counts.
```

## Programmatic Next Moves

1. Extend `sc_perspective_flip_cyclotomic_s629.py` to sample or enumerate
   `n=7` SC classes, recording anti-cycle types, `sigma^2` automorphism
   orders, and rooted-perspective fixed/transposed counts.

2. Replace any SC script that treats a self-converse class as a scalar fixed
   point by a THM-409 perspective-flip fingerprint.

3. Build a `21`-core ear-extension ledger for the five exact `u(21)=57`
   unit-distance cores from S614/S626.

4. Test whether unit-spine traceability is preserved under every gain-3 ear
   from a retained `57`-edge 21-core, and identify the first obstruction to
   a `61`-edge `n=22` extension.

5. Track `Phi_3` events as labelled carrier events:

   ```text
   scalar value,
   spine/bulk split,
   rooted/conjugation payload,
   endpoint-owner or Moser-slab side channel,
   H-gap guardrail.
   ```

## Bottom Line

The breakthrough-shaped object is not:

```text
7 equals 7, or 21 equals 21.
```

It is:

```text
conjugation fixes the carrier only after the right observer/rooted lift,
and the forbidden scalar appears exactly when the quotient forgets the lift.
```

That is a promising common language for self-converse tournaments,
unit-distance Moser slabs, and LRC observer-source quotients.
