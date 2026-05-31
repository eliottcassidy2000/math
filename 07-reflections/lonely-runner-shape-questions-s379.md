# Lonely Runner Shape Questions S379

**Session:** codex-2026-05-31-S379  
**Script:** `04-computation/lonely_runner_shape_questions_s379.py`  
**Stored output:** `05-knowledge/results/lonely_runner_shape_questions_s379.out`

This session was deliberately about reframing.  The prompt asked: what do the
things represent, what is most fundamental, and what shape determines the
problem?

The best answer I found is:

```text
The speed set is not the object.
The pulled-back circular-arc endpoint incidence system is the object.
```

Speeds generate the system, but the proof lives in the boundary structure of
the unsafe arcs.

## Translation Table

The useful dictionary is:

```text
speed v            -> character t -> v*t on R/Z
speed set V        -> family of characters pulling unsafe arcs back to time
threshold 1/(k+1)  -> Dirichlet equality radius
forbidden interval -> unsafe fiber on the time circle
endpoint           -> exact equality event
protected endpoint -> equality event hidden inside another unsafe interior
boundary-only set  -> open cover fails only at endpoints
counterexample     -> full open cover with every endpoint protected
scalar ramp        -> Dirichlet equality spine / gauge orbit
torsion defect     -> quotient leak through a proper subgroup
```

This makes the Lonely Runner Conjecture feel less like a metric optimization
problem and more like a finite circular-arc incidence theorem.

## The Small Surprise

The seven-ladder near-disproof has a famously tiny visible gap at the threshold:

```text
gap/th = 0.005411
```

But the exact max-min loneliness height is not tiny:

```text
critical/th = 1.217391
```

So the seven-ladder is close horizontally, in gap-width space, but not close
vertically, in loneliness-height space.  It also has `84` unprotected endpoints.

That is a useful warning: a tiny complement interval can be a steep valley, not
a near-counterexample.  The endpoint graph sees this immediately.

## Tight Means Endpoint Failure

In the sample, the tight examples all have:

```text
critical/th = 1
boundary-only classification
unit skeleton = true
```

The initial segments and the known sporadic `n=8` tight example are not
distinguished by having a positive interval gap.  They are distinguished by
where the open cover fails: at unit endpoints.

This reinforces the S359-S362 unit-boundary program.  A counterexample would
not merely need to close a gap.  It would need to protect every endpoint after
the exact open intervals are laid down.

## What Is Most Fundamental?

My current answer:

```text
finite circular-arc incidence + endpoint peeling
```

The interval nerve is useful, but probably not enough: interval-overlap density
does not by itself distinguish tightness from near-disproof pressure.  The
endpoint graph remembers owner multiplicity, protector degree, unprotected
leaves, and peel depth.  Those are exactly the failure modes a proof needs.

This also clarifies scalar ramps.  They are not just annoying examples.  They
are the complete/equality spine of the endpoint-incidence system.  The real
problem begins after quotienting that spine and asking how residue repairs
change the endpoint graph.

## Better Questions

The session generated questions that feel more fertile than "try more speeds":

1. Can every primitive LRC endpoint-protection graph be peeled by a potential
   depending only on cyclic order and owner/protector degrees?
2. Are abstract all-protected circular-arc endpoint systems unrealizable by
   integer speeds?
3. Can critical radius, peel depth, and repair deficit be unified as one Morse
   function on a circular-arc metagraph?
4. Can the eight `n=14` alpha stencils be characterized as endpoint leaf types,
   without mentioning alpha-cells?
5. Is product-sum target status a coordinate-level shadow of the same
   additive/multiplicative split visible in endpoint protection?
6. What is the right homology: interval nerve homology, endpoint-core peeling,
   or repair-graph cycles?

The cleanest next proof question is:

```text
Why must the endpoint-protection graph have a leaf after the scalar equality
spine is quotiented?
```

That question has the right flavor.  It points to a finite object, a local
obstruction, a descent algorithm, and a possible certificate.
