# Perspective Depth And Controlled Forgetting - S214

Follow-on to HYP-3047/T1129: that thread isolates the `n=6`
A000568/rooted-perspective defect as incident/cross-coupling payload.  This
note records the stricter exact edge/triple carrier scout around the same
first-defect window.

It is also downstream of HYP-3049/T1131, which identifies the exact
old-root/new-observer ordered-pair extension state and its directed-edge
quotient.  Here the point is the broader carrier table: edge, triple, cyclic,
and future conflict perspectives as sidecars.

It also follows the HYP-3048/T1130 matrix-atlas update: the perspective
carriers below should be encoded as sidecar-observability rows/columns, not
collapsed to spectra or scalar matrix summaries.

The old perspective coincidence fails exactly where it should if the
controlled-forgetting ladder is the right mental model.

Write `P(m)` for the number of exact rooted node perspectives:

```text
P(m) = sum over tournament classes on m vertices of vertex-orbits.
```

The familiar small values are:

```text
P(2)=2=A000568(3)
P(3)=4=A000568(4)
P(4)=12=A000568(5)
```

The first failure is:

```text
P(5)=48,  A000568(6)=56.
```

So the answer to the "what n" question is: base perspective size `m=5`, or
next A000568 index `n=6`.

The k-depth version clarifies what the failure is not.  Using directed
neighbor-refinement colors, the exact scout gives:

```text
m exact d1 d2 d3 d4
4 12    10 12 12 12
5 48    36 47 48 48
```

At `m=5`, depth `3` already recovers all exact rooted node perspectives.  The
missing `8` states in `A000568(6)` therefore are not deeper node neighborhoods.
They are observer-extension data: how the added vertex cuts across the rooted
`m=5` body.

That is exactly the controlled-forgetting lesson.  A quotient is safe only
after the hidden coordinate has a name.  Here the chain is:

```text
unrooted class
-> k-depth node view
-> exact rooted node view
-> observer-extension cut view
```

Node depth gets you to the third rung.  The first A000568 defect says the
fourth rung is genuinely new.

The non-node perspectives look useful as sidecars:

```text
m=5 exact arc perspectives        88
m=5 exact triple perspectives     88
  transitive triples              64
  directed cycles                 24
```

The edge perspective is the cleanest creative carrier.  A directed edge has a
tail and a tip, and every outside vertex lies in one of four sectors:

```text
tail->w, tip->w
tail->w, w->tip
w->tail, tip->w
w->tail, w->tip
```

That is a two-plate object, just like a residual capacitor or a pair-good
binding switch.  The cyclic-triple perspective is the local shadow of
`Omega(T)`, and the next carrier should be disjoint cycle-conflict
perspectives at `m=6`.

For LRC, the warning is precise: the source-perspective slice remains exact by
THM-381 because deleting a source is canonical.  The full observer-coupled
problem is harder because a general observer is not only a rooted node; it is a
rooted cut into the runner tournament.  Edge, cycle, and conflict perspectives
are candidate sidecars for the non-source residuals, but none of them should be
used as a raw replacement for A000568.
