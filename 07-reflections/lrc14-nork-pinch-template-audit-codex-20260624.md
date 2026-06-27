# LRC14 NORK Pinch-Template Audit

The useful reframing this session is that HYP-2956's F6 bucket should not be
handled as a vague "covered" residue.  I treated it as NORK:

```text
No Open Residual Kernel
= qdiv >= 14, no strict safe open interval, and not AP/Goddyn-Wong.
```

The computation did not find a NORK packet.  The default run generated
`705940` rows and exactly classified `141351` `qdiv>=14` rows across AP,
one-swaps through `add<=420`, two-swaps through `add<=60`, three-swaps through
`add<=34`, and a first four-swap bank through `add<=24`.  The only zero-open
rows were AP and GW.  All other hard rows had positive open intervals and split
into named positive families:

```text
F2 positive unit-petal  28762
F3 positive K33         340
F4 positive q14 front   78651
F5 positive covering    33596
```

The sharper output is the pinch-template atlas.  For each positive row, the
script keeps the shortest strict safe interval, exact endpoints, active
endpoint owners, width, slack, q-class, atom keys, and C27-normalized labels.
This is the proof object I trust more than a scalar mass.  A tiny interval can
still be perfectly rigid if its endpoints are labelled.

The recurring templates are now concrete:

```text
13L -> 12R
14L -> 13R
11L -> 16R
7L  -> 20R
5L  -> 36R
```

The fixed-margin swap-chain paper arXiv:2606.22636 was useful only as proof
architecture: local moves first, then split scalar count sectors from
non-scalar Johnson-harmonic sectors.  In LRC language, the scalar sector is
`qdiv/M/Haar`, while the non-scalar sector is endpoint owners, C27/unital
labels, K33/state-lift address, and boundary-moment route.  HYP-2966 says the
non-scalar sector must be carried until it discharges.

Two older prompts also fit here.  The recombination-over-real-factors paper
arXiv:2410.15880 suggests a healthy analogy: do not solve a polynomial by
staring at one coefficient; factor into local pieces, then recombine with the
right integer labels.  The Euler divisor-sum translation arXiv:math/0411587
has the same shape after the logarithmic derivative: product data becomes a
labelled divisor recurrence.  For LRC14, the analogous derivative is the
endpoint-owner boundary ledger.

Beurling-Selberg minorants and Paris-Harrington rank remain good guardrails:
minorants warn that positivity should be certified by a constructed lower
carrier, while PH rank warns that a finite classifier can still need a
well-founded extension coordinate.  The NORK pinch theorem is exactly such a
candidate lower carrier.  It says:

```text
Every primitive non-AP/GW AP-source-core packet creates a positive
endpoint-owner pinch template, unless it constructs HYP-2908/THM-572
or a genuinely new F7 Johnson-harmonic sector.
```

Assumption challenge: tournament vertices are not runners here.  I used proof
carriers and pinch templates as vertices.  This preserves open-vs-boundary
status, owner labels, and state-lift visibility, while destroying raw runner
identity inside already-positive families.  That is the right loss: if a row is
already positive, its job is to show which labelled pinch made it positive.

LRC14 is not proved.  But the F6 target is now smaller and nastier in the
right way: a counterexample must evade qdiv, evade AP/GW boundary rigidity,
evade unit/K33/covering pinch templates, and still avoid producing a new
state-lift sector.  That is a much better enemy.
