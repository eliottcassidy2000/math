# Spectral Dictionary Compatibility

*codex-2026-06-28.  This reflection records HYP-3205, a synthesis pass over
HYP-3202 covariance layers/traps, HYP-3203 polarized cyclotomic support,
S31al/S73d Toeplitz moment work, and HYP-3210's Joukowski bridge.*

## The Core Shift

The frontier is no longer asking for one better scalar.  The exact k=8 bank
says the AP row is isolated by a dictionary:

```text
layers + AP support + Toeplitz margin + trap/orbit sidecars.
```

That is the useful abstraction.  A decoy can be close to AP in one language,
but it leaks in another.  The closest non-AP row
`(0,2,3,4,5,6,7,8)` is close enough to be interesting, yet it has visible
deficits in all eight measured dictionary coordinates.

This makes the remaining proof look more like a finite separation theorem
than a monotone-flow theorem.

## The Creative Synthesis

The old pieces now have cleaner roles:

- Covariance layers are the local finite coordinates.
- AP-polarized support is the root-locus separating hyperplane.
- Toeplitz lambda-min is the moment-cone interiority margin.
- Perron alignment explains the ferromagnetic leading mode but is not a
  terminal inequality.
- Joukowski/Hermite-Biehler is the transport bridge to the tournament
  real-axis language.
- Exchange and compression traps are not failures; they are test stalks for
  a sheaf of certificates.
- Incoming HYP-3204 ordered-tail exchange is the coefficient chart for
  `L_y`: it prices `q3` central mass against `q0+q6` bimodality, while this
  packet asks that the chart stay compatible with layer/support/Toeplitz
  certificates.

The circuit-complexity analogy becomes concrete: a proposed proof must name
its input basis.  If the basis is only `p0` or raw root spread, the circuit is
missing inputs.  If the basis is the dictionary vector, the traps become
addressable.

## What This Suggests

Try to prove a small certificate-Helly lemma:

```text
AP-tight(layer bundle) cap AP-tight(support) cap AP-tight(Toeplitz)
= AP/dilation.
```

Then prove primitive uniqueness by removing dilation.  This is compatible
with the `14 -> 7 -> 2` descent: each lift should transport the certificate
vector, not just the witness-time scalar.

The wild but useful picture is a moment-lattice body.  The q-vector lives in a
rational slab arrangement; the seventh-cyclotomic ideal is a cubic irrational
direction; AP is the rational support point exposed by that direction.  The
Minkowski/Bravais version would prove every non-AP rational relation vector
exits at least one certificate slab.

## Guardrail

Do not promote raw cyclotomic norm, entropy, exact `1/7`, or Perron alignment
to terminal proof targets.  The scout explicitly demotes them.  They remain
diagnostics or sidecars.  The terminal candidate is compatibility across the
dictionary.
