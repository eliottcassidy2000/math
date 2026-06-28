# LRC14 C6 Residue-Magnitude Factorization

codex-2026-06-28 reflection for HYP-3310.

## What Changed

The user's formulation gives a cleaner way to name the two pieces that were
being mixed:

```text
binding runners = units = C6/C3 cyclotomic contact skeleton
covering runners = evens + apex 7 = nonunit CRT layer
12->24 = 2-adic magnitude hinge inside the even cover branch
```

The exact scout confirms the residue split:

```text
units = (1,3,5,9,11,13)
even cover = (2,4,6,8,10,12)
apex7 = 7
C3 binder orbit = (1,13) -> (3,11) -> (5,9)
```

It also records the important non-overclaim:

```text
12 -> 24 is not residue-preserving mod 14.
12 mod 14 = 12, 24 mod 14 = 10.
```

So this is a 2-adic/magnitude hinge, not a pure CRT quotient.  That makes it
compatible with S256b's warning that the census is magnitude-level even though
`Q(sqrt(-7))` organizes the residue skeleton.

## Proof-Route Use

The split suggests a proof program with four load-bearing lemmas.

1. Binding/contact lemma: prove one complement-pair obligation, then transport
   by `C3` to the other two slots.
2. Covering-floor lemma: split nonunits into even cover and ramified apex-7,
   then prove strict floor/off-grid positivity outside named equality cells.
3. Magnitude-hinge lemma: isolate `12->24` as the only integer equality flex
   in the covering manifold.
4. Observability lemma: glue these packets through HYP-3300, forbidding any
   quotient that forgets residue, v2 magnitude, apex ramification, contact
   graph, endpoint owner, or off-grid floor route.

This is not a new scalar proof.  It is a coordinate system for preventing the
existing proof routes from forgetting the wrong layer.

Incoming HYP-3266 makes the handoff more concrete: this packet is mostly an
O15 tight-locus rigidity coordinate dictionary, with an O16 `Q(sqrt(-7))`
signed-floor sidecar guarded by the 2-adic magnitude layer, and an O12
off-grid bulk route whenever the unit skeleton is killed or no longer global.

## Assumption Challenge

I explicitly considered runners, gaps, fixed sections, section boundaries,
wall crossings, residue fibers, cover arcs, Fourier/cyclotomic modes, matroid
contact circuits, and proof obligations.  The chosen tournament vertices are
proof obligations / sidecar columns because raw runners preserve too much
irrelevant identity and too little proof status.

The quotient preserves:

```text
unit-contact certificate
covering-floor route
2-adic magnitude sidecar
apex ramification flag
finite chamber / named-debt exit
```

It destroys exact heights and off-grid safety profiles unless those sidecars
are retained.  This is why `Q(sqrt(-7))` is an organizer, not the terminal
proof object.

## Tournament Readout

Tournament Analysis uses proof carriers as vertices.  The fingerprint is
transitive:

```text
vertices = 9
directed_3cycles = 0
hamiltonian_path_count = 1
path =
  observability_morse_glue
  -> c6_unit_group
  -> seven_adic_residue_skeleton
  -> two_adic_magnitude_layer
  -> jacobsthal_doubling_hinge_12_24
  -> ramified_apex7_cover
  -> c3_binding_slot_orbit
  -> c2_qr_nqr_conjugation
  -> raw_runner_partition
```

The ordering is useful mainly as a guardrail.  The `C6` skeleton is strong
only after the global observability/Morse glue names the magnitude and
covering sidecars.  Raw runner partition is last because it has exact names
but almost no legal forgetting rule.

## Next Pull

Build a small finite chamber interface where rows are classified by:

```text
surviving C3 unit-contact graph
even-cover / apex-7 nonunit branch
v2 magnitude and hinge status
strict off-grid floor witness
AP/GW equality or named residual debt
```

Then test whether the S256b covering-flex manifold can emit the same fields.
If yes, HYP-3310 becomes the coordinate dictionary joining the local
equioscillation proof to the global covering census.
