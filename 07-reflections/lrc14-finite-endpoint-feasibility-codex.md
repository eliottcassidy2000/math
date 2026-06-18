# LRC14 Finite Endpoint Feasibility

Codex 2026-06-18.

The useful result of this session is a negative one: the finite endpoint left by
the colored-resonance program should not be treated as a naive bounded
enumeration.

The raw count below `V=711` is already around `1.675e27` S3 coordinate shapes
before covering and primitivity.  That does not mean the endpoint is hopeless,
but it means the endpoint must inherit structure from the proof, the way the
proved `k=2` slice did.  In that slice, the computation became small only after
the drop-max safe-width theorem bounded the hard core.

The session also corrected a tempting shortcut.  The `q=14V` colored CRT grid
is the right large-`V` placement object, but not the whole low-`V` endpoint.
The MILP obstruction at `V=15`,

`(1,2,3,4,5,8,9,10,11,12,13,14,15)`,

has no `q=210` witness and yet exact `M=2/23`.  So the finite endpoint cannot
just ask for a `q=14V` residue.  Low denominators and exact THM-524 crossings
still matter.

The path forward feels more sharply shaped now:

- large `V`: colored reservoir plus resonance discrepancy;
- small `V`: classify q-grid obstruction families, then run exact-M only on the
  structurally reduced cores;
- bridge: generalize the k=2 drop-one-large-speed safe-width floor to k>=3, or
  replace it with endpoint-protection cycles / denominator-family pruning.

The tournament quotient over proof strategies was transitive, with
`k2_finite_core` on top and naive enumeration at the bottom.  That is the real
message: the successful endpoint object is not "all bounded rows"; it is "the
bounded residue after a theorem has already done most of the work."
