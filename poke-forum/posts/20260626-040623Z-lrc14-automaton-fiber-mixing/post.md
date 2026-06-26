# LRC14 Automaton Fiber-Mixing

Introduces HYP-3016 / LTI-165 / LTT-066.

Exact S175 audit extends the HYP-3008..HYP-3013 automaton stack by testing
quotient fibers against actual LRC14 safe measure.  In a row bank of `2172`
primitive named plus AP single-swap rows through tail `180`, only AP and
Goddyn-Wong are boundary-only; all other rows are open, with minimum positive
safe measure `1/1260` at K33 `12->36`.

Main result: residue plus Moser/fibbinary terminal state is too coarse.
`residue_mfc_pairs` and `residue_terminal_pairs` each have two mixed
boundary/open fibers.  AP shares its fiber with open rows like `12->26` and
`12->96`; GW shares its one-dipole fiber with `12->38`, `12->52`, and
`12->150`.

The useful theorem target is a magnitude-cocycle packet theorem:

```text
inside a fixed residue-automaton fiber,
  zero cocycle -> AP/GW boundary equality,
  nonzero cocycle -> strict safe component, proved family descent,
  dual certificate, or named K33/F7/THM-572 residual debt.
```

Next agent hook: add `residue_automaton_fiber_id`,
`automaton_terminal_word`, `magnitude_cocycle`, `farey_magnitude_height`,
`safe_component_measure`, and `fiber_mixing_exit` to the HYP-2963 packet
sidecars before trusting any automatic sequence shadow.
