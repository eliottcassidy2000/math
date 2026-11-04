# The complete assembly signature still loses a later named flip

2026-11-03. **FINITE-EXACT explicit counterexample, with independent
permutation and subset-path checks.** No minimal-order claim is made.

The physical six-component return produces a useful internal hostile to
using its seven-coordinate evaluator as an arbitrary-operation state. At
steps8 and9 of [the actual schedule](SIX_PROBE.md), take the literal suffixes

```text
B = (66,68,121,99,107,67,0),
C = (66,68,121,67,107,75,0),
W = {2,4,6} (mask84).
```

They differ only by reversal of the named arc `3->5`, yet their complete
ordered-component signature is identical:

```text
sigma(B,W)=sigma(C,W)=(17,448,259,217,112,17,9).
```

The same available named action, reversal of `4->5`, produces

```text
sigma(B',W)=(27,576,349,223,114,27,9),
sigma(C',W)=(17,448,269,167,84,17,5).
```

Every displayed suffix has terminal singleton6 marked. Consequently all
four full attachments are strong at every nonempty prefix clock `u>=1`.
Initially both have `H_u=-95+259*2^u` and physical root defect112. After
the identical named action, their counts differ by `80*2^u-20` and their
defects differ by30. At `u=1` the common input count423 becomes611 versus471;
the defects become114 versus84.

The first failed implication is therefore explicit: equality of all seven
assembly observables does **not** imply equality after an active-edge flip.
The lost coordinate is a mixed named response: arc `3->5` is invisible to
the current observer but becomes visible after changing `4->5`. The exact
ordered-sum law survives unchanged; it was never a law for arbitrary flips.

The repair is to retain the [paired endpoint forest](../transport/PROOF.md)
with these names active, together with their actual orientation. Its named
join operators keep the mixed response before either flip. Equal current
Hamilton counts, equal root defects, and even the full seven-coordinate
assembly signature cannot replace that state when proposing a reusable
physical grammar.

The remaining compression question is the smallest additional named-response
state closed under the intended preparation and return words. The current
seven-coordinate equality alone cannot answer it.

Reproduce with

```text
python3 04-computation/multiport_grammar_20261103/observer_hostile.py
python3 -O 04-computation/multiport_grammar_20261103/observer_hostile.py
```

The checker independently evaluates all four suffix signatures by subset DP
and the auditor's permutation/cut routine. It also checks all eight full
attachments at clocks1 and2, including actual root readouts and reachability.
The all-clock statement follows from the inherited attachment identity and
the invariant terminal port, rather than extrapolating those eight cases.
