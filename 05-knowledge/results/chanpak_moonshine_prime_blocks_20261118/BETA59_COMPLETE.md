# One entire signed lane, with the genuine forbidden destination exposed

For every integer m≥1, let x=2^m+59 and choose the sign of 3x±1 that permits division by at least four, then take its odd part y. There is a **literal merged-metagraph route of at most two black lines** from a tournament with count x to one with count y for every m except 4. At m=4 the source count 75 is attainable, but its signed destination is 7, genuinely forbidden. The theorem does not assert that intermediate line counts always decrease or that any chosen representative of x works.

The uniform large-width proof is the [rim exchange in the physical inverse atlas](INVERSE_METAGRAPH.md): put w=m−6≥2. Its endpoint relation is (3x−1)/16=12·2^w+11, and two explicit black lines suffice. The only remaining cases use different source graphs, certified by literal frames:

| m | Actual path counts along merged lines |
|---:|---|
| 1 | 61→23 |
| 2 | 63→171→47 |
| 3 | 67→105→25 |
| 4 | 75 attained, but prescribed terminal count 7 is forbidden |
| 5 | 91→65→17 |
| 6 | 123→23 |
| 7 | 187→35 |

The second row fills the uncapped hole from the low-height atlas without contradicting its all-order exclusion below height 74. The intermediate count 171 permits the change of graph needed to re-enter height 47. The fourth row cannot be repaired by allowing overshoot, since no endpoint exists. In fact neither signed shortcut can ever have 21 as its destination, because 3x±1 and its odd part are prime to 3; 21 can only obstruct as a source in this selector. These facts distinguish a real missing endpoint from a missing one-line realization.

Replay: `python3 beta59_boundary.py`.
