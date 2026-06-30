# The bridge splits by level: the deeper decorrelations are genuine, the top level IS the existence question

*klein-2026-06-29-S18. Working the bridge from the apex skeleton (THM-590 cyclotomic gap) to the full per-level rho_j (THM-580). The honest answer: there is no single measure-inequality bridge; the bridge splits by descent level, and the obstruction localizes entirely to the top level, where it equals the original LRC question.*

## What I computed

For real descent chains (consecutive prefixes, the tightest coverings, 3000+ random coverings) I computed
the ACTUAL per-level `rho_j = meas(lonely S^{(j)}) / [meas(lonely O_j) . meas(lonely S^{(j+1)})]` (THM-580)
and correlated it with the apex cyclotomic gap `g(O_j mod 7)` (THM-590) and the 2-sheet `CV(N2_{O_j})`.

## Finding 1 — the full measure-`rho_j` is NOT bounded by the apex skeleton

`rho_j` does NOT track the apex gap. Conditioning min `rho_j` on the gap value:
```
g = 0.000 (full Z_7) : min rho = 0.00     g = 1.000 (singleton/co-singleton) : min rho = 0.19
g = 0.198 (DOUBLET)  : min rho = 0.44     g = 2.000                          : min rho = 0.71
g = 0.308            : min rho = 0.68
```
The doublet (the skeleton's binding core, `g = 4cos^2(3pi/7) = 0.198`) does NOT bind `rho_j` (its `rho >= 0.44`).
The binding is the OTHER end: conditioning on cusp-proximity `|O_j mod 7|`, min `rho_0` falls monotonically
`0.91, 0.65, 0.53, 0.33, 0.14, 0.10` as `|O_j mod 7| -> 7`. **The full `rho_j` binds at the cusp end
(`|O7| -> 7`), opposite the skeleton (the doublet, `|O7| = 2`).** So the apex cyclotomic gap is NOT a lower
bound for the measure `rho_j`: the measure goes below `0.198` (down to `~0.10` and `-> 0`) where the gap is
large. The skeleton governs the discrete/cyclotomic content; it does not control the measure.

## Finding 2 — the obstruction localizes to LEVEL 0, and level 0 IS the existence question

min `rho_j` BY descent level:
```
level 0 : 0.05  (-> 0 at the boundary)      level 2 : 0.77
level 1 : 0.56                              level 3 : 1.00
                                            level 4 : 1.07
```
The DEEPER levels (`j >= 1`) are bounded away from 0, and the bound INCREASES with depth (toward `rho -> 1`,
independence): deeper cores are smaller, lonely-richer, and overlap more. **THM-580's per-level floor
`rho_j >= c` holds for `j >= 1` (verified `>= 0.56`).** Only the TOP level vanishes -- and it vanishes for a
structural reason:
> `rho_0 = meas(lonely S) / [meas(lonely O_0) . meas(lonely S')]` -- the numerator IS `m_S`.
Verified along coverings approaching the boundary: `rho_0 = 0.34, 0.23, 0.19, 0.17, 0.16, ... -> 0` exactly
as `m_S -> 0`. So `rho_0 -> 0` is not a new decorrelation failure -- it is the tautology that the top-level
odd/even overlap is the original lonely set. `rho_0 > 0 <=> m_S > 0 <=> LRC(S)`.

## The honest bridge

There is no single inequality `rho_j >= skeleton`. The bridge splits:

- **Deeper levels (`j >= 1`): genuine bounded measure-decorrelations.** `rho_j >= 0.56` (increasing with
  depth). Here the descent has real content -- the odd part and the descended set decorrelate, with a
  uniform floor. The apex skeleton is comfortably satisfied (these `rho_j` exceed the cyclotomic content).

- **Top level (`j = 0`): the existence question itself.** `rho_0 = m_S/(...)` carries the original `m_S`;
  it vanishes at the boundary (measure), so no `m_S`-independent measure bound can hold. The bridge here is
  the MEASURE -> EXISTENCE passage (HYP-3597): one needs `rho_0 > 0` (the overlap is NONEMPTY), not
  `rho_0 >= c` (a measure bound, which is false). And `rho_0 > 0 <=> LRC(S)` -- the descent restates, does
  not reduce, the top level.

So the descent's value is precise and limited: **it peels off the easy deeper levels (genuinely bounded
decorrelations) and leaves the top-level odd/even overlap, which is the original existence question.** The
apex skeleton (THM-590) is the right object for that top level, but ONLY in the existence sense -- it must
certify that the overlap is nonempty, not that its measure is bounded. The measure `rho_0` provably vanishes
at the cusp.

## What this corrects

- **HYP-3600 consequence 3** ("the floor as a bounded product `meas(lonely S) >= (4cos^2(3pi/7))^d . cap^d`")
  is FALSE as a measure statement: `rho_0` is not `>= 4cos^2(3pi/7)`; it `-> 0` at the cusp. The deeper
  factors are bounded, but the top factor carries `m_S`. The bounded-product picture holds only for the
  discrete/existence skeleton, not the measure.
- **THM-580 reduction (a)** (`rho_j >= c` uniform) holds for `j >= 1` but FAILS at `j = 0` (`rho_0 -> 0`).
  The per-level measure floor is genuine below the top and tautologically the original question on top.

## The picture

The bridge from skeleton to full `rho_j` is not one inequality but a split: the descent's deeper levels are
genuine, uniformly-bounded measure decorrelations (the skeleton is satisfied with room to spare), and the
top level is `m_S` itself -- the original existence question, where the measure `rho_0` vanishes at the cusp
and only the existence statement (certified, if at all, by the discrete cyclotomic skeleton) can carry the
floor. This is the same measure/existence dichotomy as the global floor (HYP-3597), now localized to the top
of the descent. The descent does not dissolve the hard core; it isolates it.

See HYP-3599 (this), THM-580 (the descent), THM-590 (the skeleton), HYP-3597 (measure vs existence),
HYP-3600 (the finite families; consequence 3 corrected), HYP-3575/3576 (mac-mini's rho_j reduction).
