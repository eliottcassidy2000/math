# Actual one-arc controls on the same filler order

**FINITE-EXACT bank, PROVED lifting.** This sidecar treats local edits, not merely ordinal recipes that replace a filler by a different size. Keep the factor-seven exterior unchanged. An ordinally embedded block can be edited internally, leaving all its arcs to every other block and to the exterior fixed.

A **pure width** block has a marked edge whose reversal changes \(D\) but not \(Q\); a **pure translation** block changes \(Q\) but not \(D\). The two kinds of control therefore act on the three readouts by
\[
(\delta H_-,\delta C,\delta H_+)
   = \delta D (19,21,23)
       \quad\text{or}\quad \delta Q(840,840,840).
\]
The derivatives here are exact finite differences. For odd \(p\nmid840\) they span the plane of arithmetic progressions modulo \(p\). More importantly, each edit is a physical arc reversal and preserves both marked fork arcs. Two edits in different blocks commute on the actual labelled host; their square of orientations does not rely on projecting down to scalar counts. For a pure translation, the four vertices obtained by toggling that arc and either one of the marked exit arcs form a literal commuting square with the same signed exit gap on both sides.

Take for example these masks, with the least significant bit for vertex zero:
- Pure width: \((198,76,184,81,35,75,148,58)\), arc \(\{0,1\}\): \((D,Q)=(653,8816)\longleftrightarrow(659,8816)\).
- Pure translation: \((42,12,9,16,39,14)\), arc \(\{0,3\}\): \((37,352)\longleftrightarrow(37,376)\).

All reads are units at 5. The respective ratios for the changing coordinates are \(3,13\bmod25\); each generates the units at every power of 5. The latter rotation changes \(Q\bmod5\) but the physical center only sees that change at the next 5-adic digit, because \(5\mid840\) exactly once.

In general take \(L=\varphi(p^a)-1\) copies of each pure block, all initially in its first state. If the two changing ratios generate \((\mathbb Z/p^a)^\times\), independently toggling \(i\) blocks of the first type and \(j\) of the second realizes every unit pair modulo \(p^a\), as \(0\le i,j\le L\). The ordinal size, the physical fork exterior and the two possible marked exit arcs all stay fixed throughout this experiment. Finding \(i,j\) is just two discrete logarithms relative to the initial pair. The standard lifting test for an odd prime is that a ratio generates modulo \(p\) and its \((p-1)\)st power is not one modulo \(p^2\); binomial expansion then multiplies its order by \(p\) at every step. For principal units alone, the analogous test is value \(1\bmod p\) but not \(1\bmod p^2\).

At \(p=3\) use \(3^{a-1}-1\) copies per type for \(a\ge2\), and obtain the full product of the two principal-unit cosets based at the all-first filler. The two absolute residues modulo 3 cannot be changed independently at fixed order because \(Q\equiv(-1)^nD\bmod3\).

The included checked masks give the two generators for **every prime \(5\le p<100\)**, and both principal-unit directions at \(p=3\). Thus for \(11\le p<100\) the same-order controller can target any center and gap residues for which the solved \(D=g/2,\ Q=(c-21D)/840\) are units. This is more restrictive than arbitrary residues, but stronger physically than existence of separately synthesized fillers: a single bank of marked arcs and a fixed order suffice at each prime power. The distinction from universal persistence is deliberate: a state centered at zero modulo an odd \(p\) with \(D\) a unit exits to opposite nonzero residues \(\pm2D\).

This does not provide an all-prime same-order bank, and does not identify the exact integer count with a modular target. All toggles have been certified on their actual small modules; the exterior conclusion uses the existing exact two-read formula rather than subset enumeration at the much larger composed orders.

Replay: `python3 seven-two-read-check.py`.
