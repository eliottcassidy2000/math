---
id: HYP-1782
status: EXPLORATORY
source: codex-2026-05-30
related:
  - THM-343
  - THM-344
  - HYP-1781
---

# HYP-1782: Single-Core Target Spectrum

## Statement

For a single-core tournament obtained by inserting one core vertex into a
transitive tournament, the complete-Omega target value is governed by the
weighted signature statistic

```text
r_core(s) = sum_{1...0 pairs} 2^{max(gap-1,0)}
H = 1 + 2*r_core(s)
```

The image of `r_core` appears to exclude the small target residues needed for
complete-core H=7 and H=21:

```text
r_core = 3   absent through m=40
r_core = 10  absent through m=40
r_core = 31  first appears at m=7
```

Thus the complete-core mechanism cannot realize H=7 or H=21, while it does
realize H=63 exactly when the n=8 complete-Omega classes appear.

## Evidence

`projection_defect_bridge_s12.py` reports the target-signature search through
`m=40`:

- `r=3`: absent;
- `r=10`: absent;
- `r=31`: first appears at `m=7`, examples `1001100` and `1100110`;
- `r=42` and `r=63`: first appear at `m=8`.

The two THM-344 n=8 H=63 classes are single-core, with signatures:

```text
1001100
1100110
```

Both have `r_core=31`, `Omega=K31`, and `H=63`.

## Why This Matters

The older H=63 obstruction targeted a disconnected Omega shape.  The n=8
counterexamples show that complete Omega unlocks 63.  This hypothesis asks for
the complementary theorem:

```text
complete-core Omega still cannot unlock 7 or 21
```

If proved, H=7 and H=21 become sharper: not only are some multiplicative
Omega-shapes impossible, but the most concentrated complete-core shape also
misses their required target residues.

## First Proof Target

Characterize the image of `r_core` modulo small targets.  A minimal useful
lemma would prove:

```text
r_core(s) notin {3,10}
```

for every binary signature `s`.

## Sources

- `04-computation/projection_defect_bridge_s12.py`
- `04-computation/omega_extreme_fingerprints_s11.py`
- `07-reflections/h63-unlocks-as-complete-omega.md`
- `07-reflections/residue-calculus-feedback-loop.md`
