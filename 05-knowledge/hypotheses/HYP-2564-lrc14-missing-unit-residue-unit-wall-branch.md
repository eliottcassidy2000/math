# HYP-2564 - Every one-14-multiple row that misses a unit residue has some unit-wall branch witness

**Status:** OPEN, with strong bounded exhaustive evidence.

**Source:** codex-2026-06-17-S1. Extends THM-540's exact first-branch lemma.

**Computation:** `04-computation/lrc14_perturbed_unit_wall_codex.py`; stored output
`05-knowledge/results/lrc14_perturbed_unit_wall_codex.out`.

## Statement

Let

```text
S = A U {14m}
```

be a primitive 13-speed row with exactly one multiple of `14`, where `A` misses at
least one unit residue modulo `14`. Then:

> **Conjecture.** `S` has a strict witness on some unit wall branch
>
> ```text
> t = a/14 + delta,
> ```
>
> with `a in {1,3,5,9,11,13}` and `|delta| < 1/14`.

This is the global version of the user's perturbation idea. THM-540 proves the
exact **first branch** whenever

```text
1/(196m) < min_v (13 - a v mod 14)/(14v)
```

or its negative-side analogue. HYP-2564 says later branches should finish the job
for the remaining unit-incomplete rows.

## Why this matters

If HYP-2564 is true, then the one-multiple `C'(14)` search peels sharply:

```text
unit-incomplete A  -> unit-wall branch witness
unit-complete A    -> genuine residual
```

This is exactly where the new exact low-L rows already live:

```text
{1,2,...,11,13,36}
```

is the THM-522 / HYP-2561 extremizer and it is **unit-complete**. So the bounded
exact extremal story and the unit-wall branch story point to the same hard core.

## Evidence

1. **First-branch theorem (THM-540).**
   The first branch is elementary and exact. It already certifies most bounded
   rows with missing unit residue.

2. **Bounded exhaustive all-branches sweep.**
   `lrc14_perturbed_unit_wall_codex.py` does a full breakpoint search near the six
   unit walls and finds:

   - all `1372/1372` rows `S = A U {14m}` with `A subset [1,16]\\{14}`, `|A|=12`,
     `A` missing a unit residue, and `m <= 4`, have a unit-wall branch witness;
   - all `15666/15666` rows with `A subset [1,18]\\{14}`, `|A|=12`, `A` missing a
     unit residue, and `m <= 3`, have a unit-wall branch witness.

   Zero misses in either census.

3. **AP one-drop families.**
   For `{1,...,13}\\{j} U {14m}`, THM-540 already discharges every odd-unit drop
   `j in {1,3,5,9,11,13}` at `m = 1`. The only residual AP one-drop families are
   exactly the unit-complete ones (`j` even or `7`).

4. **A useful counterexample to overclaiming the theorem.**
   The row

   ```text
   (1,2,3,4,5,6,7,8,9,10,11,14,26)
   ```

   misses unit residue `13`, but THM-540's first branch fails. Nevertheless the
   exact breakpoint search finds a **later** branch witness

   ```text
   t = 657/8008 = 1/14 + 85/8008.
   ```

   So the strong conjecture survives, but only after allowing later wraps.

## Tests / next steps

1. Prove a second-branch / finite-branch extension of THM-540. The computational
   counterexample above shows exactly where the first branch stops being enough.
2. Push the exhaustive all-branches census past `[1,18]` and `m <= 3` using cached
   breakpoint sets by residue class.
3. Compare the branch count needed for a witness against residue structure:
   which missing unit residues force first branch, which force later branch, and
   which correlate with the new low-L extremal families?
4. Connect this to THM-522's exact extremizers. If every unit-incomplete row opens
   on some unit wall, then the singular-series/measure hard core really does lie
   in the unit-complete rows.
