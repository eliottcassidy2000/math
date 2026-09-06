# Clean-room referee audit: THM-4451

**Verdict:** the two sharp strict-component theorems are correct after two
minor proof-text repairs below.  This audit used a self-contained exact wall
decomposition and did not import either author program.

## Verified statements

- The lane capacity is
  `cap_n(L)=[2m/7+min(r,2/7)]/(2n)` for `2nL=m+r`.  Subtracting the mean
  contribution `2L/7` gives a surplus in `[0,5/(49n)]`.
- On a physical failure interval two distinct odd-tail lanes are active, so
  `sum s_n(L) >= 8L/7`.
- Successive use of `s_n<=5/(49n)` gives the asserted `a`, `b`, and `c`
  reductions.  The resulting exact candidate sets contain 123 all-odd triples
  and 209 odd-3-unit triples.
- The six capacity-unbounded all-odd pairs and three capacity-unbounded
  odd-3-unit pairs are exactly those stated.  The independently recomputed
  essential carrier `(w,g)` table agrees entry-for-entry.
- For a nonempty fixed carrier, `ell<g` makes the carrier/tooth incidence graph
  a star forest and gives component width at most `max(ell,w+2ell)`.  For an
  empty carrier the sharper and necessary clause is simply width at most
  `ell`.
- Exact strict-wall topology gives the unique maxima `17/693` at `(1,9,11)`
  and `19/1001` at `(1,11,13)`, with runners-up `1/42` and `29/1547`.
- The endpoint control `(1,7,13)` has four components of width `1/91` and
  eight of width `1/98`; `x=1/14` is absent.  Thus the apparent a.e. width
  `1/49` is not a strict component.
- Direct quotient-coordinate recomputation on every one of the 332 residual
  rows gives exactly twice the physical maximum.  Abstractly, the failure set
  is antipodally invariant, is a proper open subset of measure at most `3/7`,
  and its components therefore occur in antipodal pairs on which `x -> 2x`
  is injective.  This proves the factor two generally, not only in the boxes.

## Required text repairs

1. The all-odd bounded pair `(7,11)` has formal cutoff `165/14`, which is
   greater than `11`; it is omitted correctly because the next admissible odd
   tail is `13`.  Replace “all other cutoffs are at most b” by “all other
   pairs admit no further admissible c.”
2. Formula (9) is a valid uniform bound for **nonempty** `P*_{ab}`.  For the
   empty carriers `(1,3)` and `(1,5)`, the term `w+2ell=2ell` describes a
   nonexistent star.  State separately that `F subset C_c`, hence the width
   is at most `ell`.  With this repair, the advertised safe thresholds `7`
   (all odd) and `11` (odd 3-unit) are strict exactly as required.

Neither repair changes a finite box, a maximum, an equality case, or the final
theorem.

## Reproduction

```powershell
python -B 04-computation/lrc14_dyadic_strict_component_width_thm4451_independent.py
python -O -B 04-computation/lrc14_dyadic_strict_component_width_thm4451_independent.py
```

The two runs are line-identical and end in `PASS`.
