# HYP-1806: Packet compatibility transfer

**Status:** EXPLORATORY cross-problem hypothesis.

## Statement

The OCF lesson should transfer beyond tournaments:

```text
counting local packets is weaker than understanding packet compatibility.
```

For a problem with local witnesses or obstructions, build the packet conflict
graph first, then study its independence polynomial, deletion residues, and
incidence pivots.

## Evidence

In tournaments, `Omega(T)` turns odd cycles into packets and

```text
H(T) = I(Omega(T), 2).
```

Paley/Interval shows why packet compatibility matters: Paley can have more odd
cycles while the interval has more packable odd-cycle structure.

The same pattern appears in:

- disjoint-cycle and matching problems, where witnesses must be compatible;
- Lonely Runner endpoint covers, where protected endpoints form an incidence
  conflict system;
- Caccetta-Haggkvist/rainbow variants, where return paths or colors must be
  chosen without collision;
- active ranking, where ambiguity packets should guide the next comparison.

## Prediction

Feature blocks of the form

```text
packets -> conflict graph -> independence polynomial -> deletion residue
```

will separate examples that scalar counts merge.  In particular, they should
distinguish abundance-dominant examples from packability-dominant examples in
circulant `H` maximization and in active ranking benchmarks.

## See

- `07-reflections/kernel-residue-trick-atlas-2026-05-30.md`
- `07-reflections/applied-residue-phase-incidence-programs.md`
- HYP-1800 and HYP-1801

