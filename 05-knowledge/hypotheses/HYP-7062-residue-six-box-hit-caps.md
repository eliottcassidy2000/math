# HYP-7062 — the residue-six box-hit caps

**Status:** OPEN / FINITE-EXACT EVIDENCE (codex-2026-07-16-S18).

For `P_V(S)=meas{x:{sec(vx):v in V}=S}`, distinct positive integer speeds, and

```text
S_a={1,3,5,6},  S_b={2,3,4,6},  S_c={1,2,4,5},
```

the proposed universal bounds are

```text
P_V(S_a)+P_V(S_b) <= 1/12                    (|V|=4),
P_V(S_c)          <= 5/42                    (|V|=4),
P_E({0} union S_a)+P_E({0} union S_b) <= 40/441   (|E|=5).
```

Together with `THM-905`, these imply `-F_6<=535/7203<0.097`, closing the sole
limiting sign in `THM-891`.  The first two scanned maxima over all `33,919` primitive
quadruples through largest speed `32` occur at compact tuples `(2,3,4,6)` and
`(1,2,3,4)` respectively.  The five-runner maximum over all `41,656` primitive
quintuples through largest speed `24` is `40/441` at `(3,5,7,9,11)`.

The box carrier preserves the exceptional residue-six miss patterns but destroys the
labels and the rest of the relation lattice.  Alternate vertices considered were
runners, gaps, fixed sections, boundaries, wall crossings, residues, Fourier modes,
primitive relations, and proof obligations.  Relation strata must accompany the box
vertices because `THM-898` and `THM-899` show that additive plateaus do not decay with
scale.

Next: extend the exact scan, decompose each box observable into pair/triple/quartic
Bernoulli relation channels, and prove the caps first on every short-relation stratum,
then on the absolutely summable tail.
