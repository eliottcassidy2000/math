---
id: THM-905
title: THE RESIDUE-SIX BOX-CERTIFICATE REDUCTION — the remaining negative K6 kernel has the sharp pointwise majorant 6H+6J+H6, reducing the limiting sign to two four-runner and one five-runner sector-box inequalities
status: CLAIMED / CHECKPOINT STUB — pointwise reduction proved exactly; finite scans identify sharp candidate box constants; universal box bounds remain open
source: codex-2026-07-16-S18
depends_on: [THM-891, THM-903-reflection-frame-residue6]
related: [THM-898-fourrunner-boxhit-relation-stratified, THM-899-lattice-law-boxhit-constants, THM-904, HYP-7062]
verification: verifier and exact scan in progress
---

# THM-905 — the residue-six box-certificate reduction

This file claims the theorem number and records the proved algebra before the longer
box-hit referee.  Let `n_s` be the number of the five moving runners in sector `s`, and
put

```text
H  = n_1 n_3 n_5 n_6 + n_2 n_3 n_4 n_6,
J  = n_0 H,
H6 = n_1 n_2 n_4 n_5.
```

Extend `L_6(n)=-K_6(M(n))` by zero unless the inner missed set has size one or two.
Then every sector-count state satisfies the exact pointwise inequality

```text
L_6(n) <= 6 H + 6 J + H6.                                      (1)
```

The coefficients of `H` and `J` are sharp.  For an exceptional missed pair
`{1,5}` or `{2,4}`, the four complementary inner sectors are occupied.  If the fifth
mover duplicates one of them then `H=2,J=0`; if it occupies pinned sector zero then
`H=J=1`.  Both cases have `L_6=12` and force the two coefficients to be at least six.
The only positive singleton state not already caught by `H` misses sector six, where
`H6=1=L_6`.  All remaining kernel rows are nonpositive or zero.  This proves (1)
without computation; the verifier will also enumerate all 462 states.

For a set `S` of distinct sectors and a tuple `V` of the same size, write

```text
P_V(S) = meas{x in [0,1): {sec(vx): v in V} = S}.
```

Set

```text
S_a={1,3,5,6},   S_b={2,3,4,6},   S_c={1,2,4,5}.
```

Expanding the occupation monomials over mover subsets turns (1) into

```text
49(-F_6) <=
  6 sum_{|V|=4} (P_V(S_a)+P_V(S_b))
  + 6 (P_E({0} union S_a)+P_E({0} union S_b))
  + sum_{|V|=4} P_V(S_c).                                      (2)
```

The exact primitive scans currently support the three diameter-free bounds

```text
P_V(S_a)+P_V(S_b)                         <= 1/12,
P_V(S_c)                                  <= 5/42,
P_E({0} union S_a)+P_E({0} union S_b)     <= 40/441.            (3)
```

If (3) holds universally, (2) gives

```text
49(-F_6) <= 30/12 + 6(40/441) + 25/42 = 535/147,
-F_6 <= 535/7203 = 0.07427... < 0.097.                         (4)
```

This route has much more slack than `THM-904`'s triple bound.  The honest remaining
task is to prove the three box inequalities uniformly, relation stratum by relation
stratum.  `THM-899` explains why a naive scale-decaying remainder cannot do that: fixed
additive relations contribute nondecaying Bernoulli constants.

## Assumption challenge and tournament carrier

The quotient vertices here are target sector boxes, not runners or arcs.  It preserves
the exact exceptional-miss obligation in (1), but destroys labels within a hit box,
wall chronology, and the full relation lattice.  The planned tournament observable is
the exact excess of one primitive speed tuple over another for a fixed target box;
the gauge orients toward the larger hit mass, and ties use lexicographic order.  Its
Hamiltonian path ranks finite obstructions but cannot itself prove the universal caps.

## Honest boundary

- [x] prove the sharp pointwise box certificate;
- [x] derive the exact expectation reduction;
- [x] identify finite-exact candidate constants;
- [ ] bank the dependency-free verifier and extended scans;
- [ ] prove the three diameter-free box caps;
- [ ] compose the finite-`t` wall remainder.
