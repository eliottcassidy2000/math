# LRC(14) `z1=216`: third ruler prefix and laminar status calculus

**Status: INDEPENDENT VERIFIED-EXACT CROSS-AUDIT OF THM-3308/THM-3313;
PROJECTED-ATLAS SCOPE ONLY.**  Concurrent exact companions now canonize the
threshold compiler in
[THM-3308](../01-canon/theorems/THM-3308-threshold-chain-modular-multicovers-and-three-layer-status-circuit.md)
and the four-row closure in
[THM-3313](../01-canon/theorems/THM-3313-projected-k3-z216-third-ruler-cost-prefix-multicover-closure.md).
This audit independently freezes the same counts and laminar inequalities from
the inherited ray/status engine.  It is not a separate physical-cover theorem
and does not prove LRC(14).

**Subsequent continuation.**
[THM-3320](../01-canon/theorems/THM-3320-projected-k3-z216-fourth-ruler-prefix-and-affine-multicover-closure.md)
closes the next singleton and three-row family identified here, moving the
historical endpoint `373157 / 353 / 31` to ledger `373153`, `349` wall rows,
and `29` families.  This note's exact `373161->373157`, `357->353`, and
`33->31` audit remains the preceding step, not the current frontier.

The exact companion is
[`lrc14_j7_k3_z216_third_complete_ruler_cost_prefix_laminar_status_audit_20260803.py`](../04-computation/lrc14_j7_k3_z216_third_complete_ruler_cost_prefix_laminar_status_audit_20260803.py),
with [frozen output](../05-knowledge/results/lrc14_j7_k3_z216_third_complete_ruler_cost_prefix_laminar_status_audit_20260803.out).

## Inheritance and typed operation

The closest mechanism is the exact ray/status upper relaxation used by
THM-3139, THM-3270, and THM-3281.  The preceding
[cost-prefix reflection](lrc14-z216-intrinsic-ruler-cost-prefix-and-status-cut-taxonomy-codex-20260803.md)
left

```text
projected ledger 373161,
357 live z1=216 wall rows in 33 complete (gcd(216,L),L) families,
projected cap 216.                                           (1)
```

Its first unaudited intrinsic-cost prefix was the singleton
`gcd72/L72072`, followed by the three-row family `gcd24/L25872`.  The lawful
operation is to screen both **complete** families, retaining every row label
and inherited filter.  The scalar invoice `L*r` only orders work; it proves no
safety statement.

Canonical hostile controls are retained in the same run:

- index `14` still closes as `59=24 crude+35 status+0 residual`;
- index `64` still leaves eight residual masks, all of quotient type `8`.

Thus an empty selected screen is not an engine convention that kills every
input.

## Exact four-row descent

The complete result is

| family | indices | states | crude | exact status | residual |
|---|---|---:|---:|---:|---:|
| `gcd72/L72072` | `157` | 111 | 93 | 18 | 0 |
| `gcd24/L25872` | `53,293,357` | 312 | 196 | 116 | 0 |
| **total** | four rows | **423** | **289** | **134** | **0** |

The one-way consequence is therefore

```text
ledger       373161 -> 373157,
wall rows        357 -> 353,
families          33 -> 31,
cap                         216 (unchanged).               (2)
```

The `134` exact status states split as

```text
128 coordinate-union cuts,
  2 zero-reduced coordinate unions,
  3 two-fan cuts,
  1 laminar modular cut,
  0 unclassified.                                          (3)
```

Every inherited Farkas certificate is independently rebuilt over exact
rationals, but `(3)` records only deterministic solver-free mechanisms.

## Laminar modular tail lemma

Let a hypothetical common status table put mass on the sixteen patterns
`P subseteq {0,1,2,3}`.  Write `q` for total mass, `m_i` for its four
marginals, `c(P)` for capacity, `H_t` for the screen's required tail demand,
and

```text
E_t={P:c(P)>=t},       T_t=mass(E_t),       T_t>=H_t.       (4)
```

The `E_t` are nested.  Suppose a selected threshold set `A` and nonnegative
integer coefficients satisfy, pointwise on all sixteen patterns,

```text
sum_(t in A) 1_Et(P)
  <= c_0 + sum_(i=0)^3 c_(i+1) 1_(i in P).                 (5)
```

Summing `(5)` against the table and using `(4)` gives

```text
sum_(t in A) H_t <= sum_(t in A) T_t
  <= c_0 q + sum_(i=0)^3 c_(i+1) m_i.                      (6)
```

This is a modular upper bound on a sum of nested threshold events.  Since the
left side of `(5)` is at most `|A|`, every coefficient may be truncated to
`0<=c_j<=|A|` without weakening a pointwise majorant.  The companion exhausts
one-, then two-, then three-threshold sets `A` and that bounded nonnegative
integer language.  In all twelve displayed certificates `c_0=0`.  Required
tail demand violates `(6)` for the one new laminar state and for all eleven
states left as “weighted cores” by the preceding reflection.

Ten old cores need two thresholds and one needs three.  Together with the new
two-threshold state, the six proof shapes are

```text
(1,2;   0,2,2,1,1) x5
(1,4,5; 0,1,3,1,1) x1
(2,5;   0,2,1,1,0) x1
(3,4;   0,2,1,1,2) x1
(3,4;   0,2,1,2,1) x2
(3,5;   0,2,1,1,0) x2.                                   (7)
```

Here the entries before the semicolon are thresholds and the five entries
after it are `(c_0,...,c_4)`.  The lone three-threshold certificate is minimal
only inside this equal-weight nonnegative-integer language; no universal proof
complexity claim is made.

## Connection and loss ledger

The useful transfer is a finite Boolean-capacity compiler, not an object
identification:

```text
source:    labelled projected wall row and exact ray states;
map:       status capacities -> nested threshold events -> modular bound (6);
preserved: infeasibility of the maintained necessary upper relaxation;
lost:      physical speed entry, endpoint origin, owner, phase, current,
           chronology, and all cover geometry;
sidecar needed for LRC: a lawful physical lift retaining those coordinates.
```

The result closes the selected rows only in the direction

```text
empty necessary upper screen => projected wall row empty.   (8)
```

No converse from screen feasibility is used.

## Reproduction, boundary, and next moves

Run

```text
python 04-computation/lrc14_j7_k3_z216_third_complete_ruler_cost_prefix_laminar_status_audit_20260803.py --processes 3
python -O 04-computation/lrc14_j7_k3_z216_third_complete_ruler_cost_prefix_laminar_status_audit_20260803.py --processes 2
```

Normal and optimized frozen runs are byte-identical.  LF-normalized hashes are

```text
script  cd7f7e529ebf009af055cd92f4449f8ba1eaa768a5493157505cbb9db462dc9f
output  5861334b1f9ebfd2b2c612b1646b45f13a97dc71eeaca68301cc5070732be25c
semantic c97d088f75478c11c0b82d36c6f0dffe0bf99ae50d2bad9a803630609440424e.
```

Exactly `353` wall rows in `31` families remain.  The next cost head is the
singleton `gcd24/L76440` at index `141`, followed by the three-row
`gcd24/L30576` family.  Ranked next operations are:

1. use this artifact as a concurrent replay of THM-3308/3313 and compare any
   future compiler changes against both routes;
2. extend the four-bit monotone/laminar event catalogue and identify
   which inequalities persist beyond the equal-weight language;
3. screen the next complete prefix, with the eight-mask hostile retained; and
4. keep this projected compiler separate from the two-axis physical lift,
   current, arbitrary `k<=1`, rung, and LRC exit.
