---
id: HYP-2061
status: PROGRESS; THM-396 proved, full blocker-cover contradiction open
source: codex-2026-06-02-S558
related:
  - HYP-2059
  - HYP-2060
  - HYP-2058
  - THM-396
  - THM-397
---

# HYP-2061: HYP-2060 reduces the small-pinch route to sum-multiple shields or blocker covers

## Claim

The direct HYP-2059 target

```text
every 13-set has a reduced-sum <= 14 pair whose pinch clears all runners
```

is false.  The correct proof route from HYP-2060 is sharper:

```text
small-pinch failure = sum-multiple shields OR a finite non-shield blocker cover.
```

THM-396 proves the first branch exactly.  If a single speed blocks every
n=14-safe pinch residue of a pair `(a,b)` with reduced sum
`(a+b)/gcd(a,b) <= 14`, then that speed is divisible by `a+b`.

THM-397 adds a necessary condition for the collective branch.  If several
non-shield blockers cover every safe pinch residue of a small pair with
denominator `D=a+b`, then at least one blocker must lie in the endpoint window
modulo `D`:

```text
||c/D|| < 1/14.
```

Consequently, when the actual pair sum satisfies `D <= 14`, collective
non-shield covers are impossible.  Any failed small pair with `a+b <= 14` must
have a genuine sum-multiple shield.

## Evidence and obstruction

The script `lrc_n14_small_pinch_shield_s558.py` gives three facts.

1. **N3 is false.**  A sieve-covered, primitive, checkable-HYP-2060-style set
   can have no reduced-sum-`<=14` clearing pinch:

   ```text
   (1, 2, 9, 26, 110, 153, 166, 170, 178, 190, 192, 196, 201).
   ```

   It is still lonely; its exact optimum is not a counterexample.

2. **The HYP-2060-like N3 failures seen in bounded random probes are shielded.**
   In the example above, each small pair is killed by sum-multiple shields:
   `(1,2)` by multiples of `3`, `(1,9)` by multiples of `10`, `(2,9)` by a
   multiple of `11`, and `(2,26)` by a multiple of `28`.

3. **The universal-shield theorem is not enough.**  The pair `(3,12)` with
   denominator `15` has safe residues covered by the five non-shield residues
   `{14,8,10,11,13}`.  Thus multiple runners can collectively kill a small
   pinch without any one of them being a sum-multiple shield.

## Tournament Analysis declaration

Vertices are small pair-cells `(a,b)` with

```text
(a+b)/gcd(a,b) <= 14.
```

Pairwise observable:

```text
debt(a,b) = (
  best n=14 collar over m/(a+b),
  number of universal blockers,
  number of sum-multiple shields,
  reduced sum,
  pair sum
)
```

Switch/gauge:

```text
(a,b) -> (c,d)
```

when `(a,b)` has larger best collar; ties prefer fewer universal blockers, then
fewer sum-multiple shields, then smaller reduced sum.  The tie Hamiltonian path
is lexicographic by `(reduced sum, pair sum, a, b)`.

Fingerprints to report: score histogram, directed 3-cycles, SCCs,
Hamiltonian-path count, and edge flips against the plain reduced-sum order.

## Proof route now

To turn HYP-2060 into an n=14 proof, prove one of these sharper statements:

1. A HYP-2060 counterexample cannot have all short-resonance pair-cells
   sum-multiple shielded without violating primitivity, distinctness, or the
   mod-7 singleton coupling.
2. A collective non-shield cover of every short-resonance pair-cell forces
   endpoint blockers modulo the corresponding pair sums; the resulting
   endpoint ledger is incompatible with the endpoint-core/pressure-owner
   constraints, or creates a large-pair pinch with radius `>=1/14`.
3. Under A1-A4 plus B1, the shield graph must contain an additive cycle
   `c_i | a_i+b_i` whose gcd descends, contradicting F1 minimality.

The session did not prove those residual statements.  It did prove the
universal-blocker branch and refuted the overly broad N3 formulation.

## Files

- `01-canon/theorems/THM-396-lrc-n14-small-pinch-universal-blocker.md`
- `01-canon/theorems/THM-397-lrc-n14-collective-pinch-endpoint-blocker.md`
- `04-computation/lrc_n14_small_pinch_shield_s558.py`
- `05-knowledge/results/lrc_n14_small_pinch_shield_s558.out`
- `07-reflections/lrc-n14-small-pinch-shield-residual-s558.md`
