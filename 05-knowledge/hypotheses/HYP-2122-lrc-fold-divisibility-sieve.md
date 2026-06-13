---
id: HYP-2122
status: OPEN proof route; S587 proves the denominator-shield algebra and audits n=14 calibrations
source: codex-2026-06-03-S587
related: [HYP-2121, HYP-2120, THM-400, HYP-2119, HYP-2118, HYP-2117, HYP-2116, HYP-2115, HYP-2114, HYP-2113, HYP-2112, HYP-2110, HYP-2108, HYP-2107, HYP-2101, HYP-2100, HYP-2096]
---

# HYP-2122: Lemma B is a fold-divisibility sieve over pair denominators

## Claim

The structural unit behind Lemma B is not only the visible triple
`a+b=c`.  It is the denominator gate

```text
D = a+b,        t = m/D.
```

At every such pair-pinch clock,

```text
a*t + b*t = m.
```

Therefore any speed `v` with `D | v` satisfies `v*t in Z` and shields that
entire `D`-pinch family.  The usual visible fold `c=a+b` is only the first
case `v=D`.  A multiple such as `v=24` can shield the `D=12` pinches even
when the direct speed `12` is absent.

In THM-400 language, the shield is observer-coupled: if `v=qD`, then
`q*a + q*b = v` has nonzero augmentation `2q-1`.  Balanced 4-term energy is
still background; denominator divisibility is the part that references the
observer.

After the incoming HYP-2120 and HYP-2121 perspective results, the same point
can be said in tournament language: denominator shields are rooted/source
gates with incident threshold payload.  They distinguish the observer-coupled
slice from the source-less regular/worry fiber; balanced energy is the
unrooted projection that forgets which clocks are actually killed.

## Evidence From S587

S587 computes exact maximin values, unbalanced/balanced relation counts, a
low-denominator shield ledger, visible-fold deletion shadows, deterministic
sample controls, and Tournament Analysis.

At total `n=14`:

- AP has `M=1/14`, kills every appearing pair denominator `D<n`, and leaves
  `D=n` unshielded.
- `V*=(1,2,3,4,5,6,7,8,9,10,11,13,24)` also has `M=1/14`, kills every
  appearing `D<n`, and leaves `D=n` unshielded.  The missing direct fold
  `D=12` is explained by divisibility: `24` shields all `m/12` pinches.
- Unit-shift AP has the same raw 4-term energy and many folds, but the speed
  `14` shields `D=n`; the delta-clock survivor is killed and the row is loose
  with `M-delta=+0.05357`.
- Far-shift AP has no low pair-denominator scaffold and is very loose:
  `M-delta=+0.28571`.
- The doubled-apex stress row kills all low denominators and is close but
  positive: `M-delta=+0.00265`, matching the expectation that endpoint/Phi
  residuals begin once the clean clock scaffold is perturbed.

Deletion shadows show why the 2025 prime/divisibility machinery belongs after
the fold.  Removing a visible fold and testing the shadow witness is blocked
in all AP and `V*` fold shadows (`36/36` and `31/31`), but the residual
denominators include nontrivial prime/composite values such as `17,19,23,25,
31,35`.  Thus the fold gate creates a blocked-pinch residue ledger, and the
prime/Phi/endpoint machinery consumes that ledger.

The alternating Lemma A/B sample audit keeps the control honest.  Rows with no
unbalanced relations (`u=0`) have large positive margins in the tested
regime, while hardness correlates positively with unbalanced count and weakly
or negatively with balanced energy.  At `n=14`, the sampled fold-rich rows
still have positive margin unless a terminal clock or residual gate is present.
So S587 supports the route, not a standalone proof.

## Proof Route

1. Prove the denominator-shield lemma globally: for every pair `(a,b)` and
   every clock `t=m/(a+b)`, every speed divisible by `a+b` lies at the
   observer and blocks that clock.
2. Upgrade visible Lemma B to a low-denominator sieve.  A row is structurally
   fold-controlled when every appearing `D<n` pair-pinch is shielded by some
   speed divisible by `D`.
3. Split the remaining case by `D=n`.
   If `D=n` is unshielded, the ordinary `j/n` delta-clock is the survivor
   candidate, matching AP and `V*`.
   If `D=n` is shielded by a multiple of `n`, the direct clock is gone and the
   row must route to Cprime/Phi/endpoint positivity.
4. Route rows with no coherent denominator-shield scaffold to Lemma A:
   circuit-free or balanced-only structure should give a positive discrepancy
   margin because it lacks observer-coupled denominator gates.
5. For `n=14`, use AP and `V*` as the two calibration floor rows, unit-shift
   and far-shift AP as controls, and doubled-apex rows as the stress interface
   to HYP-2112/HYP-2108/HYP-2107.

The missing theorem is a compression statement: killed low denominators plus
an unshielded `D=n` clock must deliver `1/n`, while killed `D=n` or endpoint
residuals must deliver positive `Phi`.

## Tournament Analysis

S587 uses proof lenses and denominator gates as tournament vertices, not
runners.

The pair observable is:

```text
(Lemma-B delivery, Lemma-A control, divisibility power,
 observer coupling, maturity)
```

The transitive ranking is:

```text
D_denominator_divisibility_shield
> visible_fold_a_plus_b_eq_c
> Phi_endpoint_prime_residual
> D_eq_n_unshielded_delta_clock
> augmentation_nonzero_count
> circuit_free_A_margin
> balanced_energy_background.
```

Score histogram is `{0:1,1:1,2:1,3:1,4:1,5:1,6:1}`, directed 3-cycles `0`,
and Hamiltonian path count `1`.

This tournament encodes the hidden no-shortcut fact: balanced energy cannot
jump over observer-coupled denominator gates.  It must first become a labelled
virtual fold, endpoint component, residue state, or proof obligation.

## Assumption Challenge

Candidate vertices considered in this session included runners, pair
denominators `D=a+b`, shield speeds `v` with `D|v`, visible fold clauses,
deleted-fold shadows, residues, endpoint components, `Phi` terms, balanced
relations, and proof obligations.

The chosen quotient preserves the LRC predicate "which low pair-pinches are
killed before they can become lonely clocks, and whether `D=n` survives as the
delta-clock."

It destroys full endpoint-owner geometry and full prime-fibre languages.
Those are deliberately delegated to HYP-2112, HYP-2108, and HYP-2107.

The challenged assumption is that Lemma B is only the direct relation
`a+b=c`.  The useful structural unit is `D|v` at the pair denominator
`D=a+b`.

## See

`04-computation/lrc_fold_divisibility_sieve_s587.py`,
`05-knowledge/results/lrc_fold_divisibility_sieve_s587.out`,
`07-reflections/lrc-fold-divisibility-sieve-s587.md`,
HYP-2121, HYP-2120, THM-400, HYP-2119, HYP-2118, HYP-2117, HYP-2116,
HYP-2115, HYP-2114, HYP-2113, HYP-2112, HYP-2110, HYP-2108, HYP-2107,
HYP-2101, HYP-2100, HYP-2096.
