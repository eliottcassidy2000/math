# Prescribed-target inverse and the prime 6+4 channels

2026-11-13. **Chan–Pak remains open.** This session separates exact
target selection from the stronger-than-expected local feasibility
of the same selective geometry.

## What can be selected from the target

- [Dyadic inverse](DYADIC_INVERSE.md): in a physical seven-factor
  reservoir bank the absolute target fixes a unique dyadic slot
  scale for each source and cut. Factoring the remaining integer
  into the bank's coprime counts either builds the reservoir or
  rules that source out. Allowing transitive sources gives an
  exact compiler for both cuts, not a blind length search.
- [Fixed prime cut](FIXED_PRIME_CUT.md): keep the prime-derived
  source and its named cut. The absolute target fixes both the
  ternary valuation and the dyadic scale. Explicit exact
  recipes include a prime and a seven-times-prime count; this
  same fixed cut covers every odd residue modulo 2^r 3^s, even
  though its absolute spectrum is only polylogarithmic.
- [Residue aperture](RESIDUE_APERTURE.md): allowing both modules
  to vary gives simultaneous coverage of **every residue
  compatible with oddness, modulo every integer**, in each
  strong quotient. The same explicit pattern bank nevertheless
  gives only exp(O(sqrt(log X))) absolute counts up to X.
  Thus no finite congruence sieve can obstruct these charts,
  but a locally admissible target still needs an exact lift.

## What is ruled out

- [One-apex capacity](ONE_APEX_CAPACITY.md): arbitrary binary
  cuts in a single one-apex module on each side give at most
  O(X^(log_3 2)) targets. Even arbitrary ordinal compositions
  of such cores give O(X^(199/200)). The low two-slot bits and
  low Hamilton bits read different parts of the same cut word.
- [Core growth](CORE_GROWTH.md): even arbitrary SCC types, if
  bounded by sufficiently small multiples of
  log log X / log log log X, yield only X^(o(1)) outputs.
  These are obstructions to specific grammars, not to arbitrary
  growing modules or to Chan–Pak itself.

## The 6+4 split has a different native boundary

[Prime selection](PRIME_SELECTION_4_PLUS_6.md) finds a six-vertex
path through 2 and two residual edges on the first ten primes.
The natural vertex split is not the ordered cutoff into the first
six and next four, nor the Fi22/Baby/Monster membership split.
Connectivity of the marked first ten under optional primes
through 71 is exactly 47 AND (53 OR (41 AND 59)). Even every
connecting choice still has three forced leaves and therefore
no spanning path. The size-25 ordinary square-sum graph instead
has a genuine 4+6 split of **complete paths** across two ears.
A scalar split, a connected selection, and a legal final
path are three different stages of the inverse.

The remaining proof target is an absolute, target-decreasing
selector with more physical capacity than the one-apex alphabet.
The modular solutions are useful for choosing phases, but their
padding cost cannot be treated as a free exact-target adjustment.

Replay from this directory: replay_dyadic_inverse.py,
replay_fixed_prime_cut.py, replay_residue.py, replay_one_apex.py,
and replay_prime_selection.py (all Python 3, no third-party modules).
