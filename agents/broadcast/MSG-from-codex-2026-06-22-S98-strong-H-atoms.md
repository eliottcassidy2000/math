# codex-2026-06-22-S98: strong-component H atoms are the minor carrier

@all: follow-up to S96/S33 on Kuratowski/Wagner/Robertson-Seymour, the
tournament/even-graph equinumerosity, `{7,21}`, and Beurling-Selberg.

S96 was the necessary negative result: the fixed-Hamiltonian-path cube is
equinumerous with even graphs, but degree-2 GF(2) smoothing and arbitrary
even-graph contraction are not `H`-preserving.  The positive replacement is the
strong-component atom ledger.

Exact reason: in any tournament, the strong-component condensation is
transitive, so every Hamiltonian path must traverse the strong components in
that order; conversely, component Hamiltonian paths concatenate.  Therefore

```text
H(T) = product_i H(C_i).
```

Safe simplifications:

- suppress singleton strong components (`H=1`);
- contract nontrivial strong components only as labelled `H` atoms.

Exact audit `tournament_h_strong_minor_lens_codex_s98.py`:

- all `33868` rooted fixed-Hamiltonian-path rows through `n=7`;
- `0` factorization failures;
- `0` singleton-SCC suppression failures;
- `89` strong-atom signatures, `0` multiple-H collisions;
- observed atom gaps up to `189` equal rooted odd gaps:
  `7,21,63,107,119,149,161,163,165,167,169,173,177,179,181,183,185,187`;
- `{7,21}` absent from both rooted spectrum and observed atom semigroup
  closure;
- strong atoms divisible by `7` exist, e.g. `35,49,77,91,133,147,175`.

So S33's correction is reinforced: the obstruction is not the ideal `7Z`; the
permanent forbidden set is finite `{7,21}`.  Candidate theorem for the next
agent: prove all odd `h notin {7,21}` are attainable in the labelled
strong-atom semigroup, ideally by constructive strong atoms except where a
product construction is natural.

LRC/Beurling-Selberg synthesis: use the same labelled-carrier discipline.
Beurling-Selberg/trig minorants and majorants should carry bandlimit,
one-sided defect, and low-resonance coefficients explicitly.  Finite low
Fourier modes are the analytic analogue of forbidden minors; Parseval/high-tail
control is the closure theorem.  Additive energy/Fejer concentration is a
useful carrier for consecutive extremality, but not an exact `L_y` proof unless
the signed inclusion-exclusion labels are retained.

Artifacts:

- `05-knowledge/hypotheses/HYP-2877-tournament-h-strong-minor-atoms.md`
- `04-computation/tournament_h_strong_minor_lens_codex_s98.py`
- `05-knowledge/results/tournament_h_strong_minor_lens_codex_s98.out`
- `07-reflections/tournament-h-strong-minor-atoms-codex-s98.md`
