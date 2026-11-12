# Visible symmetry and protected 3-adic mass

[Proof](PROOF.md) · [Exact verifier](verify.py) · [Receipt](REPLAY.json)

The all-order q-color orbit inventory and endpoint floor are inherited from
Session284. New results describe their stochastic and arithmetic interaction:

- a certified limiting mean and variance, LLN and CLT for the3-primary
  color-stabilizer exponent of balanced ternary tournament towers;
- a Poisson symmetry-erasure window at k≈(log_2 n)/2 named random attachments;
- a quantified separation between this visible symmetry and Hamilton-path
  quotient divisibility, retaining the full-Aut tuple-orbit correction;
- the inherited floor's sharp necessary contact budget k>=ceil(n/6) for
  any3-adic unit path count within the extension family.

This is a complete large-scale law for a specified family/model. It is not
an exact formula for all path-count valuations or a Chan–Pak proof.

Run from repository root:

```bash
python3 05-knowledge/results/session287_padic_cores_20261112/symmetry/verify.py
python3 -O 05-knowledge/results/session287_padic_cores_20261112/symmetry/verify.py
```

The two saved receipts agree byte-for-byte. Their exact finite controls
exhaust all512 one-vertex extensions of the nine-tower and its explicit81
parent automorphisms. The limiting constants are enclosed by rational tail
bounds; the listed larger Poisson parameters are formula evaluations, not
an enumeration of their enormous tournament layers. The all-order results
are proved in the note, not inferred from these controls.

Inheritance search: targeted searches for p-adic symmetry, tower valuations,
random stabilizers and the exact Session284 orbit recurrence found
THM-4755's finite-exact tower valuations and Session284's full orbit/floor
laws. These are explicitly retained. No external novelty claim is made for
the general probabilistic methods or for related colored-tree limit laws.
