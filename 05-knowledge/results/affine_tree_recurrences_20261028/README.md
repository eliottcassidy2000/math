# Affine stopping-tree recurrences and a source-refinement budget

**2026-10-28. Status:** PROVED statements in the three proof notes;
FINITE-EXACT coefficient/guard/hostile checks; VERIFIED decimal limit
probes. Universal Collatz and native all-order wiring remain OPEN.
This packet makes no external novelty or priority claim.

The useful addition is a recurrence which keeps the affine carry, followed
by a cheaper stopping grammar whose carry tail can actually be paid.

- [Ordered recurrences](RECURRENCES.md): all integer-carry moments follow
  an explicit triangular recurrence with the argument shifts `u -> 3^i u`
  and `z -> 2^j z`. There is a port-indexed matrix version. The complete
  carry distribution retains the chronological word and source guard;
  finitely many moments need not.
- [Moment laws](MOMENTS.md): critical ternary stopping has mean slope one
  but infinite mean carry and infinite slope second moment. A family of
  `d`-child rules branching at valuations divisible by `r` gives exact
  clock, slope and carry recurrences. Its critical arities are
  `3,7,15,31,...` for an explicit probability reason.
- [Odd-arity extension](ODD_ARITIES.md): the stronger incoming clock
  compiler extends the carry law to every odd arity. The slope tilt reads
  the binary selector in base four. A five-child rule has finite slope
  variance and infinite mean carry, separating these two obstructions.
- [Refinement budget](REFINEMENT.md): the `(d,r)=(3,3)` rule has
  `EJ=7/4`, `ED=7/2`, `EB=7/12`, and an exponential carry-tail bound.
  For any odd source cylinder `I_K` of mass `2^{-K}`,

  \[
  E[B;m\ge3K+h\mid I_K]
  \le(245/193)(27/32)^K(3/4)^h.
  \]

  This bound preserves the conditioned source while paying for refinement
  with an increasing size budget. It does not eliminate a measure-zero
  exceptional integer or certify descent on a first stop.
  Three is the smallest integer branch-per-bit budget in this particular
  exponential-moment conditioning argument; actual conditional laws may
  permit a better bound.
- [Audit and recovery](AUDIT.md): independent ordered-shape sums, direct
  native integer checks, exact historical repairs, and the limits of the
  analogy with older continuant and moment work.

## Inheritance and concept board

**Anchor:** carry-sensitive Collatz recurrences and paid refinement.
**Niche:** August continuants and reciprocal coefficient jets, including
their hidden-coordinate return. **Wildcard:** the September Catalan/graph
carrier, audited rather than imported as a tournament identity.

Closest proved mechanisms are **THM-4604**'s affine carry cocycle
([exact file](../../../01-canon/theorems/THM-4604-collatz-words-are-the-reducible-locus-carries-are-the-extension-cocycle-and-traces-only-hear-the-tuning.md))
and the [Session94 valuation-tree law](../session94_recursive_ports_crt_clocks_20261029/dynamics/VALUATION_TREE_CLOCKS.md).
The canonical hostile is `41 ->107` on `(2,1,1,1)`: auxiliary extinction
does not mean descent. The corrected near miss is **THM-1770**
([file](../../../01-canon/theorems/THM-1770-the-localisation-lemma-as-first-return-renewal-and-the-pair-only-closure.md)):
primitive return atoms cannot be declared separately zero after their
sum cancels. The less-used sidecar recovered here is the full exit
coordinate of **THM-3186**, refined by **THM-3192**; see the audit links.

| Live concept | Retained information | Decisive test / boundary |
|---|---|---|
| Prefix-free first returns | chronological word and valuation clock | first stop need not descend |
| Affine carry jets | ordered extension coordinate | equal `(J,D)` with carries 103 and 151 |
| Critical versus subcritical trees | actual probability law | mean slope one with finite or infinite carry |
| Source cylinders | same marked residue under refinement | three extra branches per source bit pay the bound |
| Shared ports and reciprocal jets | joint endpoints, hidden coordinates | separate port tests do not imply a common transporter |
| Orbit comparisons | time and value orders | `(1,3,2)` refutes the old tournament reading |

## Incoming work and what was not redone

Startup incorporated `origin/main` through `d3741fe6ed02`, especially
`complement_port_crt_20261028/PORT_LAWS.md` and `CLOCKS.md`. Their two
nine-vertex subgroups have the same fixed-set data but different joint
transporters. Their shared-sign valuation covariance is a first-step
identity, not an all-generation independence law. The corrected distinction
between global converse and path-native complement is retained.

The earlier [weighted CRT packet](../crt_weighted_ports_20261028/README.md)
already proves the LCA energy obstruction for arbitrary source weights.
We do not reprove a weaker energy statement: this packet changes the
stopping grammar, keeps the affine extension coordinate, and states the
new source-measure budget with its exact pointwise limitation.

The mid-session pull through `d92d6c0d74` supplied
[Session95](../session95_niche_recurrences_20261029/README.md). Its odd-arity
compiler supersedes a singleton-only view of critical trees. Its sharp
finite-support moment reconstruction strengthens our known-support jet
statement, and its joint Schur boundaries explain why a scalar circulation
recurrence misses a real degree-12 attachment response. The new
ODD_ARITIES.md uses those stronger mechanisms explicitly. It adds the
base-four tilt and the finite-variance/infinite-carry separation; it does
not claim the incoming clock compiler as our result.

## Reproduction and next tests

From this directory run `python3 -B reproduce.py`. It compares normal and
optimized replays against the saved receipts. Proofs are in the notes;
computational scope is in [AUDIT.md](AUDIT.md).

The next productive target is a source-preserving return controller for
expanding first stops of the `(3,3)` grammar. Test it against `3 ->5`,
equal-clock/different-carry words, and a fixed shrinking dyadic cylinder
before inferring any uniform ROOT inequality. A finite-state port or
residue quotient must prove that its retained information suffices for
the endpoint inequality and for legal composition.
