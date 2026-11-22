# Chan–Pak: matching coefficients, signed folds, and branching selectors

**2026-11-22 · PROVED in the scopes below. Chan–Pak remains OPEN.**
The question is about the *total* number \(H\) of Hamiltonian paths
in an ordinary tournament, not a fixed-path count or a decoded number.
See [Chan–Pak, Conjecture 6.2](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf).

| Lane | Result | Proof |
| --- | --- | --- |
| Anchor: constructive interfaces | A single transitive twin block has an exact Stirling expansion. The largest exponential base is one plus a bipartite matching number; a permanent is a leading coefficient. Unit dyadic slope forces the single offset \(+1\). | [Euler gluing and twin blocks](EULER-GLUING.md) |
| Niche: reversal symmetry | The signed fold has explicit Walsh coefficients indexed by alternating linear forests. They identify the first mod-four response. The inherited result that every odd number, **including 7 and 21**, occurs as a fixed-path count is an essential guard. | [Signed fold](SIGNED-FOLD.md) |
| Anchor: feedback-one exact scope | Short-range congruences do not rescue the one-apex binary compiler: the image has upper relative density at most \(\log 2\) in every fixed progression. Its self-converse subfamily is locally uniform among compatible residues but has only order \(\sqrt X\) distinct total counts below \(X\). | [One-apex boundary](ONE-APEX.md) |
| Wildcard: formal languages | On a fixed core with transitive spine, each long repeat is either erasable or exponentially expensive. Consequently a polynomial-growth regular attachment language has a polylogarithmic count image; a complete regular compiler would need recurrent branching. | [Periodic erasure](ERASABLE-CLOCKS.md) |

A further result of the twin-block construction is **#P-completeness** under a
one-call polynomial-time metric reduction. It applies to exact Hamiltonian
path counting in tournaments, even strongly connected self-converse ones,
and also to palindromic counting with a supplied anti-involution. These
statements concern extracting a count from a structured query; they do not
make the decoded integer a tournament path count. No literature priority
is claimed.

## Inheritance and mechanism

- **Closest canon:** [THM-5271, metagraph/seam selector](../../../01-canon/theorems/THM-5271-the-count-is-a-potential-on-the-flip-metagraph-the-seam-ladder-implies-chan-pak-the-microcosm-3-1-6-1-21-the-square-sum-cycle-at-32-and-the-order-10-near-join-reach.md), [THM-4878, graded spine grammar](../../../01-canon/theorems/THM-4878-graded-spine-grammars-cut-state-automaton-hankel-and-seam-ranks-cutwidth-censuses.md),
  [core-changing descent](../core_descent_endpoints_20261112/README.md), and
  [periodic apex recurrences](../chanpak_cross_problem_theorems_20261114/APEX_RECURRENCES.md).
  The proofs build on these coordinates rather than re-running their scans.
- **Hostile examples:** the \(7,21\) holes pertain to total paths, yet both
  occur in the fixed channel; the all-forward seam makes its clock width
  invisible, yet one backward outer arc makes that width exponential;
  a coefficient can encode a permanent and still say nothing about hitting
  that permanent as an absolute target.
- **Corrected near miss:** positive-density input words, arbitrary residue
  mixing, and long exact arithmetic rails are not output surjectivity.
  Conversely, a branched input language cannot hide all its growth in
  neutral padding; a sparse image would need collisions on active spines.
- **Least-used sidecars:** ordered cover profile, last-exit directed tree,
  owner signs and a terminal sign, and boundary activation. Each restores
  information lost by the previously scalar view.

**A common principle.** A path count is a sum over assignments subject to
endpoint compatibility. The assignment inside a transitive block is a
surjection; its boundary is a matching problem. Under reflection, the same
compatibility appears as a doubled sign test at a turn. On a long repeated
spine, it appears as a choice between freely erasable padding and an
exponentially large family of odd cycles. The relevant symmetry is the
symmetry of the **endpoint constraints**, not just the symmetry of the
unmarked tournament.

## Next proof obligation

Use the two-port/seam route to build a genuinely branching family in which
the two choices retain an actual endpoint interface and have disjoint
arithmetic images that can be inverted on target counts. For a symmetric
route, keep the ordered pair \((F,R)\), where \(F\) counts reversal-fixed
paths and \(R\) counts two-element path orbits, so that \(H=F+2R\).
Controlling \(F\) alone is now known not to control the holes. The task is not
to find another finite residue clock, but to make the reflected-pair cost
compatible with a target descent.

Exact finite checks of the new results (not needed for the
proofs): `python3 check_matching_shell.py` from this directory.
For the folded Walsh and total-profile formulas, run `python3 check_signed_fold.py`.
