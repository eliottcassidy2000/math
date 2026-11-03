# Independent cross-audit: toothpick geometry and the self-converse order image

**AUDIT PASS**, 2026-11-02. This review is by the parity-lane agent, who did
not author `toothpick.md`, `toothpick_probe.py`, `spectrum.md`, or
`spectrum_probe.py`. All claims remain in their declared scopes. The
unrestricted all-order Hamilton spectrum is still OPEN.

## Toothpick physical rule and reflection proof

The length-two integer-grid representation is faithful: any crossing of
horizontal and vertical sticks occurs at an integer endpoint or midpoint.
Thus retaining endpoint multiplicity and midpoint occupancy suffices for
the stated exposure predicate. The code counts endpoint incidences rather
than distinct neighboring source orbits; it excludes all occupied centers.
An exposed endpoint has exactly one parent, so newborn orientation is
well-defined, simultaneous generation has no same-center ambiguity, and
the parent transport checks have the stated meaning.

The reflection induction is valid. If an off-axis stick has an endpoint on
an axis, its distinct reflection in that axis supplies another endpoint
incidence at the same physical point. That point is blocked. The three
initial on-axis centers exhaust the only remaining source of on-axis
endpoints: the seed's endpoints already carry the two stage-two centers,
while the endpoints of those horizontal sticks are off both axes. Hence
there are no later on-axis centers. The stabilizer conclusions follow
without assuming generic orbit size: the seed has orbit one; the two
stage-two sticks have orbit two; every other stick has orbit four. The
claimed fixed-set sizes and the congruence `T(n)=3 mod4` for `n>=2` follow.

At stage three, four endpoint incidences over the two ports `(1,0),(-1,0)`
really give two sources per port. The orbit quotient's one neighboring
source orbit is not exposure. The exact weighted incidence formula follows
from constant degree on each port orbit; no division by a variable
stabilizer is hidden in it. Rooted tournament insertion similarly retains
the needed vertex-orbit and groupoid weights: deleting the distinguished
root gives exactly the old automorphism action, and summing reciprocal
root stabilizers gives `(n+1)/|Aut(U)|`. The `3z+z^3` versus `z+z^3`
hostile is correct.

The planar population obstruction has the correct target: one separate
active lattice object per ordinary order-n class at the same clock n. It
does not exclude encoded weights, histories, or changes of clock. A
reflection-equivariant obstruction, if added from incoming work, has its
own different quantifiers and is not established merely by this growth
comparison.

## Spectrum formulas and all-order scope

Ordered SCC decomposition gives a unique word of strong classes. Each
Hamiltonian path visits those SCCs in that order, so H multiplies and the
star-product identity is correct. Its unit is `z`, not polynomial `1`.

An ordinary self-converse SCC word is uniquely determined by its first
half and its central SCC (or empty center). Converse reverses the word
and reverses each letter. Thus the map `U -> C -> U^op` has no hidden
half-swap quotient, and its costs really are `2|U|+|C|` and `H(U)^2 H(C)`.
The square-central generating formula is therefore both necessary and
sufficient with its stated realization guard.

For each parity epsilon, every witness has order at least
`min_(u^2|h)(2 mu(u)+theta_epsilon(h/u^2))`; a minimizing actual U and C
construct a witness attaining every finite minimum. Adding a source and
sink in converse positions raises order by two and leaves H fixed.
This proves exactly one tail or the empty set in each parity. The formula
is an exact relative image theorem, not an explicit classification of
which strong-center inventories are nonempty.

The primary [Busch paper, Corollary 1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/)
was re-opened for this audit. Its piecewise exponential minimum gives
`f(7)=25` and increases thereafter, so all strong tournaments with
`H<=21` have order at most six. This justifies the paper's finite-to-global
cutoff for the complete small-H table. In particular H19 can occur in
ordinary tournaments but cannot occur in a self-converse tournament:
19 has no nontrivial square divisor, and both strong H19 classes in the
complete inventory form a converse pair. This all-order exclusion is not
an inference from lack of witnesses at an arbitrary finite cutoff.

## Independent checks actually performed

Both producer scripts were rerun in isolated Python mode, normally and
with `-O`; their respective JSON outputs agreed byte for byte. The parity
lane's full-permutation orbit census and the spectrum lane's degree-cell
augmentation census produced identical `F_n` and self-converse `F_n`
for every `n=0,...,6`. These are different generation/canonicalization
paths, not two calls to the same implementation.

An additional toothpick replay used **closed endpoint segments** rather
than center/orientation triples or endpoint counters. For each candidate
endpoint p, it scanned all existing closed segments and grew a perpendicular
stick exactly when only one segment contained p. Through stage 32 it
independently verified every birth count against the dyadic recursion,
every reflection, and the exact three on-axis centers. Selected totals:

| Stage | 3 | 4 | 8 | 16 | 32 |
|---:|---:|---:|---:|---:|---:|
| Segment replay | 7 | 11 | 43 | 171 | 683 |

The audit's independent growth predicate was:

```python
def contains(segment, p):
    a, b = segment
    return (min(a[0], b[0]) <= p[0] <= max(a[0], b[0])
            and min(a[1], b[1]) <= p[1] <= max(a[1], b[1]))

for p in {p for segment in sticks for p in segment}:
    sources = [segment for segment in sticks if contains(segment, p)]
    if len(sources) == 1:
        # Place the unique perpendicular length-two segment centered at p.
        ...
```

All segments are axis-aligned and their endpoint pairs are sorted, so this
test detects both midpoint occupancy and multiple endpoint incidence by
geometric containment. No producer growth function was imported.

Targeted inheritance searches for `alternating trail`, `odd fibre/fiber
involution`, `least presentation ... pair`, and `involution principle` did
not recover this packet's particular ordered-fibre class pairing. They did
recover the already credited THM-643 odd-fibre mechanism and prior warnings
that complement lines are a correspondence rather than a class function.
This search is not a proof of novelty. The elementary cancellation lemma
and its explicit implementation carry no historical-priority claim.
