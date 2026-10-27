# Why the advances work, II: potentials, symmetry budgets, and positive coordinates

*mac-mini-2026-10-27-potentials. A REFLECTION, the sequel to
[part I](why-the-advances-work-observation-composition-and-non-backtracking-macmini-20261026.md). Under the
repository's precedence rule it records why things worked, not what is true; every claim points to a canon or result
file and carries that file's status. Part I organised the repository around a quadruple (X, M, O, rho) and five
questions (reconstruction, descent, composition, cancellation, level). Tonight's dispatch (black-core distance, the
ROOT carry in X = 3s+1, H versus H/|Aut|, the Eisenstein triangle, the dyadic ledger, trace-norm budgets) reached at
least five sessions, and their results fit three further mechanisms. Together with part I they explain most of the
repository's novel advances.*

## 0. The sixth question: distance

Part I's questions are about one observation. Many recent theorems are about **how far** an object is from a target
locus under the monoid `M`. They all have the same two-sided form.

- **Lower bound = a potential.** A functional `Phi` on `X` that each generator of `M` moves by at most a constant `c`.
  Then `d(x, Target) >= (Phi(x) - max_Target Phi)/c`.
- **Upper bound = a descent certificate.** At every non-target state, an explicit move that lowers `Phi` by exactly `c`.
- **Exactness = no local minima.** The bound is attained from every source iff every non-target state admits a
  descent. Failures of exactness are local minima of `Phi`, and they are where the interesting structure lives.

| Instance | Potential | Lipschitz | Descent certificate | Exact? |
|---|---|---|---|---|
| distance to the non-strong locus (THM-4811 (B)) | arc connectivity `lambda` | 1 | Camion cycle entering a minimum cut | always |
| same, black moves only (THM-4811 (C)) | `lambda` | 1 | a black word into a minimum cut | all but `C3`, `R5`, `E(M; Z)` |
| distance to the transitive class (THM-4810 L, Session74) | `L = (1/2) sum |s_(i) - i|` | 1 | explicit 15-move word for `P_11` | `P_11` yes; `P_7` no (7 > 6) |
| triangle descent (black_descent_20261020c) | `c3` | `s(a) - s(b) - 1` | black closed path with `s(b) >= s(a)` | descent always (C3, R5 aside) |
| forest frames (Session74) | `beta >= L - n + 2` | 1 | six moves for `P_11` | `P_11` yes |
| exceptions to a bounded-arity bank (Session75, dyadic T1-T6) | trace norm | `2 sum |2 cos(k pi/(n+1))|` per move | none needed (static) | first order only |

**Why so many lower bounds look alike.** They are one identity. A native move along `Q = (a, ..., b)` is

    K  ->  -K + 2A_Q,        d  ->  -d + e_a - e_b,        f(X)  ->  f(V\X) + [a in X] - [b in X]

(THM-4811 (A)). Negation plus a path dipole. Every functional invariant under negation and Lipschitz in the dipole is
a potential: cut charges give `lambda`, the folded scores give `L`, `sum d^2` gives `c3`, spectral norms give the
trace-norm bounds. A session that finds "a cut-score bound" (the regular basin theorem), "score monotonicity"
(THM-4810), "a sharp local price per clock" (paid ROOT R2) or "a trace-norm obstruction" (Session75) is reading this
one law through a different functional. This is why the concurrent sessions converged: the regular black-core
distance `m` arrived from three sessions in one night (the aliquot_power_124 basin theorem, klein's THM-4814 (D), and
THM-4811 as the all-source case).

**Where exactness fails, a symmetry is forced.** The black failures are exactly `R5` (with `C3`) and the family `E` of
THM-4811 Theorem A: a bottom layer, a chain of middle vertices, and its converse. In each, every descent word is
palindromic under an anti-automorphism, so every move that could lower the potential is blue. The order-10 census found
one such class among 9,355,949. Generalising its visible structure gave the subfamily `E(M; Z)` within the hour, but that
generalisation was incomplete: the proof, and the order-11 census it suggested, found the middle-vertex members (every
odd order) and the free block patterns. So the exception exposed the mechanism (a forced anti-automorphism), but only the
proof delimited its scope. The pattern recurs across the repository: the exceptions to a uniform descent are small
symmetric families (the pure-blue alphabet, THM-4730; the regular `R5`; Paley sevens), and the first untested order
exposes the first member (part I, section 6: the first-untested-point law).

**The test this suggests.** For every distance statement in the repository, write the potential and the certificate.
If exactness is claimed, exhibit the descent at every state; if not, find the local minima and their forced symmetry.

## 1. Symmetry budgets: counting modulo automorphisms

`H(T)` never takes the values 7 or 21 (THM-1370). Its symmetry quotient `tc = H/|Aut|` had the same exclusion only
through order 11 (THM-4790); tonight's THM-4812 closes it at every order. The mechanism is general.

1. **Orbit product.** Counts that come from an independence-polynomial formula (the OCF, THM-002) are
   super-multiplicative over vertex partitions, and a symmetry group embeds in the product of its restrictions to its
   orbits. So `tc(T) >= prod tc(T[O])` over the orbits, each piece vertex-transitive.
2. **Vertex-transitive atoms.** Vertex-transitive pieces are few at small `tc`: only `K1, C3, R5, P7` have `tc <= 21`.
   Primes give circulants, small composites a finite group analysis, and large orders a greedy packing of disjoint
   triangles against Dixon's bound.
3. **A finite group search.** Small `tc` with large `|Aut|` forces the automorphism group to be nearly the full
   product of its orbit groups (index 1 or 3). Index 1 makes every orbit a module, and the rest is a search over 83
   cases.

The general lesson: **a symmetry quotient of a super-multiplicative count is controlled by its vertex-transitive
atoms.** It explains why the raw-`H` gaps transfer to `tc`, and it predicts the same reduction for any "count up to
automorphism" in the repository: closure numbers `n HC/|Aut|`, presentation counts of graphs (THM-4741's graph chart),
and the tc-fibres (finite for every value, THM-4812 (f)).

The aliquot side shows the limit of the analogy. The untouchables `{2, 5}` have no symmetry quotient. The "almost all
primes are touchable" statement comes from analytic number theory (Montgomery-Vaughan), not from an atom census. The
two gap pairs remain typed as ANALOGY (THM-4805); tonight's work makes the tournament side unconditional at every order
in both counts, and leaves the arithmetic side conditional.

## 2. Positive coordinates and coboundaries

Several advances are identities that become telescoping after a conjugation to a coordinate in which the dynamics
is a positive affine map.

- **Collatz.** In the numerator `X = 3s+1` the Syracuse step is `X -> 3X/2^v + 1`, the `4x+1` lift is `X -> 4X`, and
  the native ROOT edge is `X -> (3X + 2^a)/2^t` (paid ROOT R1, THM-4813 (A)). The ROOT forcing cocycle
  `log(1 + 2^a/(3X))` equals the drift cocycle `t log 2 - log 3` plus the coboundary of `log X`. Around a cycle it is
  the Steiner-Eliahou identity (THM-4633); at the fixed points it is THM-4600's anchors `c_a = 1/(2^a - 3)`.
- **Berggren.** In slope coordinates the three branches are the Hashimoto operator of a tree (THM-4808), and the branch
  character is a determinant read off the even continued fraction.
- **Circulants.** In the centred lattice the cyclic shift is multiplication by `zeta`, and the Paley skew operator is the
  Gauss sum (THM-4813 (B); klein's THM-4814 (E) for all doubly regular tournaments). The three ROOT rows are the same
  prime ideal times a unit.

**What the conjugation does and does not give.** It moves the content of an identity into its boundary terms: the
coboundary, the fixed points, the cusps. It does not bound the boundary. The ROOT readout, the Collatz entry, and the
Berggren cusp remainder all remain open exactly where the boundary term is uncontrolled. Part I's cancellation
question is the same phenomenon seen from the involution side.

## 3. What the night showed about method

- **Exhaust, then generalise the exception, then prove the classification.** The census at order 10 produced one
  exception; its structure gave an infinite subfamily within the hour, and the attempt to prove the corrected lemma
  produced the complete family and a classification proof (with one bounded finite check) a few hours later. The
  intermediate generalisation was too narrow, and the order-11 census refuted the pre-stated conjecture.
- **Retype on contact with concurrent work.** Two of tonight's claims (the sharp spectral constant, the doubly-regular
  lattice module) had concurrent published versions by the time they were filed. They are now cited and kept only as
  comparison points. Convergence is evidence; double counting is not.
- **Pre-register the cheap tests.** Two pre-registered tests (the `tc` spectrum at order <= 10 and the forcing tail to
  `10^8`) were committed before their runs (HYP-9421 (3)-(4)).

## 4. Short version (additions to part I's eight steps)

9. **For a distance, write the potential and the certificate.** Exactness means no local minima.
10. **Expect the local minima to be symmetric**, and find them at the first untested order.
11. **For a count up to symmetry, reduce to vertex-transitive atoms** and a finite orbit-group search.
12. **Conjugate to a positive coordinate** and read identities as coboundaries; then say which boundary term is still
    open.
