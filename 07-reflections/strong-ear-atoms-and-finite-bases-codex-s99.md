# Strong ear atoms and finite bases

The S98 strong-component atom ledger answered which contraction preserves
`H`.  The S99 extension answers what the inverse move should remember.

Adding one vertex `x` to a tournament `T` is controlled by the cut
`sig[v]=1` iff `x -> v`.  The Hamiltonian-path count is not a mysterious new
global invariant after the insertion.  It is an exact boundary cut polynomial:

```text
H(T+x) = start(sig=1) + end(sig=0) + Q(sig=0, sig=1).
```

Here `Q[a,b]` counts old vertex orders with every adjacent edge valid except
possibly the slot from `a` to `b`.  So the new vertex either starts a path,
ends a path, or repairs/occupies one exposed slot.  This is the tournament
version of the LRC finite-witness character sum: a main/boundary term plus a
finite resonance ledger.

Two facts make this the right atom calculus.

First, if `T` is strong, every nonconstant cut gives a strong child.  The new
vertex has at least one exit into `T` and at least one entrance from `T`; the
strong parent supplies all remaining reachability.  Strongness is therefore a
cut condition, not a graph-minor shadow.

Second, through `n=8`, every strong tournament can be reduced by deleting a
vertex while staying strong, except the base `C3`.  The exact deletion-count
distribution has minimum `2` for every `n=4..8`.  The resulting ear generation
is complete through `n=8`: nonconstant ears from strong parents generate
exactly the full strong H-spectrum in each transition `3->4` through `7->8`.

The finite-basis phenomenon is the striking part.  In the largest audited
transition `7->8`, balanced ears of weight `3` already cover `295` of the
`297` strong H-values.  They miss only `49` and `75`; adding weight `1` ears
covers those.  On the parent side, a greedy basis of `12` parent atoms covers
the whole spectrum, with the first five parent atoms covering most of it.

S36 clarifies what this boundary exception means.  The S99 `w` is cut weight,
not atom count, but the missed values `49` and `75` are themselves verified
single irreducible strong atoms at `n=7`.  So the `w=1` boundary ear should be
read as a way to expose a single-atom obstruction unit, not as a three-atom
product branch.

This is the same shape as the incoming S34 rational-witness basis.  One
denominator can resonate to zero, so S34 uses a finite denominator basis.  One
balanced ear weight can miss two boundary values, so S99 uses a finite cut
basis.  In both cases, the proof state is not the raw set of runners or
vertices.  It is the labelled residue/cut boundary plus the resonance ledger.

After the S35 pull, the analogy becomes an atom-to-atom dictionary.  HYP-+2878
uses residue-covering atoms: the unsafe APs cover `Z/p` only when a finite
resonance component is over-aligned at that denominator.  HYP-2879 uses
strong-ear atoms: a cut is evaluated by boundary vectors and the exposed-slot
matrix `Q`.  A single denominator or a single cut weight can fail by resonance,
but a basis breaks the alignment.  S35's correction matters: prime-covering is
not an exponentially rare per-prime event, it is a structured interval-covering
event measured around `10%`, so the LRC proof should expect a many-prime finite
basis plus CRT over-determination.  The tournament cut-weight basis is much
smaller in the audited range, but the structural reason is the same rather
than cosmetic.

This gives a more concrete program for the strong atom conjecture:

1. Prove strong-ear reducibility for all strong tournaments on `n>=4`.
2. Formalize the insertion formula.
3. Use balanced ears as the generic branch and boundary ears for exceptional
   low values.
4. Reinterpret `{7,21}` as no-solution cut-polynomial values before the
   strong-min boundary.  The even-graph shadow can still index states, but it
   is not the proof carrier.

KPS S31e makes the E7/odd-hole thread concrete: `E_7` has exactly `1496`
chordless `C_5` holes and `196 = 14^2` chordless `C_7` holes.  The S99 route
predicts that if E7 non-chordality, winding tournaments, and `{7,21}` share an
obstruction, it should appear as an odd cycle in the exposed-slot/cut-resonance
graph built from `Q`, not as an ordinary minor of the runner set.  The next
test is to build that `Q` graph for the boundary parents that produce or miss
`H=49,75` and ask whether its odd holes align with the E7
pentagons/heptagons.  The S31e enrichment says to look for complement-paired,
high-centrality heptagon classes, not a naive 7-edge cycle model.

Tournament Analysis: the vertices were proof carriers:

```text
ear-boundary-Q-cut -> strong-component-H-atoms -> OCF-packets
-> finite-denominator-residues -> even-graph-shadow -> raw-H-scalar.
```

The challenged assumption is that atoms are vertices or even-graph edges.  For
growth, the atom is the boundary state `(start,end,Q)` together with a
nonconstant cut.
