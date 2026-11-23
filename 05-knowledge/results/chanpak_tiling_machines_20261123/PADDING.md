# The two singleton gaps are an unbounded tape at fixed Hamilton height

**PROVED scoped.** This is a construction in the variable-order
marked model with singleton insertion, singleton deletion and a
finite instruction table. It is not a claim that the raw finite
complement-line graph, or testing whether \(h\) occurs, is universal.

Fix any attained \(h>1\). By deleting all singleton SCCs, choose
one residual tournament \(C\) of height \(h\) with at least one
nonsingleton component and no singleton components. Define
\[
 T(a,b)=\mathrm{TT}_a\ \triangleright\ C\
          \triangleright\ \mathrm{TT}_b,\qquad a,b\in\mathbf N.
\]
Here \(\triangleright\) is the ordinal join and
\(\mathrm{TT}_0\) is empty. Every Hamilton path is the unique
transitive left path, a Hamilton path through \(C\), and the
unique transitive right path. Thus \(H(T(a,b))=h\) exactly.
The SCC order recovers \(a,b\) separately from an ordinary
isomorphism class: count leading and trailing singleton
components. In a converse-merged model the register ordering
must be retained as a mark unless already recoverable from the
oriented core.

Use these two gaps for a [two-counter machine](https://www.cs.ox.ac.uk/people/taolue.chen/pub-papers/Tamc09.pdf).
The finite program counter is a mark; the program has increment
and conditional decrement instructions with their exact zero
branches. Add one singleton next to the marked core at the
chosen side to increment. If there is a singleton on that side,
remove the nearest one to decrement; otherwise take the
zero branch. These operations leave all retained arcs unchanged
and preserve \(h\) at each step. Take any fixed Hamilton frame of
\(C\) and extend it by the two transitive frames; the insertion
or deletion has a consistent actual tiling description at the
new size. Halting is the marked terminal instruction, not a
gap in the Hamilton spectrum. The map from counters to classes
is injective, both update directions are exactly transported,
and overflow does not occur in the unbounded-size model.
Consequently **uniform marked halting is undecidable inside
every attained height fiber \(h>1\)**. It is already possible
at \(h=3\) with \(C=C_3\), the smallest nontrivial height.

There is a sharp low-end contrast in this same scope. At \(h=1\)
every tournament is transitive (a directed triangle already
contributes three by the odd-cycle packing formula, or by
extending its paths). With graph order as the only unbounded
storage, the model is a one-counter finite controller. Its
deterministic halting is decidable. At positive counter values,
control goes through a finite fixed transition graph until it
halts, reaches zero, or repeats a control state. A repeated
positive segment with nonnegative net displacement will
continue forever; with negative displacement it reaches zero
in finite time, since unit steps cannot jump over zero.
There are finitely many zero-control states; a repeat of one
of those is an actual configuration repeat. This gives an
algorithm. Permitting separate movable markers inside a
transitive chain would add memory and is **outside** this
one-counter model.

## What the standard observables forget

The entire ordinary slot polynomial sees the sum, not the split:
\[
 Q_r(T(a,b))=r^{a+b}Q_r(C).
\]
In particular \(T(0,1)\) and \(T(1,0)\) have the same order,
height and every slot count. An instruction testing the left
gap can halt from the first and take a permanent loop from
the second. This gives a single explicit obstruction to
predicting the next controller branch even with all homogeneous
slot counts. The relevant information is register position.
It is accessible to the gap marker but not to these scalar
observations.

At a fixed finite order there are only finitely many
configurations in any specified finite-controller model.
For the fixed core and order \(n\) there are just
\(n-|C|+1\) choices of the two gaps. Thus this universality
does not arise at a fixed order; it uses a coherent family of
increasing allowed orders. This is the dynamical counterpart
of the [inherited fixed-height series](../chanpak_support_birth_grammar_20261122/FIXED_COUNT-SERIES.md):
an ordered strong core with \(k\) components exposes \(k+1\)
independent singleton gaps when order is free.
