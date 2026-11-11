# Independent audit of the 105 reset

**Verdict: the explicit reset and the stated finite bottlenecks are valid.**
Audit by the source-class lane, independently of the reset author. The
all-order conclusion concerns only neutral single-vertex operations; the
bottleneck9 assertion is correctly restricted to fixed orders8 and9.

The repaired proof of
[THM-1900(B)](../../../../01-canon/theorems/THM-1900-the-insertion-response-calculus.md)
was read directly. Its terminal-SCC cut criterion is equivalent to the new
vertex being a singleton SCC. Thus neutral births/deletions preserve the
ordered word of nonsingleton SCC isomorphism classes, while singleton gap
lengths can vary. The converse in R1 is correct for unmarked tournaments;
its caveat about ports on removed singletons is necessary and present.
The disjoint105 branch kernels therefore do prohibit this operation menu
at every order. Count-dipping replacement steps are outside that menu.

The displayed matrices A,B,C were reconstructed independently using
`../targets/check.py`. Their counts are105. Each stated deletion card,
new-vertex cut and exact surviving-name map was checked entry by entry:

* Delete A's vertex5: card count19; cut123 reconstructs B under the
  printed mapping, with the new vertex mapped to B's vertex6.
* Delete B's vertex4: card count9; cut34 reconstructs C under the
  printed mapping, with the new vertex mapped to C's vertex4.
* Cut186 on A gives317. Cuts26 and31 on C both give319.
* Composing the two maps preserves exactly the displayed six original
  names. Every arc between any two surviving names is unchanged at
  its replacement step. This does not imply preservation of their
  Hamiltonian deletion or attachment statistics.

The complete source universes were rebuilt separately from the inherited
strong bank. Ordered exterior padding gives19 classes at order8 and37
at order9, with distinct Hamilton-order canonical codes. No3·35 source
fits at either order. For each source every vertex deletion was computed
and cards were grouped by their canonical code. Rather than rerunning
the author's maximum-bottleneck algorithm, this audit tested ordinary
connectivity after deleting edges below a threshold:

| Source order | Classes | Branch families connected at threshold9? | At threshold11? |
|---|---:|---|---|
|8|19|yes|no|
|9|37|yes|no|

All residual counts are odd, so these two tests independently establish
the optimum9 in both fixed-order graphs. An isomorphism between deletion
cards is exactly the information needed for one legal delete/insert
replacement with all surviving arcs fixed. No path computation or
threshold test here addresses larger orders, multiple-vertex moves, or
old-edge reversals. R3 explicitly preserves those limitations.

One checker-maintenance suggestion was sent to the author: derive the
exported count trace from the selected witnesses, or assert the two
card counts19,9, rather than relying only on a hardcoded trace. The
currently displayed witnesses and counts were checked independently and
are correct. There was no mathematical correction to the proof.
