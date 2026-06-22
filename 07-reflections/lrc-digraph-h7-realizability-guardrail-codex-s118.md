# LRC digraph H=7 realizability guardrail

The useful correction is simple: binary is not enough.

Tournaments have binary orientation states per unordered pair, and this is the
setting where `H=7` is forbidden.  Ordinary digraphs also have binary states
per ordered arc, but they realize `H=7` already on four vertices.  Incomplete
oriented graphs realize `H=7` on five vertices.  So the obstruction is not
"two states plus seven"; it is "two tournament states plus OCF
realizability."

This makes the LRC14 target sharper rather than weaker.  A counterexample
would have to induce the forbidden `K_3` in the tournament conflict-graph
image.  If the apex-7 over-cover only induces a generic binary relation, the
forbidden theorem cannot be imported.  If it induces an OCF packet with
`Omega=K_3`, THM-200 blocks it immediately.

The winding-tournament check has the same lesson.  Sampling on wall times
creates ties and is no longer a tournament, so evenness or `H=0` in a grid
output is not meaningful evidence against Redei.  On tie-free AP7 cells, `H`
is odd and avoids `7`; at the exact lonely point `x=1/7`, the winding
tournament has `H=175`.

So the next proof move is not another scalar analogy around `14=2*7`.  It is a
realizability theorem:

```text
LRC14 apex over-cover -> labelled tournament OCF conflict graph Omega=K_3.
```

That is the precise "impossible to disprove" route.  The repo should keep the
digraph guardrail visible, because otherwise it is too easy to prove a false
statement about binary shadows that ordinary digraphs refute immediately.

Post-rebase S48 gives the right packet to try next: the AP14 boundary has
seven tied diameter comparisons under the antipodal order-2 symmetry.  Those
ties are not tournament arcs until resolved, which is exactly the point.  The
proof should show that any resolution capable of pushing below `1/14` realizes
the forbidden `Omega=K_3`.  If the argument stops at "there are seven tied
binary choices," it is still in the digraph shadow; if it reaches labelled OCF
conflict data, THM-200 can do the work.

Post-rebase S31z also cleans up the phrase "impossible to disprove."  For
LRC14, a disproof is a finite checkable integer counterexample; ruling out all
such counterexamples is the conjecture.  So the forbidden-`K_3` route is not a
weaker metamathematical claim.  If completed, it is the proof.
