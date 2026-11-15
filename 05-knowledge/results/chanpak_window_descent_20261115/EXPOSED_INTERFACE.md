# The bounded-ratio exposure theorem

A proposed descent starts from a witnessed smaller count \(h\) and builds
a literal tournament for \(N\le C h\). A tempting black box is to insert
the whole donor as one uniform module. The following obstruction is
pointwise: the host may change with the target, and only this one donor
needs to work. No neutral padding or universal quantifier over donors
is needed.

**Theorem.** Let \(R\) be an \(n\)-vertex tournament, and substitute it
for a vertex \(q\) of a host \(G\), with uniform arcs between \(R\) and
each other vertex. Write \(h=H(R)\).
* If \(q\) is in a singleton strong component of \(G\), the new count is
  exactly \(H(G)h\).
* Otherwise the output has at least \((n-1+H(G))h\ge(n+2)h\)
  Hamiltonian paths.

**Proof.** A path in the output gives a host word by collapsing each
maximal run in \(R\) to a symbol \(q\). It contains each other host
vertex once and \(q\) some positive number \(j\) of times with no
consecutive copies. Between successive symbols it respects the arcs
of \(G\). For each such word, its lifts are in bijection with ordered
nonempty \(j\)-path covers of \(R\). In particular \(j=1\) gives
\(H(G)h\), and one \(j=2\) word adds at least \((n-1)h\): every split
of every Hamilton path of \(R\) into two nonempty sections gives a
different ordered cover.

A word using \(q\) twice gives a directed return through \(q\), so is
impossible if \(q\) is a singleton strong component. Conversely, if
\(q\) is in a nontrivial component, it belongs to a directed triangle.
For instance, start with a shortest directed cycle through \(q\) and
shorten it along a chord from \(q\) to its second successor, or use
the resulting triangle if that chord goes into \(q\). Write the
triangle as \(q,a,b,q\). Insert the other vertices one at a time in
this word by ordinary Rédei insertion. The insertion proof works
unchanged when \(q\) is repeated: if the new vertex does not go at an
end, scan for a neighboring pair whose arcs permit insertion. It
preserves both copies, hence gives a two-visit word. Finally
\(H(G)\ge3\), by inserting into each of the three paths of the
triangle, with different initial paths staying distinguishable on
restriction. This proves the claim. ∎

## What this says about descent

* For fixed \(C\), whenever \(n>C-2\), a bounded-ratio construction
  reusing the full donor as a uniform module is necessarily a
  product. This holds for target-dependent hosts of unbounded size.
* For a prime target, no nontrivial product donor can work. For
  \(N=7p\) with \(p\) prime and \(p\ne7\), no proper nontrivial
  product donor works: its factors would include the excluded count
  7. In particular a proposed \(p\mapsto7p\) lift as a single
  uniform module works only if the donor has at most five vertices,
  so it cannot work for unbounded \(p\). This is consistent with
  63: \(p=9\) is not prime, and hitting 63 is not uniform
  multiplication by 7 in a large module.
* A near-identity widening step has it worse: if \(N<3h\), a cyclic
  host is impossible even without a growing donor. A singleton
  host SCC gives an odd product ratio, so \(N=h\) is the only
  possibility with \(h\le N<3h\).
* This does not exclude grafts that expose selected endpoints,
  overlapping substitutions that do not leave the donor a module,
  two smaller donors at very different scales, or moves that
  reconstruct a witness instead of extending it.

The slot expansion under arbitrary uniform hosts is inherited
([Outlet width](../chanpak_carry_admission_20261115/OUTLET_WIDTH.md)).
The earlier [fixed-template obstruction](../chanpak_least_missing_jacobian_20261114/BLACK_BOX_AND_SYMMETRY.md)
uses arbitrary padding to reject universal scalar templates.
The pointwise consequence here instead forces exposure of a large
donor even after picking a host for each particular target. This
makes the positive supply problem sharper: find small **selective**
interfaces that expose the required carry without reopening the
whole donor.
