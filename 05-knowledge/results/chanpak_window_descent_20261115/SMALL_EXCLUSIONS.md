# The two exclusions are small seam obstructions

This is a short structural certificate of the **known** exclusions,
not a novelty claim. The key source is [Busch, Lemma 1 and Corollary 1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf);
the repository's original bounds-and-enumeration route is
[THM-1370](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).

Hamilton counts multiply over strongly connected components. Thus
unless a component has 7 or 21 paths, neither total is possible:
7 is prime, while a nontrivial factorization of 21 needs 7. Busch's
floor is at least 9 at order five and at least 25 at order seven
and above. The only strong tournaments of orders three and four
are the directed triangle and the one extreme reversal of a
transitive four, with counts 3 and 5. This settles 7, and leaves
only orders five and six at 21.

Busch's Lemma 1 says every strong tournament on at least five
vertices is either nearly transitive (with $2^{n-2}+1$ paths)
or consists of two strong induced subtournaments sharing one vertex.
At order five the overlap is triangle/triangle and leaves four
unknown seam arcs. At order six it is triangle/strong four and leaves
six. The triangle has one rooted type; the strong four with outgoing
rows $(6,8,10,1)$ has two rooted types **up to converse**, since
$(0,1,2,3)\mapsto(3,2,1,0)$ is an anti-isomorphism. That reduces the
remaining check to 16 plus 128 seam orientations. The replay counts
paths by direct permutations and shows neither union contains 21.
Conversing both components at once preserves the full count.

What this shares with the growth problem is the need to keep the
**actual seam**. Busch uses a shared vertex to inject a product of
path sets into the glued count; the injection provides a *floor*, not
a recipe for the excess. At 21 the small excess is a finite forbidden
choice of seam response. For cofinal realization the goal is the
dual kind of theorem: some large family of legal responsive seams
contains the exact necessary excess, including the carry. An injection
or a congruence can prove neither of those selection statements.
