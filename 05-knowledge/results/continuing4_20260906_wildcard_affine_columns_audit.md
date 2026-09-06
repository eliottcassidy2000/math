# Independent audit of the affine-column decoder

**Status: PASS; PROVED elementary statements + FINITE-EXACT controls.**
Root read the complete candidate proof and standalone producer and rebuilt
the relevant finite objects without importing it. No mathematical correction
was required. This audit supports promotion of the linked primary report;
it does not claim a solution of no-three-in-line extremal asymptotics.

The exact triple event is correctly typed. Three distinct rows and columns
ensure a nonzero modular slope, with all inverses taken in F_p. Dividing the
row gaps by their gcd makes an integer collinearity lift have differences
km,kn with integral k. Their total span is |k|n. The representative of the
first transformed column must be in the stated closed interval; this is
both necessary and sufficient, including its two endpoints. Since n>=2,
the allowed signed k values cannot collide modulo p. A nonzero modular
determinant cannot become zero as an integer determinant after the affine
column map. Thus the event probability and forbidden-interval formula hold
on their full stated universe, not only on triples that were initially
integer collinear.

The full-board iff follows by a union of these exact events, after disposing
of three-in-a-row or three-in-a-column cases. If exactly two points share a
row or column, a collinear triple would require its third point to share it
too; therefore no discarded mixed triple is bad. An empty affine orbit is
not an empty unrestricted column orbit. The producer retains that scope.

Sharp two-transitivity proves the one-/two-column cylinder equality.
This extends to quadratic cell-indicator observables because each product
imposes images of at most two source columns. It supplies no equality of
triple laws or of the zero event. The independent code counts actual
collinear triples by canonical integer line equations and their complete
point sets. Across all120 S5 column permutations it recovers four successes,
the affine histogram {2:8,3:4,5:4,14:4}, and the full mean64/15. It directly
checks every one-/two-column cylinder and the one-swap repair(1,3,0,4,2).
The minimal odd prime assertion follows from AGL(1,3)=S3; no statement
about non-prime board sizes is implied.

The midpoint argument is a bijection for equal-step triples. An endpoint
pair of equal coordinate parities has an integer midpoint in the box and
on the modular line; distinct endpoint pairs cannot give the same unordered
three-point AP because its middle point is unique. Conversely every AP
has such an endpoint pair. Convexity of binomial(n,2) gives the balanced
four-class minimum. The slope-two line has exactly these balanced classes:
split x at floor(p/2), where parity of [2x]_p changes. Every non-axis modular
line is sent to it by an affine column map. Hence the minimum is sharp for
the AP substatistic; unequal-step collinear triples can add more. The two
parallel lines are disjoint, so their internal AP contributions can be added.

The short-vector proof is also complete. The image of the index-p tangent
lattice under (x,y)->(x+y,x-y) has index2p. There are more than2p points in
the stated (h+1)-square, so a same-coset difference lies in that image and
has sup norm at most h. Its inverse is integral by image membership, and
has L1 norm equal to that sup norm. A shortest nonzero vector has coordinate
gcd d<p; division is allowed modulo p and would shorten it unless d=1.
The determinant functional vx-uy is a nonzero linear functional on F_p
with the prescribed line as an affine level set. Its values in an N-square
span an interval of width(N-1)r; a single residue class occupies at most
floor((N-1)r/p)+1 integer levels. Each level contains at most two selected
points. This proves the necessary capacity bound also for N!=p and for
axis lines. The orthogonal-basis sharpness example checks out, without an
unproved claim that infinitely many primes have its special form.

The independent exact source covers all7956 distinct-row/distinct-column
triples at p=3,5,7; the entire S5 example using canonical line sets; all
1023 modular lines of finite slope at p=3,5,7,11,13,17,19; the omitted
vertical direction; midpoint bijections; box sizes2 and p+2; and the four
declared sharp tangent lattices. The proof supplies unbounded quantifiers.
All88570 always-active gates pass. The producer's117652 gates are separate
evidence; both paths must also replay unchanged with optimization before
filing. Reproduction from the repository root:

    python -B 04-computation/continuing4_20260906_wildcard_affine_columns_audit.py
    python -B -O 04-computation/continuing4_20260906_wildcard_affine_columns_audit.py

The retained missing coordinate is the triple affine ratio together with
its integer lift interval. No random-event independence, heuristic entropy
constant, or large-k theorem is a proved dependency of these results.
