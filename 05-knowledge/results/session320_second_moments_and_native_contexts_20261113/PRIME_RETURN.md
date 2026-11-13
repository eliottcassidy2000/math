# A common strong context returns to a prime count

The [three-way signed context](SIGNED_TOURNAMENT_CONTEXTS.md) has a
second exact instance, at a prime target, with a different cancellation
ratio. Use the two strong ten-vertex source masks
```
X = (766,248,1018,832,264,280,944,312,3,434)
Y = (32,909,1,821,295,454,415,541,133,373)
```
with (H,Q(2),Q(3)) equal to (1915,56008,710883) and
(1915,55288,711603). Thus ΔF_2=720 and ΔF_3=−2880.
The common strong quotient with row masks (44,33,2,6,15,28), marked
at 0, has indegree two there and literal one-, two-, and three-run
counts (25,24,6). It evaluates a source to
```
25H+24F_2+6F_3 = -5H+6Q(2)+6Q(3).
```
Therefore both order-15 strong outputs have **4,591,771** Hamiltonian
paths. Trial division through the integer square root certifies that
this count is prime. Direct subset DP independently verifies both
counts and gives Q(2)=222226400 versus 224279840. The resulting
tournaments are neither isomorphic nor converses.

This is not a new coverage range: its point is a genuine prime
endpoint obtained by common-context cancellation of **opposing**
responses, outside the permanent factor-of-nine regime of the older
equal-count switch. It does not claim the prime endpoint is ready
to repeat the same context.

## An exact obstruction to repeating the first return with only three visits

For the order-13 outputs at H=532623 in the companion note, direct
DP gives Q(3)=216103869 and 219387789. Consequently the difference
in their nonempty two- and three-run cover counts is respectively
−134400 and −2880720. For any common strong homogeneous context that
can visit the marked factor at most three times, the Hamiltonian
difference is strictly negative: the itinerary coefficients are
nonnegative, and the two-run coefficient is positive (start from a
directed triangle through the marked block, repeat it once, and
insert the other quotient vertices). This argument applies to any
strong quotient of at least three vertices with the stated visit
bound; other substitution factors may be arbitrary nonempty
tournaments, since each can fill its one run in at least one way.
A further exact return therefore requires a new response at four
or more visits, an asymmetric context, or a selective contact.
