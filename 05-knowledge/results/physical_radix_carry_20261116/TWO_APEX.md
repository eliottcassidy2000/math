# A literal two-move carry square with a binary residual tape

Use a transitive spine, two apices a,b, and a symbol z=2[a beats the
spine vertex]+[b beats it] at each spine site. The two changes below
are physical: swap two b-row arcs in fixed locations and reverse the
single mutual apex arc. Both changes are independent as edge edits,
but their count increments depend on their order.

For any t>=1 and any binary word w using symbols 0 and 1, consider
the four state tournaments with spine
```
A=1230 w 0^t,       B=1320 w 0^t
```
and mutual arc either a→b (write +) or b→a (write −).
The leading tapes differ at only two b-row entries; their a-rows
are identical. All four tournaments are strong: b beats the first
spine site, a reaches a later site, and every site reaches the last
spine site, which beats both apices. From a, reaching that last
site gives a route back through b to the first.

Let p=2^|w|, rho=sum_j w_j 2^j, and beta=sum_j w_j 2^(|w|−1−j),
with j running left-to-right from zero. The physical carry square is
```
A+ → B+: +(6p+beta) 2^t − 2rho
B+ → B−: −(6p+beta) 2^t + 2−6rho

A+ → A−: −(4p+beta) 2^t + 2−6rho
A− → B−: +(4p+beta) 2^t − 2rho .
```
Hence H(B−)−H(A+)=2−8rho, independent of the guard height t.
This realizes every nonpositive offset congruent to 2 modulo eight
as well as +2, in a genuine exact cancellation of growing terms.
The same guard can be grown arbitrarily after choosing rho.
The literal changes never exceed three arc flips, independent of t
and |w|.

With empty w, put L=36·3^t+6·2^t−11. The four counts are
```
A+ : L                      B+ : L+6·2^t
A− : L−4·2^t+2              B− : L+2 .
```
In this case the overall +2 is physical at unbounded height via either
order of the two moves. The four states do not form an adjacent
four-count interval; only the diagonal closes to +2. For general w,
the inequality 2−8rho=2 forces rho=0. Choosing an arbitrary residual
offset does not choose a baseline independently: both words and
output order change with rho. No interval theorem or proof of
Chan–Pak follows.

## Exact proof by a finite linear-state certificate

Use the inherited odd-cycle packing formula H=1 plus the sum of
packings weighted by 2 to the number of cycles. Each cycle includes
an apex and uses one or two increasing runs through the spine.
There are exactly three contribution types: one cycle through
one apex (weight 2), two disjoint cycles, one through each apex
(weight 4), and one cycle through both apices (weight 2). In the
third type the two spine runs go from a to b and b to a, one possibly
empty if its mutual arc exists. A loop on one apex needs a nonempty
even run; a mixed cycle needs odd total run size.

Each run needs just five states: unused, open odd/even, closed
odd/even. At a spine position either skip it, or assign it into
exactly one run. A new run may start only if its source apex beats
that vertex; an open/new run may close there only if that vertex
beats the target apex. Once a run closes it has no more assignments.
The Python replay specifies these rules literally, without
optimized or randomized linear algebra. Its four block dimensions
5,5,25,25 track the two individual loops, the pair of loops and the
pair of cross runs.

For any fixed active prefix, appending t zero symbols with t>=1
leaves only completed and already-open runs available. Each open
run must receive at least one additional vertex, then close at its
last assignment. The possible new assignments are counted by
inclusion-exclusion for nonempty runs and projections onto their
required parities. With at most two runs, the bases of the resulting
powers are 1, 1±1 and 1±1±1. For positive t these reduce to
1, 2^t, 3^t, and (−1)^t. The four sequences evaluated at t=1,2,3,4
are independent. This supplies the all-height proof from a bounded
number of readouts, without relying on interpolation in t beyond
the proved basis.

For the binary finishing tape, take the two coupled prefix states,
along with virtual coordinates (rho,p,beta,1) initialized at
(0,1,0,1). Reading b=0 or 1 sends those four coordinates to
(rho+b p,2p,2beta+b,1). Build the exact rational span of the
initialized vector under the two transitions, admitting a new
basis vector only when independent. This span has dimension 11;
closure under the two transitions is verified by pushing and
reducing each basis candidate. Thus every word w is in the span.
The three odd-cycle contribution differences B− minus A+ are,
on every basis element at t=1,2,3,4,
```
2^(t+1) p − 2rho,       0,       2−6rho−2^(t+1) p.
```
The two same-orientation sum differences B minus A are respectively
(6p+beta)2^t−2rho and (4p+beta)2^t−2rho. All the observations
and claimed answer expressions are linear in the coupled state.
Adding the three diagonal observations and subtracting the
same-orientation observations gives every edge of the displayed
carry square. The suffix basis argument extends these certificate
identities to every t>=1.

Run `python3 check_two_apex.py`. A separate direct enumeration of spine
allocations into runs checks the four individual empty-w states'
polynomial coefficients, so the reported (L) is an actual
Hamilton path count, not just a predicted difference.
