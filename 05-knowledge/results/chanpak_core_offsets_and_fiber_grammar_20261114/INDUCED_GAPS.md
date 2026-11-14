# A fixed positive difference at unbounded height forces non-nesting

This strengthens the one-vertex obstruction in [METAGRAPH](METAGRAPH.md).
Let A be **any induced subtournament** of T, by deletion of an arbitrary
number of vertices. If H(T)-H(A)=J>0, then the ordered pair of counts
belongs to a finite, effectively enumerable set depending only on J,
not on the number of deleted vertices or on the tournament order.

Use the [odd-cycle-packing identity](https://arxiv.org/html/2412.10572v3#S3.SS4).
Factor out the product m of H over untouched nontrivial SCCs of T
(singletons, whether deleted or not, contribute 1). Denote the remaining
nontrivial SCCs meeting the removed set by C_i, i=1,…,r. Put
a_i=H(C_i∩A) with H(empty)=1, b_i=H(C_i), and d_i=b_i−a_i.
Then J=m(∏b_i−∏a_i). Taking any missing vertex in C_i, its
odd cycles of lengths 3,5,… furnish disjoint choices in the packing
identity for H(C_i)−H(C_i−v). Hence d_i>=2 and |C_i|<=d_i+2.
Expanding ∏(a_i+d_i) and using a_i>=1 gives both
sum d_i<=J/m and 3^r−1<=J/m. The difference of the odd products
is even, so m is an odd divisor of J/2. For each candidate m enumerate
the finite list of pointed-subset strong classes of order <=J/m+2,
with r<=floor(log_3(1+J/m)). Pair (m∏a_i,m∏b_i) occurs exactly
when some such product has gap J and m has a physical witness;
for sufficiency join that witness and the listed C_i ordinally,
then delete their specified subsets. Singleton padding is arbitrary.

For example J=2 forces r=m=1 and a strong block of order at most four,
so the endpoints are still only (1,3) and (3,5), even for arbitrarily
many deletions. A construction that lifts a pair of consecutive odd
counts to unbounded absolute height must do something beyond just
making one member an induced subtournament of the other. The
changed-core windows and the simultaneous return in this packet do
exactly that: their endpoint tournaments need not be nested.
