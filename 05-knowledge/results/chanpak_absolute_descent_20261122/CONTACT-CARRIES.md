# When separated contact switches cannot encode independent summands

Start with the transitive spine \(0<1<\cdots<n-1\) and one apex. Write \(x_i=1\) when the apex beats spine vertex \(i\). Its exact ordinary path count is
\[
 H(x)=1+\sum_{i<j} w(j-i)x_i(1-x_j),\qquad
 w(1)=2,\quad w(d)=2^{d-1}\ (d\ge2).
\]
Indeed every directed cycle uses the apex, so the odd-cycle packing identity reduces to one plus twice the number of odd cycles. If the first and last spine vertices on such a cycle are \(i,j\), intermediate vertices form an arbitrary even subset of the open interval. The empty interval has one even subset; any nonempty interval has half its subsets even.

A binary switch on a contact block is a pair of fixed masks on that block; its change is a vector with entries \(0,\pm1\). Say two switches are **additively independent** if the mixed difference of the four count readings is zero, so each switch's contribution is independent of the other's choice. Consider nontrivial switches on disjoint ordered blocks, one at indices \(\le b\) and one at indices \(\ge b+1\). Let their differences be \(\delta_i,\epsilon_j\), respectively, and define
\[
 A=\sum_{i\le b}\delta_i2^{b-i},\qquad
 B=\sum_{j>b}\epsilon_j2^{j-b-1}.
\]
Then the mixed difference is **exactly**
\[
 -AB-\delta_b\epsilon_{b+1}.
\]
This follows by taking the cross term of the cut formula: the bulk weights factor, and the adjacent boundary has one extra unit. Both \(A\) and \(B\) are nonzero because the largest changed dyadic place outweighs all smaller ones.

**Classification.** If separated by even one unchanging spine vertex, two nontrivial switches never commute additively. For touching blocks, they commute additively exactly when both boundary positions change, one switch changes only its boundary position, and the other is a nontrivial carry: a consecutive run from its boundary with sign \(\eta\), plus the next position with sign \(-\eta\). To see this, cancellation forces \(AB=-\delta_b\epsilon_{b+1}\) and hence \(|A|=|B|=1\). If a signed binary sum of this type has more than one nonzero place and magnitude one, every place below the largest must be present with the opposite sign, by equality in the geometric sum bound. Its resulting sign is then opposite its boundary sign. For a one-place switch the signs agree. The desired opposite product therefore requires exactly one carry.

Consequently a one-apex transitive completion with three or more disjoint ordered nontrivial switches cannot implement an additive digit bank: the first and third always have a nonzero mixed coefficient, regardless of fixed surrounding contacts. This statement is about independent switches on different blocks, not about feedback, interleaved supports, or other skeletons. For a possible circuit, the carry exception shows that eliminating the interaction is an exact boundary phenomenon rather than a modular one.
