# Quantitative selective-apex sparsity and a two-contact inverse

Let B be a nonempty base tournament, and let a new vertex v beat a
set U⊆V(B). Write H for the full Hamilton path count. The proofs
below use the Hamilton odd-cycle packing formula H(T)=∑Packing 2^c
and the corresponding ordered two-cover count Q(T)=∑Packing 2^(c+b),
where c counts odd cycles and b counts cycles plus uncovered vertices.

## There is no bounded-complexity selective decoder

Set L=log₂X and C=3/log₂3. Consider target counts at most X obtained
from |U|≤r and H(B)≤H₀. The same results hold if at most r arcs have
the opposite polarity, by taking converse.

Order the SCCs of B from sources to sinks. Let a be the earliest
component containing an out-neighbor of v and b the latest component
containing an in-neighbor of v. If either does not exist or a>b, v is
inserted as a singleton SCC and H(B+v)=H(B). Otherwise the new SCC
containing v consists exactly of v and the components in [a,b].
Indeed v can reach each via the first SCC, and each reaches v through
the last SCC, whereas those outside the interval have only the
indicated ordinal orientation toward the new SCC.

All singleton SCCs strictly outside [a,b] can therefore be discarded
without changing H(B+v), since H multiplies along ordinal SCCs and the
discarded factors are one. If a<b, choose a fixed u in component a
beaten by v and w in component b beating v. Let s count the singleton
SCCs strictly between them. For each even-sized subset of those
vertices, use the cycle v→u→(the selected vertices in order)→w→v.
It is directed and odd, and contributes two to the single-cycle
terms in the packing formula. Consequently H(B+v)≥2^s for s≥1;
the s=0 bound is automatic. If a=b, there are no interior singleton
components. In either case at most L+2 singleton SCCs are retained.

Let m count all vertices in nonsingleton SCCs of B, including those
lying outside [a,b]. The Busch bound 5^((n−1)/3) for a strong
tournament, together with the order-three and order-four minima 3
and 5, gives n≤C log₂H for each such component. (For n≥5 the Busch
bound alone suffices; check n/(n−1)≤5/4.) Multiply the component
counts to obtain m≤C log₂H₀. This uses only a lower bound, not any
classification conjecture.

Now choose a labeled tournament on m vertices to encode the
nonsingleton ordinal word. Its SCC segmentation is already encoded,
since labels may always be chosen to respect the word. There are at
most 2^binom(m,2) such choices. Relabeling the encoded components in
word order, insert the retained t≤L+2 singleton SCCs in at most
binom(m+t,m) ways, and choose the marked subset of at most r
vertices. Thus a coarse bound on the number of outputs is
H₀ (for the noncycling case) plus
[
 sum_{0le mle Clog_2H_0, 0le tle L+2}
 2^{inom m2}inom{m+t}{m}sum_{j=0}^{r}inom{m+t}{j}.
]
Out-of-range binomial coefficients are zero. This construction is an
overcount and ignores symmetries, not a claim that every coded choice
satisfies the base count cutoff.

In particular, for fixed H₀ and r the count is polylogarithmic in X.
For H₀=2^(α√L), r≤ρL, with 0≤ρ<1/2 fixed, the count is at most
[
 X^{C^2alpha^2/2+mathcal H_2(
ho)+o(1)},
]
where the binary entropy is continuous at 0. Indeed the position
factor has log₂ O(√L log L), the graph factor has the stated leading
exponent, and the mark factor is bounded by the usual binomial
entropy estimate for at most (1+o(1))L slots. For bounded r the mark
factor instead has log₂ O(log L) and the threshold in α is √2/C≈.747.
For α=.5, ρ=.1, the exponent is below .918. Thus allowing even a
subexponential-in-log-X source path count and logarithmically many
marked contacts still does not fill the odd targets.

Reference for the strong bound: Arthur H. Busch, A Note on the Number
of Hamiltonian Paths in Strong Tournaments, DOI 10.37236/1141.

## A physical inverse in a four-block chart

Take B=p→A→q→D, with A,D nonempty and p,q singleton, and let v beat
exactly p,q. Put h_X=H(X), Q_X=Q(X). Then
$$
 N=H(B+v)=(2Q_A+h_A)Q_D-(Q_A-h_A)h_D .
$$
An odd cycle through v leaves an even-order path from p or q to an
unmarked endpoint. Those paths from p ending in A contribute h_DQ_A
to H, including the odd-cycle factor two. Those from q ending in D
contribute h_AQ_D. Here signed parity cancels in the full two-cover
convolution of the nonempty A or D: each odd cycle or singleton
independently chooses a slot with signed weight ±1. An odd subset is
automatically nonempty. Lastly the paths starting p and ending in D
contribute 2Q_A(Q_D−h_D), as the optional q fixes parity uniquely
whatever the two chosen subsets in A and D. These classes exhaust the
cycles through v. Paths counted without new cycles contribute h_Ah_D.

For A transitive of order a≥1, setting t=2^a gives
N=t(2Q_D−h_D)+(Q_D+h_D). Given D and N this is an exact inverse:
(N−Q_D−h_D)/(2Q_D−h_D) must be a positive power of two of exponent
at least 1. Given A and N it is instead a linear Diophantine relation
on h_D,Q_D; solving it arithmetically does not supply a physical D.

The fixed-A coefficients have gcd gcd(3h_A,Q_A−h_A), and outputs
also carry gcd(h_D,Q_D). More sharply, reducing the packing expansion
mod 3 gives Q_X≡(−1)^|X|h_X because every packing block has odd size.
The formula therefore gives 3∤N exactly when |A| is odd, |D| is
even, and neither h_A nor h_D is divisible by 3. Thus one does not
need to search the other parity branches for prime or 7p targets
other than the one allowed multiple of 3 in each target category.

For empty A and nonempty D the count is 3Q_D−h_D; for empty D and
nonempty A the count is Q_A+h_A. Both empty gives 1. These boundary
formulas have a different residue test and are not included in the
nonempty-block assertion.
