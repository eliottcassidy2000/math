# Closing the even-exponent growth question on the 7 times powers of 3 ray

In the [Session 24 multiplicity study](../../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md),
section 3 left a precise obligation on the even ray H=7*3^(2m):
there can be at most degree m-2 growth in n for the self-converse
class multiplicity, because the center 63 does not exist. Equality
was not assumed. A literal center of count 567=7*3^4 closes it:

    470  508  408  145  480  205  396  256  40

The row indexed i is the bitmask of vertices beaten by i,
0-based. Reflection i -> 8-i reverses all arcs. The tournament has
the increasing Hamilton frame and passes the strong-prefix-cut
criterion. Direct subset DP gives H=567 (see find567.cpp and the
independent direct mask check below). It is not an arbitrary new
7p prime selector: its quotient is 81.

For every m>=2 the count at issue is now obtained by a middle copy
of this tournament and two opposite banks with m-2 copies of C3
each. The smallest literal construction in this family is on
9+6(m-2) vertices; no global minimum claim is made. For any larger
odd n, distribute (n-9)/2-3(m-2) singleton components over the
m-1 gaps in one bank and mirror that bank. These are distinct classes,
because their ordered strong component words are distinct. This
gives a lower bound, as n tends to infinity through odd orders, of

    binom((n-9)/2-2(m-2),m-2),

with leading coefficient 1/(2^(m-2)(m-2)!) in n. The known upper
bound has the same degree: the square half must consume at most
m-2 triangles' worth of factors of three, since centers 7,21,63
are impossible for the indicated choices; finitely many bounded
nontrivial cores and their singleton gaps then give O(n^(m-2)).
At even orders there are zero classes by the even blue mod-4
law (THM-791), since 7*3^(2m)=3 mod 4. At m=1 all orders are
empty. Thus the even ray has exactly the previously only bounded
degree m-2, on its odd-order sheet. The earlier odd ray was already
closed by the H=189 center; these are complementary mechanisms on
the same 7 times powers of 3 chain.
