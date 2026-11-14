# The old third vertex is not a free clock

The same seven-vertex cores have a sharp boundary when one inflates
vertex 6 rather than inserting a new painted bank. Allow arbitrary
nonempty tournaments X,Y,Z replacing vertices 0,1,6 synchronously
in the 67 and 69 arms. Write h_U=H(U) and c_U=C_2(U), the number
of ordered pairs of nonempty paths covering U. The exact difference
69-arm minus 67-arm is

    2 h_X h_Y h_Z - (2 h_X+c_X)(h_Y+c_Y)c_Z.             (1)

The complete projected path-word response only differs in five
visit entries, at (1,1,1), (1,1,2), (1,2,2), (2,1,2), (2,2,2),
with signed differences 2,−2,−2,−1,−1 respectively. The marked
subtournament on 0,1,6 is transitive in both arms. Between any two
successive occurrences of one marked symbol an unmarked symbol
must occur (otherwise it would close a cycle inside a transitive
set), so the maximum multiplicity is five. This makes the listed
table complete by finite exact recursion, not by fitting H on a
few module orders. [Domino proof and independent controls](TRIANGLE_6769.md).
For any such projected word, the lifts multiply the ordered
nonempty cover numbers of its corresponding modules. This gives (1).

## Unit-gap equality cases

For any nonempty U on n vertices, splitting each Hamilton path at
all n−1 possible cuts gives c_U>=(n−1)h_U. The first and second path
together recover the original path and its cut, so there is no
overcount. For n=2, c_U=2 and h_U=1; consequently c_U>=2h_U for
every nonsingleton U. If Z is singleton, (1) is +2 h_X h_Y:
the gap is exactly +2 iff both X,Y are transitive.
If Z is nonsingleton, (1) is at most −2 h_X h_Y h_Z. Equality
at −2 requires h_X=h_Y=h_Z=1 and c_X=c_Y=0, c_Z=2.
Thus X,Y must be singletons and Z the transitive pair.
These are all cases with absolute gap two. A new painted gate
is fundamentally different from inflating this old guard.

At X=Y=singleton and Z=TT_r the actual two counts are

    34·2^r−1 (67 arm),       32·2^r+5 (69 arm).

Indeed the vertex has exactly one in-neighbor in each old arm,
so there are at most two visits, with coefficients (67,34) and
(69,32). The boundary at r=2 is therefore 135 vs 133,
where 133=7*19 falls on a critical 7p lane. More generally the
69 arm in this family is divisible by 7 iff r=2 modulo 3 and
the other arm never is. This follows from the order-three
cycle of 2 modulo 7; no infinitude of prime cofactors is asserted.

## Ternary phase of the original two-scale pair

Without inflating Z, the full two-module tensor of the 67 arm has
rows (67,81,24), (32,86,42), (0,22,22); the 69 arm differs by +2
only in (1,1). Thus for arbitrary X,Y, the difference is 2h_X h_Y.
Modulo 3 all C_3(U) vanish because the ordered-cover count is divisible
by 3!. Also c_U ≡ ((−1)^|U|+1)h_U modulo 3: the unrestricted
two-bin convolution counts each odd-cycle packing with a factor
2 per independently assigned odd block. The number of those
blocks has parity |U|, so the convolution is (−1)^|U| times
the ordinary path count mod 3. The two empty-part choices contribute
2h_U and are removed to leave c_U.

Reducing the actual two-module tensor therefore gives, when 3 does
not divide h_X h_Y:

- The 67 arm is h_X h_Y mod 3, except at |X| even and |Y| odd,
  when it is 2h_X h_Y mod 3. In particular it never has a factor 3.
- The 69 arm is divisible by 3 in all the nonexceptional parities,
  and is not divisible by 3 when |X| is even and |Y| odd.

If 3 divides either source count, both arms are divisible by 3.
Thus the base prime versus 3*23 split persists across three
of the four transitive module-size parities, but not the fourth.
It is a local arithmetic distinction, not a supply of all primes.

The literal B core is self-converse by reversing the order of its
seven names, whereas A has no anti-isomorphism. This explains the
alias of B's two marked positions (0,1) with (6,5) for transitive
modules under global path reversal. The two base cores are prime
with respect to homogeneous substitution. These two finite graph
facts are checked by pair_structural.py; they are not claims that
minimum-order witnesses at their whole Hamiltonian fibers are unique.
