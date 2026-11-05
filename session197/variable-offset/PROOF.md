# A closed physical family with a changing ternary offset

## Construction and exact count

Start with the negative-clock core

`T=[2,12,145,149,99,143,175,19]`, root7,

and retain the same five old ports `[0,1,4,5,6]` throughout. First add one fixed gate B by `L(mask12,drop5)`: B loses only to old4,5 and beats every other old vertex. The resulting nine-vertex seed is

`S=[2,12,145,149,355,399,175,19,207]`.

From S allow two operations, always forgetting only the newborn exposure:

- R: `R(mask7,drop5)`, beating old0,1,4 and losing to every other old vertex;
- L: `L(mask24,drop5)`, losing to old5,6 and beating every other old vertex.

Let m be one plus the number of R births, and t the number of L births after the fixed B gate. Then m>=1,t>=0. The root vertices form one transitive module of size m. The L vertices form another transitive module of size t; every L vertex beats B and every root vertex. Opposite-side births commute, so every interleaving with the same m,t gives the same graph up to the explicit labeling below.

Its Hamilton count is

`H(m,t)=A(t)*3^m+C(t)`,

where

`A(t)=40*3^t+160*2^t-20`,

`C(t)=2^(t+1)-5`.

This is an all-order physical identity for every m>=1,t>=0. It is not a universal scalar-spectrum theorem.

At each fixed t, repeating R returns

`H_new=3H+10-2^(t+2)`.

The first offsets, at t=0,1,2,3, are `6,2,-6,-22`. A subsequent L birth changes the offset parameter; it must not be treated as preserving one fixed ternary digit.

## Exact graph equality and strongness

A convenient canonical build starts from T, performs m-1 R births, then adds B, then performs t L births. All these births retain the five original ports. An arbitrary word from S has a canonical relabeling as follows:

- Old0,...,6 are fixed.
- The canonical root list maps to original7 followed by the R newborn labels in chronological order.
- The canonical B maps to its original physical label8.
- The canonical L list, indexed in birth order, maps to the L newborn labels in chronological order.

Within the L list, later vertices beat earlier ones in both graphs. All L vertices beat all R vertices regardless of chronology. Their incidences to the old vertices are determined by the same masks. Thus every edge agrees under this mapping, not just the scalar count.

T contains the directed Hamilton cycle `(0,1,2,4,6,5,3,7,0)`. B has an inneighbor and an outneighbor in T; every nonempty L block has both as well. Equivalently, for t>=1 the graph is a nonempty transitive blowup of the strong ten-vertex core with one root and one L vertex. Strongness is preserved under these blowups. For t=0 it is a blowup of the strong seed S. Hence all graphs here are strong.

At m=t=1, the graph is exactly the negative420 core from the reverse-sign bridge, with its two added gate labels exchanged. This alternative construction keeps old6 exposed, unlike the five-port route which first retained the positive bridge vertex by forgetting old6. The two histories give the same unmarked graph but different available port exposures.

## Uniform run bounds and all-order proof

For any marked root replaced by s transitive clones, the exact run formula is

`H_s=sum_j a_j j!S(s,j)`.

Each internal fragment between root occurrences goes from the root's outneighbor side to its inneighbor side. Choosing one crossing edge in each fragment gives a matching, because fragments have disjoint vertex sets. A crossing-graph vertex cover of size d therefore bounds j by d+1.

For the root module in this construction, all crossing edges are among

`1->2,1->3,4->5,4->6,4->B`.

The fixed singleton cover `{1,4}` bounds its run index by3, regardless of t.

For the L module, its inneighbors are exactly old5,6. The only crossing edges from its outneighbors into that set come from old4 or B. Thus the fixed singleton cover `{4,B}` bounds its run index by3, regardless of m.

For positive m,t, these two uniform bounds put the count in the nine-dimensional tensor span `{1,2^m,3^m} tensor {1,2^t,3^t}`. This is justified by successive coefficient recovery with the invertible three-by-three matrices at exponents1,2,3, each with determinant12. A finite grid is used only after that all-order span is proved.

For t=1,2,3, the exact root-run vectors are respectively

`(a1,a2,a3)=(1259,1260,420)`,

`(2943,2940,980)`,

`(7031,7020,2340)`.

The verifier checks these full finite vectors and all nine direct Hamilton counts at m,t in{1,2,3}. They give zero coefficient of2^m, and the remaining two coefficients are exactly A(t),C(t) above. Tensor uniqueness proves the formula for every positive m,t.

The boundary t=0 is checked separately, not inferred by substituting an empty module into the positive-multiplicity argument. The seed S has complete positive root-run vector `(537,540,180)`. Hence its root clone count is `180*3^m-3`, precisely the displayed formula at t=0.

The initial B gate therefore multiplies the original negative-clock count by3 at every fixed root multiplicity. This is a one-generation identity. Repeating that B birth does not return a tripler: the next core has H1483 rather than3*537=1611.

## A closed four-coordinate physical interface

The changing offset can be retained explicitly. Define

`X=3^(m+t)`, `Z=3^m`, `D=2^(t+1)`.

On the actual family, the four-coordinate state `(H,X,Z,D)` has the exact updates

`L: (2H+40X+20Z+5, 3X, Z, 2D)`,

`R: (3H-2D+10, 3X, 3Z, D)`.

The initial state at S is `(537,3,3,2)`. These formulas follow directly from the proved count formula and the changes in m,t. They commute, as the underlying physical births do. They are not claimed valid on arbitrary numerical tuples or arbitrary five-port tournaments.

Thus this is a closed family carrying both graph growth directions and its variable ternary offset. The negative term in the R formula is an exact readout on physical states, not a claim of a nonnegative transfer on these compressed coordinates.

## Exact inverse test within the family

For any integer target N, enumerate t>=0 while `3*A(t)+C(t)<=N`. This minimum value is strictly increasing in t. For each t, check whether `N-C(t)` is divisible by A(t) and the quotient is a power3^m with m>=1. Every successful pair is a physical witness, obtained by the stated construction. Conversely every graph in this family has such a pair.

The search terminates: `H(1,t)>=120*3^t-3`. Its purpose is exact inversion of this sparse family, not a universal construction of all odd values. No claim is made that the representation is unique.

## Verification

Run `python3 verify.py` and `python3 -O verify.py`; their receipts agree. The checks include literal directed-path enumeration of S, all nine tensor-grid counts, full run vectors, the crossing covers, strongness, explicit interleaving graph identities, the four-coordinate updates, the inverse test, the reverse-core label swap, and the failed repeated-B tripler.
