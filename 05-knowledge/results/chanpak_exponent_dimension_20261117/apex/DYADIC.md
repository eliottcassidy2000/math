# A fixed-order decoder for every dyadic phase, with an exact +2 mate

For every k>=2 and odd residue N modulo 2^k there is an explicitly
decodable pair from PROOF.md with lower endpoint L congruent to N, upper
endpoint exactly L+2, and both strong and on O(k^2) vertices. The order
depends only on k, not on the residue. This theorem does not give L=N
as an absolute integer.

Take K=max(k,6). Use odd a=2^K x+r with r=1 or 9; its partner parameter
is b=a+1=2^K x+r+1. We will construct x ending in K ones. Once x mod
2^K is fixed, appending a specified K-bit suffix to x adds a constant
to F(x) modulo 2^k. This is because at each zero bit the addition
q+[q odd] only uses the preceding prefix q modulo 2^k. Likewise the
polynomial 4a^2+7a modulo 2^k is fixed by the suffix. Thus changing
x with this ending changes the lower endpoint modulo 2^k by
exactly 2 times the change of F(x).

The suffixes r=1 and r=9 give L=3 and L=1 modulo 4, respectively,
for any such x. For clarity, a new zero contributes 2 mod 4 iff
the two preceding bits differ, and 0 otherwise. Each initial zero
run in the four suffixes (r=1,2,9,10) has at least two zeros and
follows two ones; it contributes 2. Thereafter the words
1, 10, 1001, 1010 have additional increments 0,2,0,0 mod 4 from
the zero-zero state. The common polynomial term 4a^2+7a is 3 mod 4
and 2F(x)=2 mod 4. This gives the stated phases.

Start x with K ones. For each j=1,...,k-2 append a block

    0^(j+epsilon_j) 1^(K+1-epsilon_j)        epsilon_j in {0,1}.

Each block has fixed length j+K+1, begins after at least K ones and
ends with at least K ones. A zero run of length l<=K following K ones
deposits 2-2^l in F modulo 2^k: its first zero contributes 0 and the
remaining zeros contribute -2,-4,...,-2^(l-1), from the current prefix
modulo 2^k. Ones contribute nothing. The blocks therefore contribute
additively and independently at this precision. Changing epsilon_j
from 0 to 1 subtracts 2^j from F(x) and hence 2^(j+1) from L.
These are exactly the binary place weights 4,8,...,2^(k-1). Choose r
by the target modulo 4, then use binary expansion of (L_0-N)/4 mod
2^(k-2), where L_0 is the endpoint with all epsilon_j=0. The whole
construction uses one fixed address width at each k, including k=2
when no blocks are needed. It gives 2^(k-1) distinct lower counts at
one order (one in each odd phase), all with their exact +2 companion.

This is a controlled instance of an exponent becoming a coordinate:
k-2 binary variables set k-2 correction digits. What it costs is
large absolute height, not just O(k^2) arc positions. Raising the
precision increases that height again. It is therefore an explicit
example of why having every residue of every modulus, even while
maintaining an exact physical +2 move, can coexist with only
Theta(sqrt X) endpoints through X.
