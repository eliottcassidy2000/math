# The zero slot changes the grammar, and Fibonacci blocks retain ordered transport

23 October 2026. **PROVED, scoped:** the zero-slot fibre formula,
compatible-row family, and prime-power affine clock identities. **FINITE-EXACT:**
the declared controls. The preferred word below matches the user's current
prefix; it is not asserted to be the only infinite continuation.

## 1. Keep the empty position even when its length is ignored

Use the indexed weights

    h_0=1_A, h_1=0, h_2=1_B, h_3=1_C,
    h_k=h_(k-1)+h_(k-2), k>=2.

The initial values at2,3 agree with the recurrence. Thus the weights are
1,0,1,1,2,3,5,8,... . A representation is a finite subset of indices,
with no consecutive selected indices. Never select index1: it is an empty
placeholder. Retain its position when testing adjacency.

This is a genuine change from the older alphabet1_A,1_B,1_C,2,3,... .
Previously A and B occupied consecutive positions. Now they occupy0 and2,
so A+B=2 is legal. The zero has no geometric length but deleting its INDEX
would undo this permission. This distinguishes contraction of an empty
geometric block from erasure of its combinatorial boundary information.

Let Z(n) be the ordinary Zeckendorf support in indices3,4,5,..., using
weights1,2,3,5,... . For n>=0 the complete zero-slot fibre is the following
set of valid forms, omitting expressions with negative arguments:

    Z(n),
    {0} union Z(n-1),
    {2} union Z(n-1),       if3 is absent from Z(n-1),
    {0,2} union Z(n-2),     if3 is absent from Z(n-2).       (1)

Proof. Separate whether indices0 and2 occur. Neither gives ordinary
uniqueness. Index0 has no conflict with any ordinary index. Index2 forbids
ordinary index3 and no other ordinary index. They may both occur because
their distance is2. These four cases are exhaustive and disjoint.

There are three or four representations for every positive n. Indeed
Z(m) and Z(m+1) cannot both contain the ordinary unit: adding a unit to a
representation already containing it replaces the pair by weight2, and
subsequent upward merges never recreate weight1. Hence the last two
conditions in(1) cannot both fail, apart from the directly checked n=1
boundary. The third/fourth cases can both occur; n=4 is the first example.
Zero has one representation under the convention of never selecting its
empty atom. This extends, rather than contradicts, the old two-or-three
fibre theorem for the DIFFERENT unspaced alphabet.

## 2. The preferred coloured blocks and leftmost growth

Keep A,B,C terminal. For k>=4 use the proper nonconsecutive splitting

    k odd:  h_k=h_(k-1)+h_(k-3)+...+h_2,
    k even: h_k=h_(k-1)+h_(k-3)+...+h_3+h_0.             (2)

The sums telescope by the Fibonacci recurrence; selected indices differ
by2, except the final gap3 in the even case. They use only smaller-index
atoms and terminate recursively. Listing the pieces from largest index to
smallest and expanding gives words

    W_0=A, W_1=empty, W_2=B, W_3=C,
    W_4=CA, W_5=CAB,
    W_k=W_(k-1)W_(k-2), k>=6.                          (3)

Equation(3) follows by comparing the alternating lists in(2), or by
induction on their recursive expansions. Its adjacent-index concatenation
is a WORD identity, not permission to select adjacent atoms in(2).

The nested limit from W_3 onward begins

    CABCACABCABCACABCACABCABCACABCABCACABCAC...

It agrees with the current supplied prefix

    C,A,B,C,A,C,A,B,C,A,B,C,A,C,A,B.

For every n, expanding the ordinary Z(n) in descending piece order gives
exactly the first n letters U_n of this limit. Proof: if h_j<=n<h_(j+1),
greedy subtraction leaves less than h_(j-1); the concatenation
W_(j+1)=W_j W_(j-1) supplies precisely that next prefix. The small cases
j=3,4 are direct, and induction treats the remainder.

For the user's smallest-to-largest, left-aligned convention, reverse this
word. Then

    R_n = c_n R_(n-1),                                  (4)

where c_n is the nth letter of the displayed sequence. Thus exactly one
new unit square is introduced at the left. The expanded colours have this
property; the unexpanded atom cuts may change when carries merge blocks.
A row grid can place R_n in row n, giving n(n+1)/2 unit squares through
row n. The equality of inherited subwords makes the diagonal stripes
precise, without identifying a spatial row index with a dynamical clock.

## 3. A new compatible family repairs the historical 35th-colour obstruction

Under the old alphabet, the authoritative historical prefix had B at35,
where(3) has C. The old note proved that this isolated row35 could not
extend to row36 in its stated grammar. That proof remains correct for that
alphabet. The new zero slot changes its critical step:

    old preferred2: C+A  -> word CA;
    new alternative2: B+A -> word BA.                   (5)

Hence34+B extends to34+(B+A)=36, keeping the same34-prefix. This is a
literal admissible expansion, not a renaming of the historical35th colour.

**Theorem (compatible recolouring family).** Let S be ANY subset of the
positions at which the preferred infinite word has C. Change precisely
those occurrences of C to B. Every finite prefix of the resulting word
is the expansion of an admissible zero-slot representation, using only
proper recursive splittings into smaller nonconsecutive indices. The
choices agree on every position shared by two rows.

Proof. In the expansion tree of any W_k, k>=4, every C leaf is the final
h_3 immediately followed by h_0 in the splitting of some even-index node.
At that one node replace h_3 by h_2. Both have length1; the preceding
piece, if any, has index at least5, and the following piece has index0.
The new selection remains nonconsecutive. It changes exactly that leaf
from C to B, leaving every other leaf, its order and its length unchanged.
Different C occurrences have distinct parent decisions and can be changed
independently.

For a general ordinary Zeckendorf row, apply this rule inside every larger
atom. Its only possible remaining C is a terminal ordinary-unit atom of
the row itself. Replace that index3 by index2 if its absolute position is
in S. The preceding ordinary index, if present, is at least5, so this is
also admissible. These operations construct precisely the desired prefix
for each n. Retaining the SAME S of absolute positions proves agreement
between all rows. Reversal gives the same leftmost-growth law(4). QED.

Thus S={35} yields an infinite compatible continuation of the older
35-symbol row, now under the updated zero-slot convention. The theorem
also gives many other continuations, even after fixing any finite supplied
prefix. The construction treats equal-valued atoms as typed occurrences;
it does not force every occurrence of the numerical weight2 to choose the
same expansion. Imposing that additional context-free rule would be a
stronger, different requirement.

This is a useful repair rather than a claim of unique extrapolation. The
historical obstruction located exactly the missing adjacency permission.
The new parameter S is the explicit sidecar carrying shared choices.
There is no fourth neutral colour here: the neutral state belongs to the
separate XOR-charge convention, not this word grammar.

### The zero is a precise hidden charge, not a fourth visible colour

The formal two-register lift makes the missing information explicit. Use

    V_0=(1,0)=A, V_1=(0,1)=z,
    V_k=V_(k-1)+V_(k-2).

Then V_2=(1,1)=B and V_3=(1,2)=C. The first coordinate is precisely h_k;
the second is F_k, with F_0=0,F_1=1. The empty geometric slot has nonzero
charge z, the generator of the kernel of the length projection. Requiring
that zero slot to have neutral charge would identify A=B=C and lose the
three typed units.

Every preferred proper split(2) satisfies the lifted identity

    V_k = sum_(i in proper(k)) V_i + z.                 (5a)

Telescoping the recurrence leaves V_3-V_2=z in the odd case and
V_2-V_0=z in the even case. Replacing its terminal C by B loses one more
z. Consequently a fully expanded row obeys exactly

    original second register
       = weighted leaf colours + number of expansion nodes
                                + number of C-to-B replacements,       (5b)

where leaf weights are A:0, B:1, C:2. The root row is the ORIGINAL ordinary
Zeckendorf support; a replaced terminal unit is included in the final
replacement count. Every recursive node decreases its atom index, so this
is a finite source-preserving accounting identity, not a dynamical
termination assumption.

For k>=4 the number of preferred internal nodes is F_(k-2), proved by(2)
or induction. If the ordinary support of n is Z(n), its second register is

    Y(n)=sum_(i in Z(n)) F_i = 2n-floor((n+1)/phi^2).      (5c)

This is exactly the shifted-value register recovered in the earlier
[Zeckendorf guard reader](../zeckendorf_guard_automaton_20261003.md), with
our indices shifted to include the empty slot. Formula(5b) gives it a new
interpretation through the user's actual coloured decomposition tree.
It is checked against the inherited integer-square-root floor formula,
not floating-point approximations.

The shared mechanism with the chart problem is specific: a projection
keeps the visible object while omitting a nontrivial kernel coordinate.
Here that coordinate is the zero-length refinement debt; for chart gluing
it is the simultaneous attachment choice. Their algebraic groups differ,
and no map making this refinement debt a Collatz descent bound is supplied.

## 4. Ordered block transport: a useful bridge to odd cells

Suppose each letter is assigned an actual bijection of one finite state
set. Reading a word w_1...w_m means applying its letters in that order;
write P(w)=p_(w_m)...p_(w_1), with composition read right to left. Then

    P(uv)=P(v)P(u),
    (length(uv),P(uv))=(length(u)+length(v),P(v)P(u)).      (6)

A block summary with its source/target state sets and exact integer length
therefore composes without expanding the word. Partial bijections work
as well if the domain of the composite records every intermediate guard.
The main packet's private-path theorem gives the boundary-section proof.
If an internal state is shared with another attachment, retain that state
or its joint boundary relation before contracting.

For the preferred blocks(3), (6) gives

    P_k=P_(k-2)P_(k-1),
    (P_(k-1),P_k) -> (P_k,P_(k-1)P_k).                  (7)

The group-pair update is invertible: from(y,z) recover(x,y)=(z*y^-1,y).
Every orbit is therefore periodic from its start. The integer lengths
still grow by the Fibonacci recurrence. A periodic phase is not a
periodic word or clock.

The order is indispensable. Give A the rotation x->x+1 and C the
reflection x->-x on F3. Then CA and AC have equal colour counts and length,
but products x->1-x and x->-1-x. Commutative counts cannot replace the
shared chart transport. Conversely, a zero-length macro may be discarded
only when its retained transport is the identity between the identified
ports: zero LENGTH alone does not imply that equation.

These observations do not assert that the tournament's actual attachment
phases follow the Fibonacci word. They give a correct compressed verifier
for paths whose construction really does follow a stated block grammar.

## 5. The exact 3/7 phase clocks, with the integer clock retained

For S3, the affine permutation group on three points, the36 pairs in(7)
have exact periods

| Period | Number of pairs |
|---:|---:|
|1|1|
|3|9|
|6|18|
|8|8|

The universal period is24. Here is an elementary explanation. The identity
pair is fixed. A nonzero pair in one of the three C2 subgroups follows the
Fibonacci matrix modulo2 and has period3. The eight nonzero pairs in the
unique C3 subgroup follow that matrix modulo3 and have period8. The other
18 pairs generate S3; substitution using a^2=1, r^3=1, ara=r^-1 gives
period6. All36 pairs are also checked directly; all216 assignments of
A,B,C to S3 are compared against literal word evaluation.

The actual internal automorphisms of P7 are

    H={x->a*x+b : a in{1,2,4}, b in F7}, |H|=21.

They are relevant to the six large-symmetry order31 cases. On H^2 the
441 pairs in(7) have exact periods

| Period | Number of pairs |
|---:|---:|
|1|1|
|8|56|
|16|48|
|56|336|

The universal period is112. This 56=7*8 has an exact mechanism, not just
an occurrence of the same two integers. Write group multiplication as

    (a,b)(c,d)=(ac,b+ad) over F7.                        (8)

The two slope exponents lie in F3^2. If both vanish, the translation pair
follows the Fibonacci matrix over F7. Its eighth power is -I and its
sixteenth is I; its48 nonzero vectors have exact period16.

Otherwise the nonzero slope pair has exact period8. After eight steps the
translation pair is acted on by I+N, where N is a nonzero rank-one2x2
matrix over F7 with N^2=0. The eight matrices, indexed by the initial
slopes(a,c), are:

|(a,c)|I+N|
|---|---|
|(1,2)|[[1,0],[6,1]]|
|(1,4)|[[1,0],[5,1]]|
|(2,1)|[[1,6],[0,1]]|
|(2,2)|[[2,6],[1,0]]|
|(2,4)|[[5,1],[5,4]]|
|(4,1)|[[1,5],[0,1]]|
|(4,2)|[[5,2],[6,4]]|
|(4,4)|[[3,5],[2,6]]|

Each matrix follows by eight substitutions in(8). Its seven fixed vectors
have period8 in the full pair state. Its other42 vectors have period56,
since(I+N)^j=I+jN and the characteristic is7. Thus8*7=56 pairs have
period8 and8*42=336 have period56. Together with1+48 this exhausts441.
The direct group enumeration and the eight-step matrix derivation are
separate checks of the same result.

This identifies a precise bridge: Fibonacci word composition can carry
actual permutation memory, whose slope and translation return on different
clocks. It does not identify every seven-cell block relation with an H
translation, or equate the group-pair recurrence with Collatz. Source,
ports, domains, word order, and the unbounded integer length remain fields
of the certificate.

## 6. A prime-power lift: the eight-step return retains arbitrarily deep digits

The56-step example extends to a proved family. Let p be a prime other
than3, k>=1, R=Z/(p^k), and choose omega in R satisfying
omega^2+omega+1=0 with omega!=1 modulo p. Let H_R consist of the affine
maps x->a*x+b with a in{1,omega,omega^2}, b in R. Such a root modulo p
lifts uniquely at every depth because2omega+1 is a unit; the digit lift
solves the quadratic equation one power of p at a time.

Apply the same group-pair update (x,y)->(y,xy). For every nonidentity
slope pair, the slopes return after exactly eight updates. If v=(b,d)
is the translation pair, that return acts as

    v -> (I+3Q)v,       Q^2=0.                           (8a)

There are only four formulas, with r either omega or omega^2:

| Initial slopes | Q |
|---|---|
|(1,r)|[[0,0],[r,0]]|
|(r,1)|[[0,r],[0,0]]|
|(r,r)|r*[[-1,1],[-1,1]]|
|(r,r^2)|[[-1,-r],[r^2,1]]|

Eight substitutions of the affine law(8) prove this table already over
Z[omega]/(omega^2+omega+1), before any modular reduction. Every Q has
square zero, determinant zero and a unit entry. Since3 is a unit, N=3Q
has one unit Smith factor and one zero factor over R. In particular its
image and kernel each have p^k elements.

Write t=min(k,v_p((Nv)_1),v_p((Nv)_2)), with v_p(0) capped at k.
Then the exact return period of the FULL affine pair is

    8 p^(k-t).                                          (8b)

Indeed(I+N)^j v=v+jNv. The least positive j killing Nv is p^(k-t),
and the slope pair forbids a return before a multiple of eight. Thus
retaining slopes alone loses every deeper translation-valuation layer.
For all eight nonzero slope pairs together, the exact period counts are

    period8:       8p^k pairs;
    period8p^r:    8p^k(p^r-p^(r-1)) pairs, 1<=r<=k.     (8c)

Each image element has p^k preimages; the cyclic image has exactly
p^r-p^(r-1) elements of additive order p^r. This proves the count and
exhausts8p^(2k) pairs. The all-identity slope pair is a different branch:
its translations follow the ordinary Fibonacci matrix over R.

For p=7 the order of that matrix is16*7^(k-1). At depth one its eighth
power is -I. Over the integers its sixteenth power is

    [[610,987],[987,1597]]
      = I + 7*[[87,141],[141,228]].

The correction is nonzero modulo7. The binomial identity lifts its
valuation by exactly one on each seventh power, giving the stated
matrix order at every depth. Hence the exact order of the pair update
on the ENTIRE group H_(Z/7^k)^2 is

    lcm(16*7^(k-1),8*7^k)=16*7^k.                       (8d)

It is112,784,5488,... at successive depths. This is a genuine unbounded
memory tower within the explicitly defined affine recurrence. It is not
a period of the Fibonacci word, a tournament Hamiltonian path, or a
Collatz orbit. The prime3 is excluded for a reason: the factor3 in(8a)
would no longer be a unit. The literal P7 automorphism interpretation is
used at p=7,k=1 only; other rings here are algebraic controls unless a
separate graph realization is supplied.

The checker derives all eight matrices over the integer quadratic ring,
then tests every nontrivial-slope pair modulo7,13,19 and49 by direct
iteration until its first return. It also tests512 stated samples
modulo343 and the Fibonacci matrix order modulo7^k for1<=k<=12.
These corroborate the general proof, not substitute for it.

## Reproduction and scope

`check_fibonacci.py` verifies all zero-slot representation fibres through
2000 by independent subset enumeration and formula(1). It checks all
preferred prefixes, three full recolouring patterns through2000 (including
S={35}), and all256 recolour subsets on the first21 positions. It also
checks every group pair and the direct/recursive ordered-product identity
for all216 S3 letter assignments.

`check_clock_lifts.py` independently computes the eight-step matrices in
the integer Eisenstein ring and replays the prime-power controls above.
Both checkers compare frozen JSON receipts in normal and optimized Python;
their failure checks do not depend on Python assertions.

All statements concern the explicit grammar above. Other conventions for
selecting the zero atom, globally sharing the expansion of equal-valued
atoms, or indexing the triangular display must retain their own rules.
