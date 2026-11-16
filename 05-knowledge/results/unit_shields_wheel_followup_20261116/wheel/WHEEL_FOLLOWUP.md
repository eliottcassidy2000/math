# Zsigmondy wheels: slice necessity, prime levels, and the two-sheet test

Source definitions: [THM-5254 packet](../../two_three_wheels_synthesis_20261116/WHEELS.md)
and its HYP-9823. A clock label is a distinct prime p dividing 2^N−1;
its class has modulus d=ord_p(2). Multiplicity in the factorization of
2^N−1 is not a second clock. A wheel is an ambient period: any set with
least common multiple dividing N can be used there.

## Exact corrected nesting criterion

Let N=2^a m, m odd. Either the purely dyadic labels cover Z/2^a, or a
cover exists if and only if the non-dyadic labels, with their orders
replaced by their odd parts, cover Z/m. The dyadic-only alternative is
independent of m.

Proof. In any full wheel assignment whose dyadic labels do not cover
all dyadic residues, take one uncovered t in Z/2^a. Restrict the cover
to the entire CRT fiber {t}×Z/m. A class mod 2^i d with d odd contributes
nothing or one class mod d on this fiber. Thus the active nondyadic labels
cover Z/m, and the others can be placed arbitrarily. Conversely, place
one dyadic label at each level 2^i, i=1,...,a, choosing nested phases
0 mod 2, 1 mod 4, 3 mod 8, and so on. This leaves the single fiber
2^a−1 mod 2^a. Use CRT to lift the odd assignment into that fiber.
The first five levels have exactly one dyadic prime each; their summed
density is 1−2^(−a), which is below one for a≤5. At level six two primes,
641 and 6700417, both have order 64 and can each take one of the last
two slots after the first five levels. This is the exact W(64) cover.
Hence HYP-9823(W1) holds after explicitly including the dyadic-only
alternative. The literal unqualified iff is false already at N=64:
there are no nondyadic labels, but the wheel covers.

For m=pq with distinct odd primes, let r,c,s count labels with odd
parts p,q,pq respectively. The odd wheel covers iff
(p−min(p,r))(q−min(q,c))≤s. Interpret the three types as whole rows,
whole columns and individual cells in the p-by-q CRT grid. Optimize
by choosing distinct rows and columns; every leftover cell requires
its own individual identity. This is a count at selected coordinates,
not a density bound on the average whole wheel.

## Exact least prime layers, and an eventual uniform answer

For prime q, one label with odd part q covers only one point of the
odd CRT fiber. Thus the least exponent is the first a≤5 for which
sum_{i=0}^a ν(2^i q)≥q, or six if none. The unconditional bound a(q)≤6
follows from W(64), so the true least exponent is automatically
O(sqrt(q)). A problem where all dyadic-only covers are excluded would
be different. If the selected classes must have lcm exactly 64q rather
than dividing it, append an order-q (possibly redundant) clock:
2^q−1 has a prime of order q for q prime.

Exact layer counts and periods:

| q | ν(q), ν(2q), … through layer 5 if needed | least exponent | least period |
|---:|:---|---:|---:|
|11|2,1,2,2,3, at least 1|5|352|
|13|1,1,3,2,1,5|5|416|
|17|1,1,3,2,3,3|6|1088|
|19|1,1,3,3,3,6|6|1216|
|23|2,1,4,1,3,3|6|1472|
|29|3,2,2,2,5,6|6|1856|
|31|1,1,4,3,2,6|6|1984|

The rows 11,13,17 have independently replayable proofs using
the exact cyclotomic value and recursive primality certificates in
[the verifier](verify_wheel.py) and [its data](prime_certs.json).
The remaining entries use the complete factorizations in the
[Cunningham Project Main Tables, 2026-11-15](https://homes.cerias.purdue.edu/~ssw/cun/pmain1126.txt)
(Tables 2− and 2+, primitive portions of 2^q−1 and 2^(2^(i−1)q)+1).
The verifier verifies orders by testing each prime divisor of 2^i q,
so a cyclotomic factor which also has smaller order is not
mistaken for a primitive label: 17 at d=136 has order 8, for example.

For every odd prime q≥77000 the exact least exponent is six. Indeed if
a≤5 worked, there would be q distinct primes of order divisible by q
inside 2^(32q)−1. These are all 1 mod 2q. Their product is at least
(2q)^q q! ≥ (2q²/e)^q > 2^(32q), using log(q!)≥∫_1^q log x dx and
e<11/4. That is impossible for a divisor of 2^(32q)−1. Thus the
unresolved classification is finite; this does not classify all of
the intervening primes from 37 to the cutoff.

## Period 120 is not self-dual, even with single-sheet odd-order labels

Let one chosen coefficient determine both sheets k·2^e−1 and
k·2^e+1. A prime of order d even has paired phases separated by d/2;
an odd order can serve at most one sheet. This model is at least
as permissive as an even-order-only self-dual wheel.

The exact orders, one per prime, are
2,3,4,5,8,10,12,15,20,24,30,40,60,60,120.
Take the order-2 label to cover one parity of each sheet, opposite
parities. If it is inactive, activating it only enlarges coverage.
Each sheet now has four live dyadic fibers, each is Z/15. A label
with dyadic part at least four can serve at most one sheet, since its
phase and opposite phase have the same parity. Orders 4 and 8 can
fully delete one parity-pair of fibers or a single fiber. A parity
pair here means a congruence class of the divided address modulo 2:
within each sheet the pairs are (0,2) and (1,3).
Orders 10 and 30 at best give every live fiber a mod-5 column and a
single cell respectively, regardless of which phases would really be
compatible. Orders 3,5,15 give their own single-sheet rows, columns or
cells; orders 12,20 can give a row or column in a parity-pair;
orders 24,40 can give one row or column on one fiber. All inconsistent
odd-phase choices are now relaxed, never ruled out.

**Short obstruction without enumeration.** Select five of the still-live
fibers after orders 4 and 8 delete at most three. The shared column
and singleton credits leave at least 11 uncovered cells per fiber,
55 in all. The mixed rows (12,24) can remove at most 12 cells and
the mixed columns (20,40) at most nine, for 21 total. The three odd
labels 3,5,15 remove at most 4,3,1 cells respectively per fiber on
their chosen sheet. If the five fibers split (3,2), they can remove
at most 24, leaving at least ten before the last five credits.
For a (4,1) split the crude odd-label bound is 32, so initially
at least two cells remain. If either order 3 or 5 uses the small
sheet, its capacity drops by 12 or 9. Otherwise their row and
column are both on the four-fiber sheet: on each of those fibers,
one of the four row credits and three column credits must overlap.
If the column is the already covered shared column, its capacity
instead drops by all three. Either way there is an extra loss of
at least four, so at least six cells remain. Thus five final credits
never suffice.

An exact relaxed check strengthens that lower bound. For each of the 262144 possible assignments of the above support
sets, compute the remaining lower bound: on a fiber f not deleted
by orders 4 and 8, let i_f count mod-3 rows and j_f count the
mod-5 columns besides the one credited to order 10. The untreated
cell count is at least
max(0,(3−i_f)(4−j_f)−1−[order 15 serves f]).
The −1 is the optimistic shared order-30 credit. The smallest
sum is 11. Each of the two order-60 labels gives a single cell
on at most two fibers, and the order-120 label gives one more.
Five credits cannot cover eleven cells. This finite exact
certificate settles the previously timed-out N=120 case, without
attempting a general no-self-duality theorem.

## 6424 as a literal period

6424=8·11·73. The grid criterion requires row capacity r=7 from
the exact q=11 factors, column capacity c≤21, and cell capacity
s≤192. At those maxima at least (11−7)(73−21)=208 individual
cells remain, too many. Thus period 6424 does not cover.

The bounds, without expensive complete factorization, are:

| odd part | primitive capacity bounds at d=m,2m,4m,8m |
|:---|:---|
|73|3,3,6,9|
|803|24,24,48,96|

The optional replay command `python verify_wheel.py --6424` tests
every candidate p<10^9 congruent to 1 mod lcm(2,d), verifies primality
for the divisors it finds, and strips their prime powers. The unseen
distinct order-d labels must be ≥10^9 and in this same progression.
A product of k of them is at least the product of the first k eligible
integers ≥10^9; compare to the remaining cyclotomic cofactor.
No assertion that the unresolved cofactor is prime is needed, and
leaving any unwanted factors in that cofactor only weakens the bound.

## What transfers and what does not

The earlier exact covers at 24 and 36, and least prime layers 80 and
112, can now be read in the same coordinates. Period 24 at q=3 first
works at a=3, with the missing order 6 (63) excluding period 12. At
36, the odd coordinate is Z/9 rather than a prime: classes mod 3 are
full triples and classes mod 9 individual cells. W(36) has order-3
label 7, and enough order-9/18/36 clocks once the mod-3 labels
(including order 12) have removed two triples. At 6424, whole rows
and columns leave a rectangle of genuine isolated obligations. In the
tournament problem, slicing at an exposed interface is an idea worth
testing: it works here because CRT rows and columns have independently
assignable incidences and the slice is an exact finite fiber. The
physical reader imposes additional joint endpoint conditions, and no
assignment theorem for its fibers follows from this analogy.
