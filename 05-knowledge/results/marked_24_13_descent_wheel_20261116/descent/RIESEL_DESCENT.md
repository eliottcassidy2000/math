# A Riesel wheel whose divisors are actual tournament counts

This transfers a finite Riesel covering to **exact but conditional descent**
for Chan–Pak on two shifted dyadic rays. The missing row in period 24 is real:
none of its one-sign covers can avoid **7**, which is not a tournament count.
Periods 72 and 120 replace that row using composite odd parts. The labels
below are all verified by literal tournaments. Compare the finite-cover
certificate in [Cowles–Gamboa’s ACL2 paper](https://arxiv.org/pdf/1110.4671).

## From an exponent class to a physical induction step

An ordinal sum (all cross-arcs from the first tournament to the second)
multiplies Hamiltonian-path counts. Let M be a putative **least** missing
positive odd count other than 7 and 21. If M=pQ with p a certified count,
1<p<M, and Q neither 7 nor 21, then Q has a tournament witness by
minimality. Join the witness for p to the witness for Q. This contradicts
the choice of M. Equivalently a witness for Q constructively lifts to M;
the wheel itself does not construct Q.

For n≥1 and fixed positive odd K, consider N=K·2^n−1. If p is odd and
d=ord_p(2), selecting phase r mod d and imposing K≡2^(−r) mod p ensures
p divides N for every n in that class. CRT for distinct p turns any set
of covering exponent classes into actual K. Choose K with K−1>21 max p
so that every quotient is safe. The result is an all-order implication
on the **absolute targets** K·2^n−1 (also for K+2j∏p, j≥0), not a
statement that some other count has the desired residue.

## Periods 72 and 120

For period 72 take **K=66,837,283**. The triples are (p,d,r):

    (3,2,0) (5,4,1) (17,8,3)
    (13,12,7) (241,24,15) (73,9,5) (37,36,35) (19,18,11).

The first three labels cover every residue mod 8 except 7. On the residual
fiber, 13 and 241 cover exponent residues 1 and 0 mod 3. Its three
remaining points mod 9 are 2,5,8, supplied by 19,73,37.

For period 120 take **K=57,013,502,549**:

    (3,2,1) (5,4,2) (17,8,0)
    (13,12,8) (241,24,12) (31,5,1) (11,10,2) (41,20,8)
    (61,60,40) (151,15,4).

Here the dyadic fiber is n≡4 mod 8. On it the CRT coordinates mod 3 and
mod 5 produce a 3×5 grid. The two row labels 13 and 241 select 2 and 0;
the three column labels 31,11,41 select 1,2,3; and the last two labels
cover (1,0) and (1,4). All these phases also match the appropriate
dyadic coordinate of their orders. In particular the period-120
two-sign obstruction does not obstruct this one-sign wheel.

The label witnesses below are short spine words with the convention
2a_i+b_i for the two out-contact bits, followed by the mutual direction
AB or BA. Independent subset-DP verification is in `verify_descent.py`.

| count | word | mutual arc |
|---:|:---|:---|
|3|1|AB|
|5|10|AB|
|11|110|AB|
|13|120|BA|
|17|1000|AB|
|19|1010|AB|
|31|1201|AB|
|37|1200|AB|
|41|3100|AB|
|61|100010|AB|
|73|11200|BA|
|151|33100|AB|
|241|112210|AB|

No cover at period 24 can avoid prime 7: the dyadic-only labels leave a
fiber mod 8, and on it 13 and 241 cover at most two of the three odd
coordinates. The same slice principle proves **241 is the smallest
maximum prime label for a 7-free wheel of ambient period 120**. Below it,
the only available non-dyadic labels are one row (13), three columns
(11,31,41), and two singletons (61,151). Even optimizing phases leaves
(3−1)(5−3)=4 isolated cells, more than two.

**Period 72 is shortest among all 7-free wheels whose labels are ≤241.**
The verifier enumerates these small primes and computes their orders. Only
48 and 60 below 72 pass both necessary density tests: total full density
at least one and total odd-projected density at least one. Dyadic density
is always <1 in this restricted palette, so the latter necessity follows
by slicing a single uncovered dyadic fiber. For period 48 there are at
least two uncovered fibers mod 16, and only three non-dyadic labels
(13,241,97) are available. All three are needed on each fiber, but label
97 at order 48 reaches only one. For period 60, the same grid obstruction
as for the 120 wheel truncated below 241 applies.

## Why this does not prove the conjecture

A Riesel wheel varies **the exponent on a fixed ray**, not every
positive odd target. Even allowing all K in either of the CRT congruence
classes is not universal: if N≡1 modulo the product of all these odd
prime labels, no divisor label applies. The exact gap flip from
[the programmable seam](../../unit_shields_wheel_followup_20261116/apex/PROGRAM.md)
is also not an inverse row on its own: its starting count is not
assigned. Thus these two results cannot be combined by equating their
input parameters without a compatibility theorem.

There is a stronger reduction for any ordinal-join-only argument.
After 63 and 343 are verified, it suffices to prove Chan–Pak just for
odd primes and for numbers 7p with p an odd prime (apart from 7 and 21
themselves). To see this, suppose N is composite. If it has a prime
divisor p≠7, try p·(N/p). If the quotient is 7 the target is already of
type 7p. If it is 21, write N=3·(7p), which is safe unless p=3,
giving 63. If all prime factors are 7, a power 7^a with a≥4 splits
into 49·7^(a−2); the remaining cases are 7,49,343. Strong induction
then proves the stated reduction. Both isolated counts are already
in the finite census. A finite pure divisibility wheel cannot solve
the prime part; it is precisely where new affine rows become necessary.

A prospective true **descent wheel** would provide finitely many physical
rows with explicit guards: for each residue class it covers, a witness
for a strictly smaller odd q is computably sent to a witness for the
**given** N, with every bad quotient tested against a finite boundary.
An arithmetic exact-cover search may choose the rows, but cannot
certify their tournament transfer laws. Moduli 24 and 120 give natural
CRT carriers 8×3 and 8×3×5. To see a structural limit, if all candidate
row moduli are powers of individual primes, a CRT cover is possible
only if one prime's rows alone already cover that prime's whole
coordinate: otherwise choose an uncovered coordinate for each prime
and combine by CRT. In particular adding many odd primes to the
same 0 or +2 divisibility rules misses N≡1 at every prime
coordinate. The composite factorization of the **orders** in the Riesel
cover is what avoids that purely prime-power obstruction.
