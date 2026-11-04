# Exact weighted closure of prime-factor flow

Session 188. **Proved local realization and obstruction theorem.**

For one prime-power source, prescribed incoming requirements on its outgoing
valuations can be completed exactly by one LCM calculation. This extends S177's
edge-set closure to a closure on integer-valued edges. It also gives an exact
lower/upper-box feasibility test, a genuinely nondistributive arithmetic lattice,
and short valuation-overflow certificates.

This is a local theorem, not a proof about all odd perfect numbers. The rank of
apparition and LTE ingredients are classical. Historical novelty and a runtime
advantage over existing factor-chain software are not claimed.

## 1. Setup and the prime-power rank

Fix an odd prime p and a finite set R of other odd primes. Write

    U_p(n) = 1 + p + ... + p^(n-1) = sigma(p^(n-1)).

There are two source roles:

* ordinary: n is odd, and physically n >= 3;
* Euler: p = 1 mod 4 and n = 2 mod 4.

Put b=1 or b=2 respectively. Temporarily include n=1 in the ordinary role;
the zero-exponent corner is repaired explicitly below. Let Phi(n) be the vector
of v_r(U_p(n)) for r in R.

For d_r=ord_r(p), a target is role-allowed when d_r is odd in the ordinary
role, or v_2(d_r)<=1 in the Euler role. A disallowed target always has valuation
zero. A positive lower demand there is immediately impossible.

For every allowed target and integer k>=1 define

    rho_r(k) = r^k,                              if d_r=1;
    rho_r(k) = d_r * r^max(0,k-c_r),              if d_r>1,

where c_r=v_r(p^d_r-1) in the second case. Set rho_r(0)=1.

**Rank lemma.** For every positive integer n, not just role-compatible n,

    v_r(U_p(n)) >= k  if and only if  rho_r(k) divides n.       (1)

Proof. If d_r=1, odd-prime LTE and cancellation of p-1 give
v_r(U_p(n))=v_r(n). If d_r>1, the valuation is zero unless d_r divides n;
when it does, it is c_r+v_r(n), since r does not divide d_r. These are exactly
the divisibility requirements in (1). No general factorization of U_p(n) is
needed. QED.

The order-one clause is essential. Also c_r need not be one: for p=3, r=11,
d_r=5 and c_r=2. Thus valuation exactly one at 11 is impossible at every length.

## 2. The least-completion theorem

Let u=(u_r) be nonnegative integer lower demands, with u_r=0 at all disallowed
targets. Define

    L(u) = lcm(b, {rho_r(u_r): r in R}),
    C(u) = Phi(L(u)).

Every rho used here has 2-adic valuation at most that of b, so L(u) has the
correct role parity. For every role-compatible length n,

    u <= Phi(n)  if and only if  L(u) divides n.                (2)

This follows from (1), coordinate by coordinate. Consequently:

1. **Least length.** L(u) is the least positive role-compatible length satisfying
   all lower demands, both in divisibility order and in ordinary numerical order.
2. **Least valuation completion.** C(u) is the coordinatewise least realizable
   valuation vector above u. Indeed L(u)|n implies U_p(L(u))|U_p(n).
3. **Closure.** C is extensive, monotone, and idempotent.
4. **Exact realization.** A vector w is realizable on the finite target set R
   if and only if it is zero at disallowed targets and C(w)=w.

For idempotence, (2) applied to C(u)=Phi(L(u)) says L(C(u)) divides L(u).
Extensivity u<=C(u) gives the reverse divisibility, hence equality. A realizable
w=Phi(n) is fixed because C(w)<=w and extensivity gives C(w)>=w. Conversely a
fixed vector is realized at L(w), except for the physical n=1 corner below.

Equation (2) is an adjunction between valuation demands in coordinate order and
role-compatible lengths in divisibility order. The closure is its concrete,
arithmetic consequence, not an assumption about independent edges.

### The ordinary zero corner and unbounded local realizations

If L(u)=1 then all lower demands are zero and C(u)=0. Choose an odd prime ell
larger than every target prime and every target order. No nontrivial target
rank divides ell, so Phi(ell)=0, with ell>=3. Thus the fixed-point realization
criterion remains exact for positive ordinary exponents.

More generally, if w is fixed, Phi(L(w)*ell)=w for every sufficiently large
odd prime ell. No target r gains a valuation of the length, and no missing
order d_r can newly divide the length, because ell exceeds r and d_r.
Therefore **every realizable finite valuation profile occurs at unboundedly
large source exponents**. Local weighted closure alone cannot bound an exponent.

## 3. Exact box feasibility and the smallest physical witness

Suppose 0<=u<=v coordinatewise, with finite integer upper bounds v_r. The
valuation box

    u_r <= v_r(U_p(n)) <= v_r  for every r in R                 (3)

is feasible if and only if there is no positive disallowed lower demand and
C(u)<=v. When L(u)>1, its smallest physical witness is exactly n=L(u). Hence
adding a length ceiling n<=B simply requires L(u)<=B.

Necessity follows because C(u) is the least completion. Sufficiency follows
from its realization. The only subtle case is the ordinary L(u)=1 corner.
Here the smallest physical witness is the smallest odd prime ell not in the
following finite forbidden set:

* include r when d_r=1 and v_r=0;
* include d_r when d_r is an odd prime and c_r>v_r.

Why this is exact: at a prime length ell, an order-one target contributes one
exactly when r=ell; an order-greater-than-one target contributes only when
d_r=ell, in which case the contribution is c_r because r>ell. Further, if any
odd n>1 satisfies the box, a prime divisor ell of n also satisfies it, since
U_p(ell) divides U_p(n), and ell<=n. Thus a composite length cannot be a smaller
witness than every admissible prime length. The finite forbidden set also
proves existence without an exponent search. With n<=B, compare this smallest
prime witness to B.

Upper bounds of infinity can be handled by omitting that coordinate's upper
test; the same proof applies. This note's checker uses finite bounds.

## 4. The realizable profiles form a nondistributive lattice

For fixed p, role, and R, fixed points of C have the operations

    x meet y = coordinatewise min(x,y),
    x join y = C(coordinatewise max(x,y)).

The meet is fixed: monotonicity gives C(min(x,y))<=min(x,y), and extensivity
gives equality. The join is the least fixed upper bound by the least-completion
theorem. These statements assert a lattice with finite joins and meets, not
a bounded or complete lattice of all finite-valued profiles.

**The coordinatewise maximum can fail, even at tiny primes.** For ordinary
source p=3 and R=(11,23), the orders are 5 and 11, with base valuations 2 and 1.
The profiles

    Phi(5)  = (2,0),
    Phi(11) = (0,1)

are separately realizable. Their coordinatewise maximum (2,1) is not:
L(2,1)=55 and Phi(55)=(3,1). The least joint realization raises the valuation
at 11. This retains a weighted implication lost by Boolean edge closure.

**An explicit M3 sublattice.** Take ordinary source p=3 and targets

    R=(4561,368089,71).

All three targets are prime. Their orders are respectively (15,21,35), and
their base valuations are all one. The five profiles

    000, 100, 010, 001, 111

are realized at lengths 1 (or any suitable fresh prime), 15, 21, 35, and 105.
The LCM of any two of 15,21,35 is 105. Thus any two different atoms have join
111 and meet 000. These five elements form the nondistributive diamond M3.
For example, with atoms x,y,z,

    x meet (y join z) = x,
    (x meet y) join (x meet z) = 000.

All relevant primality, order, and valuation claims have short deterministic
integer checks. This is a statement about actual arithmetic profiles, not an
abstract catalogue selected without an arithmetic realization.

## 5. Odd-perfect-number consequences and the precise boundary

In a hypothetical odd perfect number N with Euler prime q, write a_r=v_r(N).
For each source p, its outgoing valuation vector must be a fixed point of the
appropriate C_p. Lower demands u^(p) therefore force at least C_p(u^(p))_r
units into each target r. Since the exact column balance is

    a_r = sum_{p != r} v_r(sigma(p^a_p)),

any branch with a valid upper bound A_r is impossible when

    sum_{p != r} C_p(u^(p))_r > A_r.                           (4)

Equality saturates that column: any additional positive contribution rejects
the branch. This is a sound weighted closure/overflow certificate before full
sigma construction or exterior-factor checking.

Two exact conditional examples:

* If 5 is Euler and 19 divides sigma(5^a_5), then 9|a_5+1 and 2|a_5+1.
  Since ord_3(5)=2 with base valuation one, closure at length 18 gives
  v_3(sigma(5^a_5))>=3. Thus this edge is impossible in a branch with a_3=2.
* If an ordinary source 3 sends edges to both 11 and 23, its contribution at
  11 is at least three. Hence such a branch requires a_11>=4 when 11 is an
  ordinary target. In particular a_11=2 is impossible. Each required edge
  separately has a witness compatible with the upper bound of two at 11.

These are conditional factor-chain implications, not new general OPN bounds.

**Local closure and global balance are distinct.** On support {3,13}, declare
3 ordinary and 13 Euler. The weighted matrix w_3,13=1, w_13,3=2 is locally
realizable: use lengths 3 and 18, hence exponents 2 and 17. The graph is strongly
connected. Its column sums prescribe exponents (a_3,a_13)=(2,1), which have the
correct Euler roles. But the prescribed Euler length is then 2 and
v_3(sigma(13))=0, not two. Local realizability, strong connectivity, and correct
parity of the column sums therefore do not imply compatible exponent selection.
There can also be prime factors outside the tested target set.

S179's primitive-divisor domains and S184's exact closed-row selection still
apply. Our theorem does not replace them or improve on their logical
completeness. It supplies an exact, compact local constraint stage and a
structural explanation of why independently chosen valuation intervals can be
inconsistent. Neither the set of support primes nor their number is bounded
here, so the global odd-perfect-number problem remains open.

## 6. Provenance and reproduction

Internal prerequisites were recovered from targeted S177, S178, S179, and S184
packets. These current GitHub files were specifically checked on 2026-11-04:

* [S178 valuation law and capacity propagation](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session178_prime_flow_honeycomb_20261104/primeflow/PROOF.md), blob 6084921f99113dfabca7e146d4037dbdcdccddbb.
* [S184 complete modular fixed-support certificates](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session184_proof_frontiers_20261104/oddperfect/PROOF.md), blob c7e032b282e721911b89b4b433d8b82f32e1e363.

Classical context: [Ray, Irmak, and Patel, *The rank of apparition of powers of
Lucas sequence*, Turkish Journal of Mathematics 42 (2018), 1566–1570](https://journals.tubitak.gov.tr/math/vol42/iss4/3/)
studies rank-of-apparition divisibility for Lucas sequences. Our U_p(n) obeys
U_p(n+2)=(p+1)U_p(n+1)-p U_p(n), with U_p(0)=0 and U_p(1)=1. The rank formula used
here is proved directly by LTE; the paper is context, not an unverified theorem
import. No exhaustive novelty survey has been performed.

Run `python verify.py` and `python -O verify.py`. The checker uses explicit
exceptions rather than disableable assertions. It compares ranks with direct
modular recurrences, closures with an independent geometric-sum calculation,
and the full interval algorithm with exhaustive bounded length search. It also
checks the weighted join, M3 witness, overflow examples, zero-length corner,
and local/global counterexample. `normal.json` and `optimized.json` contain the
deterministic receipts. Finite checks verify the implementation and examples;
the proofs above establish the unrestricted statements.

Both modes pass with identical output: 190 source/role/target configurations,
108 direct rank recurrences, 14,710 direct profiles, 1,198 closures, 76,230
adjunction comparisons, 41,040 bounded box comparisons, 3,800 meet checks,
2,532 monotonicity covers, and 1,198 fresh-prime extension checks.
