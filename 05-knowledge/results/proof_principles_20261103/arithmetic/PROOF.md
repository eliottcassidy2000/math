# Exact valuation interfaces and finite compilation of a fixed prime support

**Status: PROVED mathematical statements A1–A4 below; FINITE-EXACT controls in
`receipt.json`. Independent review is recorded separately by the parent packet.**
These are scoped extensions of Session177, with no claim of historical
novelty, a new record bound, or odd-perfect-number nonexistence.

## Inheritance and the question

Closest proved mechanism: [Session177's LCM neighborhood closure](../../session177_odd_perfect_20261104/oddperfect/lcm_neighborhood_closure.md)
characterizes possible supports of a source's divisor-sum transfers. Its
[Euler-aware graph sieve](../../session177_odd_perfect_20261104/oddperfect/odd_order_support_theorem.md)
requires strong connectivity, while the
[exact budget interface](../../session177_odd_perfect_20261104/networks/README.md)
keeps incoming valuations. The present extension combines those two coordinates.

Canonical hostile: an induced graph can be realized while omitting exterior
prime factors and failing all valuation budgets. Corrected near miss: increasing
an exponent can remove an edge (`sigma(3^2)=13`, `sigma(3^4)=121`), so a lower
exponent bound does not license an exact sigma factor. Least-used sidecar:
the prime factors of one target's multiplicative order can raise another
target's transfer valuation.

Board: local support closure; exact valuation interface; incoming budgets;
fixed versus varying support; exterior divisors. Anchor = close the fixed
support's exponent escape; niche = constructive local valuation feasibility;
wildcard = which operation-relative interfaces remain valid after composition.
No artificial tournament is introduced: the relevant object is a directed
prime-transfer graph with numerical labels and an Euler role.

## A1. Complete realization criterion for one valuation row

Fix an odd prime source p and a finite set R of distinct odd primes other than
p. A non-Euler source has a positive even exponent a; an Euler-role source
has a positive exponent a=1 mod4. Put n=a+1 and b=1 or2 respectively.
The local role terminology does not require p=1 mod4; that extra restriction
is imposed when the row is used in an odd-perfect candidate.

For r in R let d_r=ord_r(p). The target is **allowed** if d_r is odd in the
non-Euler role, or v_2(d_r)<=1 in the Euler role. For d_r>1 put

    c_r = v_r(p^(d_r)-1).

For d_r=1 put c_r=0. For every admissible positive n,

    v_r(sigma(p^(n-1))) = v_r(n),                  if d_r=1;
                       = 0,                      if d_r>1 and d_r does not divide n;
                       = c_r+v_r(n),             if d_r>1 and d_r divides n.    (1)

The geometric-sum identity gives this from odd-prime lifting. For completeness,
if r divides x-1, binomial expansion shows that multiplying the positive
exponent by r raises v_r(x^m-1) by exactly one; multiplying it by an integer
coprime to r does not change that valuation. Thus
v_r(x^m-1)=v_r(x-1)+v_r(m). Apply x=p when d_r=1 and cancel p-1;
otherwise apply x=p^(d_r), noting r does not divide either d_r or p-1.
This also explains why the order-one denominator case cannot be omitted.

Let t_r be a requested **exact nonnegative valuation** for every r in R.
For every t_r>0 reject immediately if r is disallowed, or if d_r>1 and t_r<c_r.
Otherwise set

    s_r=t_r,       k_r=r^(t_r),             if d_r=1;
    s_r=t_r-c_r,   k_r=d_r*r^(t_r-c_r),     if d_r>1.

Form L=lcm(b,{k_r:t_r>0}). The requested vector is realized by an exponent
of the chosen role **if and only if** both remaining conditions hold:

1. For every t_r>0, v_r(L)=s_r.
2. For every allowed target with t_r=0, its activation requirement does not
   divide L. That requirement is r for d_r=1 and d_r otherwise.

**Proof.** If the vector is realized, every positive requirement k_r and b
divides n, so L divides n. A positive target requires v_r(n)=s_r; since its
own requirement already puts r^(s_r) into L, this forces v_r(L)=s_r.
A zero target cannot have its activation requirement divide L, since it
would then divide n. This proves necessity.

Conversely, if L>1 choose n=L. Its 2-adic valuation is zero in the non-Euler
role and exactly one in the Euler role: all included orders satisfy that
role's test, and all target primes are odd. Consequently a=n-1 has the
correct role and is positive. Equation(1) and the two tests give exactly
the prescribed row, including its zero coordinates. If L=1, this is the
empty-positive non-Euler case. Choose an odd prime ell greater than every
target prime and order, and take n=ell. No activation requirement, each at
least3, divides ell, so the entire requested row is zero. This gives a
positive even exponent and completes sufficiency. QED.

This is a **constructive iff for the named finite target set**, not a claim
that the source has no other sigma divisors. Source rows can be realized
independently when they have independent exponent variables; a demanded
common exponent adds a further joint compatibility constraint.

## A2. A conjunction forces valuation that the support quotient loses

The exact local data are

    ord_11(3)=5,     v_11(3^5-1)=2;
    ord_23(3)=11,    v_23(3^11-1)=1.

For every positive even a,

    11 | sigma(3^a) and 23 | sigma(3^a)
        imply 11^3 | sigma(3^a).                                  (2)

Indeed the antecedents force both5 and11 to divide n=a+1. The first target
is active with valuation 2+v_11(n)>=3. Neither antecedent alone gives this:
a=4 has 11-valuation2 and no23, while a=10 has23 but no11. The bound is
sharp: a=54 gives the exact (11,23)-valuation row(3,1).
The support{11,23} is therefore locally realizable, but the exact row(2,1)
is impossible, although each single positive requirement is realizable.
This is a dependency between edge capacities caused by the target orders.

**Conditional odd-perfect consequence.** If an odd perfect N has 3-component
3^a and both displayed transfers occur, then 11^4 divides N. The source
supplies at least three11s to sigma(N)=2N, and11 cannot be the Euler prime
because11=3 mod4. Its exponent in N must therefore be even and at least4.
This clause survives arbitrary changes of a that preserve the two transfers.
It does not say that 3,11,23 merely belonging to the support forces the clause.

## A3. Every specified finite prime support has a finite exponent box

Let S be a specified finite set of odd primes, s=|S|>=2, and choose q in S
with q=1 mod4 as Euler prime. The exponent of q must be1 mod4; all other
exponents are positive even. For each r in S define

    C_r = sum over p in S\{r}, with ord_r(p)>1 and allowed for p's role,
          of v_r(p^(ord_r(p))-1).
    C=max_r C_r,       h=s-1.

Choose an integer B with

    B>=max(2,2h,C+1),       3^(B-C)>(B+1)^h.                       (3)

Such a B exists, can be found by increasing B using integer arithmetic,
and depends only on the specified primes and their Euler role.

**Theorem.** Every odd perfect number with exactly support S and Euler prime
q has all its exponents strictly less than B. More generally this holds
for every such exponent vector satisfying just the internal odd-prime
budget equations

    a_r = sum_(p in S\{r}) v_r(sigma(p^(a_p)))     for all r in S. (4)

**Proof.** Write A=max_p a_p. By(1), a transfer from an allowed source is
at most its baseline c_pr plus v_r(a_p+1), including the order-one case
with baseline zero. A disallowed source contributes zero. Therefore(4)
gives

    a_r <= C_r + h log_r(A+1) <= C + h log_3(A+1).

Apply this to a prime attaining A. Equivalently, with no numerical
logarithm needed,

    3^(A-C) <= (A+1)^h.                                          (5)

For integer k>=2h the ratio between successive values of
3^(k-C)/(k+1)^h is

    3*((k+1)/(k+2))^h >= 3*(1-h/(k+2)) > 3/2,

by the elementary Bernoulli inequality. This proves existence in(3) and
strict increase after B. If A>=B, (3) would contradict(5). Hence A<B. QED.

**Exact compilation corollary.** Enumerate every exponent below B in its
Euler parity class, then test sigma(product p^a_p)=2 product p^a_p with
integers. This is a terminating, complete decision algorithm for that
specified S and q. Testing all possible q decides the specified support.
The local order/valuation matrix and LCM implications can prune this box
before evaluating large products.

The omitted support sizes are immediate: the empty support gives N=1,
with sigma(1)=1 rather than2; a singleton odd-prime support satisfies
sigma(p^a)/p^a<p/(p-1)<=3/2<2. A support with no prime1 mod4 has no eligible
Euler role and is excluded by the inherited Euler form. Thus the decision
algorithm really does cover every specified finite set of odd primes.

The theorem is **not** a bound for all supports of a fixed cardinality:
C depends on their actual primes. It provides no finite bound on possible
supports, no efficient global search, and no odd-perfect nonexistence proof.
Nor is(4) alone asserted sufficient: sigma values could still have exterior
prime factors. The final exact integer equality retains that boundary.

## A4. A graph-admissible support eliminated by the budget bound

Take S={3,5,7,11,13}. Its only possible Euler primes are5 and13. Both
Euler-aware maximal graphs are strongly connected. The abundancy ceiling
is product p/(p-1)=1001/384>2, so those two inexpensive tests do not reject it.

For q=5, the baseline vector on(3,5,7,11,13) is(1,0,2,3,1), giving C=3
and B=13. Only3 can send a permitted transfer to13: the other orders are
4,12,12 for sources5,7,11, respectively, and fail their roles. Thus

    a_13 = v_13(sigma(3^(a_3))).

But a_3+1<=13. The13-edge is active only if3 divides a_3+1, when13 does
not divide that length, so the right side equals1; otherwise it is0.
It cannot equal the positive even a_13. This excludes the q=5 branch.

For q=13, the baseline vector is(0,0,2,4,1), giving C=4 and B=14.
Only11 can send a permitted transfer to5, because the orders from3,7,13
are all4. Thus a_5=v_5(a_11+1). The ordinary even exponent a_11 is at
most12, so its length is at most13 and the right side is at most1.
Again this contradicts the positive even target exponent. These cases
exclude this entire specified support, including every higher-exponent
branch through A3, without freezing any exponent prematurely.

This small support is a mechanism control already excluded by known general
OPN restrictions; no improved record is claimed. The useful implication is
graph feasibility -> retain valuations -> finite exponent bounds -> exact
budget contradiction. Graph feasibility alone never supplied the last step.

## Verification, losses, and continuation

Run `python3 verify.py` from this directory; stdout is the receipt. The checker
uses explicit validation active under `python3 -O verify.py` as well.

- All17920 requested exact rows, with source and three distinct targets
  among odd primes below20, both roles, and valuations0..3:2543 realized
  rows are checked by binary geometric sums;15377 have typed rejection
  certificates. This is an implementation control, not an empirical proof
  of the unbounded iff.
- An independent incremental geometric-sum implementation enumerates a
  complete44100-length period for each source11,13,17,19 and targets3,5,7.
  Its176400 lengths give exactly the predicted feasible exact0/1 triples
  in both roles, testing rejection as well as construction.
- Sharp exponent54, the separate antecedent controls, the order-one
  denominator, and forbidden Euler order4 are tested explicitly.
- Every3888 candidate vectors for q=5 and5184 for q=13 in A4's proved boxes
  fails an internal valuation budget. These finite enumerations corroborate
  the shorter hand contradictions; no OPN positive control is claimed.

Map and loss ledger: source = exponent-labelled sigma rows; target = finite
valuation constraints plus role-specific divisibility requirements. Preserved
predicate = realizability at the named targets. Lost data = exterior factors
and global perfection. Restored sidecar = exact incoming budget and, finally,
the complete integer sigma equality. Cheapest decisive tests are the joint
order overflow(2) and the supportA4, both of which pass the weaker graph view.

The concept board changes accordingly: support closure is a projection of the
valuation interface; incoming budgets bound exponents only after support is
fixed; local universality remains compatible with this finite compilation;
exterior divisors remain a separate final condition. The next substantive
problem is controlling the prime-support family itself, not extending a
cutoff on exponents that is already mathematically bounded for each support.
