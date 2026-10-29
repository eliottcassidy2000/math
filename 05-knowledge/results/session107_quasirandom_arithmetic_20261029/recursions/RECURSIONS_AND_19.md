# The Sylvester cycle at 19: an exact prime-power lift, not a numerical rhyme

**Status: PROVED elementary results; FINITE-EXACT computational certificates.**
No priority claim. No assertion of infinitely many Sylvester primes, global
squarefreeness, or random distribution of their prime divisors.

## 1. Recovered provenance and correction gate

The historical sources recovered locally are:

* `05-knowledge/hypotheses/HYP-3724-phi6-is-the-sylvester-map-egyptian-fraction-tower.md`
  (klein S31, June29): the map and Egyptian-fraction identities.
* `05-knowledge/hypotheses/HYP-3730-the-cheeky-connection-is-heegner-class-number-one.md`
  (klein S33): **MISTAKE-087 explicitly retracts the covering-minimum interpretation**.
  At n=7 the competing covering gives 2/13<7/43. Arithmetic about this particular
  construction survives; extremality does not. Its Heegner list is not used here.
* `05-knowledge/results/number_tournament_bridges_20261014.md`, §4: records
  1807=13·139 and 139=3^7−2^11 (the negative −17 Collatz-cycle denominator),
  explicitly under NUMEROLOGY. No denominator-to-Sylvester dynamics map was found.
* `05-knowledge/results/giuga_korselt_defects_20261014f.md`, §one-prime extensions:
  the paired failure 1805=5·19² and 1807=13·139.
* Incoming `atomic_tail_reset_20261029/PRIME_DEPENDENCY.md` at commit
  `201fa96604d51171909318cc482220e00df54c38`: corrected unique fixed prime bank
  {2,3,7,43}, coordinate convergence, and the explicitly unproved global
  finite-time stabilization claim. The old local March reflection still says
  {2,3} is fixed; that is false because 7−1 divides6. Incoming correction wins.

Primary literature cross-check: A. Makki Naciri,
[A note on Sylvester's sequence, Integers24 (2024), A117](https://math.colgate.edu/~integers/y117/y117.pdf),
introduction and §2, supplies the recurrence, coprimality, and congruence context.
All claims needed below are proved directly. The paper's deeper perfect-power
claims are not dependencies. The exact mod19 lift is derived here without a
literature-novelty claim.

## 2. The recurrence and exact identity of every requested small number

Put f(x)=x²−x+1, s_0=2, s_(n+1)=f(s_n). Then

    s_0,...,s_6 = 2,3,7,43,1807,3263443,10650056950807.

Induction gives s_n−1=product_(i<n)s_i. Therefore distinct terms are coprime,
and the telescoping identity

    1/s_n = 1/(s_n−1) − 1/(s_(n+1)−1)

proves sum_(i<N)1/s_i=1−1/(s_N−1), hence the infinite sum is1.
The integer 1807=13·139 is composite, but 3263443 is prime (trial division
through its square root). Thus the initial all-prime run stops at43; Sylvester
primality does NOT stop forever there.

The requested cycle is exactly the orbit modulo19:

    2 → 3 → 7 → 5 → 2.

This is verified from f(7)=43≡5 and f(5)=21≡2. Consequently19 divides no
Sylvester number. Meanwhile the complete eventual orbits are

    mod13:  2 → 3 → 7 → 4 → 0 → 1 → 1;
    mod139: 2 → 3 → 7 →43 → 0 → 1 → 1.

The same factorization therefore supplies both absorption primes13,139 and
an exclusion prime19 on the opposite neighbor: f^4(2)−2=1805=5·19².
This is a genuine dynamical role for19, not just another named prime.

## 3. Prime-divisor admissibility is not orbit accessibility

For a prime p>3 dividing f(x), x³≡−1 modp because
(x+1)f(x)=x³+1. Here x≠−1 modp (else f(x)=3), and x≠1 modp.
Thus x has exact multiplicative order6, so p≡1 mod6.
This proves the necessary restriction for every Sylvester divisor other than2,3.
The factor primes13 and139 satisfy it, but so does19, which never divides any
term. Roots of f exist modulo19 (8 and12); the seed2 orbit never reaches them.

The missing coordinate is the **pointed forward orbit**, not the splitting
condition or merely the existence of roots. A deterministic polynomial on a
finite field is eventually periodic; a hit of0 is followed by the absorbing
fixed point1. Thus a prime-divisor test is exactly reachability of0 from2,
and is decidable in at most p+1 sampled states, without constructing huge s_n.
This finite automaton is not by itself a quasirandomness theorem.

## 4. Exact periods at every power of19

**Theorem.** The orbit of2 under f modulo19^k is purely periodic, with exact
period

    L_k = 4                         (k=1,2),
    L_k = 12·19^(k−3)               (k≥3).

### Base and first-return computation

Let g=f^4 and h=g³=f^12. Exact arithmetic gives

    g(2)−2 = 1805 =5·19²,
    g'(2) =3·5·13·85=16575 ≡7 mod19.

Since7 has exact order3 modulo19 (7²≡11, 7³≡1), the first-return map has a
nontrivial linear multiplier. The orbit modulo19 has length4. The return
above implies the same length modulo19², not merely a divisor of4.

A short modular certificate, iterating f modulo19^4=130321, is

    2,3,7,43,1807,5418,27082,91376,76173,
    98195,100004,127115,116605.

Hence h(2)−2≡116603=17·19³ mod19^4, so its valuation is exactly3.
By the chain rule h'(2)≡7³≡1 mod19: each of the three g-block seeds is2 mod19.

### The lifting lemma, with proof

Let p be odd, H an integer polynomial, a an integer,
H'(a)≡1 modp, and v_p(H(a)−a)=s≥2. Then for every positive integer t,

    v_p(H^t(a)−a)=s+v_p(t).                    (1)

To prove the p-step increment, write A=H(a)−a and B=H'(a).
Taylor expansion modulo p^(s+2) gives
H(a+z)=H(a)+Bz whenever p^s divides z, since z² is divisible by p^(s+2).
Iteration therefore gives

    H^p(a)−a ≡ A(1+B+...+B^(p−1)) mod p^(s+2).

As B=1+pb, the sum is p+pb·p(p−1)/2 modulo p², hence p modulo p².
Thus v_p(H^p(a)−a)=s+1. The derivative of H^p at a is again1 modp.
Repeating with H^p,H^(p²),... handles p-powers. For u prime to p, Taylor
expansion modulo p^(s+1) gives H^u(a)−a≡uA, proving (1).

Apply this with H=h,p=19,a=2,s=3. It gives
v_19(h^t(2)−2)=3+v_19(t).

It remains to exclude shorter returns not divisible by12. Every return is
first a multiple of4 by reduction modulo19, say4r. On the disk2+19² Z,
linearization of g modulo19³ yields

    g^r(2)−2 ≡5·19²(1+7+...+7^(r−1)) mod19³.

Because7−1 is invertible mod19, this vanishes iff3 divides r. Thus any
return modulo19³ or higher has length12t, and (1) supplies the least t,
namely19^(k−3). This proves all stated periods and pure periodicity.

### What the delay does and does not mean

There is an unchanged first lift from modulus19 to19², then a factor3,
then perpetual factors19. This resembles the split/grow language of
multiplicative-order lifting, but is not the same invariant:

    v_19(2^18−1)=1.

So19 is NOT base-2 Wieferich. For multiplication by2 the analogous delay is
v_p(2^ord_p(2)−1); here it is a return defect of a nonlinear polynomial at
one specified seed, with a separate multiplier-order3. Identifying these
would erase the essential dynamical coordinate. This is a hostile control
for any proposed common explanation of19 and1093/3511.

### Exact conditional equidistribution, not independent randomness

Let C={2,3,7,43,1807,5418,6505,2209,724,2169,3978,3653}, the twelve
residues in the orbit modulo19³. For every k≥3,
the orbit C_k modulo19^k is **exactly the full inverse image of C** under
reduction modulo19³. Inclusion follows by functoriality of polynomial
iteration. The inverse image has12·19^(k−3) elements, exactly the proved
period length, so inclusion is equality.

Consequently every admissible residue modulo19^k is visited exactly once
in each complete period. Its limiting frequency is1/L_k; every inadmissible
residue has frequency0. In the19-adic integers, the orbit closure is the
union of those twelve residue balls modulo19³, of Haar measure12/19³,
and time averages on residue cylinders equal Haar measure conditioned on
that union. This closure statement means density on every finite residue
level, not an assertion that the countable orbit equals its closure.

There is also exact high-digit Fourier cancellation. For k>3 and an integer
frequency a not divisible by19^(k−3),

    sum_(x in C_k) exp(2πi a x /19^k)=0.

Indeed each fiber x=c+19³t contributes a geometric sum over
0≤t<19^(k−3), which vanishes under precisely that condition. In fact this condition is **necessary as well as sufficient**: the sum
vanishes **iff**19^(k−3) does not divide a. To prove the converse, write
a=19^(k−3)b. The sum equals19^(k−3) times
P(ζ), where P(X)=sum_(c in C)X^c and ζ=exp(2πib/19³).
If ζ=1 then P(ζ)=12. Otherwise ζ has order19^t for some1≤t≤3.
Were P(ζ)=0, the minimal polynomial Φ_(19^t)(X) would divide P(X)
in Z[X] (monic polynomial division and irreducibility). But
Φ_(19^t)(1)=19, whereas P(1)=12, a contradiction.
For completeness, Φ_(p^t)(X)=sum_(j=0)^(p−1)X^(j p^(t−1)); its translate
by1 is Eisenstein at p: modulo p the translate equals
X^((p−1)p^(t−1)), while its constant term is p and is not divisible by p².
This supplies the irreducibility used above.

More generally a sum of N roots of unity of p-power order, with repetitions
allowed and all coefficients nonnegative integers, cannot vanish unless
p divides N. The converse is NOT asserted. Here12 is prime to19, so every
low-digit frequency retains a nonzero signature of the twelve-class support. This is a
rigorous **condition first, equidistribute second** mechanism: high digits
are perfectly uniform over the allowed support, while the full residue
space is highly nonuniform. It does not make consecutive samples independent;
they satisfy the exact relation x_(n+1)=f(x_n). Global quasirandomness
and independent random-prime behavior do not follow.

In fact this measure-preserving system is **not mixing**. Let Ω be the
orbit closure and μ its conditional Haar probability measure. Polynomial
iteration preserves μ because at every finite level it cyclically permutes
the permitted residues, whose μ-weights are equal. Fix one residue ball
A={x in Ω:x≡c mod19³} with c in C. Then μ(A)=1/12, and for every integer
j≥1 the twelve-cycle modulo19³ implies f^(-12j)(A)=A, with preimages taken
inside Ω. Consequently

    μ(A ∩ f^(-12j)(A))=1/12 ≠1/144=μ(A)².

Thus correlations fail to approach the product even along an unbounded
sequence of times. Exact conditional high-digit uniformity coexists with
an explicit persistent temporal correlation; nonmixing here is proved,
not merely an absence of a mixing proof.

## 5. Fermat versus Sylvester: same product mechanism, different residues

A useful uniform construction is a_0=c+1 and a_n=c+product_(i<n)a_i,
for a positive integer c. It obeys

    a_(n+1)=a_n²−c a_n+c.

Every a_n≡1 modc. For i<j, a_j≡c moda_i; any common prime divisor would
therefore divide c and a_i, impossible. Thus all terms are pairwise coprime.
The choice c=1 gives Sylvester; c=2 gives Fermat numbers
F_n=2^(2^n)+1, with F_(n+1)=(F_n−1)²+1 and product_(i<n)F_i=F_n−2.

The distinction is structural: later Sylvester terms are1 modulo earlier
terms; later Fermat terms are2. Their quadratic discriminants are−3 and−4,
respectively, but Fermat's power-of-two formula gives much finer constraints:
if p divides F_n, then ord_p(2)=2^(n+1), hence p≡1 mod2^(n+1).
Indeed2^(2^n)≡−1 while its square is1, forcing that exact order.
No prime-divisor congruence is sufficient for actual divisibility without
retaining the exponent or pointed orbit.

If F_n itself is prime p, put d=2^(n+1). Then

    2^d−1=(p−1)²−1=p(p−2),

so v_p(2^d−1)=1. Thus every Fermat prime is non-Wieferich, including any
not currently identified; this is a theorem about the shape, not an
exhaustiveness claim about the five standard examples. This recovers the
proved mechanism in `fermat_wieferich_order_20261014e.md` rather than
relabeling its elementary theorem as new. The first composite Fermat number
F_5=641·6700417 supplies the primality hostile control.

## 6. Why the selector stops although the Sylvester orbit continues

For S(n)=product of primes p such that p−1 divides n, the complete divisor
check at1806 is essential. The values d+1 for all d dividing1806 are

    2,3,4,7,8,15,22,43,44,87,130,259,302,603,904,1807.

Exactly2,3,7,43 are prime, so S(1806)=1806. Checking only the compositeness
of1807 would not prove this. In particular13−1=12 does not divide1806,
and139−1=138 does not divide1806. Factoring the next Sylvester term and
adjoining its prime factors is NOT the same operation as this selector.

The incoming fixed-bank theorem is stronger than this calculation: after
one update every prime bank includes2; a least differing coordinate b at
one time can only cause a later discrepancy at p≥2b+1, because a differing
odd dependency q must divide p−1 together with2. This proves fixed-coordinate
convergence and uniqueness, not a single global stopping time. We retain
that distinction. The selector has an acyclic decreasing dependency order,
whereas the mod19 Sylvester system has a genuine recurrent orbit. Transferring
convergence conclusions between them is therefore invalid.

## 7. Reproduction and declared controls

Run `python3 verify_recursions.py` in this directory. The script uses no
assert statements, no third-party libraries, and raises explicit errors.
Normal and `python3 -O` runs have identical outputs. `checks.json` records
periods through19^6 as4,4,12,228,4332,82308; the theorem above, not this
finite census, establishes all k. The script independently enumerates
orbits, trial-divides the specified primes, audits every divisor of1806,
and checks exact rational Egyptian sums and both product recursions.

Controls:19 passes the local mod6 filter but cannot divide a term;19 has
polynomial prime-square delay without being base-2 Wieferich;1807 composite
is followed by prime3263443; and the selector's full divisor audit prevents
mistaking a composite next term for a general stabilization theorem.
