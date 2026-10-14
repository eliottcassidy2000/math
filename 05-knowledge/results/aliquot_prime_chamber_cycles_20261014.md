# Aliquot cycles through odd-cycle blocks and prime-cofactor chambers

Research synthesis, 14 October 2026. This note separates classical facts, direct derivations, finite exact evidence, and an open proposed classification. Novelty of the direct derivations and conjecture has not been established. No general exclusion of aliquot 3-cycles is claimed.

## 1. Result first: a narrow unbounded exclusion and a sharp next hypothesis

Let `s(n)=sigma(n)-n`. There is no aliquot cycle of three distinct integers whose members are

`2p, 2q, 2^c r`, with `c>=1` and `p,q,r` odd primes.

The proof below does not impose a bound on `c` or the primes. Its algebraic relaxation has exactly one nonconstant integral closure:

`46 -> 26 -> 16 -> 46`.

But its last cofactor is `1`, outside the prime chamber. The actual sequence is

`46 -> 26 -> 16 -> 15`.

Thus the obstruction is a precise failure of the prime-factorization hypothesis at the only surviving algebraic boundary. It is not a cancellation of a putative directed triangle.

More generally, the same calculation gives an exact divisor gate for every repeated-minimum exponent itinerary `(a,a,c)`, `c>=a`, and proves `c<=3a+1` whenever that itinerary closes integrally. This is useful before any primality testing.

**Open proposed classification.** For positive integers `a,b,c`, the unique rational closure of the dyadic prime-cofactor template is integral if and only if, up to cyclic rotation,

- `a=b=c`, with cofactors all `2^(a+1)-1`; or
- `(a,b,c)=(1,1,4)`, with cofactors `(23,13,1)`.

This is a statement about the explicitly defined affine template, not about arbitrary aliquot cycles. If proved, it would exclude nontrivial aliquot 3-cycles all of whose members have exactly one odd prime factor, occurring to exponent one. It would leave numbers with multiple odd prime factors or higher odd-prime exponents untouched.

Exact evidence covers all `256^3=16,777,216` ordered exponent triples, using their `5,592,576` cyclic representatives. There were exactly `256` constant integral closures and the one displayed nonconstant closure. This finite result does not establish the classification at arbitrary exponents.

## 2. The cycle landscape and the genuine divisor carrier

Perfect numbers are cycles of length one, and amicable pairs have length two. The [primary list maintained by David Moews](https://djm.cc/sociable.txt), dated 11 September 2026 when checked, includes lengths `4,5,6,8,9,28` beyond those. Its list contains no length-three example. [Pomerance's survey slides](https://math.dartmouth.edu/~carlp/upintconf.pdf) explicitly identify existence of a 3-cycle as unsolved. The list is a record of known cycles and specified searches, not a classification of all possible cycle lengths.

For `n>=2`, the canonical weighted divisor incidence row is

`B(n,d)=(d/n) * 1_{d|n}`.

Its row sum is `sigma(n)/n`, equivalently `sum_{d|n}1/d`, and its excess over two is exactly

`delta(n)=sigma(n)/n-2=(s(n)-n)/n`.

This gives the deficient/perfect/abundant trichotomy as a genuine signed quantity on divisor data. Its prime-power expression is

`sigma(n)/n = product_{p^e || n}(1+1/p+...+1/p^e)`.

The aliquot transition graph has one outgoing edge `n -> s(n)` when `s(n)>0`. Give this edge multiplicative weight `s(n)/n` or additive weight `log(s(n)/n)`. Every positive cycle satisfies

`product_i s(n_i)/n_i=1`, and `sum_i log(s(n_i)/n_i)=0`.

A nontrivial cycle therefore contains both abundant and deficient members and no perfect member. The five-cycle

`12496,14288,15472,14536,14264`

has sign word `AADDD`. There is no rule that the three signs require a cycle of length three: a perfect member is already trapped in a loop.

The log identity is a telescoping potential identity for any positive cycle. It cannot itself distinguish lengths three and five. Nor does a tournament obtained by sorting abundancy ratios add arithmetic information: that ordering is transitive by construction.

The natural orbifold made by using every divisor `d>1` as a cone order has Euler characteristic `2-tau(n)+sigma(n)/n`. At perfect numbers this is `4-tau(n)`, not uniformly zero: it is zero at `6` and `-2` at `28`. The reciprocal-sum resemblance therefore does not identify this trichotomy with spherical/Euclidean/hyperbolic geometry.

## 3. A canonical shared block with tournament odd-cycle methods

For an actual length-`k` cyclic itinerary, let `C` be the cyclic shift, `(C n)_i=n_(i+1)`. Then the defining equations are exactly

`(I+C)n = (sigma(n_i))_i`.

This is the shared linear block; the aliquot graph itself is not being declared a tournament.

For odd `k`,

`(I+C)(I-C+C^2-...+C^(k-1))=2I`.

Hence `det(I+C)=2`, and its Smith normal form is `diag(1,...,1,2)`. For even `k`, the determinant is zero and the Smith form is `diag(1,...,1,0)`. One direct proof uses the integer presentation with relations `e_i+e_(i+1)=0`: all generators reduce to alternating copies of one generator; closure gives `2e_0=0` for odd length and no final relation for even length. The leading `(k-1)` square minor has determinant one.

This is exactly the local `2`-per-odd-cycle factor appearing in the repository's [odd-cycle-collection expansion, THM-002](../../01-canon/theorems/THM-002-ocf.md), whose tournament cycles have length at least three. On a disjoint collection of directed cycle blocks, `det(I+C)` is the product of these local factors: an even cycle kills that product, while each odd cycle supplies two. This identifies a common canonical operator, not a transfer of Hamiltonian-path counts to divisor sums. In particular, its odd-order factor is the same for `3,5,9`; more arithmetic data are needed to separate them. The aliquot length-one block also has determinant two, but a tournament has no corresponding loop.

The repository's [THM-1926 zeta/strong-component method](../../01-canon/theorems/THM-1926-tournament-zeta-euler-product-over-strong-core.md) also transfers without changing graph type. On the finite partial functional graph with vertices `2,...,N` and an edge only when the true successor stays in that set,

`1/det(I-zA_N) = product_{cycles C wholly in the set}(1-z^|C|)^(-1)`.

Transient trees contribute one. If `t_j=tr(A_N^j)`, the number of primitive 3-cycles is `(t_3-t_1)/3`, not `t_3/3`: the perfect loops must be subtracted. A finite cutoff establishes only cycles wholly contained in that cutoff.

## 4. Odd-cycle arithmetic: gcds and prime suppliers

For any positive cyclic list of odd length, put

`g=gcd(n_i)`, `h=gcd(n_i+n_(i+1))`.

The alternating inverse shows `g | h | 2g`, so `h=g` or `2g`. More precisely, `h=2g` if all `n_i/g` are odd, and `h=g` otherwise. For an aliquot cycle these neighbour sums are the actual `sigma(n_i)`.

The resulting common-unitary-gcd obstruction is **classical, due to Dickson**, not a new theorem here: a nontrivial odd-length sociable cycle with `g>1` cannot have `g` unitary in every member. It is Theorem E in [Kobayashi, Pollack, and Pomerance, On the distribution of sociable numbers](https://www.pollack-math.net/KPP.pdf). Briefly, a common unitary gcd gives `sigma(g)|h|2g`, forcing `sigma(g)=2g`; then writing `n_i=g m_i` gives `s(n_i)=g(2 sigma(m_i)-m_i)>=n_i`, with strict inequality for `m_i>1`, contradicting a nontrivial cycle. The boundary `g=1` is not excluded.

Here is a finer usable carrier. Write `a_p=min_i v_p(n_i)`, and let

`U=product_{p: v_p(n_i)=a_p>0 for every i} p^a_p`.

Thus `U` contains exactly the common prime powers whose exponents never change. It is unitary in every member. Multiplicativity and the odd-cycle block imply the necessary condition

`sigma(U) | 2g`.                                            (Supplier gate)

This can be represented as a directed bipartite graph from a prime-power node `q^e` to a target-prime node `p`, weighted by `v_p(sigma(q^e))`. Retain the cycle position owning each prime-power occurrence. For stable common prime powers the total output to `p` is at most `a_p+1_{p=2}`. In particular, a prime outside `2g` cannot divide `sigma(U)`. Erasing the position label would mix suppliers from different integers and invalidate the test.

For odd `p`, the full exact minimum condition is

`min_i sum_{q^e || n_i} v_p(sigma(q^e)) = a_p`.

For `p=2`, it is `a_2+1` if all normalized members are odd, and `a_2` otherwise. If a three-cycle has nonconstant 2-adic exponents, exactly two of its three values `v_2(sigma(n_i))` equal `a_2`; the third is greater. This follows by considering the two edges incident to a unique minimum, or the two unequal-valuation edges when the minimum occurs twice. It is a necessary three-position supplier pattern, not a sufficient cycle criterion.

These conditions survive changes of factorization chamber. In the five-cycle above, `g=8`, `U=1`, the 2-adic exponent word is `(4,4,4,3,3)`, and the sigma-valuation word is `(5,6,3,7,3)`. The changing exponent of the common prime two is essential to escaping the common-unitary-gcd obstruction.

## 5. Exact affine closure of prime-cofactor chambers

Fix positive cores `c_i` and require

`n_i=c_i p_i`, with `p_i` prime and `gcd(c_i,p_i)=1`.

The prime/coprimality assumptions give `sigma(c_i p_i)=sigma(c_i)(p_i+1)`. Consequently,

`c_(i+1) p_(i+1)=s(c_i)p_i+sigma(c_i)`.

Equivalently, on the integer `n`, the branch is

`f_c(n)=alpha_c n+beta_c`, `alpha_c=s(c)/c`, `beta_c=sigma(c)>0`.

For a fixed cyclic core itinerary, composition gives `n_0=A n_0+B`, where `A=product_i alpha_(c_i)` and `B>0`. A positive closure therefore requires `A<1`, and then it is the unique rational point `n_0=B/(1-A)`. Integer closure, primality, coprimality, and distinctness are separate gates.

There is also an exact composition-order comparison:

`f_c(f_d(n))-f_d(f_c(n))=beta_c beta_d (t_c-t_d)`,

where `t_c=(alpha_c-1)/beta_c=1/c-2/sigma(c)`.

So comparison of these two compositions is a scalar preorder, with ties retained. This is a genuine map-level comparison, but it is transitive and cannot prohibit cycles of a switched system. Completing it to an arbitrary tournament would add no valid obstruction. A single fixed core branch is nondecreasing and cannot itself support a nontrivial finite cycle.

## 6. Dyadic three-chamber formula and the unbounded divisor gate

Set `c_0=x=2^a`, `c_1=y=2^b`, `c_2=z=2^c`. The equations are

`y q=(x-1)p+2x-1`,
`z r=(y-1)q+2y-1`,
`x p=(z-1)r+2z-1`.

Their unique rational solution starts with `p=N/D`, where

`D=xy+yz+zx-x-y-z+1>0`,
`N=(z-1)(y-1)(2x-1)+(z-1)y(2y-1)+yz(2z-1)`.

For the repeated-minimum family `a=b<=c`, put

`D=(2x-1)z+(x-1)^2`,
`Q=(2x^2-1)(3x^2-3x+1)`.

**Proposition, direct derivation.** All three cofactors are positive integers if and only if `D|Q`.

Proof. The exact polynomial identity is

`(2x-1)^2 N+xQ`
` = D(6x^3+4x^2 z-10x^2-2xz+5x-1)`.

Here `gcd(D,x)=gcd(D,2x-1)=1`: `z` is divisible by `x`, and `gcd(x-1,2x-1)=1`. Thus `D|N` if and only if `D|Q`. If `p=N/D` is integral, reduction modulo `x` gives `p=-1 mod x`, so the first equation makes `q` integral. The second gives a rational `r` with denominator dividing `z`; the third says `(z-1)r` is integral. Since `gcd(z,z-1)=1`, `r` is integral too. All three are positive by the affine formulas. Conversely, integral closure plainly requires `D|N`.

Since `D|Q` implies `D<=Q`, it also gives

`z <= (Q-(x-1)^2)/(2x-1) < 4x^3`, hence `c<=3a+1`.

For the strict bound, use `Q<6x^4` and `2x-1>=3x/2` for `x>=2`.

Finally take `a=1`, so `x=2`, `D=3z+1`, `Q=49`. Since `z>=2`, the only possibilities are `D=7` or `49`, giving `z=2` or `16`. The first yields cofactors `(3,3,3)` and the perfect loop `6`; the second yields `(23,13,1)` and the virtual loop. The latter violates primality at `r=1`, where the formula `sigma(r)=r+1` is false. This proves the unbounded exclusion stated at the beginning.

## 7. Exact checks and what they establish

The dependency-free companion script produces a JSON certificate. It checks:

- All cycles wholly in `[2,10^6]`: four loops, forty amicable pairs, one five-cycle, and one 28-cycle; no three-cycle in that finite graph.
- The cycle-block alternating identities for all orders `1,...,31`.
- The displayed known cycles of lengths `1,2,4,5,9,28` by exact factorization and divisor sums.
- All `256^3` dyadic three-chamber itineraries up to cyclic rotation, with the counts stated above.
- The rational closure formula through a separate `Fraction`-based composition for all `12^3=1,728` ordered exponent triples.
- The repeated-minimum polynomial identity and divisor equivalence over `1<=a<=64`, `a<=c<=3a+2`, including values just beyond the proved bound.
- The actual failure `s(16)=15` of the virtual closure.

Run from the repository root after staging:

`python -B 04-computation/experiments/aliquot_prime_chamber_cycles_20261014.py`

The script uses explicit error checks rather than Python assertions, so optimization does not disable verification. Normal and optimized output should agree byte-for-byte.

The next worthwhile proof target is the explicit integral-template classification in Section 1, beginning with the repeated-minimum divisor gate. A counterexample should be reported first as an exponent/cofactor tuple, then tested against the separate prime and coprimality gates. A proof of that classification would be a restricted factorization theorem; it would still not settle existence of arbitrary sociable 3-cycles.
