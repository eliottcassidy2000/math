# Exact LCM closure of sigma neighborhoods

This supplements the Euler-aware support sieve with a second, independent constraint: permitted edges cannot in general be selected independently. The local constraints form a finite closure system, analogous to implication systems in formal concept analysis.

Fix an odd prime source p, a finite set of other odd prime targets R, and a source role: non-Euler means positive even exponent a, Euler means a=1 modulo 4. Let b=1 for the first role and b=2 for the second. Keep only targets allowed by the order test in the companion theorem. For each allowed target r define

k_r = ord_r(p) when ord_r(p)>1, and k_r=r otherwise.

For a subset T of allowed targets let

L(T)=lcm(b, {k_r:r in T}),
Cl(T)={r allowed : k_r divides L(T)}.

**Theorem.** Cl is a closure operator: extensive, monotone, and idempotent. A subset T of the allowed targets is the exact sigma neighborhood of some exponent of the chosen role if and only if Cl(T)=T.

**Proof.** Extensivity and monotonicity are immediate. Adding the requirements k_r already dividing L(T) cannot change that LCM, proving idempotence.

If T is an actual neighborhood, all k_r with r in T divide a+1, as does b. Hence L(T) divides a+1 and every target in Cl(T) is also actual. Thus Cl(T)=T.

Conversely, if Cl(T)=T and L(T)>1, take a=L(T)-1. For the non-Euler role L(T) is odd and at least 3, so a is positive even. For the Euler role its 2-adic valuation is exactly one, so a=1 modulo 4. The exact divisibility criterion gives precisely T. The remaining case L(T)=1 occurs only for the empty non-Euler neighborhood. Choose an odd prime ell larger than every k_r and set a=ell-1. Every k_r is at least 3 and less than ell, so none divides ell. The neighborhood is empty as required. If no targets exist, take ell=3. QED.

## A concrete conjunction that forces a new prime

For source p=3,

ord_11(3)=5, ord_13(3)=3, ord_4561(3)=15,

and 4561 is prime. Therefore, for every positive even exponent a,

11 | sigma(3^a) and 13 | sigma(3^a) imply 4561 | sigma(3^a).

Neither antecedent alone implies that conclusion: exponents a=4 and a=2 witness this. The conjunction forces 15|(a+1), hence the order-15 edge. If a hypothetical perfect-number support contains the first two targets but excludes 4561 while requiring both edges from 3, the entire exponent branch is impossible.

The implication is a direct consequence of standard order arithmetic, not a claimed new fact about these three primes. The reusable contribution is the exact closure-system formulation. It tells a graph search to store implications between edge sets instead of treating all subgraphs of H_q as feasible.

An implementation can derive implications symbolically from the k_r labels, before selecting exponents or factoring large sigma values. This retains arithmetic information lost by the graph-only maximal support sieve. It still does not ensure closure under all prime divisors or perfect valuation balance.
