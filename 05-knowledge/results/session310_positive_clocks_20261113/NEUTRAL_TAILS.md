# Same-count refill with arbitrarily many genuine tail sites

Write Q_T(t) for labelled path-slot assignments, allowing empty slots; this is the polynomial Σ_j F_j(T) binom(t,j), where F_j counts ordered covers into j nonempty paths. For ordinal joins Q is multiplicative, by restriction and concatenation in each slot (THM-3166). For a uniformly attached module, every Hamilton path in an exterior context has a reduced itinerary recording its consecutive module runs. Each itinerary visiting the module j times has exactly F_j fillings. Hence equality of Q, not just H, permits replacement in **every uniform context**, no matter how many visits. Selective ports inside the changed module are not preserved.

## Exact graph pair and unrestricted one-shot capacity

The following strongly connected tournaments are given by outgoing bit rows, labels starting at 0:
```
A = (2,28,25,81,97,79,7)
B = (2,80,3,7,173,143,189,15)
```
Their ordered nonempty cover profiles (starting at one slot) are
```
F(A) = (95,1086,5166,13200,19200,15120,5040)
F(B) = (95,2362,18756,73464,162000,205920,141120,40320)
```
The standalone replay derives these from subset paths and disjoint subset partitions. They prove Q_B(t)=t Q_A(t)=Q_(A→TT1)(t), by the binomial multiplication identity t binom(t,j)=j binom(t,j)+(j+1) binom(t,j+1). Thus a *strong* eight-vertex core and a same-order decomposable graph cannot be distinguished by **any** such external readout, yet one has a detachable ordinal sink and the other does not.

There is a second proof explaining the match. Apply the odd-cycle packing formula separately inside each of t labelled colors, summing over vertex colorings. Each packed cycle must be monochromatic, and every uncovered vertex has independent color; a packing of cycles with lengths ℓ_1,…,ℓ_k contributes 2^k t^(n−Σ(ℓ_i−1)). A and B each have ten triangles and seven heptagons; A has eighteen pentagons and six disjoint triangle pairs, whereas B has thirty pentagons and no disjoint triangle pairs. No other nonempty packing types fit, except a disjoint triangle/pentagon at order eight; the literal replay checks none occur. A pentagon has the same exponent debt 4 as two triangles but half the weight. Thus B trades every triangle pair for two pentagons, exactly preserving the coefficient; adjoining a neutral singleton to A only raises every exponent by one. This is a cancellation of structural information in the polynomial, not a numerical coincidence at a few integer evaluations.

For r≥1 write A^r,B^r for r ordinal copies. Then Q_(B^r)=Q_(A^r→TT_r). Both have 8r vertices and count 95^r. The first has zero detachable singleton tail sites since its terminal component is strong B; the second has r. This is a neutral refill with **unbounded capacity**, preserving all exterior arcs and each outside vertex name. It does not preserve internal arcs or an inside named port. It does not contradict S309's fixed-bank theorem, which required a changed intercept; here every cover coordinate and hence the dyadic valuation are unchanged. Physical capacity is not that valuation.

In the Π frame of Session 309 put B^r on the right and any actual nonempty M on the left. Rebase to A^r→TT_r. For k=0,…,r move k trailing singleton vertices, last first, across the separator c, to the front of M. Each step is exactly one arc flip: only incidence to c differs between the old and new positions. This gives
```
H_k = 60 Q_M(3) · 1276^r · 2^(r−k) · 3^k + H(M) · 95^r.
```
The exterior and A interiors stay fixed throughout the transfers. The final right module A^r is nonempty and has no trailing singleton; no fictitious empty-module rule is used. These Π outputs are divisible by 5, so this particular Π application does not supply prime or 7p targets.

## Prime and seven-prime in a different physical frame

Use the cyclic substitution C3[B^r,TT_a,TT_b] with arcs in the displayed cyclic order. The same neutral rebase works for all a,b≥1 by full profile equality, even when visits exceed three. At r=a=b=1 the count is Q_B(2)+H(B)=2647, prime. More generally for b=1 it is
```
2·95^r + (2^a−1)(2552^r−95^r).
```
Indeed only the block of size one can be visited once, and the other two have run multiplicities 1 or 2. The three starts with multiplicities (1,1,1) and the unique starts for (2,1,1), (1,2,1), (2,2,1) give this formula. In the prime fixture the detached site is real, but the Π transfer claim does **not** automatically apply to this cyclic frame: the two cyclic slots have different arcs to the retained right core, not just a separator.

A second exact fixture is r=1,a=5,b=10: the count is 1,598,340,819,716,431 = 7 × 228,334,402,816,633, with prime cofactor. The replay enumerates cyclic itineraries: any possible run-count triple differs coordinatewise by at most one, and the valid start letters determine its multiplicity. It certifies the cofactor by the full factorization of p−1 and the multiplicative order of 5 modulo p; the largest prime divisor of p−1 needs trial division only to 254,402.

A smaller, independently replayed 7p context takes an ordinal spine with uniform incidence word 10100010 from a shared apex (the apex beats the 1-slots). Inflate the first slot by B, leaving the other seven slots singleton. Both literal versions have H=231161=7×33023, with prime cofactor. See verify_seven_splice.py.

Thus universal neutral replacement can add arbitrarily many physical sites in Π, while a different frame yields actual prime and 7p neutral releases from the same pair. These are not a repeatable prescribed-target return, do not claim infinitely many prime cofactors, and do not combine the Π moves with the cyclic-frame prime fixtures without further physical constructions.

## Source trail and structural lesson

THM-3166 supplies ordinal multiplicativity; Session 306 supplies the context/run interpretation. Session 309 separates actual tail positions from intrinsic clock valuations, and its Lonely Runner refill analogy stresses that a reserve must be located where the next move can use it. Here even *the entire path-colour polynomial* fails to locate a sink: its correct use is to certify contextual neutrality, while an independent SCC certificate measures the usable sites. The theorem deliberately keeps the two observers separate.
