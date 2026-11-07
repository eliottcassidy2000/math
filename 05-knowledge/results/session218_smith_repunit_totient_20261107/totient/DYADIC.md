# Dyadic totient repair

Status: PROVED elementary classification; historical novelty not claimed.

## E. Exact dyadic repair criterion and complete source shapes

For odd prime q and integer a>=1, 2^a q is a totient if and only if at least one of the following holds:
(i) 2^j q+1 is prime for some 1<=j<=a;
(ii) q-1=2^j for some 1<=j<=a.

Proof. Exactly one prime-power factor of n must supply the single q in phi(n). If supplied by p-1, then p=2^j q+1 for some j>=1, p occurs to exponent one, and j<=a. If supplied by p^(e-1), then p=q, e=2, and q-1 must be a power of two 2^j with j<=a. All other odd primes dividing n have p-1 a power of two and occur once. Necessity follows. For sufficiency in (i), take n=2^(a-j+1)(2^j q+1). In (ii), take n=2^(a-j+1)q^2. Both have the required totient.

This also enumerates complete fibers: choose one supplier p=2^j q+1 (once), or q^2 with q-1=2^j; choose a squarefree product of distinct other Fermat primes f with f-1=2^h; require their h sum H<=a-j. Put r=a-j-H. The 2-part is 2^(r+1) if r>0, or either 1 or 2 if r=0. A supplier q^2 excludes q from the remaining product. Each preimage appears exactly once, since the q-supplier is unique.

Thus the first dyadic repair height is the smaller of the least j for which 2^j q+1 is prime and log2(q-1) if that is integral. Existence of some such j is not asserted for every q. Examples: 14->28=phi(29), 26->52=phi(53), 34->68 remains missing then 136=phi(137). The exceptional supplier is essential to enumerate the full fiber: for q=17, a=4, 273 is composite, but 272=phi(17^2). It is not essential for membership in this example: j=3 already supplies 137, and 272=phi(548). This criterion concerns output values, whereas phi(q)=phi(2q)=q-1 concerns twin INPUTS. Their coexistence in the Smith tuple is consistent: k=2q can be a nontotient value yet be one of two inputs with totient q-1.
