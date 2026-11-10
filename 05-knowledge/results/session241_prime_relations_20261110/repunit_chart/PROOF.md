# A strong-tournament chart for the entire natural-number divisor lattice

**Status: PROVED elementary synthesis from inherited fan constructions and classical repunit identities; independently audited in the adjoining embedding audit. No historical-priority claim.** The normalized values are coordinates, not automatically Hamilton-path counts.

## 1. An actual strong tournament for every positive index

For n=1 let Q_1=C3. For n>=2, let Q_n have vertices0,...,2n+1, with1,...,2n+1 transitive, while0 loses precisely to n,n+1,2n+1. The Hamiltonian cycle0,1,...,2n+1,0 proves it is strong. THM-4748, `01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md`, gives the inherited Fibonacci-fan identity

    h_n := H(Q_n)=4^n−1.                                  (1)

It uses2n+2 vertices for n>=2 and3 for n=1. The count and construction are inherited, not rediscovered here. The proof there cancels the three apex singletons against their pair interactions. The source is pinned in `../provenance/source_pins.json`.

Define the optional normalized scalar R_n=h_n/3. This is an integer, but R_3=21 is a proved forbidden Hamilton-path count. Dividing a count by3 is not a tournament operation.

## 2. Divisibility and gcd are exact, without factoring the index

For all a,b>=1,

    gcd(h_a,h_b)=h_gcd(a,b),
    h_a divides h_b iff a divides b.                      (2)

The same identities hold for R_n, with R_1=1.

**Proof.** The Euclidean step

    (4^b−1)−4^(b−a)(4^a−1)=4^(b−a)−1

for b>=a shows that the exponent Euclidean algorithm computes the numerical gcd. At the last nonzero exponent d=gcd(a,b), the gcd is4^d−1. Because h_n is strictly increasing, gcd(h_a,h_b)=h_a iff gcd(a,b)=a. Finally gcd(3R_a,3R_b)=3gcd(R_a,R_b), and multiplying both sides of a divisibility by3 preserves it. ∎

Thus the ordered set {h_n:n>=1}, with **numerical divisibility restricted to this image**, is isomorphic to the positive-integer divisor lattice. Its bottom is3, its meet is the ordinary numerical gcd, and its join *inside the chart* is h_lcm(a,b). A chart atom above3 corresponds exactly to an ordinary prime index, although its numerical value is composite: the prime index3 maps to63.

Every positive integer is represented, including the index7, whose image is16383. This does not realize the forbidden target7; the encoded number and the actual path count are different coordinates.

## 3. Totient and arithmetic functions have an explicit preserved predicate

For n>=1, Euler's totient satisfies

    phi(n)=#{1<=a<=n : gcd(H(Q_a),H(Q_n))=3}.             (3)

Indeed (2) makes gcd(h_a,h_n)=3 equivalent to gcd(a,n)=1, including n=1 with phi(1)=1. Consequently

    sum_{a: h_a divides h_n} phi(a)=n.                   (4)

The finite divisor interval below h_n recovers the divisor count tau(n), and its incidence-poset Möbius function from3 to h_n equals the ordinary arithmetic mu(n). This is transport along the explicit poset isomorphism, not a new formula for the totient or Möbius function. The source prime index is retained. No claim that phi(h_n)=phi(n), or that an adjacency score equals phi(n), is made.

The image-value totient is a **different** function. The inherited initial-segment theorem THM-4801 (`01-canon/theorems/THM-4801-collatz-meets-the-constructible-polygons-trunk-five-initial-segment-principle-collatz-closed-family.md`) gives the exact corollary

    phi(h_n) is a power of2 iff n in {1,2,4,8,16}.        (5)

This is unconditional and does not assume only five Fermat primes exist. For completeness, if n has an odd prime divisor r, any prime q dividing2^r−1 has ord_q(2)=r, so r divides q−1. The prime q divides h_n, preventing phi(h_n) from being a power of2. Thus n=2^j. Then h_n=F_0...F_j. The known prime values F_0,...,F_4 suffice for j<=4. For j>=5, Euler's factor641 of F_5 divides h_n and640 divides phi(h_n), so5 divides it. This is precisely the inherited initial-run mechanism specialized to the tournament chart.

The identity R_n=(4^n−1)/3 also identifies the **same integer coordinate** with the familiar Collatz trunk3R_n+1=4^n. Only this equality is transferred; no Collatz trajectory or conjecture is proved here. The forbidden R_3=21 is the cheapest test against mistaking normalization for a legal tournament constructor.

## 4. Prime visibility reads an order, not primality

For any odd prime p, let d_p=ord_p(4). Then

    p divides h_n iff d_p divides n.                    (6)

If n=d_p t, the lifting-the-exponent identity gives

    v_p(h_n)=v_p(4^d_p−1)+v_p(t).                        (7)

For p=3, d_p=1 and this says v_3(h_n)=1+v_3(n). These are classical multiplicative-order/LTE identities applied to actual path counts, not new primality tests. They distinguish prime-factor visibility from divisibility between whole chart values. Equal prime orders still need not determine richer continuations, as Session219's641/6700417 escape demonstrates.

## 5. Exact addition law, and the relations that do not transfer

The chart satisfies the shifted multiplication law

    h_(a+b)=h_a h_b+h_a+h_b.                             (8)

It is just (h_a+1)(h_b+1)=4^(a+b). It gives a numerical description of adding *indices*, not an already-built local gluing operation on arbitrary tournaments.

The following failures are essential:

- **Ordinary multiplication:** h_2 h_3=15*63=945, but h_6=4095. Ordinal joining Q_2 and Q_3 realizes945, not the image of2*3.
- **Ordinary lcm:** lcm(15,63)=315, whereas the least upper bound inside the chart is h_6=4095. The image is gcd-closed, not lcm-closed in the ambient integers.
- **Ordinary addition:** h_1+h_1=6, while h_2=15. As all tournament H-values are odd, no map into actual H-values preserves unrestricted ordinary addition.
- **Primality:** all h_n are divisible by3, and h_n>3 for n>1. Prime indices are lattice atoms in the image, not numerical prime targets.
- **Induced containment:** Q_2 is an induced subtournament of Q_3, although2 does not divide3 and15 does not divide63. `exact.json` gives a full vertex map, checked arc by arc. Thus the arithmetic order is the relation on counts, not the ordinary tournament-containment order.
- **Coverage:** exactly floor(log_4(X+1)) chart values lie at most X for real X>=0, so this all-index representation is sparse. It does not address arbitrary prescribed H.

The complementary Fermat-gadget map in `../constructions/PROOF.md` preserves multiplication and ambient gcd/lcm; this chart instead provides a factorization-free construction from n, linear vertex order, and a transparent totient predicate. Their different operation laws explain why both are useful and why neither supplies Chan–Pak surjectivity.

## 6. Connection contract and checks

Source: positive integer n. Target: explicitly labelled strong tournament Q_n. Pairwise observable: forward chain order, with three specified apex reversals. Gauge and ties: fixed labels, no ties. Preserved predicate: divisibility and gcd of actual path counts as in(2). Lost predicate: ordinary addition/multiplication/lcm, numerical prime status and arbitrary target placement. Sidecar: the specified chart membership, with n decoded from H+1=4^n. It is not an intrinsic orientation of the partial divisibility relation between integers.

Reproduce:

    python3 05-knowledge/results/session241_prime_relations_20261110/repunit_chart/verify.py

The standard-library verifier checks actual adjacency Hamilton paths for n=1..6; all40,000 pairs a,b<=200 for gcd, divisibility and addition identities; both totient predicates in their explicit finite ranges; prime visibility for all odd primes<=101 and n<=200; and an actual induced embedding. Finite controls supplement, rather than establish, the all-index proofs. The independent audit uses separately written code and does not import this verifier.
