# Binary precision at one exact clock in quadratic log budget

The equal-clock dyadic one-apex switches of Session 326 implement every unit residue of H modulo (2^e) at **one** exact physical (Q(2)) with (\log Q=O(e^2)), rather than making exponentially many copies of a fixed switch. They can simultaneously keep H coprime to 99. More generally, for fixed (D=3^f11^g) they realize all units modulo (2^eD) at one exact (Q(2)), with (\log Q=O_D(e^2)). These are ordinal source constructions, not assertions of strongness or exact target coverage. The orders of the coefficients are separated from the literal clock.

## 1. Physical source and ratio digit

Use the two one-apex cuts (1^k w_\pm 0^k), where
$$
w_- =10101010100,\qquad w_+=10111011100 .
$$
The chain is increasing, and the apex beats exactly its one bits. Each source is strong, of order (12+2k). Session 326's correlation formula gives their common (Q(2)=q(k)) and the two distinct Hamilton counts
$$
h_-(k)=2048\,4^k-1026\,2^k+119,\quad h_+(k)=h_-(k)+12,\quad
q(k)=354294\,9^k-242340\,6^k+50336\,4^k-2052\,2^k+2 .
$$
At any k both h's are odd. Call k **safe** when (k\not\equiv1,4\pmod{10}); it is exactly at those k that both h's are coprime to 99 (by the two quadratics in (2^k\pmod{11}), and neither vanishes mod 3). The individual plus source never has these exceptions.

Write (r_k=h_+(k)/h_-(k)) as a 2-adic unit. For (\ell>k),
$$
v_2(r_k/r_\ell-1)
 =v_2(12)+v_2(h_-(\ell)-h_-(k))=k+3.\tag{1}
$$
Indeed (h_-(k)-119=2^{k+1}(1024\,2^k-513)) has valuation (k+1). Form the ordinal pair using either signs ((+,-)) or ((-,+)) at indices k,ℓ. Both have the same order and Q, exactly (q(k)q(\ell)), and their H ratio is (r_k/r_\ell). Thus (1) is a literal switch acting first at the precise bit (k+3), not just a density or asymptotic assertion.

## 2. Independent dyadic stages with common exact Q

First take the following two literal source switches
$$
(H_-,H_+,Q)=(73,85,1088),\quad (373,379,5896).
$$
Their ratios are 5 and 7 modulo 8 and generate every odd residue there, so one copy each seeds all four classes. The first is the inherited order-seven pair; the second is a new order-eight pair with masks `(82,28,153,209,64,31,166,51)`, `(196,101,216,51,227,197,8,74)`. Within each switch both orders match, and all four H values are coprime to 99. Then for each k=0,…,e−4 use one physical opposing pair with ℓ=k+1. Its two variants have the same H residue modulo (2^{k+3}), and opposite choices of bit (k+3) up to a common unit factor. Inductively all (2^{e-1}) unit residues occur, with exactly one choice per residue; odd parts of the H values are otherwise unrestricted. The entire bank shares one order and one Q, and for these choices
$$
\log Q=(e-3)^2\log 9+O(e),
$$
because (q(k)=\Theta(9^k)) and the sum of the two indices per bit is exactly ((e-3)^2).

To keep H coprime to 99, instead use consecutive safe indices k,ℓ. Their gaps (delta=\ell-k) are at most two. Make (2^\delta-1) copies of the opposing switch and choose its second output 0,…,(2^\delta-1) times. If (g=1+2^t u) with u odd and t≥2, then (v_2(g^m-1)=t+v_2(m)) (successive squaring proves this for powers of two, and the geometric sum deals with odd m). Hence the possible choices cover all (\delta) higher bits from t, staying fixed below. A final truncated stage uses only the copies needed for the specified e. This still costs (\log Q=O(e^2)); no supposedly independent phase is hidden in the forbidden indices.

At fixed (D=3^f11^g>1), seed all units modulo (16D) using Session 327's four-switch orbit (`EQUAL_CLOCK_ORBITS.md`) (the only inherited finite group input). Raise each opposing-pair option to (P=\lambda(D)) by using P literal copies of that same option, so their H ratio is 1 mod D. As (v_2(P)=1), (1) now acts first at (k+4). Continue with the safe gaps above. The base cost only depends on D and later cost is (O_D(e^2)). Large odd-prime precision is not claimed efficient here. Unlike in the pure-dyadic construction, source orders can vary. When 3 divides D, varying over all H units mod 3 at one exact Q actually forces both parities, by (Q(2)\equiv(-1)^nH\pmod3).

All options can be selected target-first: set every future pair in its fixed first orientation, choose the initial orbit to match the target at the small base modulus after dividing out that baseline, then process the successive new bits. Every unfinished stage has ratio 1 at the already assigned digits. This avoids both a random search and an assumption that independently correct residues must coincide at the physical clock.

## 3. Absolute information content and capacity

For this alphabet (two switches at each k, finitely many fixed additional controls, and singleton padding), at order at most n there are at most (\exp(O(\sqrt n))) component multisets by the colored-partition bound. H and Q of an ordinal sum depend only on the multiset. Since (Q\ge2^n) at order n, any one Q carrying at least (2^{e-1}) distinct H residues must have (\log Q=\Omega(e^2)) *within this alphabet*. The bound and construction match on the logarithmic scale. This is not a lower bound over all tournaments.

The same source alphabet has an unconditional real separation result. Choose once at each k=0,…,m−1 rather than pairing. Because (h_-(k+1)>4h_-(k)), if (x_k=12/h_-(k)) then (\sum_{j>k}\log(1+x_j)\le x_k/3<\log(1+x_k)). Thus the log-ratios are reverse superincreasing: the first differing choice cannot be compensated by any later choice. This constructs (2^m=exp(\Theta(\sqrt{\log Q}))) **distinct exact integer H** at one exact (Q=\prod q(k)), all within one bounded multiplicative interval, and the same construction works on safe indices only. It also makes clear what is still missing for Chan–Pak: even optimal information use in this alphabet creates a sparse family, not a route to a prescribed absolute target at an arbitrary affordable budget.

Replay: `python3 replay_phase_bits.py`. The replay checks the physical cut formulas, the exact first-different-bit valuation, the full dyadic orbit and safe gaps through precision 13, and smaller mixed-modulus target inverses. The infinite precision statements follow from (1) and the odd-unit valuation law.
