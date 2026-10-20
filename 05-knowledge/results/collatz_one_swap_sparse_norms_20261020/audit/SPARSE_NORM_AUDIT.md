# Independent audit of the one-swap norm bound

Date: 20 October 2026. Scope: mathematical audit only; no publication authorization. Verdict: the proposed theorem is sound, with the harmless period convention noted below. This establishes the stated uniform content bound and residual clock window, not universal positive-cycle exclusion or historical novelty.

## 1. Root choice and cyclic cuts

The unwrapped one positions of the zero-intercept lower mechanical word are s_j=ceil(jq/p)-1. For any unequal cyclic adjacent swap, let the moved one have index t+1 and root at its predecessor one t. The predecessor is distinct and unmoved because p≥2. A legal swap cannot move the one onto this root or onto its successor. Thus rooting removes the exceptional cut without losing the cyclic boundary swap.

The rooted positions are mk+η_m, where η_m=ceil((t+m)r/p)-ceil(tr/p), with only m=1 shifted by ε. Multiplication by 2^h clears the possible coefficient denominator; gcd(D,2)=1, so it preserves content. This proves g=gcd(|D|,P(2^k)). At k=1, some formal signed types are not actual swaps. The argument holds for each admissible type; it does not presume all 2p types are legal there. At k≥2 all one positions have zero neighbors, so both directions are legal.

## 2. Ring reduction and prime powers

Put A=(-tr mod p), a=(-(t+1)r mod p). Then pη_m-rm=a_{t+m}-A. In the ring z^p=2, substitution X=3z^(-r) gives

P(X)=2^h 3^(p-1) z^(-A) [S(z)+(2^ε-1)z^a],

where S(z)=1+z+...+z^(p-1), and (z-1)S(z)=1. Consequently the plus factor is Qplus=z^(a+1)-z^a+1, while the minus factor is z^a Qminus'/2 with Qminus'=z^(p-a)-z+1, up to the unit z-1.

This is valid modulo every divisor of g, not just its prime support. If pu+rv=1 and X=2^k modulo g, define z=2^u(3/X)^v. Because gcd(g,6)=1 and X^p=3^p/2^r, this gives z^p=2 and z^r=3/X. Both z and z-1 are units. The norm/resultant divisibility therefore includes all prime-power valuations.

## 3. Exact raw resultant

For F(X)=2^r X^p-3^p, all roots are X_j=3z_j^(-r), with z_j^p=2. Coprimality of p,r makes this a bijection of the p conjugates. The exact absolute identities are

plus: |Res(F,P)| = 2^[r(p-1)-A] 3^[p(p-1)] |Norm(Qplus)|;

minus: |Res(F,P)| = 2^[r(p-2)] 3^[p(p-1)] |Norm(Qminus')|.

For minus, η_1=(r+a-A)/p is 0 or 1 and h=1-η_1, explaining the cancellation of the apparent type-dependent exponent. Therefore the prime-to-6 part of the raw resultant is exactly M. In particular, g divides M.

## 4. Nonvanishing and norm bounds

The polynomial z^p-2 is Eisenstein at 2. Each Q, after reducing z^p=2, is a nonzero polynomial of degree below p, so its norm is nonzero.

For 1≤a≤p-2, the three exponents in each Q are distinct modulo p. With ρ=2^(1/p), root-of-unity averaging cancels the cross terms and gives the claimed averages of |Q|². Each average is strictly below 9; AM-GM yields |Norm(Q)|<3^p.

The endpoints are exactly:

- plus a=0: Q=z, absolute norm 2;
- plus a=p-1: Q=3-z^(p-1), absolute norm 3^p-2^(p-1);
- minus a=0: Q=3-z, absolute norm 3^p-2;
- minus a=p-1: Q=1, norm 1.

After stripping factors 2 and 3, M is positive odd and strictly below odd 3^p. Hence g≤M≤3^p-2.

## 5. Endpoint exactness and sharpness

For minus a=0, reducing Q sets z=3. Conversely, if d divides both 3^p-2 and 2^k3^(r-1)-1, then z=3 satisfies z^p=2, X=3z^(-r), F(X)=0 and P(X)=0 modulo d. Therefore

g=gcd(3^p-2,2^k3^(r-1)-1).

For plus a=p-1, the equations z^(p-1)=3 and z^p=2 force z=2/3. The same converse argument proves

g=gcd(3^p-2^(p-1),2^(k+r)-3^(r+1)).

Taking r=1, minus a=0, and k≥2 a multiple of ord_(3^p-2)(2) realizes g=3^p-2 by an actual legal swap. The bound is thus sharp separately for every p≥2. This does not make the corresponding rational fixed point integral: its clock is much larger.

## 6. Period and clock conclusion

Because g divides M,

g(k)=gcd(M,F(2^k),P(2^k)).

Thus ord_M(2) is a guaranteed k-period whenever M>1. For M=1, state the period as 1 rather than use an undefined or convention-dependent ord_1(2). Do not assert that this is always the minimal period. At k=1 the conclusion applies only to admissible physical swaps.

If the supplied one-swap word encodes a positive integer cycle, g=D>0; hence 3^p<2^q<2·3^p, forcing q=ceil(p log_2 3). This is a necessary residual window, not a proof that its candidates are impossible.

## 7. Independent finite checks

Two standalone scripts in this audit directory use SymPy only for exact resultants and direct positional carries for the word computations.

- audit_sparse_gate.py: 78,938 actual directed cyclic unequal-adjacent swaps, for 2≤p≤23, all coprime 1≤r<p, and 1≤k≤15. All polynomial content identities and norm divisibilities passed, including 2,722 legal k=1 swaps. Both endpoint gcd identities passed 2,394 cases each. All 550 signed type norms passed 0<M≤3^p-2.
- audit_resultant_identity.py: 748 exact raw resultant identities, for 2≤p≤12, every coprime residue r, all t and both directions. All passed.

The finite checks supplement the arguments above and are not the proof of the unbounded claims.
