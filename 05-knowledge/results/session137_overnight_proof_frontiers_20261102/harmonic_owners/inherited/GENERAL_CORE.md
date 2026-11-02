# The least missing prime controls an all-q resonant-tail spectrum

Session114, 2026-10-30. **PROVED** by the elementary Farey-tent argument below. This extends the fixed q=12 certificate in `RESONANT_TAIL.md`; it does not prove unrestricted Lonely Runner Conjecture or claim literature priority.

## The theorem

Let q>=3 and let

 C_q={1,2,...,q-1,q+1},
 s(q)=min{u>=2:gcd(u,q)=1}.

Thus s(q) is the least prime not dividing q. Append an arbitrary positive integer speed v. Write M for maximum common clearance.

1. If q does not divide v, then M(C_q union {v})=1/q.
2. If v=qm with m>=1, then

 **M(C_q union {qm})=max{1/(q+2), m/(qm+s(q))}.**

There are q+1 distinct relative speeds unless v is already in C_q, so this family satisfies the (q+2)-runner bound, uniformly in q and v. For q=12 the least absent prime is5, giving the complete spectrum proved separately by a small fixed-core certificate.

For 2m>s(q), the two and only two maximizing points are

 t=a/q - (1/q-r)/s(q), and 1-t,

where r=m/(qm+s(q)) and a is the inverse of s(q) modulo q, represented in {1,...,q-1}. Equivalently, writing h=(a s(q)-1)/q,

 t=(a m+h)/(qm+s(q)).

The last-cell mass slope at r is exactly

 **B=2(1/s(q)+1/(qm)).**

The conditional resonant-tail theorem in `RESONANT_TAIL.md` supplies this once the complete top chamber is proved. The crucial input is the minimum over *binding directional owners*, after resolving each tie by its largest speed. It is not the minimum speed of the source: speed1 is present but its old direction is replaced by q+1.

## Tight infinite subfamilies

The multiple-tail family has exactly the conjectured lower value 1/(q+2) if and only if 2m<=s(q). In particular:

- m=1 is always tight: it is simply the full AP {1,...,q+1}.
- For each m>=2 it is tight exactly when every prime <=2m divides q.

Consequently, if P_m is the product of the primes at most 2m, every q divisible by P_m supplies a tight set

 {1,...,q-1,q+1,qm}.

The m=2 criterion is 6|q; the fixed q=12,m=2 example is the familiar GW-type tight row. Higher m give the same precise primorial divisibility criterion. These are proved consequences for this family, not a classification of all tight sets. The number of initial tight multiples can be arbitrarily large as q acquires more small prime factors.

## 1. The exact Farey tent, with proof

Set N=q-1. Consecutive fractions a/b<c/d in the Farey sequence of order N satisfy

 bc-ad=1, and b+d>N.

These elementary neighbor facts follow from the determinant-one mediant construction of Farey sequences. On their closed interval put

 x=bt-a, y=c-dt.

Then

 **min_{1<=j<=N} ||jt|| = min(x,y).**

Here is a direct proof, also clarifying why this is an exact envelope rather than a heuristic. For any integer k and 1<=j<=N, determinant one gives unique integers A,B such that

 (k,j)=A(a,b)+B(c,d).

Both A,B cannot be positive because j>=b+d>N; both cannot be negative because j>0. If they have opposite signs or one is zero, then

 |jt-k|=|A x-B y|=|A|x+|B|y>=min(x,y).

The last inequality uses that (A,B) is a nonzero integer pair. Thus every j has distance at least min(x,y). Conversely the terms j=b and j=d have distance at most x and y respectively. This proves equality, including chamber endpoints.

The tent has peak 1/(b+d), at the reduced mediant (a+c)/(b+d). Indeed d x+b y=1, and the peak is where x=y. Since b+d>=q, above radius 1/(q+2) only tents with b+d=q or b+d=q+1 can possibly survive. Tents with b+d>=q+2 have peak at most the excluded boundary.

## 2. The extra speed q+1 removes all denominator-(q+1) tents

Consider a tent with b+d=q+1 and center c0=(a+c)/(q+1). On either side write t=c0+h. The AP clearance equals

 1/(q+1)-e|h|,

where e is the corresponding positive integer side slope b or d. The additional speed q+1 has integer phase at c0, hence

 ||(q+1)t|| <= (q+1)|h|.

If both constraints were at least r, they would imply

 r <= 1/(q+1+e) <= 1/(q+2).

Thus no point in any such tent survives above 1/(q+2). This argument is valid even if the added speed wraps: ||z||<=|z| is an upper bound for all real z.

## 3. What remains at denominator q

The denominator-q tents are indexed by primitive grid points a/q, gcd(a,q)=1. Let p be the inverse of a modulo q in {1,...,q-1}, and n=q-p. Their AP side slopes are p and n. One way to see the exact indices is to take left-parent numerator (ap-1)/q and right-parent numerator (an+1)/q. Their denominators are p,n, their determinant is1, they lie on opposite sides of a/q, and their mediant is a/q; their denominator sum q>N makes them Farey neighbors of order N.

At a radius 1/(q+2)<r<=1/q, write epsilon=1/q-r. The unmodified AP window is

 [a/q-epsilon/p, a/q+epsilon/n].

At center a/q, the added speed q+1 has phase a/q.

- If a is neither1 norq-1, then p,n>=2. Throughout the AP window,

  |(q+1)(t-a/q)| <= (q+1)epsilon/2
                  < (q+1)/(q(q+2)) < 1/q.

  The central phase lies between2/q and1-2/q, so the added speed remains at distance greater than1/q, and hence greater than or equal to r. It makes no cut.

- If a=1, the original AP window is [r,(1-r)/(q-1)]. Its leftmost raw added-speed phase is (q+1)r>1-r, because r>1/(q+2). Therefore any portion with raw phase below1 is unsafe; it cannot supply a second wrapped component. The raw phase never reaches2 on this window for q>=3. On the remaining branch, safety requires (q+1)t>=1+r. Thus the positive owner1 is replaced by q+1. The new left endpoint is 1/q-epsilon/(q+1)=(1+r)/(q+1); the right endpoint remains 1/q+epsilon/(q-1). The reduced phase at the left endpoint is exactly r. The added-speed phase at the right endpoint is at most1-r throughout this radius interval: both sides are affine in r, so it suffices to check the endpoints. At r=1/q this says1/q<=1-1/q. At r=1/(q+2), the left side is

   (q+3)/((q+2)(q-1)) <= (q+1)/(q+2),

  equivalent to q^2-q-4>=0, true for q>=3. There is no phase wrap. Thus no further cut occurs.

- The case a=q-1 is the reflection of a=1.

We have therefore proved, for every q>=3, the **complete top-window description** above1/(q+2): at each primitive a/q the owner pair is (p,q-p), except that every occurrence of owner1 is replaced by q+1. At r=1/q this gives all and only the core maximizers. The minimum binding owner is exactly s(q), occurring on precisely two reflected sides. Every unit u in {1,...,q-1} occurs once as a positive owner and once as a negative owner; replacing u=1 leaves the least remaining one s(q).

## 4. Finish the tail spectrum

If q does not divide v, its phase at each primitive a/q is a nonzero multiple of1/q. Thus the maximum core clearance1/q is retained.

If v=qm and 2m>s(q), the candidate r=m/(qm+s(q)) lies strictly above1/(q+2). Apply the zero-phase tail theorem to the complete windows just proved. This gives the exact maximum, both maximizers, and the terminal derivative stated above. For the derivative, a maximizing point has exactly one core endpoint owner s(q) and the tail on the opposite side. No other core wall can tie strictly inside this chamber: its affine constraint is strictly positive at interior radii unless it is the same wall, whose slope would force the same speed. Reflection doubles the reciprocal-owner sum.

If 2m<=s(q), the same necessary inequality inside the top chamber would force r<=m/(qm+s(q))<=1/(q+2). Hence no point survives above1/(q+2). To attain1/(q+2), choose any primitive grid point a/(q+2). Every core speed is between1 andq+1, so none vanishes there. Neither does qm: divisibility (q+2)|qm would require

 m >= (q+2)/gcd(q,q+2) >= (q+2)/2,

whereas m<=s(q)/2<=(q-1)/2. Thus all distances are at least1/(q+2), completing the proof. This includes the boundary case2m=s(q), possible only for s(q)=2,m=1.

Finally s(q) must be prime: if it had a proper prime factor, that smaller prime would also be coprime to q, contradicting minimality. For m>=2 the even integer2m is not prime, so2m<=s(q) is equivalent to there being no prime at most2m absent from q. This proves the primorial criterion.

## Verification and limits

`check_general.py` checks the exact top-window formula against literal safe-interval intersections, compares complete maxima against all opposite-slope rational clocks for a finite q,m bank, and checks the all-m formula's explicit witnesses at large parameters. The proof, not the finite bank, covers every q and v.

The Farey construction works because the core includes every speed1,...,q-1. Arbitrary sources have no such full Farey tent envelope. Thus these results close an infinite structured family and explain its prime-controlled tightness, while leaving the general LRC entry and universal mass-bound problems open.

## 6. Inheritance: what is recovered and what is closed

The family and its prime arithmetic are not new discoveries in this session.

- [THM-709](../../../../01-canon/theorems/THM-709-doubling-tight-locus-singleton.md) already records the Goddyn–Wong acceleration criterion. In its notation, replacing q in {1,...,q+1} by mq is tight exactly when q has a common factor with every integer in [2,2m-1]. For m>=2 this is exactly the primorial tightness condition above. We recover this inherited criterion rather than claim it as new.
- [THM-1230](../../../../01-canon/theorems/THM-1230-n14-gap-is-populated.md) records the fixed-q12 ladder and leaves its all-m proof open beyond the verified bank.
- The stored `lrc14_ladder_realization_crossN_kps_S128c86.out` contains exact F-ladder tables and the binder values2,3,5 across N. This session explains those binder values uniformly as s(q).
- [THM-1286](../../../../01-canon/theorems/THM-1286-D-graded-gate-tower-primorial-cascade.md) already proposes the all-D primorial gate. [THM-1271](../../../../01-canon/theorems/THM-1271-rung-floor-seals-and-D7-rung.md) proves several reductions and per-instance certificates, but explicitly leaves the uniform e-channel maximum/gate completeness open. [THM-1276](../../../../01-canon/theorems/THM-1276-D10-rung-picowindow-and-D9-off-diagonal.md) gives large exact checks through D=10. In particular6/1271,7/16183,9/270287, and10/5105119 are inherited values, not new numerical discoveries here.

The all-q Farey proof supplies the missing uniform upper bound and exact formula, together with a complete top-window description, exact strict-branch maximizing times, and its terminal mass derivative. It closes those scoped uniform family obligations without claiming a classification of arbitrary sources. No Lean build or formal-kernel verification of this new uniform proof was performed in this session.

### The uniform slack-one gate, now a corollary

In the inherited notation N=q+1 and D=m, put p=2D-1. For D>=2,

 M(F_D(N))=D/((N+1)D-1)

if and only if p is prime, every prime smaller than p divides N-1, and p does not divide N-1. Indeed, the desired value lies strictly above the floor, so the exact spectrum forces s(N-1)=2D-1. Since s is the least missing prime, this is precisely the stated gate. Conversely the gate sets s=p and the formula gives the desired value. This proves the fitted THM-1286 gate without a per-N candidate scan; composite p cannot be the least missing prime.

### Arbitrarily small normalized positive excess

For every odd prime p>=5 choose q as the product of all primes below p, and m=(p+1)/2. Then s(q)=p=2m-1, and the proved exact clearance is

 M=m/(m(q+2)-1)>1/(q+2).

Writing k=q+1 for the number of relative speeds, the excess satisfies

 M-1/(k+1)=1/((k+1)(m(k+1)-1)),
 k^2 (M-1/(k+1))<1/m ->0 as p increases.

Thus there is no universal positive constant c such that **every positive excess above the conjectured floor**, over all runner counts, is at least c/k^2. This conclusion uses explicit above-floor families and does not assume unrestricted LRC. If a positive-side spectral gap is defined as an infimum over such families, it likewise has no uniform positive k^2 lower bound. We do not assert that this infimum is attained, nor use the ambiguous phrase 'non-tight' to include hypothetical below-floor counterexamples.

The old saturation language in THM-539 was already superseded by the later exact rungs. The advance here is the uniform proof for arbitrarily high primes and the exact consequence, not rediscovery of a previously verified rung.
