# Short physical slope transport and an exact response-fiber normal form

Session 200. **PROOF CANDIDATE, pending independent audit.** This addendum
gives an independent short proof of odd-numerator response coverage, using
only the previously proved even-numerator theorem. The transfer costs
O(log b) extra vertices and preserves the denominator. It also completely
identifies the physical 5U-9D=24 fiber with the old 8U-9D=29 fiber. It does
not establish prescribed positive Hamilton-path count coverage.

## 1. A short block moves any interior slope to an even numerator

Use actual nontrivial start-one source profiles (D,E,U,d), and chronological
O/Z operations as in `dependencies/POSITIVE_SOURCE_GRAMMAR.md`. A unit
source has d=1; a strong unit source has E=1. Direct iteration gives, for
ell,j>=1 and arbitrary input,

    U(s.O^ell Z^j)=2^(ell+j) U-(2^j-1)d,
    D(s.O^ell Z^j)=D+(2^j-1)2^ell U-(2^j-2)d,
    E(s.O^ell Z^j)=d.

The initial O erases the input E. Hence for positive coprime a,b,

    aU(s.O^ell Z^j)-bD(s.O^ell Z^j)
      =a' U(s)-bD(s)+c d(s),                     (1)
    a'=2^ell[b-2^j(b-a)],
    c=2^j(b-a)+a-2b.

This is an exact physical source transformation, not rational rescaling of
an abstract profile.

**Theorem 1.** Suppose b is odd and b/2<a<b. There are positive ell,j with
ell+j<=2ceil(log2 b) such that a' in (1) is even, coprime to b, and
b/2<a'<b. The same fixed physical block therefore transfers the entire
even-numerator response image to the response image at a/b.

Choose j minimally with

    t=b-2^j(b-a)<b/2.

The initial a is above b/2, so j>=1. The previous value was above b/2,
because equality would equate an integer and half an odd integer; applying
t_next=2t_old-b shows t>0. Also t is odd. Choose ell minimally with
2^ell t>b/2. Then ell>=1 and b/2<a'=2^ell t<b, with strictness again from
odd b. Furthermore gcd(t,b)=gcd(2^j a,b)=1, so gcd(a',b)=1.
Each of j and ell is at most ceil(log2 b), giving the claimed vertex cost.

For unit inputs c is an ordinary integer. Since j>=1, c has parity a.
Therefore if r has parity a-b, r-c is odd. The old even-numerator theorem
constructs infinitely many strong unit sources with a'U-bD=r-c.
Apply the fixed block and use (1) to obtain aU-bD=r. Both positive
coordinates increase and the output retains d=E=1 and strongness, so
infinitely many distinct inputs give infinitely many distinct outputs.
The map on unit (D,U) is injective because U has a nonzero scaling factor
and D has coefficient one. This independently proves the new parity-class
surjectivity without extending the old cyclic-clock formula.

Examples are

    slope 3/5:  O^2 Z,  effective numerator 4, c=-3;
    slope 5/7:  O Z,    effective numerator 6, c=-5;
    slope 5/9:  O^3 Z,  effective numerator 8, c=-5;
    slope 7/9:  O^3 Z^2,effective numerator 8, c=-3.

The theorem concerns the signed response. D and U are not individually
preserved, and the positive graft readout generally changes. The additive
shift includes the genuine deletion sidecar d; arbitrary denominators,
fractional sources, and free assignment of integer coordinates are not used.

## 2. Every member of one odd-numerator fiber has this suffix

Write g(1)=1, g(2h)=g(h), and g(2h-1)=g(h)+h-(h mod 2). The physical unit
source profiles are exactly (g(n),n mod 2,n,1), n>=2. To justify completeness
from the full source grammar: any J retaining d=1 must occur when D=1,
at a transitive source (1,0,2^r,1). Its output (2^r+1,1,2^r+1,1) is also
the profile of t.Z^r. After D>1 no J can retain d=1. Replacing this possible
J segment leaves an O/Z recipe. Odd n>=3 give exactly the strong unit
profiles.

**Theorem 2 (exact fiber transport).** The physical strong unit profiles
satisfying 5U-9D=24 are exactly the profiles q.O^3 Z for strong unit
sources q satisfying 8U_q-9D_q=29. The predecessor profile is unique.

It suffices to prove the arithmetic equivalence

    n>=3 odd, 5n-9g(n)=24
    iff n=16q-1, q>=3 odd, 8q-9g(q)=29.           (2)

Write n=2m-1. The g recurrence turns the left equation into

    m-9g(m)+9(m mod 2)=29.

If m is odd, then m>=3 and the strong cone gives g(m)>=m/2+1,
making the left side at most -7m/2<0. Hence m is even. Write m=2^s q,
s>=1 and q odd. Since g(m)=g(q), the equation is

    2^s q-9g(q)=29.

The case q=1 would require 2^s=38, impossible. Thus q>=3. If s<=2,
g(q)>=q/2+1 makes the left side negative. If s>=4, using g(q)<=q yields
29>=(2^s-9)q>=7q, so q=3. But g(3)=3 then makes the left side divisible
by three, unlike 29. Therefore s=3, proving the forward direction of (2).

Conversely g(16q-1)=g(q)+8q, so the right side of (2) gives
5(16q-1)-9g(16q-1)=8q-9g(q)-5=24. The physical operation O^3 Z has exactly
this profile map. Its inverse q=(n+1)/16 is unique, proving the theorem
at the profile level, without claiming a unique graph isomorphism class.

The old even-numerator theorem at 8/9 supplies infinitely many actual q,
so this is a complete nonempty fiber description, not only a necessary
congruence. For example

    q=(600339,1,675385,1), 8U_q-9D_q=29,
    q.O^3 Z=(6003419,1,10806159,1).

Both are certified by explicit O/Z words and the integer g recurrence.
Taking left source v=t.ZJOZO=(17,0,42,3), the two C2 presentations
N(v,J(q.O^3 Z)) and N(Jv,q.O^3 Z) share the positive target 2972894005,
as verified in the preceding packet. No new spectrum gap is claimed filled.

## 3. Exactly what remains for a prescribed positive target

Theorem 2 refines the previous parameterization. Write

    U_q=18j+7, D_q=16j+3, j>=0.

Every physical member of the 8/9 response-29 fiber has this form by
integrality and parity. Its physical membership is exactly

    g(18j+7)=16j+3.

For the displayed left source v, the common C2 target becomes

    n=79232j+30133.                              (3)

Thus this whole repaired family occupies one explicit residue class,
and for a prescribed n the index j in (3) is fixed. The response theorem
proves infinitely many accepted j but not every j. For example j=0 gives
g(7)=5, not 3. A fiber's being nonempty or Diophantine-definable cannot
replace this exact physical membership condition.

The useful transfer is from one witnessed arithmetic fiber to a second,
with a literal short graph operation, an integral inverse on this specific
fiber, and the deletion sidecar kept. It is a concrete constructive use of
the proposed arithmetic themes, not an appeal to a conjecture about totient
multiplicities or integer rings.

## 4. Verification

Run `python3 verify.py` and `python3 -O verify.py`; receipts must agree.
The checker verifies block identities on the full four-coordinate basis,
all coprime interior slopes with odd b<=101, parity and word cost, literal
small examples, both profile directions on a declared bounded universe,
the large exact witness independently through g, and the excluded j=0
control. Dependencies are frozen and hashed. The all-order conclusions
come from the proofs, not the bounded checks.
