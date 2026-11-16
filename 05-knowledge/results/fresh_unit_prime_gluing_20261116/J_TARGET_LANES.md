# An exact target-lane change inside literal prime gluing

Use the [two-port gluing](PRIME_GLUE.md), with shared modules
M,N, u=H(M), v=H(N), and R=P_M(2), T=P_N(2) counting optional
ordered two-slot covers. Its three counts satisfy

    H=p+q−uv,
    p−2uv=(4R−u)(2T−v),
    q−2uv=(8R−u)(T−v).

## Why this particular lane cannot supply the prime part of Chan–Pak

For any tournament M, the odd-cycle slot polynomial gives
R≡(−1)^|M|u modulo 3, since 2≡−1 and all exponents have the
same parity as the order. Likewise for N. Write e,f for the
two parity signs. Then

    p/(uv) ≡ 2+(e−1)(−f−1)  (mod 3, when 3 does not divide uv),
    q/(uv) ≡ 2+(−e−1)(f−1)  (mod 3).

Both are automatically divisible by 3 if uv is. Hence if p,q
are primes at least 11, then u,v are units mod3 and the two
orders have the **same parity**: odd/even makes p divisible by
3, and even/odd makes q divisible by 3. With the same parity,
p≡q≡2uv mod3, so H≡0 mod3. In particular even the enlarged
richer-module version does not give a prime host (nor 7p
for p≠3) using prime inputs. This is a proved limitation
of that template, not a primality heuristic or a claim about
J in general.

For its singleton-first target line H=29T−19v, the same
calculation shows an odd-order right module always gives a
multiple of 3. Also T≥(n+1)v at order n, by splitting each
Hamilton path at its n−1 possible cuts and counting the two
empty-slot choices. Thus the inverse must pick v≤H/(29n+10)
in the correct congruence. The equation is a specific
descent requirement on a realized **pair** (v,T), not a
licence to choose the two coordinates independently.

## An asymmetric template genuinely escapes the mod-3 lock

Using the complete small-order atlas, take the order-eight
tournament with outgoing rows

    (2,8,11,97,15,87,151,63).

Let A delete vertex 2 and B delete vertex 7. Inflate their
shared vertex 4 into a transitive tournament of length t≥1,
identically in all three graphs. The numbers of projected
Hamilton words using the shared port once, twice, three
times, or four times are

| graph | one | two | three | four |
|--|--:|--:|--:|--:|
| host |79|36|0|0|
| A |43|18|0|0|
| B |37|18|0|0|

A t-site transitive port fills two distinct ordered nonempty
runs in 2^t−2 ways. Thus for all t≥1, with X=2^t,

    H=36X+7,  p=18X+7,  q=18X+1,

and exact J holds. For t=2, this is the prime H=151
from primes 79 and 73; for t=10 it is the prime H=36871
from primes 18439 and 18433. Trial division certifies the
six summand/host claims. All are literal induced-gluing
witnesses, however this exponential ray has O(log X)
outputs up to X and there is no assertion that its prime
triples continue indefinitely.

More generally the unchanged skeleton certificate reads an
arbitrary port U of count v and optional two-slot count T
as (H,p,q)=(36T+7v,18T+7v,18T+v), so the offset changes to
−v. Both summands and host are units mod3 whenever v is.
This gives a different physical target line. With v=1 the
offset is exactly −1, but the module is transitive and
T is only dyadic. The next exact selector must retain the
unit offset while changing the available factors or must
compensate for a new v.

The [optional census](search_j_ports.py) scans both small atlases:
29 prime-input deletion pairs at order seven and 93 at order
eight satisfy J. There are respectively 52 and 174 choices
of an inflation port preserving J to all sizes. Among them
4 and 16 admit an additional prime-input prime or 7p host
at inflation size 2 through 12. The count is by representative
and choice of marked vertices, not up to marked isomorphism.
The displayed asymmetric witness and its coefficients have
a separate replay in check_j_escape.py.
