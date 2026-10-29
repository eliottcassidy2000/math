# A second full recursion: Fricke coefficients with their level retained

**CITED** eta expressions for the modular functions; **PROVED** elementary
coefficient recurrences; **FINITE-EXACT** verification through q^40 at each
level by [check_eta.py](check_eta.py), with [JSON](check_eta.json).

For p in{2,3,5,7,13}, put h=24/(p-1) and

    P_p(q)=product_(d>=1,p does not divide d)(1-q^d)^h,
    t_p=q^-1 P_p(q),   J_p=t_p+p^(h/2)/t_p+h.           (E1)

These J_p are the classical Fricke functions; the normalization is given in
[Matsusaka, introduction](https://link.springer.com/article/10.1007/s40993-017-0090-x).
At p=2 this is the Zagier/Ramanujan function already used in the prior packet.
The following recursions do not require a numerical evaluation of eta.

Write P_p=sum c_n q^n and P_p^-1=sum d_n q^n, with c_0=d_0=1. Define

    B_p(j)=sum_(a|j,p does not divide a) a
          =sigma_1(j)-p sigma_1(j/p),

where sigma_1 of a noninteger is0. Logarithmic differentiation gives

    n c_n = -h sum_(j=1..n) B_p(j)c_(n-j),
    n d_n = +h sum_(j=1..n) B_p(j)d_(n-j).              (E2)

Therefore J_p=q^-1+sum_(n>=1) j_p(n)q^n, with the all-order recursion

    j_p(n)=c_(n+1)+p^(h/2)d_(n-1).                     (E3)

The constant term is c_1+h=0. Equations(E2) follow by equating coefficients
in qP'/P=-h sum_j B_p(j)q^j; the inverse has the opposite sign. They keep
the prime level in both the divisor filter and reciprocal scale.

| p | h | reciprocal scale | first positive coefficient |
|---|---:|---:|---:|
|2|24|4096|4372|
|3|12|729|783|
|5|6|125|134|
|7|4|49|51|
|13|2|13|12|

Thus at level13 the actual recurrence starts

    J_13=q^-1+12q+28q^2+66q^3+132q^4+... .

Here the12 is exactly `c_2+13=-1+13`; its source is distinct from the
order12 unit clock modulo13 in [CLOCKS](CLOCKS.md). Retaining level and
normalization lets the two appearances be compared without conflation.
The supplied list includes p=13 but not11; substituting11 into this table
is not justified by the stated modular-function formula. Nothing here
excludes other modular descriptions at level11.

This gives a useful comparison with the A000568 recursion(W6). Both use
an inventory of prior pieces tested by a new piece: gcd interaction weights
for tournament cycles, filtered divisor weights for eta products. The
logarithmic derivative closes linearly on coefficients, whereas(W6) must
update a whole vector of cycle weights. These are explicitly different
composition laws; a common recurrence framework is not an isomorphism
of their counted objects.
