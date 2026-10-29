# The same unit sees 11+1=13-1=12, with different signs

**PROVED** U1-U2; finite regressions through prime exponent4. This addresses
the annotated g58 claim. [Exact script](check_crt.py), [JSON](check_crt.json).

## U1. Split versus nonsplit clocks

Use the actual Ramanujan unit u=(5+sqrt(29))/2 in O=Z[u],

    u^2=5u+1,  Norm(u)=-1,
    u^12=13741001+71351280u.                             (U1)

At p=13, discriminant29 is a square (3=4^2 mod13). The two roots are11
and7; both have order12 in F13*. Thus u^12=1 in O/13, and the primitive
unit has order12. At p=11, discriminant29=7 is nonsquare, so O/11=F121.
Frobenius sends u to its conjugate -u^-1, hence u^12=-1. Checking proper
prime divisors of24 gives exact order24.

So the **same** unit has the two clocks

    split p=13:   p-1=12,   u^12=+1;
    inert p=11:   p+1=12,   u^12=-1.                    (U2)

This is a field/norm mechanism, not an inference from neighbouring integers.
It explains why the twelfth-power trace is constant at13 and alternating
at11. The inverse and conjugation operations stay in their named quadratic
ring; these cannot be replaced by arbitrary bit permutations.

The factorization in the inherited
[trace-carrier theorem](../cursor_cm_phase_20261028/TRACE_CARRIER.md) also has
11^2, which is significant. In the basis(1,u),

    u^12=-1 mod121,  but u^12 != -1 mod1331,
    u^12=+1 mod13,   but u^12 != +1 mod169.

The exact all-h>=1 lifting laws are

    ord_(O/11^h)(u)=24*11^max(0,h-2),
    ord_(O/13^h)(u)=12*13^(h-1).                       (U3)

Proof: the two coefficients of u^24-1 have common 11-adic valuation2;
those of u^12-1 have common13-adic valuation1, directly from(U1). For an
odd prime p and integral ring element v=1+p^s z with z nonzero modp,
the binomial expansion gives v^p=1+p^(s+1) z' with z'=z modp. Induction
raises the valuation exactly once at each p-power. Reduction modulo p
forces every other period to be a multiple of the base period. This proves
(U3), rather than extrapolating it from the four checked exponents.

## U2. What the sampled trace forgets

Source states are **all** x in O/m, not just powers of u. The named operations
are multiplying by u^12 and its inverse, and conjugation; observation is Tr(x).
The inherited exact Smith theorem has first invariant1 and determinant

    29*71351280 = 2^4*3^4*5*7*11^2*13*29.

Consequently each sampled-trace history has fiber gcd(m,29*71351280).
At m=11^h this fiber has size11^min(h,2), and at m=13^h it has size13.
Adding the primitive clock u changes the determinant to29. At either prime
and every h it recovers the entire m^2-state ring.

For m=143, CRT combines the two signs: u^12=(-1,+1). The sampled carrier
has143 states from20449 ring states; the primitive carrier has20449. This
is a concrete example of a common sampled clock concealing different local
phases. It is not the Lucas golden ring from(C1), whose defining polynomial
is X^2-X-1, whereas this one is X^2-5X-1.

## Two different Zagier contexts, with their parameters retained

The annotated Ramanujan setting connects directly to the classical identity

    t=(eta(tau)/eta(2tau))^24,
    T_2A=t+4096/t+24=q^-1+4372q+96256q^2+...,
    64(g_n^24+g_n^-24)+24=T_2A(q).

This is the level-two Fricke formula treated in Zagier's *Traces of singular
moduli*, section8; [Matsusaka's primary paper](https://link.springer.com/article/10.1007/s40993-017-0090-x)
also gives the eta expression. The [prior Ramanujan note](../cursor_cm_phase_20261028/RAMANUJAN.md)
fixes q=e^(-pi sqrt(n)) and all normalizations. Here g58^2=u, so g58^24=u^12:
the sampled exponent is actually supplied by that formula. The finite-field
signs(U2) explain the resulting residue losses. They do not explain or prove
the modular transformation exponent24 itself.

A second search recovers the quite different
[Zagier six-triples note](../zagier_six_triples_20261018c.md). Its recurrence is

    (n+1)^2 v_(n+1)=[a n(n+1)+b]v_n-c n^2 v_(n-1).

The Apéry row is(a,b,c)=(11,3,-1); another sporadic row is(12,4,32).
The accessory b is lost by the singular-location polynomial1-ax+cx^2.
That note gives explicit nonintegral hostile substitutions. The connection
to this packet is the **need to retain the coefficient that future recursion
tests**, just as(W6) retains cross-cycle weights and(C3) retains relative
phase. No isomorphism between these recurrences or their moduli is asserted.

The prime-Lucas reflection rule supplies yet another typed12: n=11 mod12
selects the odd reflection and S_(L_n), versus A_(L_n) in the other eligible
residues. Its mechanism is Fibonacci/Lucas parity modulo4, not(U2). These
are compatible examples of clock-sensitive quotients, with separate sources.
