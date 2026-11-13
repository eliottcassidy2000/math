# The physical norm chamber has no Pell return

**PROVED:** the cone, unit-orbit exclusion, prime-target uniqueness and
inert-prime obstruction below. **FINITE-EXACT:** the accompanying small
controls. The tournament norm constructor itself is inherited from
[MIRROR_MODULE](../../sixty_three_bridge_20261113/MIRROR_MODULE.md).
The polarized replacement is developed by the root agent in the current
packet. Chan–Pak and ordinary Collatz remain **OPEN**.

## 1. Inheritance and the precise source

The closest older mechanism is the
[Gold/Pell half-step,§§1–2](../../gold_pell_collatz_halfstep_20261019.md):
an integral square root of a positive parameter recurrence leaves the
positive-dimensional chamber; its square returns. Its explicit norm
coordinate is X²−5Y², but the positive-dimension entry condition cannot
be discarded. The hostile here is even stronger: no nonidentity Pell
unit returns to the tournament's physical parameter cone at all.

The corrected near miss is to treat an infinite Pell orbit of one target
norm as infinitely many alternative module realizations. The missing
sidecar is an actual nonempty tournament M with BOTH observations

    h=H(M), p=P_M(2),

where P_M(2) counts ordered two-slot spanning path covers with empty slots
allowed. The inherited mirrored construction, followed by NR, yields

    T=12p²+8hp−h²,
    Y=3T+2h²=(6p+2h)²−5h²,                             (1)

and balanced rooted state(Y,T,T,Y). Its readiness defect before NR is
h², not necessarily1. Put x=6p+2h.

Every actual nonempty M has

    h positive odd, p positive even,
    p≥(|M|+1)h≥2h,
    x≥14h, x≡2h (mod12).                              (2)

Oddness of h is inherited Rédei parity. Slot exchange is a fixed-point
free involution on the two-slot covers, proving p even. Cutting each HP
at one of its |M|+1 positions injects into those covers: concatenation
recovers the HP and the first-slot size recovers the cut. This proves
the inequality.

The empty module is a separate boundary: h=p=1 gives Y=59 and x/h=8.
It is not included in the nonempty cone(2). No order or cone bound below
silently counts it as a singleton.

## 2. Every nontrivial norm-preserving Pell unit leaves the cone

Let u,v be integers with u²−5v²=1, v≠0. Multiplication by u+v√5 acts by

    (x,h) ↦ (ux+5vh, vx+uh).                           (3)

Allow also an overall sign and conjugation, so a transformed pair may be
returned to positive coordinate signs when possible. If x≥14h>0, every
such nontrivial transformation, after taking these absolute signs, has

    |x'|/|h'| < 3.                                    (4)

In particular it cannot represent a physical pair in(2).

**Proof.** Change the overall sign so u>0 and write w=|v|. Testing
w=1,2,3 in u²=5w²+1 gives no square, hence w≥4 and

    √5 < u/w ≤9/4.

Write r=x/h≥14. For v>0 the new slope is
(ur+5w)/(wr+u), which is less than3 because
(3w−u)r+3u−5w>0. For v<0, ux−5wh>0 and −wx+uh<0. The positive absolute
slope is(ur−5w)/(wr−u), again less than3 because

    (3w−u)r+5w−3u
      ≥(3/4)wr−(7/4)w >0.

The denominator is positive since r≥14>u/w. These inequalities prove(4).
The same argument works for every integral norm-one unit directly; it
does not assume a finite list of unit powers. ∎

The larger quadratic integer ring does not evade this parity guard.
A half-integral norm-one unit has the form(a+b√5)/2 with a,b odd. Since
x is even and h odd, both transformed numerators ax+5bh and bx+ah are
odd, so the coordinates are not integers. If a,b are even the preceding
integral case applies. Norm-minus-one units change Y to−Y and do not
preserve the positive target.

Thus each norm-unit orbit, including sign and conjugation, contains at
most one positive physical pair. Different unit orbits may share a norm;
this theorem does not assert global uniqueness for composite targets.
Unlike the old Gold half-step, taking a second or later unit step cannot
repair the physical chamber violation.

## 3. A prime target has at most one arithmetic candidate

If Y is prime, there is at most one positive integral pair(x,h) satisfying

    x²−5h²=Y, x≥14h.                                  (5)

**Proof.** The cone gives h²≤Y/191 and x²≤196Y/191. Suppose two pairs
(x1,h1),(x2,h2) satisfy(5). Then

    0<x1h2+x2h1≤28Y/191<Y.

Modulo Y,

    (x1h2−x2h1)(x1h2+x2h1)=0.

Primality and the strict sum bound force x1h2−x2h1≡0 modY.
Its absolute value is also less thanY, so it is zero. The positive pairs
are proportional; substituting into their common nonzero norm makes the
proportionality factor1. ∎

This proof uses neither a class-number theorem nor a claim that every
norm representation belongs to one Pell orbit. It turns a prime-target
question into one physical interface test, if a candidate exists.

Primality cannot be dropped from this arithmetic claim: Y=604679=29²·719
has both(h,p)=(11,126) and(37,118), with x=778 and782 respectively.
Both obey the arithmetic cone and congruences. This does not assert that
both pairs come from tournaments: the latter has order bound2 and is
therefore impossible at H=37. The final physical test carries real content.

For any fixed positive targetY, a finite arithmetic candidate generator is:

1. Enumerate positive odd h≤floor(sqrt(Y/191)).
2. Require Y+5h² to be a square x².
3. Require x−2h divisible by12, and put p=(x−2h)/6.
4. Require p≥2h.

These are necessary arithmetic conditions, not sufficient module
realizability. Any actual module for a surviving pair has

    |M|≤floor(p/h−1).                                 (6)

It can therefore be checked in a finite tournament bank, with its actual
(H,P2) pair retained. For primeY, at most one pair reaches that check.
The boundary Y=191 has the singleton pair(h,p)=(1,2), so the numerical
cone constant is attained by an actual module.

The inherited second bound p≥2^|M| sharpens(6) to
|M|≤min(floor(p/h−1),floor(log2 p)); see
[the complete inverse audit,§1](../../target_directed_cores_20261112/ports/COMPLETE_INVERSE_AUDIT.md).
Its elementary proof is that every one of the2^|M| vertex subsets has
at least one HP on itself and on its complement, so every summand in
P_M(2)=Σ_S H(M[S])H(M[V\S]) is at least1.

## 4. Inert primes obstruct the norm chart, including the7p lane

Let ℓ be an odd prime for which5 is a quadratic nonresidue moduloℓ.
For any nonzero integers x,h with Y=x²−5h²≠0,

    v_ℓ(Y)=2 min(v_ℓ(x),v_ℓ(h)).                       (7)

Indeed, divide out the common power ofℓ. If the remaining norm vanished
moduloℓ, its two coordinates would both vanish: when the second is a
unit this would make5 a square, and when it is zero the first is forced
zero as well. This contradicts the normalization and proves(7).

Forℓ=7,5 is a nonresidue. Hence

    7|Y ⇒49|Y.                                         (8)

No value7q with q a prime other than7 is in this norm chart, regardless
of how many module types are supplied. This is a constructor obstruction,
not a tournament nonattainment claim.

There is a stronger conclusion for a primitive balanced output in(1).
If gcd(Y,T)=1, NO inert prime can divideY. Forℓ≠3, equation(7) would
giveℓ|h and thenℓ|T=(Y−2h²)/3. Forℓ=3, equation(1) shows directly
that3|h forces3|T. The same contradiction applies. This includes the
prime3 because5 is a nonresidue modulo3.

Finally, p even and h odd giveY≡7(mod8). Primitivity gives3∤h and thus
Y≡2(mod3). If Y itself is prime, it exceeds5; the identity Y≡x²(mod5)
then givesY≡1 or4(mod5). Combining these congruences gives the narrow necessary condition

    Y≡71 or119 (mod120).                              (9)

It is only necessary: actual (H,P2) realizability remains after every
congruence and norm test.

## 5. The repaired construction releases the equality of the two modules

The root agent's polarized version uses independent nonempty modules M,N
with data(h,p),(j,q). In the same prepared seven-block skeleton its count
and NR target are

    Z=12pq+2hq+6pj−hj,
    Y=3Z+2hj=36pq+6hq+18pj−hj
      =(6p+h)(6q+3j)−4hj.                              (10)

This is an ACTUAL change of the module interface, not an application of a
formal Pell unit to(1). M=TT2,N=TT1 has(h,p,j,q)=(1,4,1,2), so

    Z=123, Y=371=7·53,

with primitive balanced state(371,123,123,371). Both7 and53 occur to
odd valuation, making the inert obstruction visibly inapplicable to the
new bilinear form. The incoming companion-word/reflection construction
already realizes this particular count; the interpretation here explains
which equality restriction had prevented it in the symmetric chart.

The three-apex-arc change between the unprepared and prepared skeletons
has exact count drift

    4(pj−hq).                                          (11)

Thus independent modules do not usually give a count-neutral surgery.
The determinant of the two observation vectors is the missing coordinate.
Neutrality holds exactly when p/h=q/j, not only when N is M's converse.
For instance TT2 and C3 have proportional observations(1,4),(3,12),
giving a same-count669 surgery with readiness defect0→3. This does not
become unit readiness; the defect is hj. Equation(11) connects to the
ordered-flank determinant in the old Farey/Stern theory, but coordinate
addition alone is not a tournament construction.

The resulting exact-target task is now a bilinear factor/physical-module
selector, developed in the parent packet. It escapes a proven obstruction
of(1); it does not prove coverage of all targets or ordinary Collatz.
