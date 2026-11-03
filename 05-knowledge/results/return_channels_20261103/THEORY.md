# Error transport, a common ternary family, and the irreversible parity boundary

Status: **PROVED** for the stated65-symbol operations; finite certificates
are labelled separately. No uniform positive answer for every required scalar
is asserted. Words are in processing order, reverse source-SCC order.

## T1. Inherited objects

Use row states v=(h,x,y), literal graphs and matrices P_l from
[inverse65/INVERSE I1-I2](../inverse65_zeckendorf_20261103/INVERSE.md).
Appending l sends v to vP_l. Equal observer states certify all future same-end
H responses, not equality of graphs or original vertex order. For a signed
row delta, a common append transports a difference by delta->delta P_l.

The inherited actual processing seeds are

    U=r5 a a,          v(U)=(31,14,25),
    U2=a a q1 a a,     v(U2)=(95,41,75),
    U4=q0 a a a a a,   v(U4)=(97,40,75).

With t=(2,-1,0), these satisfy v(U2)=3v(U)+t and v(U4)=3v(U)+2t.
Their construction and failure of naive repeated affine lifting are credited
to [golden_inverse_defects/TRANSPORT T1-T5](../golden_inverse_defects_20261103/TRANSPORT.md).
They are source replacements, not retained-edge vertex attachments.

## T2. A concrete bridge into the norm13 channel

Let W=r43 a r13, with actual vertex cost11. Direct multiplication gives

    P_W = [[1165,692,965],[2324,1386,1932],[2016,1204,1678]],
    t P_W = 2r,            r=(3,-1,-1).                  (T1)

Both r21 and r38 have the same observer matrix R and satisfy rR=r. The
terminal triangle b has `(rP_b)_h=-5`. Therefore, for EVERY word w over
{r21,r38},

    H(U2 W w b)=3H(U W w b)-10,
    H(U4 W w b)=3H(U W w b)-20.                        (T2)

This is the missing common seed for the two negative-offset families of the
previous session. It lies outside that session's H<=2000 experiment, so it
does not contradict the recorded finite absence.

Before the terminal b the common base state is

    v0=(119051,70956,98913),  J=2x+2y-3h=-17415.

Write M_k=H(U W r21^k b). Vertex clocks for the base and companions are
21+5k,22+5k,23+5k. The literal r21/r38 choices must remain in the source
word even though they have the same state matrix. The entire source has a
preserved Hamilton spine, and its endpoint cross edge reverses under K;
this closes the spine to a Hamilton cycle and proves every target strong.

The inherited conserved-charge law gives

    M_(k+2)=40M_(k+1)-13M_k+174150,
    M0=798527, M1=31810907,
    sum M_k z^k=(798527-928700z+304323z^2)
                   /((1-z)(1-40z+13z^2)).             (T3)

The generating function follows by multiplying the recurrence and summing;
its coefficients are H values of this family, not class counts.

Since M0=M1=2 mod7 and the recurrence preserves that value,

    Q_k=(3M_k-20)/7

is integral for all k and obeys

    Q_(k+2)=40Q_(k+1)-13Q_k+74710,
    Q0=342223, Q1=13633243,
    sum Q_k z^k=(342223-397900z+130387z^2)
                   /((1-z)(1-40z+13z^2)).             (T4)

Equivalently the second companion's conserved charge is
`3J-52=-52297=-7*7471`, so its b attachment has count7(h-14942).
The first quotients are342223,13633243,540955531,21461063791.
The first is7*48889, a composite hostile. Trial division through146495
certifies21461063791 prime; through893 it certifies798527 prime. These are
individual certificates, not a prime-producing recurrence theorem.

The new source phase is also certified. For Q put r0=-7471 and
`Z_k=13Q_k-5r0`; then `Z_(k+2)=40Z_(k+1)-13Z_k`. Z0 is a unit modulo13
and `(27Z0-Z1)/13` is a unit modulo13. These are exactly the hypotheses of
[inverse65/CHARGES C4](../inverse65_zeckendorf_20261103/CHARGES.md): the pair
clock has tail e and period13^e modulo13^e, for every e>=1, with compatible
CRT with all coprime moduli. The source phase is recomputed rather than
inherited from the old seeds. Exact checks through13^3 and at91,1183 agree.

Modulo7 the Q clock has period6 and is zero exactly at k=0 mod6. Modulo13
it has tail1 and is zero exactly at k=1 mod13 (k>=1). Since all Q_k exceed13,
these are certified composite terms. Their union excludes exactly18 of the
78 index classes; the other60 classes are only candidates, not certified
primes. k=3 passes these guards and has the separate trial-division certificate.
k=2 also passes, but Q2=29*1327*14057 is composite; this is the hostile to
treating the residue selector as a primality test.

**Independent graph check:** a frontier DP reads only the literal adjacency
rows, retaining endpoints of partial path components. It agrees with the
matrix counts for k=0,1,2,3 on all three families, up to38 vertices. Backward
width is4 and the observed peak is338 frontier states. The algorithm is
inherited from the earlier path-cover transfer, not a new claim of generic
polynomial-time Hamilton counting. All65 atoms are independently checked by
permutation enumeration, and the original seeds also by subset DP.

## T3. Primitive errors cross a one-way parity boundary

Every65 matrix has first column congruent to(1,0,0) modulo2. Let delta be a
primitive integer row with odd first coordinate. For every word w, delta P_w
still has odd first coordinate. Its coordinate gcd is odd, so dividing by
that gcd leaves a primitive row whose first coordinate is odd.

Consequently no word sends r=(3,-1,-1) to a nonzero integer multiple of
t=(2,-1,0), or of q=(2,-2,1). The same statement holds for the odd-first
directions(3,-2,1) and(3,-2,0). This is an all-length result, not a failed
bounded search. Scaling an actual correction by2 does not evade the proof:
divide its initial content first, and apply the same integer equation.

More generally, content(gcd of coordinates) can only gain divisors under an
integer matrix. Once the2-adic valuation of the first coordinate equals that
of the content, the equality persists for this bank. The bridge(T1) changes
t into2r and crosses precisely this absorbing boundary. Any return needs a
different operation, such as an explicitly justified source replacement or
re-marking; these are not excluded by the append-only theorem.

The four next affine states from the inherited U2/U4 constructions remain
absent after multiplying by P_W:

    (1071483,638596,890209), (1071489,638594,890207),
    (1071501,638590,890203), (1071507,638588,890201).

The complete state inverse rejects each. Thus(T2) provides an infinite common
family but not closure under another round of the same ternary state updates.

## T4. Native loops with nonunit error gains

There is another route that stays in the even-first channel:

    t P_(r24 a a)=2q,      q P_a=-t.

Four explicit loop words return t with the following gains and characteristic
polynomials; all symbols are in the original bank:

| Loop | Word | Vertex cost | Gain g | Characteristic polynomial |
|---|---|---:|---:|---|
| L2 |r24 a a a|8|-2|(X+2)(X^2-367X-30)|
| L6 |r24 q4|9|6|(X-6)(X^2-1133X+90)|
| L45 |r19 r32 r48 q4|19|45|(X-45)(X^2-1581580X+5985)|
| L13 |r19 r19 a a r18|17|13|(X-13)(X^2-348170X+81225)|

For any source v, H(v L^k) satisfies the displayed cubic recurrence by
Cayley-Hamilton. For the three inherited seeds, the companion errors are
exactly2g^k and4g^k. More generally a mixed sequence multiplies these gains
in order, while its full3x3 matrix product controls the base H.

**FINITE-EXACT:** every return to the line Zt using at most five bank letters
has gain in

    {-45,-36,-6,-2,-1,1,2,6,13,36,45}.

The exhaustive search merges identical signed error rows at each length.
It discards one-sided rows because all matrices are nonnegative: such a row
can never return to a mixed-sign direction. This proves completeness of the
finite census. Longer words may introduce additional gains and primes.

The gain and vertex clock do not determine the base: U L2 L6 and U L6 L2
both add17 vertices and have error gain-12, but their H values are7450383
and12892021. The missing data are the noncommuting bulk state and actual word.
