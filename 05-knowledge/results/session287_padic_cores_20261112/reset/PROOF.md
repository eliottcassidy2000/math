# A single physical reset, two non-Archimedean boundary geometries

**Status: PROVED scoped compactification and finite physical realization; independently reviewed.**
The physical fusion and scalar response formula are inherited from
[Session286, reset theorem](../../session286_precise_bridges_20261112/reset/PROOF.md).
The new result is its exact two-ended inverse-limit geometry at two primes,
including the finite physical precision bill. No completion point is asserted
to be a tournament count, and Chan–Pak remains OPEN.

## 1. Inherited physical source and the two observables

The deletion core has outgoing masks `(18,8,3,5,14)`. Attach a root with cuts
`A=9, B=3, E=18, C=24`. For n>=2 let s=(s_0,...,s_(n-1)) be binary. Prepend
and append zeros to make the edge word of the state chain
`0,s_0,...,s_(n-1),0,0`, interpreting `00=A,01=B,11=E,10=C`.
Thus k=n+2 blocks fuse at the common root, their deletion blocks ordered
left-to-right. Let H_n(s) be this tournament's path count. Reroot at the
specified first vertex of the final A block and fuse a rooted triangle on
the left. Let F_n(s) be the new count, on 5n+13 vertices. The inherited
physical identity is F_n(s)=H_n(s)+2V_n(s), where V_n is the corresponding
nonnegative two-path ending-root observable. Both H_n and F_n are actual
Hamiltonian-path counts, of two explicitly related tournaments.

Set q=104/471 and u=13/471, and the rational constants

    A=-34104960/367, B=3945600/367, C=-140160/367,
    D=-27648000/471^2, E=-829440/471.

The inherited response formula, merely rewritten with n=k-2, gives

    g_n(s) = [V_n(s)-V_n(0)]/471^(n-1)
           = sum_i [A q^i+B q^(n-1-i)+C q^(n-1)]s_i
             + sum_(i<j) d_(j-i) s_i s_j,
    d_1=E q^(n-2),  d_d=D q^(n-1-d) for d>=2,
    f_n(s) = [F_n(s)-F_n(0)]/[2*471^(n-1)]
           =g_n(s)+13u^(n-1)j(s),
    j(s)=sum_i s_i-9sum_i s_i s_(i+1).                       (1)

The exact central difference is H_n(s)-H_n(0)=2*13^n*j(s).
Nothing in this note rebrands these finite response identities as new.

## 2. Boundary surface and the finite precision theorem

For infinite binary sequences x,y define

    X=sum_(i>=0) q^i x_i,   Y=sum_(i>=0) q^i y_i,
    Phi(x,y)=A X+B Y+D X Y.                                 (2)

These series converge in both Q_13 and Q_2, because v_13(q)=1 and v_2(q)=3.
The two evaluations of the same formal address will be denoted Phi_13 and
Phi_2. Set Psi_2=Phi_2/128. All coefficients of Psi_2 are2-adic integral.

**Finite physical precision theorem.** Let r=floor(n/2). Split s into its
left ceil(n/2) bits and its right floor(n/2) bits read backwards, and extend
both sequences by zeros. Then

    f_n(s) = Phi_13(x,y) mod 13^r,
    g_n(s)/128 = Psi_2(x,y) mod 2^(3r).                     (3)

In particular a depth-r address pair is realized at precision13^r and8^r
by an explicit tournament on10r+13 vertices (choose n=2r). The second
observable uses the two physical counts H_n,F_n and their zero-spin
baselines; it is not F_n alone and is not itself claimed to be a path count.

Proof. Cross-pairs i in the left half and j in the right half contribute
exactly Dq^(i+n-1-j) to (2), matching (1), except an adjacent cross-pair.
Same-half nonadjacent interactions have exponent n-1-(j-i)>=r; hence their omission
has valuation at least r at13, and at least3r after division by128 at2
(the constant D/128 has2-valuation6). Same-half adjacent interactions have coefficient E q^(n-2), also within the
claimed bound whenever such a pair exists (n>=3). Linear terms directed toward the
opposite boundary and the C term have exponents at least r. For the
adjacent discrepancy,

    E-D=-40335360/24649,
    v_13(E-D)=1,  v_2((E-D)/128)=4.

Its exponent is n-2, so even n=2 obeys the required precision. Finally
v_13(13u^(n-1)j)>=n. This last term must be subtracted at2, where u is a
unit: it does not vanish uniformly. These estimates prove (3). The
interior leftover bit for odd n has depth r and is invisible at that
precision. ∎

The same proof allows arbitrary infinite completions of x,y after their
known prefixes: their tails have valuations at least r and3r respectively.
Thus (3) is a uniform approximation theorem, not a subsequence-only limit.

## 3. Exact13-adic tree: the first different boundary layer determines distance

**Theorem.** If t is the first index at which (x_t,y_t)!=(x'_t,y'_t), then

    v_13(Phi_13(x,y)-Phi_13(x',y'))=t.                       (4)

Consequently the image is a compact four-branch13-adic Cantor set. Modulo
13^r it occupies exactly4^r residue classes, and each occupied class has
exactly four occupied children at the next level. Its Haar measure is0,
it is nowhere dense in Z_13, and its Hausdorff dimension in the usual
13-adic metric is log(4)/log(13).

Proof. The coefficients (A,B,D) reduce to(11,3,4) modulo13. At t=0 the four
outputs are0,11,3,5, distinct. For t>=1, put DeltaX=X-X', DeltaY=Y-Y'.
The bilinear difference is

    (A+D Y')DeltaX+(B+D X')DeltaY+D DeltaX DeltaY.

After dividing by q^t, the last term vanishes modulo13. The leading pair
is (A+D y_0,B+D x_0). For outermost states00,01,10,11 these are respectively
(11,3),(2,3),(11,7),(2,7), whose four binary subset sums modulo13 are all
distinct. This proves (4). The argument also gives an effective inverse:
try four next bit-pairs and retain the unique one matching the next digit.

Compactness follows from continuity of the product of finite binary spaces.
The level-r cylinder cover has4^r balls of radius13^(-r); the uniform address
measure gives each such cylinder mass4^(-r). The standard upper-cover and
lower mass estimates give the stated dimension, with no probabilistic model
for tournaments assumed. Haar mass is the limit(4/13)^r=0; compact Haar-zero
sets contain no open ball. ∎

## 4. Base8 geometry at the prime2: bounded distortion, not a false isometry

**Theorem.** With t as above and

    delta=5(x_t-x'_t)+7(y_t-y'_t),

one has

    v_2(Psi_2(x,y)-Psi_2(x',y'))=3t+v_2(delta),             (5)

where v_2(delta) is0,1 or2. The image occupies exactly4^r residue classes
modulo8^r, has Haar measure0, is nowhere dense in Z_2, and has Hausdorff
dimension2/3.

Proof. The coefficients (A/128,B/128,D/128) are(5,7,0) modulo8; the last
coefficient actually has valuation6. At the first differing pair, all later
digits contribute multiples of8 after division by q^t. The leading
expression is5Delta x+7Delta y modulo8. The four single-pair values are
0,5,7,4, distinct. Since their nonzero differences have valuation at most2,
higher terms cannot cancel their leading power of2. This proves (5),
including t=0. The exact factor(q/8)^t is a2-adic unit and changes no
valuation. The metric is thus within constants1 and1/4 of8^(-t), not
identically8^(-t). Covers by4^r depth-r balls give upper dimension log4/log8;
uniform cylinder mass and the bounded distortion give the matching lower
bound. Haar and interior statements follow as above. ∎

The common address induces a bijection between the13-adic and2-adic Cantor
sets. It is a homeomorphism satisfying the precise two-sided Hölder estimate

    (1/4)|z-z'|_13^alpha <= |Theta(z)-Theta(z')|_2
                               <= |z-z'|_13^alpha,
    alpha=log(8)/log(13).                                  (6)

This is a concrete relation between prime geometries on one physical
construction, not an identification of Q_2 with Q_13.

## 5. Finite realization and the crucial correlation sidecar

For n=2r, all4^r physical words produce distinct normalized residues both
mod13^r and mod8^r by(3)–(5). The paired residue set also has exactly4^r
elements, not16^r: both coordinates encode the **same** bit address. The
ambient residue product has104^r elements. CRT cannot select an independent
13-adic address and2-adic address in this family. Equality of the recovered
addresses is a complete compatibility test at depth r.

The inverse is a genuine finite algorithm. For either prime, start with empty
prefixes, try the four next pairs, evaluate the rational polynomial(2) on
the finite prefix modulo the next prime power, and retain the unique matching
pair. Every denominator is a unit at the relevant prime. A target that fails
at a digit is outside the specified cylinder set; a target that passes r
digits determines a physical length2r word. It only guarantees the stated
normalized congruences, not equality of its full integer count with a chosen
target. One must finally evaluate the actual count for that claim.

## 6. Exact additive saturation in the completions

Let K_13 be the image of Phi_13 and K_2 the image of Psi_2.

**Theorem.** K_2 is an additive basis of Z_2 of minimum order2, and K_13
is an additive basis of Z_13 of minimum order3:

    K_2+K_2=Z_2,     K_13+K_13+K_13=Z_13.                  (7)

This is a proved local additive-surjectivity result, despite both summand
sets having Haar measure zero. It is not a physical closure-under-addition
statement for tournament path counts.

Proof at2. After any selected prefixes, extension by one address pair
changes the next base8 digit by a common unit times one of{0,5,7,4}.
The twofold sum of this set is all of Z/8Z. Given a target2-adic integer,
choose two root pairs to match its first base8 digit. Whatever choices have
already been made, the same twofold sum property chooses the next pairs.
Induction and continuity give an exact sum representation. One summand is
not enough because it occupies only four of eight root residues.

At13, the first derivative alphabet depends on the two outermost bits,
but never on later bits. Enumerate outermost states as

    0=(0,0), 1=(0,1), 2=(1,0), 3=(1,1).

Their root values are0,3,11,5. Their later derivative alphabets are

    D_0={0,11,3,1}, D_1={0,2,3,5},
    D_2={0,11,7,5}, D_3={0,2,7,9}.

The following table selects three outermost states for every root target;
each row has D_i+D_j+D_k=F_13 (directly check the at most64 sums).

    target:  0       1       2       3       4       5       6
    states: (1,3,3) (0,1,2) (3,3,3) (0,2,3) (1,1,2) (0,0,3) (1,2,3)
    target:  7       8       9       10      11      12
    states: (2,2,2) (0,1,3) (0,2,2) (0,3,3) (0,0,2) (1,2,2)

Once this root choice is fixed, arbitrary existing prefixes contribute a
known carry; varying their next pairs contributes a common nonzero factor
(q/13)^t times that same full derivative sumset. Every next target digit can
therefore be chosen, forever. Continuity proves the claimed exact sum in
Z_13. Two summands cannot suffice: their root sumset is
{0,1,3,5,6,8,9,10,11}, missing2,4,7,12. This proves minimality. ∎

**Finite physical corollary.** For every r>=1, sums of two length2r
normalized physical V-responses cover all residues modulo8^r, and sums of
three length2r normalized physical F-responses cover all residues modulo
13^r. The digit-lift algorithm constructs those words, and(3) transfers their
residues to actual tournament observables. The two prime conclusions are
separate; they do not assert that the same two or three words solve
independently prescribed residues at both primes.

This pinpoints a potentially useful missing operation: a tournament
construction implementing these correlated sums with a controlled absolute
baseline and reusable receiver. The usual ordinal-sum operation multiplies
path counts, not adds them. Formula(7) therefore does not fill the integer
spectrum. It removes the local additive obstruction **in the stated
completion model** and leaves a concrete physical realization/height bill.

## 7. What this contributes and what it does not

Session286 supplied an injective13-adic inverse at each fixed k. This note
adds a common inverse-limit surface, exact finite precision uniform in all
words, the prime2 code in a separately observable response, two dimensions,
and an explicit cross-prime compatibility law. It isolates a reusable
principle: **retain the address when changing arithmetic completion**.

The compactification shows a real regularity emerging from long tournament
chains, but its four children per level are an exact local digit shortage,
not a route around the shortage. Neither single code is locally surjective, although the additive convolutions
in(7) are. The
rational normalizations include a length-dependent baseline and unit factor;
forgetting them would manufacture a nonexistent integer selector. The2-adic
observable additionally subtracts a physical central count. There is no
claim that arbitrary completion points, arbitrary p-adic integers, or
arbitrary positive odd integers are actual path counts. Chan–Pak remains
open. A physical construction enlarging these compatible digit fibers while
retaining absolute count and a renewable receiver is still needed.

## 8. Reproduction and provenance

`python3 verify.py` and `python3 -O verify.py` use only the standard library.
`core.py` and `transfer.py` are verbatim copies of the reviewed Session286
literal definitions, retained here for self-contained reproduction, not new
algorithms. The new verifier checks their physical responses against(1), all
finite words through n=8, exact cylinder cardinalities through depth6,
leading digit inverses, exceptional adjacent coefficients, seeded metric
controls, the finite root-sum certificate, additive digit inverses, and literal
normalized physical sum controls. These finite checks support the all-order proofs; they do not
replace them. The formulas, dimensions and address-correlation conclusions
have the explicit infinite scope proved above.
