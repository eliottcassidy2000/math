# The square-subgroup layer: a complete Boolean-quotient classification

**PROVED uniform family classification with bounded exact checks.** This is the
final enlargement of the specifically named families in this packet, not a
classification of all tournaments or a new Ramsey bound. Earlier character and
valuation-frame results remain unchanged. The mechanism is elementary: retain
the quotient by squares, then impose reflection oddness on that entire quotient
instead of assuming every quotient function is a character.

## 1. Exact family and full-support existence

Let N>1 be odd. Write U(N)^2={u^2:u in U(N)} for the square-unit subgroup. Consider
all integer sign functions σ on the nonzero residues satisfying

    σ(-d)=-σ(d),   σ(qd)=σ(d) for every q in U(N)^2.

Thus σ gives a full circulant tournament and is invariant under square
multipliers. We do not require a character describing its response to other
units, either globally or separately on each stratum.

Fix an additive-order stratum S_m={(N/m)t:t in U(m)}. Reduction of U(N) onto
U(m) is surjective, so the image of U(N)^2 is exactly U(m)^2. The square orbits
in S_m are therefore identified with

    U(m)/U(m)^2 ≅ {±1}^{ω(m)},

using the Legendre signs ((t/p)) at the distinct primes p|m. CRT and the
prime-power square classification prove the identification and surjectivity.
Every square-invariant stratum symbol is an arbitrary sign-valued function
F_m of these ω(m) signs. Negation acts on its p-coordinate by multiplication
with (-1/p).

If a prime p≡1 mod4 divides N, consider the stratum m=p. Negation is then
itself a square action, so invariance would require σ(-d)=σ(d), contradicting
oddness for a nonzero sign. Conversely, if every p|N is three modulo four,
negation flips all quotient coordinates. It has no fixed point on each
nontrivial stratum quotient, and an odd sign function can always be chosen.

Hence full support exists **if and only if every prime divisor of N is three
modulo four**, the same existence condition as the local-character family.
The resulting families need not have the same size.

## 2. Complete parametrization and uniform count

Assume N=∏p^{a_p} and every p≡3 mod4. On a stratum with r=ω(m) distinct prime
factors, choose arbitrarily

    F_m:{±1}^r→{±1},    F_m(-s)=-F_m(s).

Then set σ((N/m)t)=F_m(((t/p))_{p|m}). This is a complete, unique parametrization:
all quotient points occur, and every square-invariant odd symbol factors in
exactly this way. The 2^r quotient points form 2^{r-1} negation pairs, each
with one freely chosen sign. Thus this stratum has exactly 2^{2^{r-1}} choices.

The strata are independent, giving the exact total

    # square-invariant odd symbols = 2^C,
    C = Σ_{m|N,m>1}2^{ω(m)-1}
      = ( ∏_{p|N}(1+2a_p) - 1 )/2.

The product identity follows by choosing exponent zero at p with weight one,
or any of its a_p positive exponents with weight two. The excluded divisor one
contributes one before division by two. This is a labelled-symbol count, not
an isomorphism-class count.

## 3. A precise dimension threshold: characters cease to suffice at three primes

The local-character family is a subfamily: on each quotient it uses the functions

    F_m(s)=ε∏_{p in E}s_p,  |E| odd.

These are the signed odd Boolean characters. There are 2^r of them on an
r-coordinate quotient, whereas the full odd-function family has 2^{2^{r-1}}.
The corresponding exponents agree exactly for r=1,2 and differ strictly for
r>=3. Consequently the local-character and square-invariant full families are
equal **if and only if N has at most two distinct prime factors**. This follows
both from the exact counts and from the explicit inclusion. If N has at least
three distinct prime factors, one of its divisors has three, and its quotient
already supplies strictly more choices.

In the notation of the local-character theorem,

    L=Σ_{m|N,m>1}ω(m)=Σ_{p|N} a_p∏_{q≠p}(a_q+1),

so its size is 2^L and the larger family's size is 2^C. This identifies the
exact missing information: a three-prime stratum can carry a nonlinear odd
Boolean response, not merely a sign character.

## 4. Complete three-coordinate quotient classification

For r=3 there are sixteen odd sign functions. Eight are signed odd characters:
±s1, ±s2, ±s3, ±s1s2s3. The other eight are exactly

    ±χ_even(s) M(s),
    M(s)=majority(s1,s2,s3)
        =(s1+s2+s3-s1s2s3)/2,
    χ_even in {1,s1s2,s1s3,s2s3}.

Proof: M is sign-valued and odd. Multiplying by any even character and either
sign preserves oddness. These eight functions are distinct because M never
vanishes, so equality would force equality of the signed even characters.
They cannot be characters: M has four nonzero Fourier coefficients, namely
1/2,1/2,1/2,-1/2 on the four odd characters. Character multiplication merely
permutes that four-element support, and a global sign changes no support.
The sixteen total odd functions therefore consist precisely of these two
eight-member families. This proves completeness, not only existence of a
majority counterexample.

These functions are reflection-odd under simultaneous sign reversal. Their
Boolean truth tables are balanced, not odd-weight. This must not be conflated
with the separate extreme-Walsh/odd-weight Boolean families studied earlier.
No circuit lower bound follows from this quotient classification.

## 5. The first new full-support modulus:231

The three smallest primes congruent to three modulo four are 3, 7, 11. Thus the
least eligible modulus with three distinct prime factors is N=231. Here

    C=(3^3-1)/2=13,    L=3·2^2=12,

so there are 8,192 square-invariant odd symbols and 4,096 local-character symbols.
Only the order231 stratum has three prime factors; the proper strata still
have at most two. This accounts for the exact factor-two enlargement.

One explicit new symbol is obtained as follows:

- On units modulo231, use M((d/3),(d/7),(d/11)).
- On every proper additive-order stratum S_m, write d=(231/m)t and use (t/p),
  where p is the smallest prime divisor of m.

Every value is a sign. Each formula is square-invariant and odd, so the whole
symbol gives a full circulant tournament. It is not in the local-character
family: on units its value at 1 is 1, while its values at 2, 5, 10 are all -1.
Since 2·5=10, multiplicativity fails. A sign twist cannot repair this because
normalization at 1 has already fixed that twist.

The checker supplies all 115 positive differences of this explicit symbol,
checks all 3,465 equations from its fifteen square-unit multipliers, and verifies
the displayed failure of multiplicativity. It performs no 231-vertex isomorphism,
Hamiltonian-path, or extremal Ramsey search.

## 6. Reproduction and boundary

`check_square_symbols.py` enumerates all odd quotient tables for r=1,2,3,4,
checks the signed-character subfamilies and the complete majority classification
at r=3, including their exact Walsh supports. It verifies the arithmetic exponent
identities for every eligible odd N<=500 and the literal231-residue construction.
Normal and optimized executions agree with `square_receipt.json`; all failure
checks are explicit exceptions.

This unifies an arithmetic threshold (three distinct odd primes), a quotient
dimension (three square-class signs), and a new family of odd functions
(majority-type rather than character-type). The classification provides exact
objects and counts for later extremal investigations. It does not prove that
the extra objects improve a Ramsey bound, Hamiltonian maximum, or circuit record.
