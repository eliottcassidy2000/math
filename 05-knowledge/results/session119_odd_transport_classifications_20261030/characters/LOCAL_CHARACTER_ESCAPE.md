# Escape from the single-character obstruction: independent characters by stratum

**PROVED constructive family and exact counts; FINITE-EXACT spectral/cycle controls.**
This extends the stratum/stabilizer mechanism of [CHARACTER_STRATA](CHARACTER_STRATA.md),
which in turn credits the inherited prime-power valuation-frame theorem. It does
not claim a new Ramsey number or classify all circulant tournaments. A single
unit character on all strata is an extra constraint, not a consequence of oddness.
Dropping precisely that constraint gives the following constructive extension.

## 1. The named family and its complete all-modulus count

Let N>1 be odd. For each additive-order stratum

    S_m={ (N/m)t mod N : t in U(m) },  m|N, m>1,

choose its own odd sign character χ_m:U(m)→{±1} and its own sign ε_m, and put

    σ((N/m)t)=ε_m χ_m(t).

Different strata may use different characters. This means that multiplication
by a unit acts on each stratum with its prescribed sign, not necessarily with
one common global sign. Oddness holds because every χ_m(-1)=-1. Consequently
x→y iff σ(y-x)=+1 defines a full circulant tournament.

**Theorem 1.** This full-support stratumwise-character family is nonempty if and
only if every prime divisor of N is congruent to three modulo four. When
N=∏p^{a_p} satisfies that condition, its exact number of labelled symbols is

    2^{ Σ_{m|N,m>1} ω(m) }
      = 2^{ Σ_{p|N} a_p ∏_{q|N,q≠p}(a_q+1) }.

Here ω(m) is the number of distinct prime divisors of m. No quotient by graph
isomorphism or multiplier relabelling has been taken.

**Proof.** A sign character on U(m) is a product of the Legendre characters at
some subset of the prime divisors of m, as proved in the companion classification.
It is odd precisely when that subset contains an odd number of primes congruent
to three modulo four. If m has r distinct prime factors and at least one such
prime, exactly half of the 2^r sign characters are odd: pair subsets by toggling
one fixed three-mod-four prime. Thus there are 2^{r-1} odd characters, and
multiplying by the independent sign gives 2^r distinct stratum symbols. The
sign is recovered as σ(N/m); normalizing by it recovers χ_m, so this count has
no duplicate parametrizations.

If N has a prime divisor p≡1 mod4, the required prime stratum S_p has no odd
sign character, so full support is impossible. Conversely, if every p|N is
three modulo four, every nontrivial divisor m has an odd character, so all the
independent choices exist. Multiplication of their counts gives the first
formula. Finally

    Σ_{m|N,m>1}ω(m) = Σ_{p|N} #{m|N:p|m}
                    = Σ_{p|N} a_p ∏_{q≠p}(a_q+1),

which proves the second formula. This is an exact classification of the stated
family, not merely a lower bound supplied by one construction.

At a prime power p^a the character on every stratum is forced to be the same
Legendre character, so the formula reduces to the inherited 2^a valuation
frames. At a modulus with several prime factors it provides genuinely new
cross-stratum choices that a single-character rule forbids.

## 2. Partial support for arbitrary odd N

One can allow an entire stratum to be zero, or use a nonzero odd-character
symbol on that stratum. A character chosen on a zero stratum is not counted as
additional data, since it does not change the symbol.

Call m eligible if at least one of its prime divisors is three modulo four.
The complete support polynomial for this partial stratumwise family is

    Σ_σ z^{|supp σ|}
        = ∏_{m|N, m>1, m eligible} (1+2^{ω(m)} z^{φ(m)}).

All ineligible strata are forced to zero; all eligible ones have either one
zero option or 2^{ω(m)} different nonzero symbols. Their independence proves
the formula, including its exact coefficients. The full-support theorem is the
special case where no nontrivial stratum may be zero.

For example N=15 has eligible strata m=3 and15, but m=5 is forced to zero.
The polynomial is

    (1+2z^2)(1+4z^8)=1+2z^2+4z^8+8z^10,

so there are fifteen partial symbols. By comparison the fixed global character
χ_3 admits only nine partial symbols. The extra choices come from allowing the
unit stratum to choose χ_3χ_5 independently of the χ_3 character required on
the order-three stratum. This does not manufacture a full tournament at N=15,
which still fails at the order-five stratum.

## 3. Modulus21: sixteen full symbols and exactly four isomorphism classes

At N=21, strata have orders 3,7,21. Their numbers of nonzero choices are 2,2,4,
respectively, giving exactly sixteen symbols. Both prime divisors are three
modulo four, so these are full tournaments. In contrast, the single-global-
character theorem admits no full symbol at this modulus.

The checker enumerates **all** 2^10=1,024 odd sign functions on Z/21Z and tests
local multiplicativity directly, rather than assuming the construction is
complete. Exactly sixteen pass, and they agree with the parametrized family.

Let S[x,y]=σ(y-x) be the skew adjacency matrix. For each of the sixteen
symbols, the checker computes tr(S^4) by literal integer convolution and counts
simple directed five-cycles by explicitly enumerating cycles through zero and
using translation symmetry. Cycles are counted up to cyclic rotation, with
direction respected. The observed profiles are:

| tr(S^4) | Directed five-cycles | Number of symbols |
|---:|---:|---:|
| 15,876 | 19,530 | 4 |
| 19,908 | 19,404 | 4 |
| 23,940 | 19,278 | 4 |
| 44,100 | 18,648 | 4 |

The receipt contains a positive-difference connection set witnessing each row.
Both statistics are tournament-isomorphism invariants, so different rows cannot
be isomorphic. We can also prove the matching upper bound without a general
isomorphism solver. Normalize each stratum by its representative N/m, with signs
(ε_3,ε_7,ε_21). The unit-stratum character χ_21 is either χ_3 or χ_7. A unit
multiplier u changes the first two signs independently by χ_3(u),χ_7(u), since
CRT realizes all four sign pairs. It changes ε_21 by the selected χ_j(u).
Thus for fixed χ_21=χ_j the product ε_21 ε_j is invariant, and the eight sign
choices split into exactly two multiplier orbits of size four. Across both
choices of χ_21 there are four multiplier orbits. The literal checker verifies
these orbits and assigns a different displayed trace to each one. Every member
of an orbit is explicitly isomorphic by a unit permutation; distinct orbits
cannot merge under any further isomorphism because their traces differ.

Consequently there are **exactly four tournament isomorphism classes within
this sixteen-symbol family**. This does not classify all circulant or arbitrary
21-vertex tournaments, and equal spectral data alone was not used to prove
isomorphism.

For the spectral calculation, put c(d)=Σ_t σ(t)σ(d-t). Oddness makes c even,
and direct matrix multiplication gives

    tr(S^4)=N Σ_d c(d)^2.

For the cycle count, let R be the number of distinct-vertex sequences
(0,v1,v2,v3,v4) whose consecutive arcs and final arc v4→0 all point forward.
Translation symmetry gives N R rooted sequences in total; each unrooted directed
five-cycle has five roots. Hence its exact count is NR/5. These computations
retain the actual relation, not just an arithmetic count of characters.

## 4. What is unified, and what remains open

The exact source-to-target map is

    prime factorization → local odd characters → difference symbol → tournament.

The global-character obstruction explains why the old prime-power construction
cannot simply be transplanted unchanged. The stratumwise extension restores a
complete constructive family and exposes genuinely different cycle/spectral
structures already at21 vertices. It does not prove that any member is Ramsey
optimal, maximizes Hamiltonian paths, or provides an unrestricted circuit lower
bound. Those are separate extremal or complexity questions requiring additional
structure beyond the character count.

## 5. Reproduction

`check_local_characters.py` directly enumerates every odd sign symbol for all
odd N=3,...,21 and compares the accepted set with the construction and formula.
It also checks every constructed symbol for N=27,49,63,77,99,147 against the
actual local character equations. For odd N<=15 it enumerates all odd
{0,±1}-valued functions, including zero/mixed-support hostiles, and verifies the
partial support polynomial. Its N21 profile calculation uses exact integers.
Normal and optimized runs agree with `local_receipt.json`. Explicit exception
guards survive optimization; finite tests corroborate the all-N proof above.
