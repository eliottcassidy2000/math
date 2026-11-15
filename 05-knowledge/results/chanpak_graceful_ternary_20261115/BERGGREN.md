# A Berggren address is not a tournament norm address

**Proved.** The natural integer chart from the Berggren odd roots \((s,t)\) to
the order-two path-cover profile \((q_2-H,H)\) has no physical edge, and
no bounded-length macro alphabet supports an infinite chain of physical
endpoints. There is, however, an exact dyadic membership decoder on every
fixed-\(H\) parabolic ray. This is about the literal norm carrier, not about
all tournament constructions or the validity of Chan–Pak.

## The right analogy: an invariant and an admissible inverse

The [September Berggren note](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/ternary_berggren_20260925.md)
uses coprime odd \(s>t>0\), with primitive triple
\((st,(s^2-t^2)/2,(s^2+t^2)/2)\). Its three branches are

    B₁(s,t)=(s+2t,t), B₂(s,t)=(2s+t,s), B₃(s,t)=(2s−t,s).

Their positive inverses are forced by the disjoint intervals
\(s/t>3\), \(2<s/t<3\), \(1<s/t<2\). Primitivity and parity handle
boundaries, and the parent decreases height. The earlier
[square-sum norm-one absorber forest](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/square_sum_berggren_absorber_bijection_20261019.md)
shows that this is not a feature of the 3-4-5 root alone: the same
integral \(3\times3\) action preserves the norm-one shell and has
a unique positive lower parent outside its three roots.
Its coefficient forest indexes literal square-sum absorbers because
a separate endpoint construction realizes the coefficient solution.
The operative lesson is to prove preservation of the admissibility
cone under the inverse, not merely preservation of a quadratic equation.

## A literal edge obstruction

For a nonempty tournament of order \(n\), set \(y=H\) and
\(x=q_2-H=P_2+H\). Each Hamiltonian path split at each of its \(n-1\)
internal gaps gives a different ordered two-piece cover, so \(x\ge ny\).
For \(n=2\) the only profile is \((3,1)\); for the singleton it is \((1,1)\).
In the domain \(x>y\) all profiles therefore have \(x/y\ge3\).
Complementary subsets in the two-slot count pair equal odd terms
(Rédei applied to induced subtournaments), so \(q_2\equiv0\bmod4\)
for \(n\ge2\). The natural chart \(s=x,t=y\) thus requires

    s/t≥3,  s+t≡0 (mod 4).

On an odd-root pair \(B_1\) changes \(s+t\) by \(2t\equiv2\bmod4\).
The other branches land at ratios \(2+t/s\in(2,3)\) and
\(2-t/s\in(1,2)\). Hence no single edge has two physical endpoints.

There is also a useful mod-three rejection before any graph search.
The optional-slot polynomial \(q_T(j)=\sum_i (P_i/i!)j^{\underline i}\)
has integer coefficients, and \(q_T(-j)=(-1)^n q_T(j)\).
For the latter, concatenate each ordered slot cover to a permutation
and let \(d\) count backwards arcs along the permutation. The slot
assignments contribute \(\binom{j+n-1-d}{n}\); evaluating at \(-j\)
gives \((-1)^n\binom{j+d}{n}\), the reversed permutation contribution
with sign. Thus \(q_T(2)\equiv q_T(-1)\equiv(-1)^n H\bmod3\).
In the primitive chart \(3\mid t\) is impossible (it would force
\(3\mid s\)). Otherwise \(s\equiv0\bmod3\) forces even order, and
\(s\equiv t\bmod3\) forces odd order. The remaining class
\(s\equiv-t\bmod3\) is impossible.

## Exact decoder along a fixed-target parabolic ray

Fix an odd positive \(h\). All tournaments with \(H=h\) have their
nontrivial strongly connected components in a finite effectively
enumerable catalogue. The component condensation is transitive;
Hamiltonian counts and optional two-slot counts both multiply over
ordinal sums, since each slot passes through components in their
fixed order. Each nontrivial component has at least three Hamiltonian
paths and has order at most \(1+3\log_5 h\) by
[Busch's lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3).
There are at most \(\lfloor\log_3 h\rfloor\) such components.

Enumerate all core lists whose component Hamilton counts multiply
to \(h\), and let \(\mathcal C_h\) be the finite set of their
two-slot products (for \(h=1\) use the empty core). Then

    x+h ∈ ⋃_{c in C_h} {c·2^a : a≥0}

is a necessary and sufficient physical profile test at \(y=h\);
for the empty core use \(a\ge1\) to exclude the empty tournament.
Indeed each singleton component multiplies \(q_2\) by two, and
adding the requested number of singleton components constructs the
witness. The test is exact: compare the odd part of \(x+h\) to that
of a catalogue value \(c\), then compare their dyadic valuations.
It is finite for each prescribed \(h\), without a claim of useful
complexity when \(h\) varies.

Thus in the parabolic ray \(B_1^j(s,h)=(s+2jh,h)\) the admitted
indices are exactly the nonnegative integers of the form
\((c2^a-(s+h))/(2h)\). The integrality test is periodic in \(a\)
once the power-of-two divisibility requirement holds, but the
physical indices are exponential in \(a\). In the first \(J\) ray
steps there are only \(O_{s,h}(1+\log(J+1))\) admissible endpoints
for a fixed starting pair. In particular at \(h=1\), the physical
profiles are exactly \((2^n-1,1)\) with \(n\ge2\) in this chart:
the tournament is transitive. At \(h=3\), the sole nontrivial core
is a directed triangle, and \(q_2=12\cdot2^a\); no point of the
primitive chart is admitted.

## Two stronger no-go statements

1. **Bounded macros.** Fix a finite upper bound \(L\) on the length
   of words used as graph edges. Suppose a chain has literal physical
   profiles at macro endpoints. A macro containing \(B_2\) or \(B_3\)
   ends with ratio \(<3+2L\) (after its last such letter only \(B_1\)
   remains), so its terminal order and \(H\le n!\) are bounded.
   Each occurrence of \(B_2,B_3\) strictly increases the second root
   and \(B_1\) leaves it fixed. Therefore only finitely many macros
   may contain those letters. The remaining infinite tail would be
   a bounded-step path on a fixed parabolic ray, contradicting the
   logarithmic admissible count just proved. Empty words do not
   count as steps.
2. **Rational linear recharting.** No fixed invertible rational
   \(2\times2\) linear change of the root coordinates can send the
   full Berggren tree into literal norm profiles. Iterate \(B_2\)
   from the root: \(s/t\to1+\sqrt2\) and \(t\to\infty\).
   The second row of the rational change of coordinates cannot
   vanish on \((1+\sqrt2,1)\), because that slope is irrational.
   Its output \(y\) either becomes negative, or tends to \(+\infty\)
   with bounded \(x/y\). In the latter case \(x\ge ny\) bounds
   the source order \(n\), contradicting \(y\le n!\).
   This excludes rational linear conjugacies of the two-coordinate
   systems, not nonlinear re-encodings or extra marked states.

There is a separate signature obstruction to the most literal
fixed-target proposal. The Berggren quadratic shell is indefinite and
has infinite integral orbits. A fixed diagonal target \(N=x^2+2y^2\)
is positive definite and has finitely many integral points. The
integer linear isometries of this binary norm are just independent
sign changes: the first column must have norm \(1\) and the second
norm \(2\). Only the identity preserves the positive quadrant.
A same-target ternary mechanism cannot be obtained by transporting
these two coordinates alone with integral norm isometries.

## What this leaves for exact placement

The fixed-\(H\) decoder is useful as a rejection test, but is not a
proof that every odd target is physically accessible. To get a
Berggren-style positive descent, the state needs something outside the
literal unmarked norm profile: an independent reservoir or a marked
continuation coordinate with a return rule depending on the prescribed
target. In such a proposal, check that all inverse chambers are covered,
that each inverse has a physical lift, and that the positive height
drops while the exact equality is preserved. Congruence, positivity
of an abstract norm, or a complete ternary formal language cannot
substitute for the physical lift.

## A positive payoff: the exact joint \((2,3)\)-adic primitive chart

The parity restriction at 3 has a second physical consequence. On an ordinary
diagonal cyclic norm \(D=x^2+2h^2\), odd source order gives \(x\equiv h\pmod3\),
hence \(D\equiv0\pmod3\). Even source order gives \(x\equiv0\pmod3\), hence
\(D\equiv2\pmod3\) if \(3\nmid h\), or 0 otherwise. Thus the class \(1\bmod3\)
is impossible. Since \(D\equiv3\pmod8\), every prime diagonal output \(D>3\)
is **11 mod 24**, a stronger guard than the mod-eight norm alone.

This guard is exact locally: for every \(L\ge3\) and \(a\ge1\), diagonal
profiles with \(\gcd(x,h)=1\) attain every class \(D\equiv3\pmod8\) and
\(D\equiv0\) or \(2\pmod3\) modulo \(2^L3^a\). To see this use \(k\) ordinal
copies of the fixed physical order-five block \((H,q_2)=(11,100)\) and
\(s\ge L\) singleton pads. One witness for the block is lexicographic arc
code 12, with bit 1 directed from the smaller to the larger label.
The outputs have \(h=11^k,\ q_2=100^k2^s\), and their profile is primitive.
Their norm reduces to \(3\cdot121^k\pmod{2^L}\), and 121 generates the subgroup
\(1+8\mathbf Z\bmod2^L\): its difference from 1 has valuation 3, and repeated
squaring raises this valuation by one. Select \(k\) modulo \(2^{L-3}\)
to match the dyadic target.

At the ternary precision, write \(r=q_2/h=(100/11)^k2^s\) and
\(D=h^2[(r-1)^2+2]\). The powers of 2 run through all units modulo \(3^a\)
because \(2\not\equiv1\bmod3\) and \(2^2-1\) has 3-adic valuation 1.
If the target is 0 modulo 3, its normalized value \(D/h^2\) can be solved
by lifting \(r\equiv-1\bmod3\), since the derivative \(2(r-1)\) is a unit.
Choose the corresponding \(s\equiv s_0\bmod 2\cdot3^{a-1}\), increasing
it as needed to ensure \(s\ge L\). If the target is 2 modulo 3, require \(r=1\)
and choose \(k\) to solve \(2\cdot121^k\equiv D\bmod3^a\): 121 generates
\(1+3\mathbf Z\bmod3^a\), as \(121-1=120\) has 3-adic valuation 1 and cubing
raises that valuation by one. This specifies \(k\bmod3^{a-1}\), compatible
by CRT with the dyadic exponent class. Again choose \(s\ge L\) to realize
the prescribed ratio \(r\).

This is a complete two-prime residue theorem without any common divisor
in the physical profile. It does not imply positive density of prime
diagonal outputs or exact attainment of every integer in the classes.
In fact this explicit two-parameter family only has \(O(\log^2 X)\)
outputs up to \(X\): from \(D\le X\) one gets \(h\le\sqrt{X/2}\) and
\(q_2=x+h=O(\sqrt X)\), bounding both \(k\) and \(s\) by \(O(\log X)\).

For **absolute** diagonal membership, enumerate positive odd solutions
\((x,h)\) of \(x^2+2h^2=N\) and apply the fixed-\(h\) catalogue test above.
For a prime target the candidate is unique if it exists: in
\(\mathbf Z[\sqrt{-2}]\), rounding both coefficients in a division
gives norm error at most \(1/4+2/4<1\), so the ring is Euclidean; two
prime norm representations differ only by a sign or conjugation.
This is an exact but not necessarily fast decision procedure for this
carrier, not a test for the unrestricted spectrum.

Replay with `python verify_berggren.py`.
