# Independent audit of the orbit crossing phase

**Verdict: ACCEPT for the scoped crossing identities, finite-height
Sturmian language, all-odd finite saturation, run-resource obstruction,
and conditional fixed-orbit Fourier/Benford statements in both clocks.**
No Collatz conclusion is obtained. The
reciprocal-summability input is inherited proved canon, not independently
reproved by this audit.

Auditor: `fibonacci_colours`, separately from the author of
[the crossing potential note](crossroads_crossing_20260926_potential.md).
The main dependency is
[THM-4476, thin divergent orbits and finite reciprocal sums](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md).
The lower-discrepancy irrationality corollaries of that theorem already
cover bounded and some logarithmic strips; this audit does not repackage
those as new results.

## 1. Identities and types

For successive positive odd orbit values x_j, put
a_j=v2(3x_j+1), alpha=log2(3/2), c_j=log2(1+1/(3x_j)),
phi_j={log2 x_j}, and let u_j be the height change at the first shortcut
odd step. The exact logarithmic equation is

    log2 x_(j+1)-log2 x_j = alpha+c_j-(a_j-1).

Even halvings leave the mantissa unchanged. Hence

    u_j=phi_j+alpha+c_j-phi_(j+1),
    h_(j+1)-h_j=u_j-(a_j-1).

The test `3x_j+1 >= 2^(h_j+2)` is the correct dyadic up-crossing test.
The alphabet is {0,1}, including x=1 where the first shortcut step is
1->2 and u=1. The terminal odd fixed point is not an exception erased
from the formulas: its carry is positive at every iteration and its
carry sum diverges.

On a non-eventually-periodic positive integer orbit, all states are
distinct, so they tend to infinity. THM-4476 also gives
sum 1/x_j<infinity. The elementary inequality log(1+t)<=t proves
sum c_j<infinity without an independence or average-residue assumption.
All conclusions in the following sections concern this individual orbit
when an infinite-orbit hypothesis is stated.

## 2. Finite-height cut stability: proof accepted

The cut for prefix length s in the multiplicative phase coordinate
t=2^theta is

    b_s=2^ceil(log2(3^s))/3^s.

This agrees with 2^ceil(s alpha)/(3/2)^s; subtracting the integer s does
not change the fractional logarithmic phase. These cuts are distinct,
since any equality would give an integer relation between powers of 2
and 3. Together with 1 and 2 they define L+1 intervals and therefore
L+1 binary words of length L.

For an actual block, freeze its actual carry prefix P_s=sum_(j<s)c_j.
Its prefix crossing counts are exactly
floor(phi_0+s alpha+P_s). The perturbed cut is b_s/2^P_s. If all sources
are at least m and (1+1/(3m))^L<rho_L, where rho_L is the smallest ratio
of consecutive sorted cuts including both endpoints, then each cut
remains strictly between its former predecessor and its former position.

Four details are necessary and present in the argument:

1. The carry prefixes are frozen at their actual values. The proof does
   not assert that the carries remain constant when one varies the
   starting integer or phase.
2. Including endpoint 1 ensures that no shifted cut passes the initial
   phase; thus every prefix floor at phase zero is unchanged.
3. Including endpoint 2 closes the circular seam. No separate untracked
   wrap occurs.
4. At a cut, floor is right-continuous. Endpoint values add no extra word.

The ordered changes in the prefix-floor vector are therefore identical
to those in the unperturbed rotation. Taking successive differences
preserves the complete word set. The finite-height theorem is valid.

The added equivalence log2(rho_L)=min_(1<=d<=L)||d alpha|| is also valid:
differences between the phase cuts are precisely multiples d alpha
modulo one, and the least circular pair distance is a consecutive-cut
gap. Consequently a word outside the rotation language certifies that
at least one of its sources lies below the certified minimum height.

For an escaping orbit, every fixed L eventually meets this height bound.
The valid quantifiers are `for every L, there is an index J(L)`.
They do not imply a single index after which the infinite word is exactly
one unperturbed Sturmian word. The note retains this distinction.

The local consequences were separately checked: `00` is impossible
because 2 alpha>1; `111` is impossible when its three sources are at
least 7 because (11/7)^3<4. The odd source block 3,5,1 shows why the
small core must remain in the statement.

## 3. Fixed-orbit Benford and bounded Fourier sums: proof accepted

Writing C_J=sum_(j<J)c_j gives

    phi_J={phi_0+J alpha+C_J},   C_J -> C_infinity.

For a fixed nonzero integer m, set z_j=exp(2 pi i m phi_j) and
q=exp(2 pi i m alpha), with q!=1. Direct subtraction gives

    (q-1) sum_(j<J) z_j
      = z_J-z_0-q sum_(j<J)z_j(exp(2 pi i m c_j)-1).

The bound |exp(it)-1|<=|t| yields the claimed bound

    |sum_(j<J)z_j| <= (2+2 pi |m| C_infinity)/|q-1|.

Every fixed Fourier average thus tends to zero. Continuous
trigonometric approximation, followed by approximation of interval
indicators whose endpoints have zero limiting measure, proves the
base-2 mantissa law for the odd orbit values. The use of one fixed orbit
is legitimate: there is no random starting integer in this argument.

The denominator depends on the mode. An unrestricted infinite Fourier
series requires a summability condition against those denominators.
A finite Fourier correction is bounded and cannot independently supply
an unbounded cumulative negative drift. These limitations are correctly
stated in the note.

The added base-3 clock on every shortcut iterate also checks exactly:
log3(n_t)=log3(n_0)-t log3(2)+O_t+D_t, where O_t is the integer odd-step
count and D_t is the actual carry. Reducing modulo one cancels O_t.
Reciprocal summability makes D_t converge, and the same telescoping
Fourier argument proves the conditional base-3 law for the full orbit.
The odd-event base-2 clock and full-time base-3 clock use different
sampling times. Their separate laws do not imply independence of phase
and parity, which the producer explicitly preserves as an open coupling.

## 4. Saturation and the hidden run resource

Every length-L unperturbed rotation cylinder contains an open phase
interval. Choose arbitrarily large

    n=-1 mod 2^(L+1)

with mantissa in a compact subinterval of that cylinder. Such choices
exist because the normalized arithmetic-progression mesh in dyadic
intervals tends to zero. The next L odd-to-odd steps then all have
a_j=1, and their cumulative additive carry tends to zero as n grows.
For sufficiently large n their crossing word is the desired cylinder
word. Thus **every allowed finite crossing word has arbitrarily large
actual integer realizations with no extra even-halving depth**.

The exponent L+1 in the modulus is essential for guaranteeing all L
valuations equal one. Merely imposing n=-1 modulo 2^L does not give the
last valuation. The producer used the correct exponent.

The additional resource identities in the experiment were independently
checked algebraically. For r(n)=v2(n+1), an a=1 transition satisfies
r(Tn)=r(n)-1, and

    J(n)=3^r(n) (n+1)/2^r(n)

is unchanged. But for any odd k>=3 the exact transition

    n=(2^(k+2)-5)/3 -> 2^k-1

has a=2 and changes r from 1 to k. Its J value changes from
2^(k+1)-1 to 3^k. Thus a single drop whose value ratio tends to 3/4
can recharge an arbitrarily large amount of this resource. Restricting
the source to the image of the odd map does not eliminate the family:
n is not divisible by 3 for k=1 or 5 modulo 6. This is a genuine hostile
to charging resource resets by only a_j or the one-step logarithmic loss.

Finally, the proposed potential V_c(n)=log2(n)+c v2(n+1) fails for every
constant c. Arbitrarily large a=1 transitions have increments tending to
alpha-c, so eventual nonincrease would require c>=alpha. The reset
family then has increment log2((2^k-1)/s_k)+c(k-1), which tends to positive
infinity because its first term tends to log2(3/4). This is a correct
two-family obstruction. The equality c=alpha already has a small positive
increment on each a=1 transition, but the weaker necessary inequality
used by the producer is sufficient for the contradiction.

The producer's final strengthening is also accepted. Continue the reset
through its complete a=1 run:

    s_k -> 2^k-1 -> ... -> z_k=2*3^(k-1)-1.

Both endpoints have R=1, whereas z_k/s_k is asymptotic to
(1/2)(3/2)^k. Every intermediate odd value is at least 2^k-1. Hence no
potential log2(n)+f(R(n))+g(n), with arbitrary finite-valued f and bounded
g, can be nonincreasing on every sufficiently large odd transition:
the f terms cancel between endpoints and g cannot absorb the unbounded
logarithmic increase. A positive constant times log2(n) is excluded by
the same argument.

For k=1 or 5 modulo 6, the entire odd macroblock is coprime to 30. To see
the mod-5 part, an intermediate value has the form
3^j 2^(k-j)-1, and its product term is (-1)^j 2^k modulo 5. Since k is
odd this is 2 or 3, so the value is 1 or 2 modulo 5. Modulo 3 all j>=1
values are -1, while j=0 has residue 1; the source s_k is a unit for the
specified k classes. All values are odd. Thus a restriction to unit
residues modulo 30 leaves these obstructions intact.

## 5. Independent finite checks

The independent implementation imports none of the code under review:

    python -B 04-computation/experiments/crossroads_crossing_20260926_phase_audit.py

Its [retained exact stdout](crossroads_crossing_20260926_phase_audit.out)
records the 72 language checks, 324 actual saturation blocks, and sixteen
polynomial-shadow blocks described below and in Section 7. Normal and
Python -O executions agree byte-for-byte; all checks use explicit
exceptions and integer/Fraction arithmetic.

I reconstructed the rational cut system and language separately from
the producer's script. For each L=1,...,24, three independent positive
carry profiles were used: source heights m, m+j, and 2m+3j, with m the
smallest integer satisfying the sufficient-height inequality. Enumerating
all perturbed phase intervals gave exactly the unperturbed language in
all 72 cases. The first 16 sufficient heights reproduced the producer's
list. The complete list through 24 was

    2,6,9,12,32,39,45,52,58,64,71,296,
    320,345,369,394,419,443,468,492,517,541,566,591.

I separately chose one 101-bit source in n=-1 modulo 2^(L+1) for every
one of the L+1 rotation cylinders for each L through 24. All 324 actual
source blocks had the prescribed crossing word and every a_j=1.
This finite census corroborates the general saturation proof; it does
not replace its arbitrary-height argument.

The producer's retained reproduction command is

    python -B 04-computation/experiments/crossroads_crossing_20260926_potential.py

Its controls use integer/Fraction arithmetic for all decisions. Floating
logarithms are printed diagnostics only. The experiment also checks
actual completed excursions of 27, including growing returns, separately
from the phase construction. The revised producer experiment retains
its own complete finite-saturation and run-resource controls.

## 6. Consequence and stopping boundary

The crossing phase has an exact small language and a conditional
fixed-orbit equidistribution law. Neither controls the total valuations.
Finite saturation shows that no length-L crossing-word exclusion can
alone force even depth: every allowed word survives the maximally growing
a_j=1 choice for that entire finite block. The next theorem must use
compatible continuation, a return-time condition, the exact integer
carry, or another unbounded arithmetic sidecar.

Nothing in this audit promotes the attached poset preprints to external
canon, supplies a uniform linear-extension law on actual excursions,
or proves Collatz. The existing irrationality exclusions from THM-4476
remain inherited results rather than newly claimed progress.

## 7. Root's finite polynomial-valuation obstruction: independently accepted

The root subsequently proposed the following stronger theorem. Fix a finite
prime set S, finitely many nonzero integer polynomials P_i, a constant c>0,
an arbitrary real-valued function f of the finite valuation vector

    (v_p(P_i(n)))_(p in S, i),

and a bounded real function g(n). Then

    V(n)=c log(n)+f((v_p(P_i(n)))_(p,i))+g(n)

cannot be nonincreasing at every sufficiently large shortcut Collatz step
where its polynomial features are nonzero. It also cannot be nonincreasing
at every sufficiently large accelerated odd-to-odd step. The vector has
finitely many coordinates; the coordinates and f need not be bounded.

**Verdict: ACCEPT.** The following is an independent derivation, with the
points at which hidden assumptions could otherwise enter made explicit.

Choose a positive integer B divisible by ord_p(2) for every p>=5 in S.
Choose L divisible by both ord_p(2) and ord_p(3) there. Empty sets of
orders cause no problem: B=L=1 is allowed. For sufficiently large
a=1 modulo L, 3^a>2^(a+B). Set

    D_a=3^a-2^(a+B),
    r_a=-(3^a-2^a)/D_a,
    M_a=3^a/2^(a+B)>1,       ell=a+B.

The shortcut branch word w=1^a 0^B has affine map

    F_w(x)=(3^a x+3^a-2^a)/2^(a+B),

so r_a is its fixed point. It is a **valid periodic parity word**, not
merely a formal affine fixed point: D_a is odd and

    r_a+1=-2^a(2^B-1)/D_a,

whose 2-adic valuation is exactly a. The first a sources are odd;
F_(1^a)(r_a)=2^B r_a has valuation exactly B; the next B steps are even
and return to the odd rational r_a. Equivalently, each block is a valid
accelerated odd-step word with valuations 1 repeated a-1 times, then B+1.

All denominators are units at the selected primes. At p=2 the denominator
is odd. At p=3 it equals -2^(a+B) modulo 3. At p>=5 the order choices
give D_a=3-2=1 modulo p. Consequently r_a is p-adically integral for
every selected p.

Also

    r_a=-1-(2^B-1)/((3/2)^a-2^B)

strictly increases to -1 through distinct values. A nonzero polynomial
has only finitely many roots, so one may choose a in the prescribed
progression for which every P_i(r_a) is nonzero. This step prevents a
valuation pole from being hidden in f. Fix this a from now on; it does
not change with the repeated-block length k.

Choose e_p strictly larger than every v_p(P_i(r_a)). Include p=2 as an
auxiliary congruence prime if it is absent from S. For each k, the Chinese
remainder theorem gives arbitrarily large positive integers n satisfying

    n=r_a modulo 2^(k ell+e_2),
    n=r_a modulo p^e_p for each odd p in S.

Rational congruences are legitimate because the denominators are units.
The 2-adic congruence forces the actual parity prefix w^k. Its endpoint y
therefore satisfies the exact real and rational identity

    y=r_a+M_a^k(n-r_a).

At p=2, multiplying n-r_a by M_a^k loses exactly k ell digits, leaving
at least e_2. At p=3 it gains ak digits. At every other selected odd prime
it preserves the precision. Both n and y thus equal r_a modulo p^e_p.
Polynomial evaluation with integral coefficients preserves these
congruences, and the strict inequalities defining e_p force

    v_p(P_i(n))=v_p(P_i(r_a))=v_p(P_i(y))

for every feature. The complete f term cancels between endpoints.

The ordinary-positive-integer issue is also explicit: choose n greater
than 2^(k ell) times any desired finite-core bound within its CRT class.
Since every shortcut step is at least half its positive source, all
intermediate values exceed that bound. Choose the bound above every
positive integer root of the finitely many polynomials. All intermediate
features are then defined, so every claimed local monotonicity inequality
can legitimately be summed along the block.

Finally r_a<0 gives

    y/n=M_a^k+(-r_a)(M_a^k-1)/n > M_a^k.

Thus the endpoint change in V exceeds c k log(M_a)-2||g||_infinity,
which is positive for large k. This contradicts either shortcut-step
monotonicity or odd-to-odd monotonicity. No conjecture about eventual
termination is used; only finite positive orbit prefixes are needed.
The positive source may change with k. The unique infinite periodic
shadow here is the negative rational r_a, so this construction does not
create a divergent positive integer orbit or bypass the fixed-integer
compatibility obstacle by a quantifier exchange.

Independent exact controls used nine polynomials

    n, n+1, n-1, n+5, n^2+1, n^2-1, n^3+2, 3n+1, 9

and k=1,2,4,8 for each of four prime sets. The selected parameters were

| S | B | L | a | largest endpoint bit length |
|---|---:|---:|---:|---:|
| {2} | 1 | 1 | 3 | 58 |
| {2,3} | 1 | 1 | 3 | 58 |
| {2,3,5,11} | 20 | 20 | 41 | 582 |
| {2,3,7,13} | 12 | 12 | 25 | 362 |

All sixteen actual repeated-word blocks had exactly the claimed parity
prefix, affine endpoint, identical valuation vectors at their endpoints,
ratio greater than M_a^k, and every intermediate value above one million.
The polynomial n+5 deliberately excludes the first candidate fixed point
r=-5 when B=1; selecting a=3 exercises the root-avoidance step.

The theorem excludes this precisely specified class of globally monotone
potentials. It does not exclude functions using infinitely many primes,
unbounded real phase or mantissa information, infinitely many forms,
nonlocal excursion memory, or selected-return-time monotonicity. The
independent script and retained stdout linked in Section 5 reproduce the
separate controls; this section records the hostile proof audit.
