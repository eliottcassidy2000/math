# A guarded head/cycle compiler, mixed return family, and ranked extension

2026-10-03 (America/Denver).

**PROVED:** the sufficient compiler theorem and all-height first-descent
families below. **FINITE-EXACT:** their specified replay, polynomial, and bank
comparisons. **OPEN:** universal coverage and Collatz. No novelty claim for
affine composition, negative cycles, CRT, or golden parity coordinates.

## Inheritance and the productive connection

The [recursive entry grammar](entry_20260927_recursive.md) supplies retained
source thresholds, composition of actual words, and the -5 repeated block.
The [new -17 family](collatz_minus17_return_20261003.md) supplies a second
guarded cycle block and an exact all-height return budget. The immediate
next question was whether composing these rules adds coverage while still
paying the original source. The answer below is affirmative.

[THM-4528, golden parity and holonomy](../../01-canon/theorems/THM-4528-collatz-parity-in-base-phi-golden-beta-map-and-holonomy.md),
with its [audited proof](collatz_golden_holonomy_20261001.md), reads the
ordinary-map parity itinerary as a base-phi fraction. Its original prompt
already contains phi^8+1=7phi^4, phi^10=1+11phi^5, and 11_phi=100_phi.
Those identities do not identify the parity itinerary with the literal
digits of its starting integer. Here the integration is a compiler that
retains both a legal integer transfer and its exact golden parity polynomial.
The new return proof comes from the integer guard and threshold.

Closest proved mechanism: source-preserving affine composition. Canonical
hostile: a composed child can descend below its own source while remaining
above the parent's source. Corrected near miss: the missing final oddness
bit in [THM-4512](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md).
Least-used sidecar: the forced ternary factor of a composed inverse guard,
which changes the constant in the uniform payment inequality.

The live board is head word, cycle shadow, CRT guard, retained source,
ordinary-time polynomial, and certified coverage. The anchor is new bank
coverage; the niche is rational preperiodic shadows; the wildcard is an
exact interface to the golden digit reader. No assertion identifies the
three negative cycles with the previously requested three-color grammar.

## 1. General sufficient compiler theorem

Let a=(a_1,...,a_p) be a possibly empty positive valuation word, with

    A=sum a_i, H=2^A, M=3^p,
    Q_a(n)=(M*n+S)/H.

Let w be an actual negative odd cycle based at -c, c positive odd, with
q letters, total division exponent B, P=3^q, Q=2^B, and

    Q_w(n)=(P*n+c*(P-Q))/Q.

Require every nonempty prefix of both words to have strictly expanding
coefficient: 3^i>2^(sum of its first i valuations). In particular P>Q.
The negative cycle must be checked with its actual valuations; its formal
fixed-point equation alone is not accepted by the compiler.

Choose m>=1. Define

    C=H*c+S, g=gcd(C,M), M0=M/g, C0=C/g.             (1)

This compiler also requires

    H*Q^m>C0,   g*Q^m>c.                            (2)

These are sufficient uniform positivity conditions, not necessary
conditions for all possible variants of the construction. Rejection
outside this domain does not refute existence of a route certificate.

Choose the least t>=1 for which

    D=2^t*H*Q^m-M*P^m > max(0, 2^t*C0-M0*c).        (3)

It exists by (2). Solve the coprime congruences

    H*u*Q^m=C0 mod M0,
    g*u*P^m=c mod2^t.                               (4)

All their leading coefficients are units in the indicated rings. The
second congruence forces u odd. CRT gives one positive representative
u0<M0*2^t and every positive solution u=u0+l*M0*2^t, l>=0.

**Compiler theorem.** Every such u gives a positive odd source

    n=(H*u*Q^m-C0)/M0                               (5)

whose exact first descent occurs after p+qm odd steps, with endpoint

    oddpart(z),  z=g*u*P^m-c.                        (6)

The sources form one complete dyadic cylinder

    n=n0 mod2^(A+Bm+t),
    n0=(H*u0*Q^m-C0)/M0.                            (7)

**Proof of head legality.** Its proposed target is y=g*u*Q^m-c, positive
odd by (2). Equations (1), (4), (5) give M*n+S=H*y. An integral inverse
of a positive valuation word with a positive odd target is an actual
positive odd source: reducing modulo 3 forces the last inverse edge to
be integral; divide it out and repeat. Each intermediate is positive
odd, so every indicated valuation is exact. The factor g is forced, not
optional: writing the pre-cycle shift originally as bQ^m, the full
integrality congruence forces g|b, and b=g*u yields (4).

**Cycle legality and final division.** The first m-1 repetitions of w
are exact. At a prefix of i letters of the repeated cycle, with cumulative
division b_i<Bm, its state is the corresponding negative-cycle state plus
g*u*3^i*2^(Bm-b_i). The exponent is at least one, so this state is odd;
consecutive such states certify exact earlier valuations. In the final
copy all but its last letter are therefore exact. If the last nominal valuation is d, its
actual numerator at that position is 2^d*z. Hence the last valuation is
d+tau, tau=v2(z)>=t, and (6) is the actual endpoint. This explicitly
retains the extra division instead of treating the affine quotient as odd.

**First descent and payment.** Every legal positive prefix before the
last step increases its block source: its multiplier exceeds one and
its affine carry is positive. The head grows above n, and each full
cycle and final partial cycle grows above that larger current source.
Finally

    M0*(2^t*n-z)=u*D-(2^t*C0-M0*c)>0                (8)

by u>=1 and (3). Thus oddpart(z)<=z/2^t<n, giving exactly the claimed
first-descent time. Positivity and the cylinder step in (7) are immediate
from (2), (4), (5); also 0<n0<2^(A+Bm+t).

For fixed head and cycle, v2(M0*n+C0)=A+Bm, so different compiled m
are disjoint. Moreover all sources with m>N lie in the single parent
cylinder M0*n+C0=0 mod2^(A+B(N+1)). Its density tends to zero with N.
Consequently the all-height union has an ordinary natural density equal
to the sum of its cylinder densities; countable additivity is not being
assumed for arbitrary sets with natural densities.

The input domain is deliberately explicit. For instance the head (1,1,2),
cycle (1), m=1 fails (2): H*Q=32<C0=35. The implementation records this
as outside its sufficient domain, without declaring that head impossible.

## 2. One -5 block followed by the -17 return

Take head a=(1,2), so p=2,A=3,S=5. Take

    w=(1,1,1,2,1,1,4), c=17, P=2187, Q=2048.

Then C=141, g=3, M0=3, C0=47. The theorem becomes particularly simple.
For m>=1 let t_m>=1 be least with

    2^t*(8Q^m-47)>9P^m-51.                          (9)

Its right-side debt 47*2^t-51 is positive, so this is exactly (3).
Choose positive odd u with

    uQ^m=1 mod3,   3uP^m=17 mod2^t.                 (10)

Then

    n=(8uQ^m-47)/3,
    U^2(n)=3uQ^m-17,
    first_descent(n)=7m+2,
    U^(7m+2)(n)=oddpart(3uP^m-17)<n.                (11)

The pure -17 return bound would only compare the endpoint with U^2(n).
Inequality (9), or equivalently (8), pays the smaller original n.

The first member and its actual odd route are

    27291 ->40937 ->30703 ->46055 ->69083 ->103625
          ->77719 ->116579 ->174869 ->8197.

Here m=1,t=1,u=5, and the whole cylinder is n=27291 mod32768.
The rational shadow center is -47/3: the affine head sends it to -17.
This is a finite preperiodic head followed by the inherited negative cycle,
not a newly claimed negative cycle or a positive divergent orbit.

The parameters are recoverable from the actual source:

    v2(3n+47)=11m+3,
    u=(3n+47)/2^(11m+3).                            (12)

Require a positive integer m, compute (9), then check (10). No successful
orbit search is assumed by this membership test.

## 3. Additional coverage, beyond both prior recursive families

The mixed cylinders are mutually disjoint by (12). Every m>=1 lies in
n=667 mod2048, since n=-47/3 modulo that power of two. **FINITE-EXACT
bank certificate:** all 171 old rows miss this entire parent cylinder;
the same exact intersection test holds for their 65 disjoint cylinders.
The script regenerates the rows and compares every one with the retained
[old table](reset_20260926_swaplift.out). This finite certificate proves
the all-m disjointness from that specified bank.

They are also disjoint from both previous infinite extensions:

- 3(n+5)=8*(uQ^m-4), hence v2(n+5)=5. The old -5 return family requires
  a positive multiple of 3 for this valuation.
- The mixed sources are 3 mod8; the pure -17 family is 7 mod8.

The added density among all positive integers therefore is

    delta_mix=sum_(m>=1) 2^(-11m-3-t_m)
             = approximately 3.053248656570591e-5.   (13)

For existence, the tail m>M lies in a single class
3n+47=0 mod2^(11(M+1)+3), whose upper density tends to zero. For an
explicit numerical tail, (9) gives

    2^t > (9/8)*(P/Q)^m,
    sum_(m>M) 2^(-11m-3-t_m) < 1/(9*(P-1)*P^M).     (14)

The first inequality follows by cross multiplication from
423P^m>408Q^m. The output retains the exact m=1,...,20 fraction and
strict tail bound, and independently normalizes the finite union with
the old bank. These are source-set densities, not orbit visit rates.

**Weak-budget hostile.** At m=9 the required budget is t=2. The weaker
t=1 condition accepts u=11. Its source and selected step-65 endpoint are

    18592208803347364555284980345499,
    18885261011608818665618169991037.

Every step through 65 remains above the original source. The legal word
survives; payment (8) fails. Later descent is not excluded.

## 4. A bounded compiler search retains one stronger additional head

**FINITE-EXACT search.** Enumerate heads of length 1 through 4, letters
1 through 4, retaining those whose every prefix expands. There are 13
heads. Compile each with w and m=1,2,3: 39 candidate cylinders, each with
first-descent time at most 25. Compare with the old 171-row bank and the
entire -5, -17, and mixed parameterized families. This infinite-family
comparison is exact: a candidate of first-descent time at most 25 cannot
intersect a certified family member of greater first-descent time. Thus
only -5 parameters m<=12, -17 parameters m<=3, and mixed parameters m<=3
need comparison. The inherited finite union is deduplicated before mass
subtraction. Twelve heads add coverage; head (1,2) is already retained.
The output records the full ranked table, rather than selecting examples
after observing their individual traces.

The best head is a=(1). Its three tested cylinders add exactly
4196353/34359738368; the next head (1,1) adds one quarter of that. We retain
only this strongest head as a second all-height family. With P=2187 and
Q=2048, define t_m>=1 to be least satisfying

    2^t*(2Q^m-35)>3P^m-51.

Take every positive odd u satisfying

    uQ^m=1 mod3,   uP^m=17 mod2^t,
    n=(2uQ^m-35)/3.

**PROVED:** every such n has exact first odd descent at 7m+1, with endpoint
oddpart(uP^m-17). This follows from the compiler with H=2,M=3,S=1,
C=C0=35,g=1,M0=3; both positivity conditions hold for every m>=1.
The debt 35*2^t-51 is positive for t>=1, so the displayed budget is
exactly the compiler condition. The family is one full cylinder of
modulus 2^(11m+1+t_m), and the source itself decodes its block count:

    v2(3n+35)=11m+1.

The first class is n=6815 mod8192, with first descent at step 8;
6815 reaches 5459 there. All heights lie in n=2719 mod4096, a parent
class missed by every row of the regenerated old bank. They are 7 mod8,
whereas both -5 and mixed sources are 3 mod8. Finally

    3(n+17)=2uQ^m+16

has valuation 4, excluding the pure -17 family, whose v2(n+17) is a
positive multiple of 11. These are exact all-height separations, not a
finite missing-sample claim. Distinct m are disjoint by the decoder.
The additional source-set density is

    sum_(m>=1) 2^(-11m-1-t_m) = 0.00012212994626282017... .

The script stores its exact first-20 partial sum, with strict tail bound
1/[3(P-1)P^20]. For general M the tail is less than
1/[3(P-1)P^M]: cross multiplication gives
(3P^m-51)/(2Q^m-35)>(3/2)(P/Q)^m because 105P^m>102Q^m;
each density term is consequently less than 1/(3P^m).

**Weak-budget hostile.** At m=5, the required t is 2. Weakening it to 1
admits u=5 and source 120095990063213215. Every odd step through the
selected step 36 remains above that source; its step-36 endpoint is
125078862747499259. The fixed legal word alone again fails to pay for
the original source. This family is a head into the same -17 shadow,
not a new negative cycle or a universal coverage result.

## 5. The golden polynomial certificate, with its positions retained

For a valuation word v, let L_v be its ordinary-map length, the sum of
one odd operation plus its divisions for each letter. Let H_v(z) have
coefficient one at the ordinary times occupied by odd sources. Then

    H_(uv)(z)=H_u(z)+z^L_u H_v(z).

For the compiled head and cycle, the nominal polynomial is

    H(z)=H_a(z)+z^L_a H_w(z)*(1+z^L_w+...+z^((m-1)L_w)). (15)

The actual last division adds tau=v2(z_terminal) trailing zero bits, so
the actual length is L=L_a+mL_w+tau. If r is the actual endpoint, then

    F_n(z)=H(z)+z^L F_r(z),
    Theta(n)=beta*H(beta)+beta^L*Theta(r), beta=phi^-1. (16)

This is a proved finite-head identity, independent of a claim about the
unknown future of r. If a certified suffix to 1 is later supplied, its
polynomial certificate composes through the same identity.

For the mixed family,

    L_a=5, H_a=1+z^2,
    L_w=18, H_w=1+z^2+z^4+z^6+z^9+z^11+z^13.

At source 27291, tau=2 and L=25. Exact arithmetic in Z[phi] gives

    Theta(27291)=(-9971+6163phi)
                 +(-121393+75025phi)*Theta(8197).    (17)

The second factor is phi^-25. The script compares literal ordinary-map
bits, polynomial positions, and two exact evaluations in the pair basis
(A,B)=A+Bphi. All 638 instances are additionally checked by the separate
[`golden_digit_carry_20261003.py`](../../04-computation/experiments/golden_digit_carry_20261003.py)
Laurent-polynomial reader and its signed-power function. The public
functions `parity_positions`, `golden_fraction`, and `compile_family`
provide reusable inputs for the same comparison.

**Position warning.** The algebraic equality 11_phi=100_phi preserves an
evaluated number; it need not preserve a timed parity word, its length,
or an actual Collatz edge. Ordinary Collatz parity words already forbid
adjacent ones. Literal source digits, Zeckendorf digits, and route parity
positions remain different typed inputs. The paired arithmetic and golden
certificates keep precisely the sidecars that a value-only rewrite loses.

## Verification and open boundary

Run:

    python -X utf8 -B 04-computation/experiments/collatz_mixed_return_compiler_20261003.py
    python -X utf8 -B -O 04-computation/experiments/collatz_mixed_return_compiler_20261003.py

The [saved output](collatz_mixed_return_compiler_20261003.out) records 200
mixed instances (m=1,...,50; parameter lifts 0,1,7,29), 200 winning-head
instances over the same parameter range, and 238 generic
instances from the printed five-head/three-cycle matrix, m=1,...,8 and
lifts 0,2; one parameter triple is outside (2). It also records three
compiler rejection controls, all 171 saved-table agreements, the exact
bank unions, the exact bounded head ranking, both weak-budget hostiles,
and the first complete mixed trace.

Every accepted instance is checked by literal odd and ordinary iteration,
generic affine composition, and three exact golden polynomial evaluation
paths, including the independently implemented Laurent reader. All
checks remain active under optimized Python. The theorem establishes all
positive members of its guarded families; the finite replays are controls,
not its quantifier justification. Universal successful certification and
coverage remain OPEN.
