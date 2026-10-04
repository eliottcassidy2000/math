# Golden lattices, ordered carries, and certificates over integers

2026-10-04. Owner-directed continuation of the atom, golden-ratio,
{2,3,11}, and enriched-Collatz-state investigations.

**Status:** PROVED for the elementary identities and carrier soundness;
FINITE-EXACT for the complete stated golden-lattice censuses and their
independent word-enumeration replay; OPEN for the proposed universal normalizer.
Neither positive Collatz nor the assertion that every negative integer
reaches one of the three known negative cycles is proved here.
No novelty claim or theorem ID is assigned.

## 1. Inheritance and the working board

The closest proved mechanisms are:

- [THM-4528, golden parity and holonomy](../../01-canon/theorems/THM-4528-collatz-parity-in-base-phi-golden-beta-map-and-holonomy.md):
  ordered ordinary-Collatz parity words, their golden evaluation, and the
  positive denominator-two criterion. This note retains its map convention.
- [THM-4512, exact coefficient-descent cylinders](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md):
  the final oddness bit is essential; modulus 2^(A+1), not just 2^A.
- [THM-4507, finite-valuation obstruction](../../01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md):
  arbitrary nonlinear functions of finitely many polynomial valuations,
  including prime bank {2,3,11}, cannot repair logarithmic height into a
  globally nonincreasing potential. Its scope includes fibre-bounded colours.
- [The new signed-atom / route-memory note](collatz_atom_route_memory_20261003.md):
  guarded words, shared tails, ordered affine composition, and exact
  source/internal tournament growth. Its source code is a known-basin
  certificate; universal coverage is still an independent obligation.
- [The recursive entry board](entry_20260927_board.md): parameterized
  repeat rules already certify unbounded excursions, including a family
  through 27. Reuse these rules rather than rebrand their families.
- [The crossings note, section 1.2](collatz_crossings_20260926_potential_and_seeds.md):
  the exact odd-square-window characterization of {2,3,11}.
- Incoming commit `ad236753bf`, [modular/Fourier/Fibonacci route synthesis](modular_fourier_zeckendorf_routes_20261003.md),
  [marked digit guard reader](zeckendorf_guard_automaton_20261003.md), and
  [the -17 first-descent macro](collatz_minus17_return_20261003.md): exact
  phase, shifted Fibonacci register, and a retained original-source bound.
  These give concrete input and return interfaces for the proposed carrier.

Hostiles: the actual negative cycles and their arbitrarily long positive
2-adic shadows; the noninteger cycle through 1/13; the two words (1,2)
and (2,1), with identical clocks but different carries. Corrected near
misses are THM-4512's omitted oddness bit and MISTAKE-554/556's map/code and
golden-boundary distinctions. The useful underused sidecar is the complete
denominator ideal together with the signed affine integer-realization test.

| Lane / live concept | Object and operation | Missing information / decisive probe |
|---|---|---|
| Anchor: proof-carrying numerals | Guarded word plus endpoint, value map to Z | Universal construction; empty open state must exist for every n |
| Golden parity | Ordered word evaluated in Q(phi) | Denominator quotient loses order; retain exact phase, polynomial and boundary convention |
| Niche: denominator modules | G^L-I, Smith factors, denominator ideal | Golden periodicity need not be an integer orbit; 1/13 |
| Signed arithmetic | Carry polynomial, exact dyadic cylinder | Clocks alone lose order; (1,2) versus (2,1) |
| Wildcard: rational-anchor stack | Repeated cycle charts and live residual | Counter falls but residual grows; negative-cycle shadows |
| Two recursion orders | Affine composition / mixed-radix rewriting | Translation and guard debt; compute the interchange defect |
| Prime triple | Odd-square window, cycle clocks, golden splitting | Prime ordinal is not a dynamics invariant; retain separate maps |

Cards used: search the statement before the method; type every analogy;
find the hidden second coordinate; separate local support from coverage.
No new meta-pattern card is warranted.

## 2. Resolve the 12/8 and 11/7 observation

For a star b=a+b+2ab, the all-one interchange comparison is

    (1+1) star (1+1)=12,    (1 star 1)+(1 star 1)=8.

Subtracting 1 gives (11,7), precisely

    (S_+(7),S_+(9))=(11,7),   S_+(u)=oddpart(3u+1).

The direction matters: 7 goes to 11, and 9 goes to 7. The common atom is
N=2. The failed information quotient identifies 7 and 9 at N=2 while
their targets have different atom addresses.

This equality does not generalize along the uniform interchange ray.
At a=b=c=d=t>=1 put H=8t^2+4t and L=4t^2+4t. If
S_+(4N-1)=H-1, then 6N=H. The other target is at most
3N+1=4t^2+2t+1, strictly less than L-1 for t>1. Thus t=1 is the unique
positive uniform-ray alignment with the same-atom pair. This is a precise
boundary for this particular proposed connection, not a no-go for all
other parameterizations.

The values 7 and 11 also have a structural golden interpretation. With

    G=[[0,1],[1,1]],    phi=(1+sqrt(5))/2,    psi=1-phi,

the eigenvalues are phi,psi, and tr(G^k)=L_k is the Lucas number.
Consequently tr(G^4)=7, tr(G^5)=11. The similar matrix [[1,1],[1,0]]
is the no-adjacent-ones transfer matrix, so these traces count cyclic
ordinary-Collatz parity schedules of lengths 4 and 5. This counts schedules,
not integer orbits realizing them. The identities

    phi^8+1=7 phi^4,    phi^10-1=11 phi^5

are the k=4,5 instances of phi^(2k)+(-1)^k=L_k phi^k. These supply a
shared word/matrix carrier; the isolated numeral matching above does not
by itself intertwine the two dynamics.

## 3. The prime triple has several distinct roles

The inherited theorem says that {2,3,11} are exactly the primes p with
2p below the next odd square above p. If (2j-1)^2<p<(2j+1)^2, the condition
requires 2(2j-1)^2<(2j+1)^2, hence 4j^2-12j+1<0 and j<=2.
Checking those two windows gives precisely 2,3,11. Their ordinal positions
among primes are 1,2,5. The latter indexing is an exact numerical fact;
no map preserving Collatz reachability follows from that indexing.

Two other appearances have specified maps:

- The halving counts of the ordinary cycles through 1,-5,-17 are 2,3,11,
  respectively. Their odd-step counts are 1,2,7 and their total periods
  3,5,18. The -1 cycle separately has one odd step and one halving.
- The prime 11 is a golden denominator prime for the period-five word;
  its field-theoretic role is made exact below. The prime 3 is the affine
  multiplier, while 2 is the halving base and the positive golden denominator.

The fact that three different procedures select the same small numbers
does not identify their predicates. In particular the -17 golden denominator
introduces the additional prime 19, so a fixed {2,3,11} coordinate bank
would omit part of even the known-cycle data.

## 4. A golden denominator is conserved along a certified basin

Let C(n)=n/2 on evens and 3n+1 on odds, for signed integers. Put

    F_n(z)=sum_(j>=0) b_j z^j,  b_j=C^j(n) mod2,
    Theta(n)=phi^(-1) F_n(phi^(-1)).

No '11' occurs in an ordinary parity word. The exact golden coordinate with
the boundary convention below determines the word; its denominator or ideal
alone does not. A finite word certificate also supplies an effective check
that separately supplied arithmetic and golden annotations agree.

### Boundary repair for signed integers

Use the **upper** golden map on the closed interval:

    beta_hat(x)=phi*x-d,  d=1 if phi*x>1, else 0.

For every integer n, Theta(C(n))=beta_hat(Theta(n)). At phi*x=1 the two
expansions are a terminating word and the word with tail (10)^infinity.
A nonzero integer cannot reach 0 under C, so its actual expansion uses
the latter, exactly the upper convention. At n=0 both sides are 0.
This restricts the statement to integers; it does not silently extend to
all 2-adics. For example -1/3 has the terminating expansion and would
need the other boundary choice. The upper convention keeps -1 <-> -2
as Theta-values 1 <-> phi^(-1), repairing the specific old -2 exception.

Suppose Theta(n)=(a+b phi)/q with integers a,b, q>0 and gcd(a,b,q)=1.
One forward step sends

    (a,b,q) -> (b-dq,a+b,q).

The integral linear part is unimodular, so gcd(a,b,q) and the least
positive scalar denominator q are invariant. More precisely the ideal

    D(x)={u in Z[phi]: u*x in Z[phi]}

is invariant, since phi is a unit and d is integral. This statement is
conditional on a rational-golden coordinate existing; finite orbit
certificates supply it without assuming universal convergence.

The exact coordinates of the four known cycle roots are:

| Root | Raw parity period | Theta(root) | Exact q |
|---|---|---|---:|
| 1 | 100 | phi/2 | 2 |
| -1 | 10 | 1 | 1 |
| -5 | 10100 | (-1+7 phi)/11 | 11 |
| -17 | 101010100101010000 | (9+41 phi)/76 | 76 |

These values follow by summing the indicated periodic series, not fitting
decimals. Thus their denominators label their entire certified basins.

### The cyclic modules explain the denominators

The scalar recurrence has linear part G. Periodic golden states solve
(G^L-I)v = an integral digit vector. The Smith factors are

    L=2: (1,1),  L=3: (2,2),  L=5: (1,11),  L=18: (76,76).

These give the finite modules 0, (Z/2)^2, Z/11, (Z/76)^2. In particular
G^18-I=76 G^9, with G^9 unimodular. This is why 76=4*19 appears.

At 11, phi has roots 4 and 8 modulo 11; the first has multiplicative
order 5. The denominator ideal of Theta(-5) is (11,phi-4), of norm 11:

    (phi-4)*Theta(-5)=1-2phi.

This is an actual order-five / prime-eleven connection. It is still not
enough to specify an integer basin. The legal period-five word 10000 gives

    arithmetic root 1/(16-3)=1/13,
    golden coordinate (1+4phi)/11,
    (phi-4)*(1+4phi)/11=-phi.

It has the same scalar denominator AND the same denominator ideal as -5.
The missing test is signed integer realization of the ordered word.

## 5. A finite golden gate for the signed problem

**PROVED reduction + FINITE-EXACT terminal classification.** For a nonzero
integer n, the following are equivalent:

1. Its C orbit reaches the cycle through 1, -1, -5, or -17.
2. Theta(n) belongs to Q(phi) with exact scalar denominator in {1,2,11,76}.

(1) implies (2) by the root table and invariance.
For (2) implies (1), conjugation gives
x'_next=psi*x'-d. Since |psi|<1, the conjugate eventually enters |x'|<=3
and stays there. With fixed q and 0<=x<=1, this is a finite lattice set.
Every periodic point already satisfies |x'|<=phi^2<3, so the set contains
all possible cycles. Enumerate it using exact signs of quadratic surds:

| Exact q | States in 0<=x<=1, abs(x')<=3 | Golden cycles | Period counts | Integer-cycle roots |
|---:|---:|---:|---|---|
| 1 | 4 | 2 | 1:1, 2:1 | 0, -1 |
| 2 | 8 | 1 | 3:1 | 1 |
| 11 | 321 | 13 | 5:2, 10:11 | -5 |
| 76 | 11592 | 240 | 18:240 | -17 |

This table uses exact denominators, the upper boundary convention, and
bound 3; it is not the older half-integral census with bound phi^2.
An independent enumeration of primitive cyclic no-'11' words of lengths
1..18 recovers exactly these cycle counts by evaluating their rational
golden series, without constructing the lattice graph.
The search box is exhaustive: b/q=(x-x')/sqrt(5) has absolute value below
2, and a/q=x-phi*b/q has absolute value below 5. Hence the script's ranges
abs(b)<=2q and abs(a)<=5q contain the whole stated lattice set.

For each cycle word, write its actual arithmetic map as
(3^p n+B)/2^A. Its only possible arithmetic root is

    n=B/(2^A-3^p).

Every periodic raw parity word has an odd denominator and the stated
2-adic parity; the script verifies the full rational orbit. The integer
filter leaves exactly the table's roots. Eventual periodicity of an integer's
parity forces it onto this rational cycle: matching k repeats fixes its
2-adic residue to precision kA, and letting k grow forces equality.
The zero parity orbit corresponds only to n=0, excluded here.

This gate gives an exact target for richer structures. It does not prove
the denominator-spectrum assertion for arbitrary integers; that assertion
is now the explicitly named universal obligation.

## 6. Proposed carrier: a guarded word with three synchronized descriptions

Call a finite object a **ribbon** for this note. It stores:

    binary sleeve e >= 0;
    positive valuation word w=(a1,...,ap), possibly compressed;
    signed odd endpoint m;
    exact input guard;
    ordered carry polynomial, raw parity polynomial, and their summaries;
    optionally a checked cycle seal or a pointer to a smaller certified source.

Set A_i=a1+...+ai, A=A_p, A_0=0. The ordered carry polynomial is

    B_w(X,Y)=sum_(i=0)^(p-1) X^(p-1-i) Y^(A_i),
    B=B_w(3,2),    Q_w(u)=(3^p*u+B)/2^A.

The polynomial, together with A, recovers the entire exponent word. Its
scalar evaluation is the affine carry; (A,p) alone does not retain it.
Words (1,2) and (2,1) both have (A,p)=(3,2), but B=5 and B=7.
For u followed by v, with summaries (A,p,B) and (D,q,E), composition is

    (A+D,p+q,3^q*B+2^A*E).

The exact guard is the odd residue class

    u=(2^A-B)*3^(-p) mod 2^(A+1).

The object projects to ONE integer by

    pi(ribbon)=2^e*(2^A*m-B)/3^p,

requiring the inner value to be a nonzero odd integer and the guard to
hold. Equivalently replay the inverse steps (2^a*m-1)/3, checking odd
integrality at each step. Signed sources use the same formula throughout;
the map on positive magnitudes of negative sources is 3n-1.

**Unconditional coverage of open states:** every nonzero integer has the
empty-word representation with its own odd part as endpoint and its usual
binary sleeve. Thus the ambient carrier is not defined only for numbers
already known to converge. The empty representation carries no arrival proof.

**Exact Collatz on the carrier:** halve a nonempty binary sleeve by reducing
e by one. When e=0 and w starts with a, remove that first word letter and
put a on the binary sleeve of the remaining ribbon. The represented value
becomes 3n+1. An open empty word can first be extended by one explicitly
computed guarded edge without changing its projected integer.

For a cycle-sealed ribbon stopped at its chosen odd cycle endpoint, the rank

    e + sum_i(ai+1)

falls by one at each ordinary step. It is a finite certificate rank, not
an unconditionally known function on every integer.

### The golden and binary descriptions must agree on the same word

The raw word associated to w is 10^a1 ... 10^ap, with polynomial P_w(z).
If its endpoint has a checked raw cycle polynomial P_c and length L_c,
the finite polynomial

    R(z)=z^e*((1-z^L_c)*P_w(z)+z^(A+p)*P_c(z))

certifies F(z)=R(z)/(1-z^L_c). It has three checks:

1. At z=2, read the rational value 2-adically as the parity coordinate.
2. At z=phi^(-1), obtain the golden coordinate and its denominator module.
3. From B_w(3,2), check the signed integer edge equation and the endpoint.

The common ordered word is the joining data. Independently supplied scalar
annotations are not a certificate that they encode the same legal path. This is the
concrete lesson of both the interchange defect and the same-atom split.

The prototype writes verified examples for 7,9,27,-3,-9,-27,-7,-25,151,
64,54,-18, including nonempty binary sleeves.
For instance -9 seals at -7 in the -5 cycle and has q=11; -25 is already
on the -17 cycle and has q=76. The polynomial and affine decoders agree
with an independently propagated exact golden coordinate. Evaluation at 2
also matches 64 parity bits obtained directly from each represented integer.

## 7. Creative extensions, with cheap tests and stopping boundaries

### A concrete bridge to the incoming Fibonacci reader

The incoming exact reader appends a digit d by

    (X,Y) -> (Y+d,X+Y+2d)=G*(X,Y)+d*(1,2).

Its full integer charge has an especially useful interpretation here:

    H(X,Y)=(2Y-3X)+(2X-Y)*phi=(X+phi*Y)/phi^3.

Direct substitution proves H(next)=phi*H+d. The coefficient change
[[ -3,2 ],[ 2,-1 ]] is unimodular, so H retains both integer registers.
The golden orbit update is instead Theta(next)=phi*Theta-d.
These are two actions of the same matrix, with different digit insertions.

For an ordinary parity prefix w=d_0...d_(L-1), feed THAT time word to the
Fibonacci reader, giving (X_w,Y_w). Then the exact joining identity is

    phi^L*Theta(n)=H(X_w,Y_w)+Theta(C^L(n)).

Proof: H=sum_j d_j*phi^(L-1-j); split the defining series for Theta at L.
Leading zeros are harmless. All 12 example ribbons verify this identity
through both the register reader and direct golden polynomial arithmetic.
The actual Zeckendorf digits of n are a different word. One can use them
to read source guards, while the parity-prefix reader certifies the time
word. A marked source representation keeps the extra unit marker explicitly;
neither digit sequence is inferred by identifying the two.

This suggests a two-tape carrier: the first tape stores an input numeral
and its exact modular reader; the second stores the ordered route and its
golden prefix/tail identity. Guarded affine cells verify their connection.
It extends the inherited typed certificate compiler with a destination
module, while retaining the growing precision needed for unbounded routes.

### Five compatible storage choices

**A. Golden lattice passport, attached to the carry ribbon (preferred).**
Store the denominator ideal, the numerator phase in its finite module,
and both real embeddings. The conjugate contracts and permits a finite
terminal search once q is fixed. Cheapest test: reject the 1/13 cycle
with the signed integer equation while accepting -5. Passed. Next obligation:
produce the passport from the input without computing an already-assumed
finite trajectory. The bounded golden module is a verification target,
not an oracle allowed in the initial state.

**B. A stack of rational-anchor charts.** For any word with fixed point
r=B/(2^A-3^p), use displacement n-r instead of n. The block becomes
n-r -> (3^p/2^A)(n-r). For an odd integer cycle anchor r and its word,

    n=r+2^(Ak+1)*b -> after k blocks -> r+2*3^(pk)*b.

Every indicated exponent is exact: at each unfinished block its displacement
from r is divisible by 2^(A+1). The live register update is
(k,b)->(k-1,3^p*b). This handles arbitrarily long shadows with one repeat
node. It is checked for anchors -1,-5,-17, both signs of b, and k=1..4.
It does NOT prove descent: the countdown decreases while b grows. A global
proof needs a well-founded rule across chart changes. Let the anchor library
grow to rational fixed points of newly encountered words; retaining only a
fixed prime bank would fall back inside THM-4507's obstruction.

The incoming -17 return theorem supplies a real, all-height exit rule for
some of these charts. With P=3^7, Q=2^11, choose the least t>=1 such that
2^t*(Q^m-17)>P^m-17 and require b*P^m=17 mod2^t for positive odd b.
Then n=b*Q^m-17 first descends after 7m odd steps to oddpart(b*P^m-17).
Store m, the guard, the original-source inequality and a lower-source proof
pointer. At m=2,b=1 the certificate is 4194287 -> 597869 in 14 odd steps.
Thus the stack can already compress a proved infinite family with a
source-decodable guard; it need not merely compress a trajectory after it
has been run. Universal coverage still requires further exit rules.

**C. A mixed-radix proof diagram.** Keep binary deletion and ternary append
as local actions, with carries at crossings. Adjacent base exchanges are
administrative rewrites preserving pi; an actual Collatz step changes pi.
The established precedent is Yolcu–Aaronson–Heule, *An Automated Approach
to the Collatz Conjecture*, v3, section 3.2 and Theorem 3.17:
[primary PDF](https://arxiv.org/pdf/2105.14697). Their mixed binary–ternary
termination equivalence is not a solved termination proof. Add the golden
module as a checked terminal annotation, not as an assumed ranking function.
Cheapest test: distinguish value-preserving base exchange from parity-word
shortening. An arbitrary fair scheduler or contracting conjugate alone
does not establish termination of the whole proof diagram.

This is also the useful replacement for a strict Eckmann--Hilton analogy.
Two accelerated letters have the exact order defect

    Q_b(Q_a(n))-Q_a(Q_b(n))=(2^a-2^b)/2^(a+b).

For a=1,b=2 this is -1/4, matching the carries 5 versus 7 over denominator
8. A tile crossing can record this affine defect and the guard that makes
the chosen route legal. Swapping tiles is a rewrite with a recorded carry,
not a commutative interchange law. The four-input transported-multiplication
defect and this two-letter affine defect are different formulas with the
same design lesson: the crossing needs a stored correction, not erasure.

**D. A graph or surface realization.** Use the existing ordered-SCC tournament
grammar, or a strip whose cells are the affine edge equations, with an actual
terminal boundary. Shared strips are a DAG of common suffixes; periodic ends
are sealed annuli. The map to an integer is the same guarded decoder.
Topology must retain the boundary phase and ordered carry. Counting holes,
cycles, or tiles without those labels cannot distinguish the 1/13 hostile.
No assertion that Collatz components require a tournament is made.

**E. A polynomial proof obligation rather than a guessed scalar potential.**
Keep R(z) unsealed while extending a prefix. A completed cycle supplies the
factor 1-z^L and an exact identity, making the result finitely checkable.
Parameterised repeat rules and the existing entry grammar can share proofs
over infinite families. The difficult part is finding a normalization rule
that always closes the obligation or reduces to a strictly smaller certified
integer. Finite-prefix matching or density-one coverage is insufficient.

## 8. A noncircular specification for the desired normalizer

The strongest concrete target is a terminating construction from a signed
binary numeral n of a rational pair t(n) in Q(phi), with exact bounds
0<=t(n)<=1 and the functional identities

    t(2n)=phi^(-1)*t(n),
    t(2n+1)=phi^(-1)+phi^(-2)*t(3n+2),

using the nonterminating boundary convention for nonzero integers. A bounded
solution has to equal Theta(n): iterate t(n)=phi^(-1)*(b(n)+t(Cn)) and
the remainder vanishes geometrically. If the constructor always produces
q in {1,2,11,76}, the finite gate proves signed convergence to the four
known cycles; sign preservation then gives the positive target 1.

This is a sufficient proof specification, not a claimed available algorithm.
The odd recursive call grows, so the two identities by themselves do not
define a terminating procedure. An input-bit induction would need the
unbounded carry stack or another well-founded structure to handle that call.
The candidate program is therefore a small guarded grammar with integer
registers, shared suffix proofs, and a golden-module verifier. Universal
coverage and a rank across changes of chart are its two open obligations.

## Reproduction and audit scope

Run:

    python3 04-computation/experiments/collatz_golden_carriers_20261004.py --write-certificates 05-knowledge/results/collatz_golden_carriers_20261004.json

[Source](../../04-computation/experiments/collatz_golden_carriers_20261004.py),
[output](collatz_golden_carriers_20261004.out),
[example certificates](collatz_golden_carriers_20261004.json).

No random sampling or inherited exclusion filter. Complete universes:
exact-q golden traps for q=1,2,11,76; independent primitive no-'11' cyclic
words through length 18; all integers -10000..10000; all 340 positive
valuation words of lengths 1..4 and letters 1..4, with three signed
cylinder representatives each. The 40401 surd-sign checks use a separate
rational isolating interval for sqrt(5). All assertions remain active
under Python -O. This is two computational paths by one author, not an
independent human or agent audit.

The finite integer census has 10000 positive sources reaching 1 and
negative counts 3244,3213,3543 at roots -1,-5,-17. Those counts are controls,
not a coverage theorem. The golden-lattice gate uses its complete finite
universe plus the conjugate-contraction proof; the unresolved part is
placing every signed integer's parity value in that universe's denominator
spectrum. No root, prime index, local density, or graph count pays that debt.
