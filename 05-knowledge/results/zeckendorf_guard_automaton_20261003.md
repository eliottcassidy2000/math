# Fibonacci digit registers retain arithmetic guards and the color charge

2026-10-03 (America/Denver).

**Status:** PROVED elementary register, charge, typed-representation, and
fixed-modulus automaton identities; FINITE-EXACT for the stated audits.
The owner has confirmed the red/black/blue diagonal pattern with extra
copies of 1. The XOR charge below is an auxiliary invariant, not that color
convention. The supplied finite stripe is preserved; its complete infinite
row-continuation rule is not inferred from a near-fit. No Collatz coverage
or new theorem ID is claimed.

## Recovered definitions and corrections

There are three different objects in the inherited work:

1. The owner's authoritative 35-symbol stripe, with K=black:

       RKBRKRKBRKBRKRKBRKRKBRKBRKRKBRKBRKB.

2. The Wythoff candidate: in ordinary Zeckendorf notation using weights
   1,2,3,5,..., let z be the number of final zero digits. Assign R for z=0,
   K for odd z, and B for positive even z. It matches positions 1--34,
   but gives R at position 35, where the supplied word gives B.
3. The inherited charge Q(n)=(n-2b(n),b(n)), where
   b(n)=floor((n+1)/phi^2), and its reduction c(n) in F2^2. Its three
   nonzero atom colors accompany a fourth neutral state. This is not the
   supplied stripe: 2 and 4 both have charge (0,1), but supplied colors K/R.

The primary local truth sources are
[the color interpretation audit](crossroads_crossing_20260926_colour.md),
sections 1--4;
[the marked-unit construction](reset_20260926_colours.md), sections 1--4;
and [the charge/carry note](duck_zeckendorf_20260925.md), Z1--Z14.
The marked alphabet 1_K,1_B,1_R,2,3,5,... makes the supplied row 35 legal,
but its stated descending-piece grammar cannot continue that row coherently
to row 36. That is a scoped grammar obstruction, not permission to change
the user's data. Relevant corrections in `01-canon/MISTAKES.md` are the
2026-09-26 crossing-crossroads entry and the representation-boundary entry.

Closest mechanism: Fibonacci recurrence plus unique no-adjacent-one normal
forms. Canonical hostile: 2=2=1+1 changes charge under repeated atoms.
Corrected near miss: treating a 34-term near-fit as the owner's full rule.
Least-used coordinate: the integer value after shifting every Fibonacci
digit by one place. The live comparison is between the supplied stripe,
the four-state charge, two-register arithmetic, and the route's actual
residue guard. Their projections are different.

## A faithful reader for the confirmed extra-unit construction

Use the inherited ordered atoms g_0=1_K, g_1=1_B, g_2=1_R and g_i=F_i
for i>=3. Legal representations occupy distinct nonconsecutive indices.
Let Z(n) be ordinary Zeckendorf support in indices>=2. For n>0 the full
representation fibre is exactly

    Z(n),
    {0} union Z(n-1),
    {1} union Z(n-1), if 2 is absent from Z(n-1).        (E1)

Thus numerical value together with the marker ordinary/K/B recovers the
entire support. This is a lossless map of typed representations; forgetting
the marker is a quotient with two or three elements per fibre.

For the marker polynomial, let r(n)=1 when index 2 belongs to Z(n), else 0.
Tracking the exposed unit atoms, rather than recursively expanded words,
gives exactly

    C_n(u)=u_R^r(n)+u_K*u_R^r(n-1)+(1-r(n-1))*u_B.      (E2)

The ordinary source n and the ordered unit positions are retained. This
is a polynomial on a representation fibre, not a single color of n.

There is a native arithmetic reader for these colored representations.
Read its ordinary part in indices>=2 using the registers (X,Y) defined
below, and retain its last digit r. If an extra black or blue unit is
present, then

    (X,Y) -> (X+1,Y+2-r).                              (E3)

Indeed b(a+1)-b(a)=r for the ordinary-part value a, so its shifted value
changes by 2-r. A blue marker requires r=0; a black marker permits either
r. The marker is kept separately. For the inherited full charge, assign
red/blue unit charge U=(1,0), black unit charge V=(-1,1). The raw-charge
defect relative to canonical Q(n) is (-2D,D), where

    D=1-r for a black marker, and 0 otherwise.          (E4)

Exactly two sheets preserve full charge. Red and blue remain distinct
typed units despite having equal U. Formula (E3) lets actual colored
representations feed exact modular guards, while (E1) and the marker
retain information that normalization would erase. No color is identified
with a residue or a charge without this explicit map.

For the preferred recursive piece expansion, start W_0=K,W_1=B,W_2=R,
W_3=RK,W_4=RKB and W_i=W_(i-1)W_(i-2) for i>=5. Expand support in
descending index order. The supplied row 35 is exactly W_9 B, from
34+blue1. At row 36 the two legal supports are 34+2 and 34+red1+black1;
both expand as W_9 RK. Hence they force red at position 35. Retaining
row identity, unit marker, and diagonal address faithfully records this
boundary choice. It does not manufacture a prefix-coherent continuation
that this particular grammar forbids.

The actual three labels have an ordinary C3 Fourier basis: keep, at every
retained position d, the three indicators f_d(K),f_d(B),f_d(R), and form
sum_c f_d(c) omega^(-rc), r=0,1,2, with omega^3=1 and omega!=1. The full
complex coefficients recover each indicator by inverse Fourier transform.
Cosine alone or squared magnitudes conflate blue and red in this frame.
Counts also discard positions. This is a reversible basis change on label
data, not a claim that rotating labels preserves representation legality:
black+red is legal, while black+blue occupies adjacent indices and is not.

## Exact two-register digit arithmetic

Put f_i=F_(i+2), so f_0=1 and f_1=2. For a digit prefix read high to low,
retain

    X=sum_i e_i f_i,    Y=sum_i e_i f_(i+1).

Appending a low digit d shifts every existing weight and then adds d:

    (X,Y) -> (Y+d, X+Y+2d).                            (1)

This follows directly from f_(i+2)=f_(i+1)+f_i. Starting from (0,0),
the final X is the represented integer. A previous-digit bit rejects 11;
leading zeros do not change the state.

The inherited exact shift identity is

    Y=2X-b(X),
    Q(X)=(2Y-3X,2X-Y),
    c(X)=(X mod2,Y mod2).                              (2)

Thus the familiar color action is the modulo-two shadow of (1):

    c -> M c+d(1,0),    M(a,b)=(b,a+b),    M^3=I.

This is a precise connection to integer arithmetic. It does not identify
the finite stripe with c. On arbitrary integer addition, the inherited
carry delta(a,b)=b(a)+b(b)-b(a+b), in {-1,0,1}, is still required:

    c(a) XOR c(b) XOR (0,delta mod2)=c(a+b).             (3)

For the strict diagonal 5, the pairs 1+4 and 2+3 have naive XOR colors
(1,1) and (1,0). Their corrected colors agree. Neutral occurs at n=6.

For every fixed modulus m, reduce both registers in (1) modulo m. Together
with the previous-digit bit this gives an exact finite automaton with at
most 2m^2 states, recognizing residue guards on canonical Fibonacci words.
For odd m, working modulo 2m also retains the charge by the Chinese
remainder theorem. A separate four-state color factor is equivalent data.

The second register cannot be discarded even after keeping the color:

    n=2: digits10,    (X mod3,Y mod3,previous,c)=(2,0,0,(0,1));
    n=8: digits10000, (X mod3,Y mod3,previous,c)=(2,1,0,(0,1)).

Appending 0 gives the values 3 and 13, whose residues modulo 3 differ.
The source, target, and preserved predicate of the connection are now
explicit: canonical digit words map to value/shift residue states, and
acceptance preserves the specified numerical congruence. Finite residues
forget height; neither height nor route termination follows from color.

## Connection to guarded route growth

The [route construction](collatz_atom_route_memory_20261003.md) uses
compressed choices j_i with actual valuation
k_i=2j_i+k0(m_i), where k0 is determined by the odd target modulo 3.
For a block with p preceding odd steps, its lawful fixed-head growth is

    j_p -> j_p+3^p.

The digit automaton with m=3^p recognizes the exact phase equality here.
It gives an independent numeral representation of the existing guard;
its input length is unbounded. A fixed depth has a finite automaton,
whereas arbitrary depths require growing precision or a parameterized
grammar. This construction does not give a fixed finite color clock
closed under all arithmetic operations.

The hostile boundary remains precise: a different phase can still produce
a valid route while changing preceding valuation letters. The gate protects
the retained head, not mere membership in the convergent basin. The audit
includes both rejected routes and valid changed-head routes.

## Verification

Reproduce with:

    python3 04-computation/experiments/zeckendorf_guard_automaton_20261003.py

Companion: [exact output](zeckendorf_guard_automaton_20261003.out).
The script imports the route decoder definitions and the separate guard
compiler from `04-computation/experiments/collatz_route_tournaments_20261003.py`
and `04-computation/experiments/collatz_fourier_guards_20261003.py`. It
independently reads digits, computes shifted weights, enumerates typed
supports, and replays valuation prefixes. It checks:

- Every integer 0..10000 at moduli 3,9,27,81,4,16,64: 70,007
  value/shift/charge cases, with 40,004 separate odd-modulus CRT cases.
- The 2/8 missing-register hostile; illegal 11; the supplied stripe and
  its sole Wythoff mismatch at 35; the neutral state and charge/stripe split.
- Every colored support of values 1..200, by independent exhaustive subset
  enumeration: 524 representations, 1,572 native-reader checks, the marker
  polynomial, and exactly two charge-preserving sheets. It verifies the
  legal supplied row 35, the scoped row-36 obstruction, and exact actual-label
  Fourier inversion at all 35 positions.
- All positive summand pairs a,b<=200: 40,000 carry corrections.
- All compressed words of lengths 1..5 with letters 0..4: 694 valid and
  3,211 rejected words. All 3,174 word/block macro phases agree through the
  independent digit path and preserve directly replayed valuation prefixes.
- 2,480 wrong-phase controls: 780 still encode valid routes but change
  their head. This checks the exact preservation boundary.
- 14,760 compiler/digit-automaton/strict-decoder interface cases over all
  compressed heads of lengths 1..4 with letters 0..2 and two phase periods.
  Exactly 45 arithmetic-legal index-zero cases require the isolated
  first-hit exclusion. The whole zero residue class must not be removed.

Independent theorem audit: the guard compiler has exactly 2^p arithmetic
phases modulo 3^p because each of the 2^p choices of odd/even valuation
parity yields one landing residue; distinct choices cannot give the same
residue because legal backward decoding determines every exponent parity.
The sibling chart permutes residues modulo 3^p. For target 1, index 0
creates a selected 1->1 padding block regardless of its earlier prefix;
positive indices in the same residue class remain valid. The independent
pre-implementation check covered 605 prefix/target cases and 36,905 phases
at targets 1,5,7,11,13; all counts agreed. This preliminary check was an
ephemeral audit; the 14,760-case interface check above is reproducible in
the saved script.

All checks passed. The output contains no sampled claim of color-based
descent, independence, or all-source coverage. The useful extension is a
digit-language guard product retaining exact arithmetic coordinates.
