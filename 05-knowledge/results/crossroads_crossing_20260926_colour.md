# The triangular colour prefix, two hidden unit types, and a bounded-potential obstruction

**Status: FINITE-EXACT interpretation audit; PROVED elementary constructions
and obstructions. The continuation of the user's 35-colour sequence is OPEN.**
No assertion identifies that finite word with a unique infinite sequence.
No Collatz termination theorem follows from a colour quotient.

The user explicitly chose the listed sequence as authoritative and asked us
to test possible rules. We therefore preserve the last blue rather than
silently change it to fit the simplest Fibonacci rule.

## Inheritance and live board

Closest proved mechanism: the exact integer charge and one-bit overlap
correction in [the Zeckendorf charge note](duck_zeckendorf_20260925.md),
Z1--Z14. Closest geometric mechanism:
[THM-3339, Fibonacci three-ray Berggren transplant and moving-owner obstruction](../../01-canon/theorems/THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md),
especially its period-three primitive content and period-six Farey order.
The relevant three-branch result is
[THM-3357, Berggren three-branch Walsh level collapse and parent circuit](../../01-canon/theorems/THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md):
all three branches are genuine integer maps, not one cyclic action.

Canonical hostile: two equal visible units can have distinct lifted charges;
discarding that difference can invalidate addition. Corrected near miss:
the same integer has distinct Fibonacci decompositions only after relaxing
the ordinary nonconsecutive convention. Least-used sidecar: the charge of
the distinguished boundary unit in a decomposition of a Fibonacci number.

The live board is: exact supplied prefix; substitution supertiles;
nonconsecutive representations; two boundary-unit charges; period-six
Berggren phase; orbit-coupled potentials. The anchor is a faithful rule
test, the niche is the boundary-unit construction, and the wildcard is a
precise obstruction to turning bounded phase balance into Collatz descent.

## 1. What the supplied data says exactly

Use R=red, K=black, B=blue, with indices beginning at 1. The supplied word is

    RKBRKRKBRKBRKRKBRKRKBRKBRKRKBRKBRKB.

Its blue positions are

    3, 8, 11, 16, 21, 24, 29, 32, 35.

The numbers of nonblue entries between successive blues are

    4, 2, 4, 4, 2, 4, 2, 2.

Including the initial two entries before the first blue gives exactly the
user's displayed `2,4,2,4,4,2,4,2,2`. That first 2 is a boundary gap.
Every nonblue piece is a concatenation of RK pairs. Compress RK to A and
blue to B; the exact token word becomes

    ABAABABAABAABABAABABAB.                         (C1)

The Fibonacci substitution A -> AB, B -> A gives

    ABAABABAABAABABAABABAA... .                     (C2)

Thus (C2) matches all first 34 colours but predicts R at colour 35, where
the supplied data is B. The discrepancy is token 22, not an indexing or
red/black naming issue.

There is a stronger hostile check. Word (C1) contains length-five factors
`AABAA` at token 8 (one B) and `BABAB` at token 18 (three B's). Every binary
mechanical word, including every shifted Fibonacci/Sturmian word, is
balanced: two equal-length factors differ in B-count by at most one.
Indeed sums of a mechanical word's floor increments over length m are
either floor(m alpha) or ceil(m alpha), for any fixed slope and intercept.
Consequently **no slope, intercept, shift, or exchange of its two letters
makes the entire supplied token word a binary mechanical word**.

This excludes a substantial class of simple rules. It does not exclude an
extra state, a boundary modification, a nonstationary substitution, or a
different geometric indexing convention. A finite word alone cannot choose
among them. Extending (C1) by a periodic tail and extending it by a
nonperiodic tail already give two different continuations preserving every
provided term.

## 2. The strongest near-match has an exact formula

This section defines a **candidate infinite rule**, not the established
continuation of (C1). Put phi=(1+sqrt(5))/2. In the Fibonacci A/B fixed word,
B occurs at token positions floor(k phi^2), k>=1. Expanding A to RK and B
to B therefore puts the blue entries at

    b_k = 2 floor(k phi^2) - k
        = 2 floor(k phi) + k.                       (C3)

Its first positions are

    3,8,11,16,21,24,29,32,37,42,45,50,... .

For completeness, the complementary Beatty sequences floor(k phi) and
floor(k phi^2) partition the positive integers. The counts up to N sum to
N because 1/phi+1/phi^2=1 and the relevant arguments are irrational; the
increment at each N is one. Let A occupy the first sequence and B the
second. Under A->AB, B->A, the B produced by the j-th A is at position
floor(j phi)+j=floor(j phi^2). Hence this word is fixed by the substitution,
which proves (C3) without identifying a finite experiment as a proof.

A directly coloured substitution producing the same candidate is

    R -> RK,       K -> B,       B -> RK.            (C4)

The maps commute with A->RK, B->B. The independent Beatty and substitution
implementations agree for every one of the first 100,000 positions. The
first-symbol generations have lengths 1,2,3,5,8,..., so Fibonacci-length
supertiles are real objects here. They are words containing several
colours, not a single colour assigned to each Fibonacci integer.

Write B(N) for the number of blues in the first N candidate entries and
beta=phi^3=2phi+1. Formula (C3) gives

    b_k = beta k - 2 {k phi}.

Using b_k<=N<b_(k+1), including the k=0 boundary separately, yields

    -1 < B(N)-N/beta < 2/beta.                      (C5)

Thus the candidate blue frequency is phi^-3, and its red and black
frequencies are both phi^-2. The two nonblue colours differ in count by at
most one in every prefix. These are bounded discrepancy statements in
**spatial order**. They are not yet statements about the order in which a
Collatz orbit visits the integers.

One explicit triangle convention has row r, column c (1<=c<=r), cell label
r(r-1)/2+c, and diagonal index d=r-c+1. Assigning the listed colour C(d)
then creates a left-aligned triangular display of diagonal stripes. This
is a fully specified rendering convention; the user's prose alone does
not uniquely fix diagonal orientation or labels. It does not change the
word mismatch at d=35.

## 3. A genuine repair of nonconsecutive Fibonacci splitting

Set f_0=1, f_1=2, f_(k+2)=f_(k+1)+f_k. The maximum sum of a nonconsecutive
subset of f_0,...,f_(k-1) is

    f_(k-1)+f_(k-3)+... = f_k-1.                    (C6)

This follows by the recurrence, or by comparing the choices with and
without the largest available term. Ordinary Zeckendorf representations
therefore cannot split f_k into smaller distinct nonconsecutive positive
terms. An extra unit repairs exactly the one missing unit:

    f_k = 1* + f_(k-1)+f_(k-3)+... .                (C7)

Here 1* is a distinguished boundary atom, not another unmarked member of
the ordinary Fibonacci alphabet. The smaller Fibonacci terms in (C7)
are distinct and nonconsecutive. For example

    3=1*+2,    5=1*+3+1,    8=1*+5+2,
    13=1*+8+3+1,    21=1*+13+5+2.

The inherited charge makes the boundary atom more than a drawing device.
Let alpha=phi^-2 and

    q(n)=floor(alpha(n+1)),
    Q(n)=(n-2q(n),q(n)).                            (C8)

This is the charge of every distinct Fibonacci representation of n;
its ordinary integer value is recovered by (A,B)->A+2B. Its increments
are exactly the two vectors

    U=(1,0),       V=(-1,1),                        (C9)

both of ordinary value 1. The charge of the boundary atom in (C7) is

    Q(f_k)-Q(f_k-1) = U if k even, V if k odd.      (C10)

To prove this, the Fibonacci identity
alpha f_k-F_k=alpha(-1/phi)^k implies q(f_k)=F_k and
q(f_k-1)=F_k for k even, F_k-1 for k odd. Subtract (C8).
The smaller terms of (C7) have charge Q(f_k-1), so (C10) is exactly the
correction required, not merely a compatible colour assignment.

This is a literal construction of **two hidden unit types with the same
visible integer length**. Modulo two their charges are (1,0) and (1,1).
Two independently chosen colour names may encode these states, but that
choice does not derive the supplied R/K/B diagonal sequence. In particular
the inherited four-state charge cannot equal the supplied colour under
any relabeling: integers 2 and 4 have equal charge (0,1) modulo two but
different listed colours K and R.

Ordinary single-term atom charges cycle with period three. Combining that
atom colour with the boundary type gives six different states, with exact
period six in k. The map k mod6 -> the oriented Farey channel order of
THM-3339 at index k+2 is therefore explicit via its displayed six-cycle.
It preserves this phase and Cassini parity. It forgets Fibonacci size,
Berggren word length, and the affine owner/current. Those data are not
recovered by calling the six states a tournament.

## 4. Multiplication by three exposes the missing carry

The full Q(n) is an injective integer lift, so it does preserve arithmetic
computability. Its small colour quotient does not. Put

    u={alpha(n+1)},
    d(n)=q(3n+1)-3q(n)=floor(3u-alpha).

Then exactly

    d(n) in {-1,0,1,2},
    Q(3n+1)=3Q(n)+(1-2d(n),d(n)).                  (C11)

The four carries first occur at n=7,2,1,4 respectively. Formula (C11)
follows by writing alpha(3n+2)=3alpha(n+1)-alpha. Reducing modulo two
requires the carry parity in addition to Q(n) modulo two. The smallest
retained hostile found for adding modulus 30 is n=1 and n=31: both have
charge (1,0) modulo two and residue 1 modulo 30, but their 3n+1 images
have charges (0,1) and (0,0).

For the candidate stripe rule, even colour plus residue modulo 30 does
not determine the next colour under the shortcut Collatz map:
n=1 and n=61 both have state (R,1 mod30), but their images have colours
K and B. Modulo 60 the same witness still works. These latter claims are
about the defined candidate extension, not extrapolated facts about the
user's unchosen continuation.

More generally no finite state containing the candidate colour can be
closed deterministically under n->n+1: deterministic evolution on a finite
set is eventually periodic, whereas the colour has irrational frequency
phi^-3. A fixed modulus cannot turn this irrational rotation coding into
a complete finite-state arithmetic model. Finite digit transducers with
unbounded input words are a different possibility and are not excluded.

## 5. What this contributes to an orbit-coupled potential

The bounded discrepancy D(n)=B(n)-n/phi^3 is a rigorous spatial potential.
Along any actual orbit x_i it obeys the exact telescoping identity

    sum_(i<m) [D(x_(i+1))-D(x_i)] = D(x_m)-D(x_0).  (C12)

Thus its total contribution is bounded, however many crossings occur.
It cannot itself pay an unbounded number of crossings. It can still be
a bounded correction to a potential whose unbounded part records actual
arithmetic depth, excursion height, or a crossing debt.

There is a useful general obstruction. For any bounded function H on the
positive integers and any c>0, the potential

    V(n)=c log(n)+H(n)

cannot satisfy V(Tn)<=V(n) at every n>1 for shortcut Collatz T. Indeed,
n_k=2^k-1 has k consecutive odd steps and

    T^j(n_k)=3^j 2^(k-j)-1,    0<=j<=k.

The change in c log n over the block is
c log((3^k-1)/(2^k-1)), tending to infinity, whereas the change in H is
bounded by twice its supremum norm. Summing any proposed nonincrease
inequalities gives a contradiction for large k. This includes every
fixed finite-state correction and the discrepancy (C5).

The conclusion is about a particular potential class, not about the
possibility of Collatz Lyapunov functions. A successful crossing potential
may use an unbounded sidecar or operate only at selected return times.
The construction tells us which coordinate must do the work.

Three precise next questions survive:

1. Can a spatially balanced poset/rotation statistic be paired with an
   unbounded, endpoint-stable crossing debt whose change has a sign on
   actual arithmetic excursions?
2. Does a weighted extension law retain an inequality after adjoining
   the exact carry d(n), rather than averaging away the selection of the
   starting integer?
3. Can the two boundary-unit lift define useful signed cancellation on
   a fixed orbit, with an independently bounded unmatched boundary, even
   though the visible three-colour quotient alone is not closed?

These are OPEN targets with typed coordinates, not consequences of a
numerical resemblance between phi, 1/3, 2/3, and Collatz branch ratios.

## 6. Reproduction and controls

Run

    python -B 04-computation/experiments/crossroads_crossing_20260926_colour.py

The retained [output](crossroads_crossing_20260926_colour.out) checks every
provided term, independent substitution and Beatty constructions through
100,000 entries, all four exact triple carries through 100,000 integers,
candidate transition collisions for moduli 1,2,3,5,30,60,210, boundary-unit
identities through Fibonacci index 100, and all-odd growth identities
through length 1000. All decisions use integers, including exact square
roots for the relevant Beatty floors. Checks use explicit exceptions and
remain active under Python -O.

Positive controls are the 34-symbol match, commuting substitutions, exact
carry reconstruction, and boundary-unit decomposition. Hostile controls
are the last supplied colour, unbalanced length-five factors, same-charge
different-listed-colour integers 2/4, same-charge-plus-mod30 sources 1/31,
and the arbitrarily long all-odd block. No computation assumes a random
colour or a random Collatz source, and no finite range certifies an orbit
termination theorem.
