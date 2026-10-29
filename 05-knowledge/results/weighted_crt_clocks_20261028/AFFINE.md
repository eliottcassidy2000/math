# Critical ternary branching and signed Collatz: the source cannot be erased

**PROVED** A1-A3 below. **FINITE-EXACT** replay: [check_affine.py](check_affine.py),
[JSON](check_affine.json). The classical cycle equation and sign conjugacy
are credited to the prior
[signed-cycle packet](../arithmetic_braids2_20260917_signed_cycles.md), which
also cites [Lagarias's overview](https://arxiv.org/abs/2111.02635).
Collatz(3n+1) and LRC(14) remain OPEN.

## A1. Both signs have the same complete finite valuation law

For epsilon=+1 or-1, define the accelerated odd map

    U_epsilon(x)=(3x+epsilon)/2^v2(3x+epsilon).

On odd 2-adic integers exclude the pole and its iterated preimages whenever
a word is evaluated. They form a null set. For an ordered positive valuation
word a=(a_1,...,a_r), set A_j=sum_(i<=j)a_i, A_0=0, A=A_r, and

    C(a)=sum_(j=0..r-1)3^(r-1-j) 2^A_j.

Then on its guarded cylinder the exact composite is

    U_a,epsilon(x)=(3^r x+epsilon C(a))/2^A.             (A1)

The cylinder is the **one odd residue class modulo2^(A+1)** given by

    x=(2^A-epsilon C(a)) (3^r)^-1 mod2^(A+1).           (A2)

Proof: each inverse branch y->(2^a_i y-epsilon)/3 sends odd 2-adic y to
odd x and has exactly the prescribed valuation. Composing the branches
gives(A1) in reverse and the image of all odd y is(A2). Consequently its
normalized Haar probability among odd inputs is2^-A, independently of
epsilon. This proves equality of **all finite valuation-word laws** for
the two signs, not just equality of one-step statistics. The one-step law
is P(a=j)=2^-j for j>=1.

This probability is on odd 2-adic source states. A finite starting-height
cutoff on positive integers is a different measure and need not give it
exactly. The uniform automorphism law in(W4), with offspring0 or3 of mean1,
is also a different probability space. No map identifying these probability
spaces or their dynamics has been constructed here.

The exact signed conjugacy is U_(-epsilon)(-x)=-U_epsilon(x). It preserves
forward time and the entire ordered valuation word, but exchanges positive
and negative chambers. It is **not** reversal of the iteration arrows.

## A2. Cycles test both the carry and the positive chamber

A periodic word at source x must satisfy

    x=epsilon C(a)/(2^A-3^r),                          (A3)

as well as the exact valuation guards. The carry is positive, so positivity
requires the denominator to have the sign epsilon. It is not determined
by r and A alone: (1,2) has C=5, whereas(2,1) has C=7, both with r=2,A=3.

For epsilon=-1, positive cycles checked here include

    1 -> 1                                  [word(1)],
    5 -> 7 -> 5                             [word(1,2)],
    17 -> 25 -> 37 -> 55 -> 41 -> 61 -> 91 -> 17
                                             [word(1,1,1,2,1,1,4)].

The last displayed cycle makes37 an actual state in this lane; this is
not an identification with the discriminant37 in the Ramanujan lane.
For epsilon=+1 the same words have negative cycle sources. The checked
positive fixed point1 instead has word(2). No claim of a complete positive
cycle census for either sign is made.

The script checks every ordered positive word with r<=5 and A<=12 via(A2),
and independently enumerates every odd residue modulo2^(A+1) when A<=8.
It also checks every edge and ordered valuation in the four named cycles.

### The genuine Lucas bridge changes the clock

The recovered [THM-4592, current all-period addendum](../../../01-canon/theorems/THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md)
already supplies an exact bridge: periodic rational points of the **standard**
3x+1 map, with odd denominators, have cyclic no-11 parity words. At odd
standard period n, the golden reader is a bijection onto R_n. One odd
accelerated step of valuation a costs1+a standard steps. Thus the standard
clock for our r-step word is n=r+A, not r or A separately.

For the positive3n-1 cycle through17, negate sources to obtain the negative
standard3n+1 cycle through-17. Its r=7,A=11 gives standard period18. The
inherited source places it in

    R_18=O/(phi^18-1)=O/(76),  |R_18|=76^2=5776.

Indeed phi^18-1=76 phi^9, and phi is a unit. Its chart is O/4 times O/19;
19 splits, so this is a16-state factor times two19-state factors. The
natural quotient R_18 -> R_9 has76-element fibers, since the target has76
states. It forgets actual even-period source data. Consequently the repeated
76 here has a proved ideal-factorization source, while treating R_18 and
R_9 as the same carrier would be false.

For the5<->7 accelerated3n-1 cycle, r=2,A=3 gives standard period5 after
negation. This is precisely the negative five-cycle in the older eleven-state
dictionary. The same period-five rational dictionary contains the positive
1/13 cycle. These are direct connections among5,7,11,13,76, with the map,
sign, integrality domain and changed clock all retained. They concern
periodic rational states, not a proof that every positive integer reaches a
cycle. The all-period reader theorem is inherited; we do not claim it here
as a result of the new finite checks.

## A3. Why centralizer holonomy is not a Collatz affine group

An allowed r-step word has positive affine slope3^r/2^A. For r>0 and A>0
this is never1. The group Aff_+(Q) of affine maps x->a x+b with a>0 has
no nonidentity finite-order element: finite order forces a=1, then b=0.
It follows that every homomorphism from a finite group into Aff_+(Q) is
trivial. In particular the finite ternary centralizer C_k cannot embed
faithfully into this affine group (C_k is nontrivial from k=2).

Thus the inherited complete re-marking holonomy, including its constructive
length bound, does not directly become an identity group of Collatz words.
A word can fix one source x, as in(A3), without being the identity map.
For two shared sources its difference law is

    U_a,epsilon(x)-U_a,epsilon(y)=3^r 2^-A (x-y).        (A4)

No nonempty such word fixes two distinct sources. A one-point cycle and a
universal port-preserving correction are different types of object. Finite
modular quotients may have nontrivial torsion, but they forget the rational
slope/height information used in this obstruction.

The useful surviving comparison is a **source-marked path groupoid**, with
the affine map, the valuation guard, and source/target retained on each arrow.
For asynchronous joins, also retain the two clock entries. The prior
[clock-holonomy note](../collatz_clock_holonomy_20261007.md) proves that receipts
(x,y;a,b) and(y,z;c,d) compose as

    (x,z; a+max(c-b,0), d+max(b-c,0)).                  (A5)

The difference a-b adds, but does not retain the lower time barrier. A
nonzero closed charge proves a periodic tail only after its checked entry
time. Equations(A1),(A4),(A5) explain the shared lesson of the tournament,
CRT and Collatz lanes: a composable relation needs its seam and clock;
an unlabelled loop count does not supply them.
