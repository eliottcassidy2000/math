# Critical growth bands contain exponentially enough Collatz words

**Status: PROVED, with FINITE-EXACT controls; independently audited by the
geometry and automata lanes of this session.** This is a finite-word theorem.
It does not assert that any single positive integer follows an infinite
non-descending word. Its use is the sharp finite-horizon editing-price theorem
in the companion flow note.

## 1. Inheritance and connection contract

The source mechanisms are the cyclic-word lower bound in
`collatz_procgen_20260922_choice_ladder.md`, the retained pathwise deadline
versus aggregate capacity distinction in
`continuing8_20260906_coin_rational_handoff.md`, and the exact affine carry in
`collatz_guards_20260921_discrepancy.md`. The new operation is to concatenate
cycle-rotated, almost-critical blocks, retaining a uniform prefix growth cap.
The old ballot/Newton ratio note concerns a different observable; its
quadratic identities are not imported here.

Source: binary words with a specified number of ones. Target: Collatz parity
cylinders whose every prefix slope lies in a controlled real band. The map
is cyclic rotation at a cumulative minimum, then concatenation. It preserves
length and odd count, and restores prefix control lost by an endpoint-only
count. It discards the ordinary integer source; source residue and affine
carry remain mandatory. The hostile is positive integer realization of every
finite word on both signs and at q=5.

Use the half-step map T(n)=n/2 for even n, (3n+1)/2 for odd n. For a binary
word w of length L let e_j be the number of ones before time j, and

    a_j(w)=3^e_j / 2^j,              a_0=1.

Put theta=log_3(2), h=H_2(theta), eta=1-h. Define

    B_L(K)={w: 1<=a_j(w)<=K for every 0<=j<=L}.

## 2. Cyclic-minimum construction

For b>=1 put k=ceil(theta b). Since 2 and 3 are multiplicatively independent,

    1 < 3^k/2^b < 3.

Among the binom(b,k) binary words with k ones, at least binom(b,k)/b have
every prefix slope >=1. Indeed, choose a minimum of the cumulative sums of
the increments bit*log(3)-log(2), excluding the final endpoint, and rotate
there. Prefix sums in the first portion are nonnegative; wrapped prefixes
are nonnegative because the total sum is positive. Each output has at most
b preimages (its rotations). This counts words rather than rotation orbits;
periodic words cause no problem and require no division by their own period.

Every such block has all prefix slopes <=3^b and final slope <3. Write
L=tb+r, 0<=r<b. Concatenate t independently chosen blocks and append r ones.
All resulting words are distinct, have every prefix slope >=1, and obey

    max_j a_j <= K_L(b):=3^(t+b+r).

Consequently the following is an explicit finite lower bound:

    |B_L(K_L(b))| >= (binom(b,ceil(theta b))/b)^t.       (1)

One may replace binom/b by the exact integer count of good blocks, as the
script does. No rounding of logarithms is used for acceptance decisions.

For b=floor(sqrt(L)), Stirling's bounds, k/b=theta+O(1/b), and smoothness of
H_2 near theta give

    log_2 |B_L(K_L(b))| >= h L-O(sqrt(L) log L),
    log K_L(b)=O(sqrt(L)).                             (2)

The entropy upper bound follows from e_L>=theta L:
|B_L(K)|<=sum_{e>=theta L}binom(L,e)<=2^(h L+O(log L)).
Thus the subexponential growth cap in (2) loses no exponential population.

An optional better deterministic block choice,
b of order sqrt(L log L), changes the total subexponential cost in the
editing-price application to O(sqrt(L log L)); (1) remains the exact finite
certificate. This refinement uses the same Stirling estimate, not a claimed
Brownian-meander approximation.

## 3. Why this repairs the moment failure

For each selected full word and observation time k<L, the source cylinder
has density 2^-L. Its forward affine slope a_k is the Jacobian weight in the
hub capacity. The diagonal contribution to the squared capacity is

    sum_(w,k) a_k(w) 2^-L <= L K |B_L(K)|/2^L.

Without the cap, squaring tilts the odd-step frequency to 3/4 and leaves a
nonvanishing diagonal. The cap removes that obstruction while (2) retains
the desired entropy exponent. This observation alone does not bound the
off-diagonal overlaps. The companion flow argument supplies the missing
bound using actual integer spacing and the affine carry interval.

## 4. Controls and scope

`crossroads_20260926_ballot.py` compares dynamic programming with exhaustive
word enumeration through length 12 for multipliers 3 and 5 and caps
None,2,4,8,16; independently rotates all near-critical endpoint words; tests
the block lower bound; and realizes 32 constructed words as positive
integer trajectories on q=3,5 and shifts +/-1. There are 84,918 enumeration
and rotation word-checks. All gates passed.

The q=5 construction is a hostile, not the same density theorem: log_5(2)<1/2,
so its critical-endpoint subpopulation is exponentially smaller than its
entire non-descending population. The q=3 upper-tail argument uses
theta>1/2 and must not be transferred to q=5.

Fixed narrow bands also lose exponential rate; e.g. cap 2 has no survivors
at the displayed q=3 lengths. K grows with L in (2). The construction gives
a different family of ordinary integers for every finite L, not one integer
with an infinite bounded or slowly widening itinerary.

Reproduce: `python3 04-computation/experiments/crossroads_20260926_ballot.py`.

## 5. Retaining the exact block gain improves finite certificates

The asymptotic proof used P=3^k/2^b<3. Retaining P instead gives

    K=max(P^(t-1)(3/2)^k, P^t(3/2)^r),               (3)

because a block with k ones has no internal multiplier above (3/2)^k.
This is a concrete use for rational approximations k/b to log_3(2): a small
endpoint error reduces accumulated growth. It supplies a certificate
parameter, not a cycle exclusion or an irrationality theorem.

For a completely integral capacity constant, use

    M=sum_(j<L)(j+1)(floor(j/3)+1),

which overcounts all possible odd counts. If K=3^x/2^y is the larger term
in (3), the following factored lower bound is exact:

    binom(b,k)^t 2^y / (b^t 2^L 3^x M).              (4)

The script chooses among an explicitly stated finite candidate bank using
approximate logarithms; it does not certify global optimality. Selection
does not affect the validity of (4), and comparison of the two possible
values of K uses integer powers. Descriptive decimal negative logarithms:

| L | chosen b | -log_2 of lower bound (4), approximately |
|---|---|---|
| 100 | 50 | 59.455331 |
| 1000 | 125 | 211.107883 |
| 10000 | 500 | 1020.878229 |
| 100000 | 1369 | 6821.161268 |
| 1000000 | 7812 | 55616.245172 |

These finite bounds are weak at small L even though the limiting exponent
is sharp. The output retains every exact binomial, power and capacity
factor, rather than treating the displayed decimal as a certificate.
