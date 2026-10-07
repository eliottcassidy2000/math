# Explicit absorption certificates for a bounded translation bank

**Status: PROVED** by the elementary construction below; **FINITE-EXACT**
controls have the universe in section 6. This proves positive conditional
Haar probability from every state in a finite bank. It does not prove
absorption for each individual driving integer. In fact, that universal
equal-time claim is false.

## 1. Inheritance and the retained coordinates

The inherited transition mechanism is
[THM-4569 — the Terras clock](../../01-canon/theorems/THM-4569-the-terras-clock-recurrence-of-two-collatz-orbits-is-unconditional.md).
The construction makes quantitative the finite-box accessibility step in
[THM-4581 — Haar coalescence](../../01-canon/theorems/THM-4581-haar-coalescence-affinely-related-collatz-orbits-merge-almost-surely.md).
This package does not depend on that theorem's global return drift or
absorption conclusion. The direct parity/ordinal interface was checked in
[the inverse-clock bridge](collatz_terras_inverse_clock_bridge_20261007.md).

Our concept board is: exact affine state, driver parity cylinder, shrinking
integer translation, certificate cost, and the equal-time/ROOT distinction.
The canonical hostile is the out-of-phase pair (2,1); the useful sidecar is
the driver's complete parity prefix, not merely its length.

Let T(y)=y/2 for even y and T(y)=(3y+1)/2 for odd y. At a current state
u=3^k y+e, retain the rational e, integer k, and next driver bit beta=y mod2.
Its transition is the inherited four-case table. At k=0, e is an integer.
The target of this note is the identity state (k,e)=(0,0).

| Source | Map | Preserved predicate | Lost information / required sidecar |
|---|---|---|---|
| Integer translation e at k=0 | Deterministic bit word b(e) | The word forces absorption | It does not assert that an arbitrary supplied y has those bits |
| Bit word of length L | Unique residue r mod 2^L | Exactly those sources realize the word | Actual magnitude requires a lift r+t2^L |
| Bounded bank, \|e\|<=M | Uniform deadline and probability floor | Each bank state has an absorbing cylinder | Distinct states may use distinct source cylinders |

Under full conditional Haar driver law, every specified L-bit prefix has
probability 2^-L. This premise holds at a finite parity stopping time; it
must not be replaced by an arbitrary already-fixed source condition.

## 2. A terminating integer algorithm

For positive e, repeat the following operations until e=0.

* If e is even, choose beta=0. The new state is (0,e/2).
* If e is odd, put a=v2(3e+1) and choose the word 0^a 1. The first zero
  departs to (1,(3e+1)/2), the next a-1 zeros halve the translation at
  level 1, and the final one returns to level zero with

  e'=(3e+1-2^a)/2^(a+1)=(oddpart(3e+1)-1)/2.

Here e' is a nonnegative integer and e'<e. More precisely,
e'<=(3e-1)/4<3e/4. Consequently the algorithm terminates, without any
assumption about the Collatz orbit of e. It only uses a strictly decreasing
auxiliary integer.

For e=-d<0, apply the positive algorithm to d with the two orbits swapped.
If its current state is u=3^k y+c and chosen driver bit is beta, the bit of
u is beta xor (c mod2). Use those bits as the driver word for the original
negative-translation pair. This preserves the word length. The mirrored
affine state is (-k,-3^-k c), so the last state is again (0,0).

There is no earlier identity state in either construction. There is also
no earlier literal merge on any compatible source: after a literal merge
all later parities agree, so k remains fixed. Its final value zero would
force the earlier relation to have already been (0,0).

## 3. Logarithmic cost and a uniform probability floor

Write L(e) for the constructed word length, and put

    kappa = 2/log2(4/3).

For a nonterminal odd round e -> e'>0,

    a+1 < log2(e/e')+log2(3).

Even halvings cost exactly log2(e/(e/2))=1. If there are N nonterminal odd
rounds before the last odd value e_terminal, their contraction gives
N<=log_(4/3)(e). The terminal round satisfies
3e_terminal+1=2^a, hence costs a+1<=log2(e_terminal)+3.
The logarithmic drops telescope, giving the all-integer bound

    L(e) <= log2(|e|)+N log2(3)+3
         <= kappa log2(|e|)+3                  (e != 0).

The identity kappa=1+log2(3)/log2(4/3) explains the constant. Also
3^5=243<256=4^4 implies kappa<5. Thus an exact, slightly looser deadline is

    L(e) <= 5 ceil(log2 M)+3                  (1<=|e|<=M).

For every integer M>=1, set D_M=floor(kappa log2 M+3). For each state
(0,e), 1<=|e|<=M, the conditional probability of absorption by D_M is at
least

    delta_M = (1/8) M^(-kappa).

Indeed its particular absorbing prefix has mass 2^-L(e)>=delta_M. Shorter
words may be padded after absorption to the common deadline. This is a
statewise uniform lower bound, not a claim that the same cylinder works
for every state. The boundary e=1 attains L=3: word 001, residue 4 mod8,
and cylinder mass 1/8.

The checker verifies the sharper real-log bound without floating point:
for L>=3 it is equivalent to 4^(L-3)<=e^2 3^(L-3). The exported bank
deadline uses only the bit length of M-1.

## 4. Actual sources and the exact hostile

For a length-L word, the parity decoder returns its unique residue r mod
2^L. Every sufficiently large positive y=r+t2^L has y+e>0 and realizes
the certificate by direct integer T-iteration. Thus the result supplies
arbitrarily large positive integer examples, with a checked source guard.

For e=1 and y=4, the pair (5,4) follows

    (5,4) -> (8,2) -> (4,1) -> (2,2).

For the incompatible y=1, the pair has the exact closed cycle

    (2,1) -> (1,2) -> (2,1).

These two distinct states prove that this pair never merges at equal T
time. This is an infinite obstruction proved by closure, not a finite
timeout. Both individual orbits nevertheless visit ROOT. Consequently a
universal integer theorem must retain phase or allow asynchronous common
future; equal-time coalescence cannot be its universal target.

A separate clock sidecar remains necessary even when a merge occurs. For
p=2417 and q=805, p=3q+2 and q=5 mod16, the affine identity first appears
at T-time 9, at the common even value 128. The first common odd value is
1 at time 16. The exact odd words are (2,6,8) and (4,1,1,10), each of
total cost 16. Thus the first common odd clock is tau+v2(common state),
not necessarily the identity-chain absorption time tau itself. The checker
replays both clocks and both literal odd words independently.

## 5. Scope of the quantitative improvement

The construction gives an explicit replacement for qualitative finite-box
accessibility in a separate recurrence argument. Recurrence of the
translation process into some bounded bank is not proved by this package.
Nor does its logarithmic certificate cost say that a supplied driver
integer satisfies that certificate: its source residue must be checked.
The result is independent of any global positive-integer ROOT theorem.

## 6. Reproduction and controls

Run the matching script normally and with Python -B -O:

    python -B 04-computation/experiments/collatz_finite_bank_certificate_20261007.py
    python -B -O 04-computation/experiments/collatz_finite_bank_certificate_20261007.py

The universe is all 1,024 translations e in {-512,...,-1,1,...,512}, with
three direct positive-integer source lifts for each certificate. Independent
literal parity inversion checks all 1,023 words of lengths 0 through 9.
Controls include the exact sharp length inequality, signed residue shift,
first merge, e=0's empty identity receipt, the closed nonmerge cycle, and
ten malformed-type or incompatible-phase inputs. The longest constructed
word in this finite universe has length 23; the general integer deadline
for M=512 is 48. These finite maxima are not asserted beyond that universe.

Normal, optimized, and saved outputs agree on 31,873 explicit checks.
An independent proof/code audit passed, including both clock and phase
hostiles. The LF-normalized output SHA256 is
`738bed16febf025eb7e4c49959962cfdcc277a47d6fbc8b9d4a340882ee12ebb`.
