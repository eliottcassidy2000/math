# Small recursive certificates at the insufficient-precision boundary

**PROVED scoped identities and guarded descent families / FINITE-EXACT
controls / Collatz OPEN, 2026-09-27.** No novelty claim. This continues
[the precision-reset investigation](reset_20260926_board.md). Throughout,
U(n)=oddpart(3n+1) on positive odd integers. Independent checks and exact
reproduction are recorded in [the audit](entry_20260927_audit.md).

## 1. The strongest positive result: unbounded time, three rule nodes

For each k>=0 choose the least integer t>=1 such that

    2^t*(8^(k+1)-5) > 9^(k+1)-5.                    (B1)

Let b be ANY positive odd integer satisfying

    b*9^(k+1)=5 mod2^t,
    n=b*8^(k+1)-5.                                  (B2)

Then n has **exact first-descent time 2k+2**. Its first2k steps follow
the repeated valuation word (1,2), the next step also rises, and

    U^(2k+2)(n)=oddpart(b*9^(k+1)-5)<n.              (B3)

The last inequality follows directly from (B1), because the endpoint is
at most (b*9^(k+1)-5)/2^t and b>=1. The finite prefix before it grows
above the retained original source. There is no height or k cutoff.

The proof has three parameterized nodes: repeat the guarded (1,2) block
k times, take one U step, take one U step. Parameters and integer arithmetic
are unbounded; node count is not bit complexity. Membership is decidable
from n itself: require s=v2(n+5) to be a positive multiple of3, recover
k=s/3-1 and b=(n+5)/2^s, then test (B1)--(B2). Thus a bounded atomic
lookahead obstruction does not obstruct a small recursive proof grammar.

For each k, (B2) is one complete dyadic residue class at arbitrary height.
The classes are disjoint because v2(n+5)=3k+3. The old171-row bank contains
the k=0,1,2 classes, but is disjoint from n=-5 mod4096, which contains
every new class with k>=3. These are actual additions to the inherited
descent coverage, not just different descriptions of its cylinders.

The first20 classes add exactly 2553380107527241/2^64 to the old bank's
natural density, giving union density
6990313556829423891/2^65, approximately0.18947282861673337. These densities
refer to all natural inputs satisfying local descent certificates, not
to orbit frequencies or to the entire convergent basin. The remaining
class tail after k=19 has density less than1/(8*9^20).

The [recursive note](entry_20260927_recursive.md) proves the grammar,
guards and coverage comparisons. Its exact boundary controls matter:
the wider one-extra-bit condition fails at k=5, source786427; the next
tier fails at k=11, source68719476731. They fail for the selected return
word, not for their complete Collatz orbits. Increasing the precision
threshold in (B1) repairs precisely that failure.

## 2. An actual three-type excursion controller

There is a literal three-way arithmetic decomposition, separate from any
proposed Zeckendorf colouring. For every n=3 mod8 put

    L=v2(n+5),  k=floor((L-1)/3),  r=L-3k in{1,2,3},
    b=(n+5)/2^L (positive odd).

Then

    n=2^r*b*8^k-5 -> U^(2k)(n)=2^r*b*9^k-5=:m.      (B4)

The arrow is k actual guarded (1,2) blocks. It transfers three binary
valuation units to two ternary valuation units each time. The actual
redecoded state updates (k,r,b) -> (k-1,r,9b), so k decreases while the
odd coefficient grows. The quantity 2*v2(n+5)+3*v3(n+5) is preserved
across the phase. Equivalently the j-th boundary is
2^r*b*9^j*8^(k-j)-5 using the original b and an explicit elapsed counter j.
At the phase boundary the three states have different exact behavior:

| Type r | Boundary behavior | Retained arithmetic needed |
|---|---|---|
|1|m=1 mod4, so U(m)<m|This need not pay the original n. The terminal exponent depends on b.|
|2|m=7 mod8; exactly ell=1+v2(b*9^k-1) further valuation-one steps lead to1 mod4|The new growth run can be arbitrarily long.|
|3|m=3 mod16; two further steps give oddpart(b*9^(k+1)-5)|The congruence in (B2) is sufficient to pay the original n.|

For type2, the argument of the valuation is positive: the excluded b=1,
k=0 case would have source-1. The run length follows from
v2(m+1)=2+v2(b*9^k-1). The guarded subset of type3 is exactly the successful
all-k family above; type3 without the guard need not pay the selected return.
Type1 contains the ordinary members of the closed family in section3.
The original 27-family 4*8^k-5 occupies type2, including27.

This is a small controller carrying full integer registers. It does not
replace those registers by three colours. Its operation, preserved
predicates, and failure boundaries are explicit; no identification with
the supplied red/black/blue diagonal word is asserted. A next proof attempt
can target these three return types rather than repeatedly rediscovering
their long digit-level prefixes. The controller still needs a rank across
successive excursions, not just within one of them.

## 3. A different family containing27 really is closed and convergent

The [family note](entry_20260927_families.md) constructs

    N(k,h)=2*(47*4^h+7)*8^k/3^(2k+1)-5,
    k>=1, h>=0,  47*4^h=-7 mod3^(2k+1).              (B5)

All these are actual integers with a direct membership decoder, and

    U^2 N(k,h)=N(k-1,h),   U N(0,h)=47.              (B6)

The counter k decreases at actual phase boundaries. This is the desired
kind of recursive certification of a small fact, on a proved domain.
The pair (1,0) gives27; its base is31, and31->47 is growth, so the base
case uses the explicitly checked47-to1 suffix rather than a false descent.
Every member except27 first descends at step2k+1, to47. Every member,
including27, first reaches1 at step2k+39. The special first descent for27
is at step37. Only the fixed suffix needs finite verification.

Admissible exponents h form one residue class modulo9^k. A ternary digit
lifting rule computes that class, retaining all precision. The density
of exponents permitting at least k blocks is exactly9^(-k); source values
have count O(log X) below X and therefore density zero. This is a rigorous
recursion/frequency relation, not an assumed orbit distribution.

This family intersects the older fixed-coefficient family 4*8^k-5 only
at27. The older family remains a hostile to short lookahead: its first
descent is later than2k+4+v2(k), with infinity allowed. No convergence
theorem for all of that older family has been inferred.

## 4. What a universal entry theorem would have to certify

The [entry note](entry_20260927_entry.md) proves two distinct reductions.
Every positive odd orbit already reaches1 or a locally descending value
1 mod4. Also every positive odd orbit reaches1 or a value3 mod8. The
second reduces global convergence to convergence on the fixed low-precision
interface (q,R)=(1,1). The first does not prove convergence: its landing
point can exceed the original source by an arbitrarily large factor.

We now have a precise recursive baseline. Freeze the old bank, use its
certified decreasing jumps, and allow at most d unmatched single U steps.
The interpreter always halts because each move strictly lowers

    W(n,d)=n*2^d.                                    (B7)

For an unmatched move, U(n)<2n pays for decrementing d. The interpreter
may halt with rejection; (B7) proves total checking, not universal success.
Its accepted sets C_d are nested, decidable, and their union is exactly
the basin of1. Every fixed C_d omits some *proved-convergent* sources,
as [the independent entry review](entry_20260927_entry_review.md) proves
by explicit inverse completion. This does not obstruct a parameterized
macro such as (B3) or (B6), which compresses arbitrarily many atomic steps.

On the frozen finite universe of49,999 odd inputs below100,000, this old
policy accepts15,045 with no exceptions; all49,999 with at most36 exceptions.
The all-height macro theorems, not these finite counts, justify the new
certificates. An inequality that *earns* sufficient fuel on every uncovered
source, or a recursive grammar whose success is proved for every source,
remains missing.

## 5. The Gilbreath wall gives a real local mechanism, with a sharp limit

[THM-4511, the size-four wall theorem](../../01-canon/theorems/THM-4511-gilbreath-size-four-wall-theorem.md)
uses a closed alphabet, one-way dependence, and an irreversible column
phase. The [carry-wall note](entry_20260927_wall.md) transfers an actual
coalescence mechanism: multiplication-by3 carry maps synchronize all
incoming carries exactly when the scanned bits contain00 or11. Once
synchronized, dependence stays erased for the rest of that scan.

For an odd Collatz input, the first repeated bit pair ends precisely at
a=v2(3n+1), and the exact unread-tail code is

    n=(epsilon_a*2^a-1)/3+2^(a+1)*h,
    U(n)=6h+epsilon_a,
    epsilon_a=1 for even a, otherwise5.              (B8)

The unread h is unrestricted and carries the next actual integer. It
cannot be dropped on the strength of synchronization. Moreover3=(11)_2
maps to5=(101)_2, so the internal-wall flag can disappear at the next
iteration. The spatial irreversible phase is not an orbit-time rank.
Canonical alternating words (4^j-1)/3 are exactly the direct preimages
of1; finite zero padding resolves their carry scan as well.

The signed-Hadamard doubling and Mersenne zero triangles remain valid
separate structures. HYP-9162's all-level automorphism assertion remains
a hypothesis. The transfer used here is loss of boundary dependence;
no graph-size coincidence is substituted for a source-preserving map.

## 6. Incoming connection: the general residue sieve and one missing bit

A concurrent session added
[THM-4512, coefficient-descent cylinders](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md).
The [integration audit](entry_20260927_incoming.md) independently reproduced
its finite counts and repaired its dictionary. An exact valuation word
of total division budget A needs modulus2^(A+1); the coarser modulus2^A
permits extra final divisions and provides an upper bound. The source
comparison (2^A-3^j)*n>S remains the correct sufficient certificate.

This supplies a useful comparison within the repository: the general
coefficient sieve already covers about99.73 per cent of odd integers
by41 odd steps in density, whereas the selected old bank covers about37.87
per cent. Our added cylinders improve that compact selected bank, not the
general sieve. Their useful feature is their explicit unbounded-depth
grammar and whole-class root comparison. Source4091 demonstrates the
distinction: outside the old bank, it first descends at step8 to1639.
A bank failure is not a no-descent theorem.

The concurrent statement equating moving local-region entry with Collatz
was also corrected: local entry is already unconditional; root repayment
or convergence-certified entry is the remaining target. The one-member
threshold bound is independently checked through j=5000. Extending it
to all j needs specified effective constants and a finite cutoff bridge;
an existence citation alone cannot fill that interval. These repairs
leave the useful affine threshold and the exact density computation intact.

## 7. Next proof obligations, in order

1. Enlarge the grammar beyond the type3 guard (B2), using exact type1/type2
   returns while preserving the original root threshold. The new variable
   precision requirement grows with k; failing it is an explicit branch,
   not a proof of danger or divergence.
2. Seek a well-founded rank across completed excursions. The within-phase
   k rank and the optional-fuel inequality are already rigorous; neither
   currently controls regenerated precision at the next excursion.
3. Prove coverage of a convergence-certified region, or root-paying returns
   for every unresolved original source. Positive density, local wall
   synchronization, and universal entry into local descent regions do not
   individually establish this obligation.

The inheritance pass recovers the exact reset ledger and inverse completion;
its hostile is the fixed-coefficient27 shadow; its corrected near miss is
local entry versus paid root; its least-used sidecars are ternary exponent
digits and carry-bank images. The board stayed on six concepts: exact source,
phase rank, root debt, dyadic guard, ternary index, and boundary dependence.
The productive connection is operational: compile an unbounded excursion
into a finite guarded rule, then audit its arithmetic guard and its actual
consequence separately. The outstanding step is still global coverage.
