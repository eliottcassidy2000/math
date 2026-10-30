# Independent audit: eighth-integral defect and one E12 family

## Verdict and dependencies

**PASS, conditional on the pinned E11+ theorem and its native return-kernel
criterion.** This audit independently checks the deduction and exact constants;
it does not repeat the inherited 2,509,000,000-edge certificate computation.
It closes the `(b,a)=(3,0)` family only. E12, the all-layer content bound, and
universal Collatz remain open.

Source inspected directly using `git show` at commit
`cd6d98ba58a2418fb81054659bc37962676c8dd0`:

* `05-knowledge/results/terminal_fees_phase_covers_20261030/COLLATZ.md`,
  opening E11+ statement and section 5;
* `05-knowledge/results/return_kernel_core_transplants_20261030/COLLATZ.md`,
  section 5 exact native return-kernel criterion.

The baseline worktree does not contain the former packet; the explicit commit
pin is therefore material. No producer code is imported by this audit.

## Definitions and units

For each finite positive valuation word, A is the sum of its original positive
integer letters, and S11 is the sum of integer truncated 5-adic contents over
nonempty prefixes. Set

    D = 8*5^10 = 78125000,
    gamma = 5878847/19531250,
    lambda11 = 5/8 - 1/D,
    lambda12 = 5/8 - 1/(5D).

The D used here is NOT the certificate serialization denominator 390625000
called D elsewhere in the source. They differ by a factor of five. The native
11-excursion price is called `delta11=lambda12*A-S11` in the source. Calling
it delta12 would be an alternative indexing convention and should be avoided.
The inherited E11+ applies to every finite word containing a non-one letter,
not just a terminal subclass. It states

    S11 <= lambda11*A - gamma.

Hence it applies individually to every non-one native 11-excursion with the
original clock, without altering its tangent or substituting an artificial
word. The empty tangent classes retain price infinity.

## Universal proof

Put H=5A/8-S11 and j=8H=5A-8S11. Then j is an integer and

    H >= gamma + A/D > gamma > 1/4.

Thus j>=3, including all equality/boundary cases. Rearranging this same
inequality, rather than weakening the units, gives A<=D(H-gamma). Therefore

    delta11 = H-A/(5D) >= (4H+gamma)/5
            >= (3/2+gamma)/5
            = 17587861/48828125 =: kappa.

The final inequality is inclusive. It is not necessary to assert that the
bound is attained by any native word. Consequently every triple of non-one
native 11-excursions, regardless of tangents, has total price at least

    3*kappa = 1 + 3935458/48828125 > 1.

In particular this applies to every triple whose tangent sum is zero, exactly
the source's `(b,a)=(3,0)` requirement. Taking minima within tangent classes
preserves the uniform bound. The proof retains the source's chronological
native blocks, clocks, and next-layer return payment; it does not merge
blocks, change source maps, or charge a shared terminal return repeatedly.
The inherited criterion pays one unit per actual next-layer zero-sum group.
This closes only one family and makes no assertion that every group is now
paid: the other six remain.

## Exact scalar envelope and its limitation

The congruence j=5A-8S11 implies A=5j (mod 8). For fixed integer j>=3, the
largest scalar-admissible clock is

    U_j = 5j + 8 floor((D*(j/8-gamma)-5j)/8).

Whenever that scalar class is nonempty, its minimum fee is

    F_j = j/8-U_j/(5D).

At j=3, U_3=5781487 and S11=(5U_3-3)/8=3613429. This pair saturates the
inherited scalar inequality, so F_3=kappa exactly: imposing clock integrality
and the mod-8 congruence gives no further uniform improvement in this case.
For all j>=4, the earlier universal inequality gives
F_j>=(2+gamma)/5>kappa. Hence this is the sharp scalar envelope over integer
pairs, not a claim of sharpness over actual native excursions. It also obeys
U_(j+8)=U_j+D and F_(j+8)=F_j+4/5.

The scalar pair is NOT a Collatz word witness. In particular it may violate
native return constraints or tangent conditions. Excluding it and its nearby
low-fee classes is a genuine remaining arithmetic task. The six other source
families have b*kappa+a/5<1, so uniform-price algebra alone does not close any
of them.

## Independent checks

`check_audit.py` uses only standard-library exact fractions and explicit
exceptions, with no producer imports and no `assert` statements. It checks
all stated constants, the scalar saturation pair, strict three-fee excess,
and all six residual uniform-price gaps. It additionally checks the envelope
and its exact period on 10,000 consecutive j values as regression tests;
the universal claim is proved above, not inferred from this finite test.
Both ordinary Python and `python -O` pass with byte-identical output in
`audit_receipt.json`.

## Structural addendum: low-defect excursions have exactly one heavy block

**PASS under the same inherited source assumptions.** Let a non-one primitive
11-excursion have j=5A-8S11 in {3,...,9}. Cut it at every actual return modulo
5^10. Its final endpoint is a return modulo 5^11, hence also modulo 5^10,
so there is no unfinished suffix in this decomposition. Each block is a
primitive native 10-excursion. The entering affine map is identity modulo
5^10; native and actual content counts through level 10 therefore agree.
Original clocks and S10 add over these consecutive blocks. Since the entire
word is primitive at level 11, S11-S10=1, exclusively at its final prefix.
Consequently, with Delta=lambda12-lambda11,

    delta11 = sum_blocks(delta10) - 1 + Delta*A.

Call a block heavy if it contains a non-one letter. There is at least one
heavy block because the original word contains a non-one. The inherited
terminal certificate gives every heavy block delta10>=1+gamma; each
all-ones primitive 10-block has exactly delta10=1/5. If there were two or
more heavy blocks, nonnegative prices of the other blocks would give

    delta11 >= 1+2*gamma+Delta*A > 1+2*gamma > 9/8.

But delta11=j/8-A/(5D)<j/8<=9/8, a contradiction. Thus there is exactly
one heavy block. Writing a10 for the number of all-ones 10-blocks yields

    gamma+a10/5+Delta*A <= delta11 < j/8.

Since Delta*A>0 and a10 is a nonnegative integer, the necessary limits are

    j:          3  4  5  6  7  8  9
    max a10:    0  0  1  2  2  3  4.

The bounds are inclusive integer maxima derived from a strict real
inequality. They are necessary structural restrictions, not existence
claims for words attaining them. In particular, j=3 or j=4 permits only
one total 10-block, so the original primitive 11-excursion is itself also
primitive at level 10. No alteration of the native tangent or of the
source-certificate assumption is used. No further E12 family is closed by
this addendum alone.

## Final placement and two-digit grammar audit

**PASS.** With u the unique all-ones primitive 10-block, the exactly-one-heavy
block result gives w=u^l c u^(a-l), 0<=l<=a<=4. If c has low tangent
Y=-a*v, v=(4,4), every placement is primitive11: preceding block boundaries
have tangent k*v for 0<k<=l, and boundaries after c have tangent (k-a)*v
for l<=k<=a. They are nonzero except at k=a, since 0<=a<=4. For a=0 the
only block boundary is the terminal one. Internal prefixes of each block
cannot return at level11 because they do not return at level10, even after
composition with the entering map, which is identity modulo5^10.

At q=5^10, retain actual affine coefficients
F_c=(1+q*alpha_c,q*beta_c), F_u=(1+q*alpha_u,q*beta_u)
modulo25q. Since q^2 vanishes modulo25q, composition adds the pairs
C,U in (Z/25Z)^2. Thus C+aU and the terminal native11 tangent
(C+aU)/5 modulo5 are independent of placement. The latter quotient is well
defined exactly when C+aU=0 modulo5, as required. Direct modular powering
of the all-ones affine map gives U=(4,4) modulo25, not merely modulo5.
This computation is about affine maps with invertible powers-of-two
denominators reduced modulo5^12; it does not discard a boundary digit.

The level10 native content sums and original clocks are unchanged under
these placements, and every valid placement has precisely one added
level11 return. Hence S11 and the price are unchanged. The exact all-ones
clock T=3906250 gives 5T/8-S10(u)=1/4, so

    j11(w)=j10(c)+2a-8,
    j10(c)=j11(w)+8-2a.

`check_placements.py` independently checks all 9375 combinations of a
retained pair C, a in {0,...,4}, and l in {0,...,a} by sequential modular
affine composition. All 375 combinations satisfying the endpoint constraint
have no earlier level11 return at block boundaries. It also verifies the
quotient tangent and the exact defect transfer. Internal-prefix exclusion
uses the native primitivity proof above, not the finite pair enumeration.
Normal and optimized Python outputs agree in `placements_receipt.json`.
The converse is a return-level grammar conditional on an actual native core;
it is not an existence theorem for arbitrary residue pairs, and satisfying
the grammar alone does not certify a low-defect E12 obstruction.

## Final frozen refinement: quantize the core as well

**PASS. This supersedes the earlier weaker per-j count caps**, which remain
valid necessary bounds but are not the strongest conclusions of this audit.
For the unique heavy primitive10 core c, the inherited certificate says

    H10(c)=5A_c/8-S10(c) >= 1+gamma+A_c/D > 1+gamma > 10/8.

Because 8H10 is an integer, j10(c)>=11. Combining this with the exact
j10(c)=j11+8-2a10 gives the improved caps

    j11:        3  4  5  6  7  8  9
    max a10:    0  0  1  1  2  2  3.

Thus the low-defect candidate grammar has a10<=3; the broader placement
lemma already checked remains valid through a10=4.

Transferring the same heavy-core inequality to lambda12, with its correct
fee offset 1+gamma rather than gamma, gives

    lambda12*A_c-S10(c)
      = H10(c)-A_c/(5D)
      >= (4H10(c)+1+gamma)/5
      >= (4*(11/8)+1+gamma)/5 = 1+kappa.

Each all-ones primitive10 block costs exactly 6/25 at lambda12. Subtracting
the unique terminal level11 unit therefore yields

    delta11(w) >= kappa + (6/25)*a10.

An unpaid E12 candidate with a12 all-ones11 blocks, b11 non-one11 blocks,
and total internal all-ones10 count K consequently satisfies the STRICT
necessary inequality

    a12/5 + b11*kappa + (6/25)*K < 1.

For (a12,b11)=(0,1),(1,1),(2,1),(3,1),(0,2),(1,2), its integer K maxima
are respectively 2,1,0,0,1,0. In particular every **unpaid candidate** in
(2,1) must have K=0 and be primitive10; this must not be asserted of every
paid member of that family. No whole remaining family is proved empty.
`check_final_tiers.py` verifies these exact rational endpoints and the next
excluded integer tiers, with matching ordinary/optimized receipts in
`final_tiers_receipt.json`. All inherited assumptions are unchanged. This
is the final bounded audit scope; no deeper-layer search was performed.
