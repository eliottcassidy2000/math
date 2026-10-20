# Independent audit of guarded reset-debt composition

Date: 2026-10-20. Independent audit of the companion guarded reset-debt package.

## Verdict

**PASS on the mathematical identities, guarded-cylinder theorems, concrete
examples, proofs, and finite receipts.** Two short publication clarifications
have been applied: explicitly name the hypotheses controlling the sign in (2.6),
and explicitly identify the ambient space of the repeated-A intersection.
Neither changes the main composition or guarded-bank theorem.

This is an independent computational and mathematical audit, not external
peer review. The audit does not prove Collatz convergence, arbitrary-source
bank coverage, a uniform native signed-row residual bound, or a universal
mean-current cancellation estimate. No literature-novelty claim is certified.

The producer's source files and the repository were not changed by this
audit. The published main document contains only the two approved scope
clarifications. Its exact predecessor is preserved in provenance, and
publication_overlay.diff records the complete change.

## 1. Reproduction and source identity

The original checker was rerun with ordinary and optimized Python. Both
outputs are byte-identical to both supplied receipts. All six entries in the
producer's `SHA256SUMS` verify, and all five inherited files have exactly the
hashes recorded in `SOURCE_INVENTORY.json`. The read repository HEAD is
`f3e658d7cc105d5e367e4dd565812a81e127a6fc`; its worktree was clean.

Independent receipt: **148,619 explicit checks**, with normal and optimized
outputs identical. Producer receipt: **321,690 explicit checks**, 120 valuation words,
41,370 odd residue representatives, 1,096 macro schedules, and maximum
literal shortcut clock 1,200. The normal/optimized check does not rely on
Python assertions.

Startup was read-only and used `python3 -B agents/start_session.py --topic
'Collatz reset debt composition'`. The coordination peek could not run
because this checkout has no `.machine-id`; registration was not attempted.
This is not a mathematical verification blocker. No fetch, Git mutation,
repository file change, or remote write was performed under this audit's
read-only repository scope.

## 2. Forward normalization, chronology, and slack

The exact forward map is `f_j(q)=(q+delta_j)/lambda_j`, with positive
lambda. Direct substitution gives

    f_2(f_1(q)) = (q+delta_1+lambda_1 delta_2)/(lambda_1 lambda_2).

Thus the chronological law and all-length debt sum are correct. Substitution
of inequalities requires positive lambda; algebraic identities alone do not
justify the same inequality with negative lambda. The independent checker
includes an explicit sign-reversal hostile.

Subtracting the debt sum from the telescoping product identity yields

    Lambda_m-1-Delta_m = sum_j Lambda_(j-1)(lambda_j-1-delta_j).

This proves the stated semigroup closure for every finite length, not merely
the tested lengths. It also proves the exact macroblock accounting rule:
positive later slack can pay earlier negative slack only with the displayed
chronological weights. No commutativity assumption is used, and debts do not
compose in reverse order accidentally.

The normalized debt barycentre is correct for all factors greater than one.
Its extension to lambda equal to one requires genuine identity factors,
delta zero, and at least one nonidentity factor so the denominator is not
zero. Negative weights arise if a factor falls below one, as the note warns.

Equation (2.6) is also correct. However, its sign is determined by
`delta_2 eta_1`, not by nonnegative debts alone. The exact hostile is

    (lambda_1,delta_1)=(1/2,0), (lambda_2,delta_2)=(2,1),
    r_12-r_1 r_2=-1/2.

The phrase "nonnegative debt cone" should therefore explicitly mean the
unit-budget cone `0<=delta<=lambda-1`. A precise replacement sentence is:

    If delta_2>=0 and eta_1>=0, in particular when both blocks lie in the
    nonnegative unit-budget cone, the difference is nonnegative. Nonnegative
    debts alone do not determine its sign.

## 3. Barriers, inverse order, and native norm debt

The identity `q_out-1=(q_in-1-eta)/lambda` proves invariance below one and
strict descent above one whenever lambda exceeds one. Strict eta is not
needed. The identity case provides no progress, exactly as stated.

For original resource units and barrier B, the correct certificate is
`d<=B(lambda-1)`, equivalently `e<=B(1-rho)` when `rho=1/lambda`.
For the inherited norm recurrence, which has `0<rho<1`, the exact conversion
is `delta=(lambda-1)R/B`. Its unit budget is equivalent to `R<=B`.
The weighted delta convolution divided by Lambda is the inherited
discounted residual convolution.

The original one-reset failure at scale k=23 remains a failure for the
actual signed rows, original coefficients, terminal labels and clocks. The
scalar integer-height bank does not repair it. The refinement paper carries
a uniformly priced Gram remainder under a common native map; it does not
charge a new independent refinement invoice at each reset, and it does not
pay the coarse mean-current estimate. The earlier split-product account has
a multiplicative imbalance factor and a different reset convention. The
packet correctly keeps all three carriers distinct.

The inverse real map is `G(y)=lambda y-d`; applying G1 then G2 gives
`d_12=lambda_2 d_1+d_2`. The inverse of the forward word uses the reverse
order, so the two displayed laws agree. Ray preservation follows from the
barrier inequality but does not give a legal integer preimage. For instance,
the A inverse at y=7 is 425/27, although 7 lies above the real barrier.

## 4. Complete valuation cylinders and ROOT guards

The affine carry formula and full modulus `2M`, not M, are correct. The
terminal condition says `S n+C` has exactly the required total power of two
before an odd terminal. Reducing this relation at each prefix gives the
exact odd intermediate values and consequently each prescribed valuation.
Positivity follows because the affine numerator and denominator are positive.

The unique source progression is `r+2M t`, t nonnegative. Its image is
`z_0+2S t`, with the lower bound retained; inverse integrality without the
correct image congruence and lower bound is insufficient.

Independent direct valuation enumeration covered every odd residue for all
340 words of lengths 1–4 with exponents 1–4, a total of 837,930 residues.
No imported function from the producer checker was used. For each word,
the checker also tested four source lifts and all odd target candidates in
a declared finite window, 56,414 targets altogether, against the inverse
image criterion.

ROOT stopping is separate from the affine word. The control `(4,2)` has
first lift `5 -> 1 -> 1`, which is not a valid stopped two-step trajectory;
the next lift `133 -> 25 -> 19` is valid. Thirteen of the 340 tested words
require deleting their least source. All tested positive later lifts pass
the proper-prefix ROOT guard. Generally, each prefix increases strictly
with t; as the least prefix value is a positive odd integer, only the t=0
lift can need deletion. The note's weaker "finite head" assertion is valid.

The sharp cylinder descent criterion follows immediately from
`F(n)-n=(C-(M-S)n)/M`: no positive source descends when M is at most S;
otherwise the maximum displacement occurs at the least allowed source.
Hence the stated iff is valid, including strictness at the equality boundary.

## 5. First-reset and grouped banks

For `(1^h,a)`, direct affine composition gives
`S=3^(h+1)`, `M=2^(h+a)`, `C=3^(h+1)-2^(h+1)`. This proves the stated
lambda, delta, and budget test. There is an admissible a for every h.
The budget forces slope at most
`1/[2-(2/3)^(h+1)] <= 3/4`. Consequently the arbitrary finite concatenation
theorem follows from closure plus the exact combined valuation cylinder.

The first-exit formula uses the shortcut clock correctly: h valuation-one
steps followed by j reset shortcut/halving steps, total h+j. The independent
checker audits h=0 through 150 at the least admissible final valuation.

The macro data are correct:

| Macro | Valuations | Affine map | Source guard | Image guard | Fixed point | First literal exit |
|---|---|---|---|---|---|---|
| A | (1,2,3) | (27x+23)/64 | 11 mod 128 | 5+54t | 23/37 | 5 |
| B | (1,3) | (9x+5)/16 | 19 mod 32 | 11+18t | 5/7 | 4 |
| E | (1,1,3) | (27x+19)/32 | 55 mod 64 | 47+54t | 19/5 | 5 |

The image progressions in the table use t nonnegative. A has unit slack
14/27, and B has unit slack 2/9. The first two odd steps of A grow every
source lift by `2+16t`, yet the full A macro descends. E fails the unit
budget but every actual source exceeds its fixed point; barrier four
certifies it. The common barrier-four slack for E after normalization is
1/108, positive. All three source minima exceed four, and all slopes are at
most 27/32.

The independent audit checks the uniform first exits by deriving every
literal value's affine dependence on t from two lawful lifts and checking
the signs of both coefficients. This proves no prior value is below the
source and the indicated value is below it for every t, not only finitely
many sampled lifts. The actual first-lift literal trajectories are recorded
in the independent receipt.

Every word over the three macros has a nonempty full exact valuation
cylinder. The common-barrier envelope composes by an elementary induction.
The independent finite semigroup controls cover every schedule of lengths
1–8, 9,840 schedules. These tests do not replace the all-length proof.

## 6. Finite depth is not one infinite integer orbit

The repeated-valuation-two formula is correct. Its complete integer
cylinders intersect in {1}, while their ROOT-stopped counterparts have
empty intersection. Every finite stopped depth still has arbitrarily large
positive sources.

For repeated A, `37 C_r=23(64^r-27^r)`. The exact odd-terminal cylinders
are equivalent to `37n=23 mod 2*64^r`. These nested residue classes have
the single 2-adic intersection point 23/37. In contrast, their intersection
inside the positive integers is empty: an integer in every class would
make the fixed integer `37n-23` divisible by arbitrarily large powers of
two and hence zero, which is impossible for integer n.

The rational 2-adic trajectory is explicitly

    23/37 -> 53/37 -> 49/37 -> 23/37,

with exact valuations 1,2,3. This confirms both the address and the extra
odd-terminal bit. The independent checker verifies repeated A and repeated
valuation two through depth 200.

The publication sentence should say "the intersection in the 2-adic
integers is {23/37}; the intersection among positive integer sources is
empty." Calling the positive-integer cylinders' intersection simply
"23/37" leaves the ambient space implicit and is avoidable.

There is no coverage theorem hidden in the bank. The actual source 11
follows A to 5, but 5 is in none of the A/B/E source cylinders. A uniform
contraction argument therefore rules out endless continuation within the
prescribed bank, not escape from the actual unresolved Collatz problem.
The note preserves this distinction throughout.

## 7. Independent validation universe and commands

The independent checker uses exact integers and fractions with explicit
exception checks. It does not import the producer implementation. It checks:

- 27,931 arbitrary affine parameter words, including signed debts and
  positive lambdas below, at, and above one;
- 1,728 three-stage native recurrences with varying rho and residuals;
- all 340 valuation words and 837,930 odd representatives described above;
- 56,414 candidate inverse image values and four lifts per valuation word;
- ROOT trimming and both directions of the least-source descent criterion;
- 151 first-reset words, with valuation-one run length 0–150;
- the all-lift first-exit coefficient signs of A, B and E;
- all 9,840 macro schedules of lengths 1–8;
- repeated-word identities through depth 200; and
- negative-lambda, insufficient debt-sign, inverse-ray, continued-ROOT and
  coverage hostiles.

Run from the audit directory:

    python3 -B independent_audit.py > /tmp/reset-independent-normal.json
    python3 -B -O independent_audit.py > /tmp/reset-independent-optimized.json
    cmp /tmp/reset-independent-normal.json /tmp/reset-independent-optimized.json
    cmp /tmp/reset-independent-normal.json independent.normal.json

The producer reruns were byte-identical to `receipt.json` and
`receipt.optimized.json`. The original proof supplies the infinite family
quantifiers; both computational paths supply finite-exact controls.

## 8. Publication boundary

With the two scope clarifications applied, the packet is suitable for publication
as elementary affine identities and sufficient guarded Collatz-cylinder
certificates, with the exact finite receipts attached. No new theorem ID,
global Collatz coverage claim, native 7X repair, or literature novelty claim
is warranted by this audit. Keep the inherited source inventory and the
finite-versus-all-length distinction in the published packet.
