# Eighth-defect transfer closes one of the seven E12 price families

**Status: PROVED conditional on the inherited E11 certificate and its stated
positive slack.** This packet does not rerun that large certificate. It proves
a new uniform fee consequence by elementary integer arithmetic, closes the
three-non-one E12 price family, and reduces the remaining six families to 21
unordered eighth-defect partitions. E12 itself, all-layer content, and Collatz
convergence remain OPEN.

## Provenance and definitions

Read from repository HEAD `16a4cdb898c71d54dcb2c6cf58784c8e13692d22`:

- `05-knowledge/results/terminal_fees_phase_covers_20261030/COLLATZ.md`,
  sections 1 and 5: E11+, gamma, and the seven E12 families.
- `05-knowledge/results/return_kernel_core_transplants_20261030/COLLATZ.md`,
  sections 1, 2, and 5: native excursions, tangent boundary, and the exact
  zero-sum criterion. That predecessor already observes finite cost sublevels;
  finiteness alone is **not** claimed as a new result here.

For a positive valuation word w, A is the sum of its letters,
F_w(X)=(P_w X+B_w)/2^A, and c(w)=v_5(gcd(P_w-2^A,B_w)).
S_e is the integer sum of min(c(prefix),e) over nonempty prefixes.
Put lambda_e=(5/8)(1-5^(-e)). The inherited theorem is

    S11(w) <= lambda11 A - gamma       for every non-one word,
    gamma = 5878847/19531250.

A primitive11-excursion returns to the identity modulo 5^11 for the first
time. Its retained next tangent X lies in F_5^2, and its native price is

    delta11(w) = lambda12 A - S11(w).

The exceptional all-ones excursion has clock T=2*5^10=19531250,
price 1/5, and tangent v=(4,4). All native non-one clocks are at least 3.
No entry is optimized or replaced in the argument below.

## 1. General quantized slope-transfer lemma

Let A be a positive integer, S an integer, and assume

    S <= (5/8 - 1/D) A - gamma,    D>0.

For any rho>1, set

    h0 = (floor(8 gamma)+1)/8.

Then

    (5/8 - 1/(rho D)) A - S
       >= (1-1/rho) h0 + gamma/rho.                 (Q)

**Proof.** H=5A/8-S belongs to (1/8)Z. The hypothesis gives
H>=gamma+A/D>gamma, so H>=h0. Also A/D<=H-gamma. Consequently

    H - A/(rho D) >= (1-1/rho)H+gamma/rho
                       >= (1-1/rho)h0+gamma/rho.

This is a uniform argument for every integer clock, not a sampled-word
experiment. The strict definition of h0 matters when 8 gamma is integral.
For a general rational limiting slope p/q, the identical proof uses
H in (1/q)Z and h0=(floor(q gamma)+1)/q. Thus the mechanism is rational-score
quantization, not a numerical analogy between unrelated appearances of 1/8.

## 2. An E12 native family is now paid

Take D=8*5^10=78125000 and rho=5. Then
lambda11=5/8-1/D and lambda12=5/8-1/(5D). Since 1/4<gamma<3/8,
h0=3/8. Every non-one word, and in particular every native11-excursion,
satisfies

    lambda12 A-S11 >= eta
       := (3/2+gamma)/5
        = 17587861/48828125
        = 0.36019939328.                            (F)

This is stronger than simply discarding the positive clock term in
`gamma+(lambda12-lambda11)A`. It combines that growing bound with the discrete
eighth-defect lattice, rather than assuming an upper bound on A.

Thus for EVERY three non-one native excursions, without ANY tangent condition,

    sum delta11 >= 3 eta
                 = 52763583/48828125
                 = 1 + 3935458/48828125 > 1.

This closes the previously open `(a,b)=(0,3)` E12 family. The six remaining
pairs (a all-ones blocks, b non-one blocks) are

    (0,1), (1,1), (2,1), (3,1), (0,2), (1,2).

Their tangent condition remains sum X_i=-a(4,4). The new bound alone does
not pay them: eta+3/5<1 and 2 eta+1/5<1.

**Sharpness boundary.** Eta is the exact minimum allowed by the inherited
integer-score inequality alone. The arithmetic data

    A=5781487, S=3613429, H=3/8

satisfy inherited equality and give delta=eta. No actual valuation word with
these data is supplied or asserted. Thus eta is NOT called the sharp actual
native excursion fee. Improving it requires information beyond this scalar
integer inequality.

## 3. Exact defect grammar for the six remaining obstructions

For each non-one block write

    j_i=8H_i=5A_i-8S11(w_i),

an integer at least 3. Set J=sum j_i and B=sum A_i. The inherited inequality,
clock integrality, and the price identity give

    A_i <= D(j_i/8-gamma),
    A_i = 5j_i mod8,
    total family price = a/5+J/8-B/(5D).

If the family is unpaid, its price is less than 1, so necessarily

    D(5J/8-5+a) < B <= D(J/8-b gamma),              (W)
    J < 10-2a-2b gamma.                            (J)

Conversely, for *actual native words* of these defects, the strict lower
inequality in (W) is exactly the unpaid condition. The upper bounds are
necessary consequences of the inherited theorem, not existence proofs.
Sorting the j_i loses no information for obstruction search: the tangent
boundary is additive and the criterion checks every zero-sum sequence.

The complete finite partition list is:

| a | b | possible unordered non-one defects j_i |
|---|---|---|
|0|1|3,4,5,6,7,8,9|
|1|1|3,4,5,6,7|
|2|1|3,4,5|
|3|1|3|
|0|2|(3,3), (3,4), (3,5), (4,4)|
|1|2|(3,3)|
|0|3|none|

There are 21 remaining defect partitions, with exact clock ranges in
`eighth_defects.json`. Each letter is positive, so these clock bounds also
bound every letter and word length. The finite language can still be enormous;
no exhaustive search of it is claimed or recommended here.

A particularly narrow target is `(a,b)=(3,1)`: its unique non-one word must
have j=3, clock A<=5781487 with A=7 mod8, integer score S11=(5A-3)/8,
and tangent X=-3v=(3,3). The other narrow target `(a,b)=(1,2)` requires both
words to have j=3, each A_i<=5781487 and A_i=7 mod8, and X_1+X_2=(1,1).
Ruling out these *native* low-defect/tangent combinations would close further
families; the arithmetic conditions alone do not do so.

## 4. Aggregate clock bands (weaker diagnostic)

The separate `check_clock_bands.py` and `clock_bands.json` record a weaker
necessary filter, useful for checking the eighth arithmetic independently.
For an unpaid group of TOTAL clock C=aT+B, integrality implies

    floor(lambda11 C-b gamma)=floor(lambda12 C).

Indeed its integer S11 sum lies between the right floor and the left floor,
and the left real argument is strictly smaller than the right. Write
r=C mod8 and f_r=(5r mod8)/8. The condition is equivalent to the union

    0<C<=D(f_r-b gamma),  or
    5D f_r<C<=D(1+f_r-b gamma),

with C=r mod8. Its automatic coarse bound is
C<(1-b gamma)/(lambda12-lambda11). The checker partitions every residue
progression at every exact rational floor jump and verifies all 37 constant
cells, rather than sampling a selection of clocks. These bands alone leave
some apparent three-block clocks; section 2 excludes all of them using
individual block integrality. This is why the aggregate filter must not be
mistaken for the strongest result.

## Reproduction and limits

    python3 check_eighth_defects.py > /tmp/eighth.json
    python3 -O check_eighth_defects.py > /tmp/eighth-O.json
    cmp /tmp/eighth.json /tmp/eighth-O.json
    python3 check_clock_bands.py > /tmp/clocks.json
    python3 -O check_clock_bands.py > /tmp/clocks-O.json
    cmp /tmp/clocks.json /tmp/clocks-O.json

Both checkers use exact rational arithmetic and explicit exceptions, not
assert statements. They verify derived arithmetic and partition completeness;
they do not independently verify the inherited multi-billion-edge E11
certificate, nor enumerate native11 words. Independent mathematical audit of
this packet is requested separately. No E12 or convergence claim is made.

## Structural refinement

[TEN_RETURN_GRAMMAR.md](TEN_RETURN_GRAMMAR.md) further proves that every
remaining candidate contains exactly one non-one primitive10 core and at
most three all-ones10 blocks (at most two across any unpaid E12 family). Every j=3 or j=4 candidate is itself primitive10.
The normal form retains the full two-digit tangent needed for E12, with an
independent exact modular-composition checker.
