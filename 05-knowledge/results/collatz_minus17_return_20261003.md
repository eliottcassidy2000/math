# A -17 macro adds disjoint all-height return cylinders

2026-10-03 (America/Denver).

**Status: PROVED elementary guarded-word and first-descent statements;
FINITE-EXACT inherited-bank comparison and controls; OPEN universal
coverage and Collatz.** No novelty claim for the negative cycle, affine
composition, or the proof method. The contribution here is an explicit
additional family and its measured coverage beyond the named old bank.

**Independent audit, 2026-10-03:** the family-algebra agent reviewed the
guard equivalence, prefix inequalities, terminal exponent and payment,
disjointness, and strict density-tail bound; no correction was found.
The auditor also independently parsed the retained table without importing
research code: 171 rows, 65 pairwise-disjoint listed cylinders, the stated
exact density, no intersection with 4079 mod4096, and the unique q=333
match for 2031 mod4096. This independently checks the saved-bank separation;
the present script additionally regenerates all rows from the earlier code.

Throughout, U(n)=oddpart(3n+1) on positive odd integers, and a first descent
compares every iterate with the retained original source, not its parent.

## Inheritance and scope

The closest mechanism is the -5 all-height return family in
[entry_20260927_recursive](entry_20260927_recursive.md), sections 2-4.
The fixed comparison bank is the 171 rows / 65 disjoint cylinders in
[reset_20260926_swaplift](reset_20260926_swaplift.md), each certifying a
selected source-relative descent within 41 odd steps. Neither bank is
being identified with the whole convergent basin.

The -17 cycle is inherited: [THM-4484, free and sporadic
cycles](../../01-canon/theorems/THM-4484-free-and-sporadic-cycles.md) records
its shape (11,7) and gap 139; [five-mirrors audit](collatz_five_mirrors_20260929_audit.md),
item 44, checks the odd word and carry 2363. Targeted searches for Repeat17,
2187/2048, powers of 2187/2048, and the shifted n+17 coordinate recovered
these predecessors, but no earlier instance of the return-cylinder
augmentation proved below. That search is not a priority claim.

Canonical hostile: an insufficient terminal division can leave the entire
selected word above its source (section 4). Corrected near misses are the
lost last-oddness bit and the source/bank-residual conflation recorded in
[MISTAKES](../../01-canon/MISTAKES.md), 2026-09-27 incoming THM-4512 audit.
The adjacent 2026-09-26 two-copies correction also matters: the owner's
earlier "three parameterized rules" meant three nodes of the -5 certificate,
not the three known negative cycles. This note does not reinterpret that
phrase or identify negative cycles with Zeckendorf colors.

The board here is shifted source, legal repeated word, terminal valuation,
original-source payment, disjoint cylinder, and phase ledger. The anchor
is actual additional descent coverage; the niche is the -17 shifted ledger;
the wildcard is distinguishing this binary guard from a three-color state.

## 1. The guarded repeat macro

Put P=2187=3^7, Q=2048=2^11, and

    W=(1,1,1,2,1,1,4).

The negative odd cycle and its cumulative division exponents are

    c_i: -17, -25, -37, -55, -41, -61, -91, -17,
    A_i:   0,   1,   2,   3,   5,   6,   7,  11.

Direct substitution gives Q_W(n)=(Pn+2363)/Q, hence
Q_W(n)+17=(P/Q)(n+17). The affine formula alone does not establish legality.

**Guarded repeat lemma.** For positive odd n and k>=1, the first 7k exact
valuations equal W repeated k times iff

    v2(n+17)>=11k+1.                                  (1)

On that domain the actual endpoint is

    U^(7k)(n)=(P^k*n+17*(P^k-Q^k))/Q^k.               (2)

To see the guard directly, let n=aQ^k-17 with a positive even integer.
At block l and position i, the actual candidate is

    a*Q^(k-l)*P^l * 3^i/2^A_i + c_i.                 (3)

The high summand has more binary precision than the next cycle division,
including at the last step because a is even. Thus induction proves every
valuation exactly. Conversely, the exact word implies that its final
affine endpoint is odd. Equation (2), with P odd, then forces n+17 to be
divisible by 2^(11k+1). All steps are integral and positive.

Every positive-time iterate in this guarded repeat is above n. A block
boundary increases n+17 by P/Q>1. Within a block, for i=1,...,6, the
difference from the original source is at least

    Q*(3^i/2^A_i-1)+17+c_i
      = 1016,2540,4826,3112,5684,9542, respectively.   (4)

This conservative bound already uses the smaller coefficient 1; positive
even a only strengthens it. Thus the repeat itself is an expanding phase,
not a descent certificate.

## 2. An all-height first-descent family

For each m>=1 let t_m be the least integer t>=1 satisfying

    2^t*(Q^m-17)>P^m-17.                             (5)

Choose the unique positive odd beta_m<2^t with

    beta_m*P^m=17 mod2^t.                            (6)

**Return theorem.** Every positive b=beta_m mod2^t gives

    n=bQ^m-17,
    first_descent(n)=7m,
    U^(7m)(n)=oddpart(bP^m-17)<n.                    (7)

These are complete dyadic cylinders

    n=beta_m*2^(11m)-17 mod2^(11m+t_m).              (8)

The least representative is positive and exceeds 1; no small-source
exception is needed. The chosen budget is sufficient uniformly in b,
without a claim of optimality after exploiting the actual beta_m.

**Proof.** The actual word before its last letter is W^(m-1) followed by
the first six letters of W. Formula (3), with a=b and k=m, applies there:
the last partial block uses cumulative precision at most 7<11. The block
boundaries and (4) show every preceding positive-time iterate exceeds n.
At the last position, direct substitution gives

    3*U^(7m-1)(n)+1=16*(bP^m-17).                   (9)

Thus the last valuation is exactly 4+v2(bP^m-17), and the endpoint in
(7) is the actual odd part. Congruence (6) gives v2(bP^m-17)>=t_m.
Let D=2^t*Q^m-P^m. Inequality (5) says D>17*(2^t-1), so

    2^t*n-(bP^m-17)=bD-17*(2^t-1)>0.               (10)

This proves the strict final descent and, with the preceding comparisons,
its exact first time. In particular the endpoint is not falsely treated
as the unhalved affine quotient. A certificate can use one Repeat17(m-1)
node plus a seven-step terminal block with guard (6) and payment (10).
The counters and their bit cost are unbounded even though this description
uses a bounded number of parameterized rule nodes.

Membership is source-decodable: compute s=v2(n+17); require s=11m with
m>=1, recover b=(n+17)/2^s, compute t_m, and check (6). No unbounded search
along the orbit is needed to test this family.

## 3. Exact coverage added to the inherited bank

The cylinders for distinct m are disjoint because v2(n+17)=11m exactly.
The m=1 cylinder is n=2031 mod4096 and is already an inherited bank row
(q=333). For every m>=2, (8) lies inside

    n=-17 mod4096, i.e. n=4079 mod4096.               (11)

**FINITE-EXACT, with an all-height consequence:** none of the 171 original
rows intersects (11). The script regenerates every row, compares it with
the retained table, and tests each modular intersection exactly; the same
test is repeated on the 65 disjoint cylinders. This finite certificate
therefore proves that every m>=2 cylinder is entirely new to that bank.
All the inherited -5 return cylinders have residue 3 mod8, while (8) has
residue 7 mod8, so the entire additional infinite family is disjoint from
that prior recursive extension as well. The ordinary one-step descending
region 1 mod4 is also disjoint from these source cylinders.

The first added cylinder and its first-descent certificate are

    n=4194287 mod8388608, m=2, t=1,
    U^14(4194287)=597869<4194287.

The retained output includes all 15 values; every intermediate exceeds
the source. The first budgets are t=1 for m=1,...,10 and t=2 for m=11,...,21.

**PROVED density.** The added natural density among all positive integers is

    delta_add=sum_(m>=2) 2^(-11m-t_m)
             = approximately 1.192675256472887e-7.   (12)

For existence, the tail m>M is contained in n=-17 mod2^(11(M+1)); its
upper density tends to zero. Finite disjoint unions thus approximate the
whole union. Moreover (P^m-17)/(Q^m-17)>(P/Q)^m, so (5) gives
2^(-11m-t_m)<P^(-m). The tail after m=M is strictly less than

    1/((P-1)*P^M).                                  (13)

The output stores the exact fraction from m=2 through 20 and its rigorous
tail bound, as well as an independently normalized union with the old
bank. These are source-set densities, not frequencies along one orbit.
The newly covered sources have certified descents; this alone does not
assert their entire future or every positive integer is covered.

## 4. Failure boundary and retained phase information

At m=11 the required budget first increases to t=2. Keeping only t=1
admits b=1, for which v2(P^11-17)=1. Its source and selected endpoint are

    2658455991569831745807614120560689135,
    2737200544710109691038577966784875873,

respectively. Every iterate through the selected step 77 is above the
source. The first failed implication is terminal payment (10), not word
legality. The stronger survivor is the guarded expanding phase; the repair
is precisely the extra binary budget. Later descent is not excluded.

During each complete W block, the shifted quantity n+17 loses 11 binary
valuation units and gains 7 ternary units. Consequently

    7*v2(n+17)+11*v3(n+17),

and its factor coprime to 6 are conserved across the phase. This is a
phase ledger, not a global Lyapunov function: a later word can recreate
precision. It is distinct from the old -5 ledger 2*v2(n+5)+3*v3(n+5).

The full shifted valuation now splits into a counter and one of 11
remainders, with seven positions inside the macro; no three-color state
has been shown to replace those coordinates. The Zeckendorf colors in
[duck_zeckendorf](duck_zeckendorf_20260925.md) and
[reset colors](reset_20260926_colours.md) have their own precise charge
and marker meanings. A simple exact Fourier observation is available:
exp(2*pi*i*n/8) has opposite phases on 3 mod8 and 7 mod8, whereas reducing
modulo 4 merges them. That distinguishes the old/new guards; it does not
prove payment (10) or identify their colors with Fibonacci atom colors.

## Verification and next obligation

Run from the repository root:

    python -X utf8 -B 04-computation/experiments/collatz_minus17_return_20261003.py
    python -X utf8 -B -O 04-computation/experiments/collatz_minus17_return_20261003.py

The [retained output](collatz_minus17_return_20261003.out) records:

- 150 generated legal repeats: k=1,...,30 and even coefficients 2,4,6,10,18;
- 12,288 iff-guard controls: all odd sources 1,...,8191 and k=1,2,3,
  including two accepted cases and 12,286 rejected cases;
- 200 first-descent cases: m=1,...,50 with coefficient lifts 0,1,7,29;
- all 171 regenerated old rows equal their saved table, all parent-cylinder
  intersections, the 65-cylinder normalization, and the exact union density;
- the weakened-budget hostile, and the first added source's full trace.

Two arithmetic paths are compared: literal exact-valuation replay and a
generic letter-by-letter affine composition, each checked against the
closed fixed-point formula. Bank regeneration is also compared against
the previously saved table; this is code/data cross-checking, not a claim
of independent mathematical methods for every component.

The next coverage question is whether source-preserving compositions of
the -5 and -17 rules admit useful additional guarded domains. A rule that
only reaches another locally descending source must still pay the original
threshold. Universal successful certification remains OPEN.
