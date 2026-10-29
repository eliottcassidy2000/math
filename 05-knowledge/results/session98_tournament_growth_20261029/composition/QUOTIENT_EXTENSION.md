# Exact observable quotients of finite-surplus composition jets

**PROVED:** the two observer quotients and their refinement maps below.
**FINITE-EXACT:** complete rational ranks and kernel bases through surplus 12.
**INHERITED FINITE-EXACT:** existence/counts of the order-12 prime action types
used for the actual-family hostile. No new prime representatives were generated.
This extends [RESULT.md](RESULT.md), whose sufficient-jet theorem is unchanged.

## 1. Two different meanings of 'all child palettes'

Fix surplus h. The ambient jet space V_h has monomial basis z^alpha with
sum_(d odd) d alpha_d<=h. Define

    E_b(x)=sum_(i>=1) b_i x^i,
    Phi_h(F) = F(E_b(x),E_b(x^3),...) modulo x^(h+1).           (1)

**Formal observer:** treat b_1,...,b_h as independent indeterminates.
This models arbitrary numbers of tagged child colors at each size. Its exact
linear state is V_h/ker(Phi_h). The observation matrix has columns z^alpha
and rows (r,b^beta), the coefficient of b^beta x^r in (1).

**Ordinary tournament-palette observer:** children are pairwise distinct
ordinary tournament isotypes, with exactly the unique singleton as baseline.
There are a_(i+1) available colors of size i+1, so

    b_i in {0,1,...,a_(i+1)}.                                (2)

Every vector in this finite Cartesian grid is realizable: independently choose
that many isotypes of each size. Which isotypes were chosen does not alter
inflation orbit counts, because canonical root uniqueness sees them as
pairwise distinct colors. Their internal automorphisms are irrelevant here.

Let R_h be reduction of coefficients in the polynomial ring modulo

    I_h = < product_(j=0)^(a_(i+1)) (b_i-j) : 1<=i<=h >.        (3)

Then the exact ordinary-palette linear state is

    V_h / ker(R_h Phi_h).                                    (4)

**Proof of exactness and linear minimality.** A polynomial vanishes on the
whole grid (2) iff its unique remainder with degree in b_i at most a_(i+1)
is zero. Indeed, divide successively by the monic one-variable polynomials
in (3). The remainder vanishing on the grid is zero by repeated univariate
interpolation. Thus equal states in (4) are precisely equality of every
ordinary-palette inflation series through surplus h. The formal version is
ordinary polynomial coefficient equality. Any sufficient linear map must
have kernel contained in the respective observation kernel, so cannot have
smaller rank; the observation map attains that rank. QED.

Scope: these are minimal linear quotients of the **ambient jet space**.
They need not be minimal on the span actually realizable by a fixed arity's
prime action indices. They do not describe pointed attachments or metagraph
edges. The sufficient original jets remain correct, merely nonminimal.

## 2. First formal redundancy and earlier ordinary-palette redundancy

At h=6,

    Phi_6(z3^2-z1 z5)=0,

because both terms first occur at x6 and have coefficient b1^2. This is the
first ambient formal kernel: rank equals dimension through h=5.
At h=3, however,

    (R_3 Phi_3)(z3-z1^3)= (b1-b1^3)x3 = 0,

since there is only one ordinary 2-vertex tournament, b1 in {0,1}.
Thus formal independent b_i are not synonymous with all actual ordinary
palettes. This is the first ambient ordinary-palette kernel.

## 3. Exact ranks

| h | Jet dimension | Formal rank | Ordinary-palette rank |
|---:|---:|---:|---:|
| 0 | 1 | 1 | 1 |
| 1 | 2 | 2 | 2 |
| 2 | 3 | 3 | 3 |
| 3 | 5 | 5 | 4 |
| 4 | 7 | 7 | 6 |
| 5 | 10 | 10 | 8 |
| 6 | 14 | 13 | 10 |
| 7 | 19 | 18 | 13 |
| 8 | 25 | 23 | 16 |
| 9 | 33 | 30 | 19 |
| 10 | 43 | 37 | 23 |
| 11 | 55 | 47 | 28 |
| 12 | 70 | 57 | 33 |

These are rational matrix ranks, not numerical estimates. Full kernel bases
are in `quotient_receipt.json`. Through h=12 only bounds b1<=1, b2<=2,
b3<=4 affect reduction: b4 has bound a5=12 while its exponent is at most3,
and later variables likewise never reach their bounds.

## 4. A reachable prime-family collision, not just formal redundancy

Use four actually occurring arity-12 prime action types:

    R = p1^12,
    A = (p1^12+2 p1^3 p3^3)/3,
    B = (p1^12+2 p3^4)/3,
    C = (p1^12+4 p1^2 p5^2)/5.

The inherited audited type census supplies respectively more than two rigid
classes, 14572 of type A, 2933 of type B, and 72 of type C. Select eleven
distinct roots for each family:

    L = 8 roots of type A + 3 roots of type B,
    K = 2 rigid roots + 9 roots of type C.

There is no tagging or repeated root needed. Their aggregate 3-jets satisfy

    J3(L)-J3(K) = 24(z3-z1^3),                                (5)

so are distinct. Yet every actual ordinary palette gives exactly the same
numbers of L-rooted and K-rooted inflations at sizes 12,13,14,15. The first
three fixed-point moments (M0,M1,M2) are (11,60,258) in both families;
M3 differs by -24 and the single-3-cycle moment by +24.
They are distinguishable formally by taking b1=2 (tagged size-two colors),
which is impossible for an ordinary tournament palette.

They also separate at the next actual horizon: for the ordinary palette
{TT1,TT2}, the surplus-four coefficient difference is -18. Hence this is a
true temporary observer equivalence that refines under growth. It is not
claimed to be the smallest-order or smallest-family reachable witness.

Realizability provenance: `prime_types12.inherited.txt` is copied unchanged
from `04-computation/experiments/prime_cycle_index_tower_20261027/out/t12.txt`
at source pin `5862c94555eafea4e9a0638a17f7c661f8e45d37`.
THM-4816's independent audit is the inherited basis for trusting that census.
The new algebra verifies its use, not the existence census anew.

## 5. Growth and exact refinement

Dropping weight h+1 jet coordinates and dropping output terms x^(h+1)
commute with Phi. They also commute with ordinary-grid reduction; the new
variable b_(h+1) cannot occur below degree h+1. Consequently there are
well-defined surjective linear maps

    V_(h+1)/ker(Phi_(h+1)) -> V_h/ker(Phi_h),

and analogously for the ordinary observer. Thus increasing the growth horizon
refines the classification, while lower-horizon classes are recovered by
projection. There is generally no inverse map: a discarded kernel direction
can become visible later, as (5) does. The quotient stores precisely the
observable output polynomial, not a chosen arbitrary lift back to jet space.

## 6. Reproduction

    python3 observable_quotient.py > quotient_receipt.json
    python3 -O observable_quotient.py > quotient_receipt.optimized.json
    cmp quotient_receipt.json quotient_receipt.optimized.json
    python3 search_reachable.py > reachable_receipt.json

The implementation constructs coefficient rows by exact expansion, performs
rational elimination, and multiplies every reported kernel vector back into
the original matrix to check zero. Separate modular elimination with reversed
pivot order at moduli 1000003 and 1000033 matches every reported rank,
providing independent lower-bound checks. All horizon projection identities through
12 are checked. Direct evaluations on the relevant b1/b2/b3 grid through
h=6 independently check polynomial reduction. Normal/optimized receipts
agree. The normalized integer family relation is checked again against the
ordinary-grid observation matrix, and its next-horizon difference -18 is
verified directly for TT1/TT2. The inherited core supplies only odd partitions; no symbolic algebra
package or external reference ranks are used.
