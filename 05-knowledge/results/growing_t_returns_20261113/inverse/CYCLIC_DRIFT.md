# Binary periods have a sharp, positive Hamilton-count drift

**PROVED IN NOTE:** the cyclic correlation formula, its sharp lower bound,
its equality cases, and the fixed-period obstruction below. **FINITE-EXACT:**
all32766 binary periods through length14 and four actual unit-ready graphs.
The growing-core selector, the unit-ready405 ray, and Chan--Pak remain OPEN.

## 1. Inheritance and the exact inverse interface

The inherited function is

    g(0)=g(1)=1,
    g(2n+1)=g(n),    g(2n)=g(n)+n+(n mod2).

Its exact section inverse was already developed in
[recursive seed switching S1](../../recursive_seed_switching_20261103/selector/NOTE.md),
and specialized to actual ready-root constructors in
[the preceding inverse packet](../../changing_core_returns_20261113/returns/PROOF.md).
The latter retains the extra start multiplicity: for a marked core with
H=h, E1=e, S2+P2=G, E2=f, and start count b, the literal outer funnel and
newborn have full Hamilton count

    F(n)=h g(n)+G n+e(n mod2)+(G-e-f-b).                 (1)

The newborn's actual indegree-one defect is b. A scalar solution cannot
silently change b3 to b1. The previous k17 witness is a useful hostile:
its full count is right but its defect is3.

The periodic-word formula itself is also inherited from
[paired descent, section4](../../h13_paired_descent_20261111/construction/PROOF.md).
The new result is a cyclic expression for its drift, its exact sharp
minimum, and the resulting complete test for zero drift. The live board
is exact inverse / cyclic word / positive drift / actual root state /
changing-core selector. The latter must escape the fixed-period obstruction,
rather than merely lengthen a familiar periodic word.

## 2. The drift is an explicit cyclic correlation

Fix length L>=1 and a binary word of numeric value r, allowing leading
zeros. Put Q=2^L-1. Let b_i be its bits from low to high, with indices
read modulo L. Define

    Z_d = sum_i (1-b_i)b_(i+d),       1<=d<L.

Thus Z1 counts the cyclic zero-to-one transitions. For L1 set Z1=0
and interpret the sums below as empty. The inherited suffix law is

    g(2^L n+r)=g(n)+a n+beta(n mod2)+c,
    a=Q-rev_L(r), beta=1-b_(L-1), c=g(r)-1.

Its exact repeated-word drift is

    D=beta*b_0+c-a*r/Q.

**Theorem.** It also has the phase-independent expression

    D = Z1 - (1/Q) sum_(d=1)^(L-1) 2^(d-1) Z_d.        (2)

Consequently D is unchanged by cyclic rotation, reversal, or complement
of the period. These symmetries preserve the drift, not the absolute
address or the source's exact constant term.

**Proof.** Expand the recurrence down the binary digits, putting b_L=0
temporarily at the end of the finite word:

    g(r)-1 = sum_(i:b_i=0) (floor(r/2^(i+1))+b_(i+1)).

The term beta*b0 supplies the missing cyclic transition at the boundary.
Also a=sum_(i:b_i=0)2^(L-1-i). If R_(i+1) is the cyclic right rotation
of r by i+1 bits, ordinary division gives

    Q floor(r/2^(i+1))-r*2^(L-1-i)=-R_(i+1).

Therefore QD=QZ1-sum_(i:b_i=0)R_(i+1). Expand each rotation in binary
and collect by cyclic distance d to obtain(2). Rotation invariance is
immediate. Writing Z_d as the number of ones minus their cyclic
autocorrelation proves reversal and complement invariance. QED.

## 3. Strict positivity and the exact smallest mixed drift

For a mixed word, L>=2 and Z1>=1. For each d,

    Z_d <= min(d,L-d) Z1.                              (3)

Indeed an interval of d cyclic steps from a zero to a one crosses a
zero-to-one transition; summing over the possible transition positions
gives dZ1. Traversing the complementary interval gives (L-d)Z1, since
the number of one-to-zero transitions equals Z1.

Set u=floor(L/2), v=ceil(L/2). A finite geometric sum gives

    sum_(d=1)^(L-1) 2^(d-1)min(d,L-d)
       =(2^u-1)(2^v-1).

Combining with(2) proves the sharp inequality

    D >= Z1*(2^u+2^v-2)/(2^L-1)
      >= (2^u+2^v-2)/(2^L-1) >0.                     (4)

Equality in the final minimum occurs exactly for cyclic rotations of
a single zero block and single one block whose lengths are u and v,
including complements. Such words have Z1=1 and attain equality in
every distance bound(3). Conversely equality requires Z1=1; for one
block of each symbol, the distance profile is
Z_d=min(d,L-d,number_of_zeros,number_of_ones). To attain all bounds,
the two block lengths must be as balanced as possible.

Constant periods have D0 directly. Hence **D0 if and only if the
period is constant**. The positive minimum is of order2^(-L/2); the
ratio tends to2 along even lengths and3/sqrt(2) along odd lengths.
Thus long block periods can have very small drift; smallness does not
make that drift vanish.

## 4. What this proves about an exact target selector

Let n_j be the address obtained from a fixed initial prefix p by
repeating the period j times. The inherited exact law is

    n_j=(p+r/Q)2^(Lj)-r/Q,
    g(n_j)=C_g*2^(Lj)+D*j+E_g,       j>=1.

The address identity holds also at j0. The g identity is stated from
j1 because the initial prefix may have a different parity from b0;
that first suffix can contribute a separate boundary correction.
If the prefix parity equals b0, the same formula extends to j0.

For any actual fixed-core response(1),

    F(n_j)=C*2^(Lj)+hD*j+E,       j>=1.               (5)

The coefficient h=H(core) is positive. Thus a genuinely mixed repeated
period cannot realize a pure exponential target A*2^(Lj)+B for all j,
or even eventually: exponential, linear, and constant sequences are
linearly independent. Fixed extra suffixes preserve the coefficient
of g and hence do not remove the hD*j term. A different exponential
rate is excluded already by its leading growth.

This excludes fixed-core mixed-period selectors, not arbitrary words,
not changing cores, and not two different source languages whose drifts
cancel. The latter possibility is real: the inherited13/25 construction
has a positive drift match, with its separate seed-compatibility equation
still open. The constant-period boundary is handled separately in
[REINDEXED_CLOCKS.md](REINDEXED_CLOCKS.md).

The source-to-target map here is a binary period -> its cyclic distance
profile -> its exact linear coefficient in a Hamilton-count family.
The lost information is the period's phase and prefix. The adjusted
seed and growth coefficients in the inherited periodic-return criterion
are the necessary sidecars for an actual equality.

## 5. A positive physical mixed-period family

Take the actual unique-start directed triangle as core. Its response is

    F(n)=3g(n)+15n+(n mod2)+10.

Use prefix1 and period01, so n_j=4n_(j-1)+1 and n0=1. Equation(5)
becomes the exact all-j formula

    H_j=(64*4^j+6j+23)/3,       j>=0.                  (6)

Every member is an actual strong unit-ready tournament under the
inherited funnel/newborn construction. The first four positive-j
counts are95,353,1379,5477. They are independently counted from their
literal adjacency matrices, including the predecessor's unique start
and the newborn's defect1. The linear term is2j, exactly3*(2/3)j.
This is a classified infinite family beyond pure dyadic clocks, not
the missing405 family or a source-return theorem.

Reproduce with `python -B 05-knowledge/results/growing_t_returns_20261113/inverse/cyclic_drift.py`
and then `-O -B`. Every word through length14 is checked against the
inherited drift expression, all stated symmetries, the sharp bound,
and all equality cases. Normal and optimized receipts agree.
