# Modular slope gaps sharpen the conditional Collatz bank head

**Status, 2026-10-19: PROVED elementary refinement.** The uniform excluded
head for the concrete expanding bank {K,K^4F} and reset bank {F,H} improves
from 534 to 107. Four exponent-parity classes admit bounds 42, 77, 107, 24.
These are sufficient bounds, not assertions of optimality or new source
entry. The inherited arbitrary-bank theorem and its proof remain unchanged.

## 1. What is inherited

The [arbitrary-bank result](collatz_bank_reset_20261019.md) proves the
following conditional statement. An expanding bank has maps a_i x+b_i
with a_i>1 and b_i/(a_i-1)<=C. A contraction bank has maps r_j x+d_j
with 0<r_j<1 and d_j/(1-r_j)<=P. For an arbitrary nonempty expanding word
followed by a reset word whose total multiplier contracts, the full lawful
Terras parity cylinder has

    T^N(x)=(Sx+H)/M, S=3^h, M=2^N,
    Delta=M-S>0, 0<H/M<C+P,
    x=rho+tM, 0<rho<M, t>=0.

Its exact first descending lift is

    t_min=max(0, floor((H-Delta rho)/(Delta M))+1). (1)

The previous bound used only the positive integer inequality Delta>=1.
For independently guarded adjacent-return banks, the upper source is
rho+1+tM with the same lift index, and its descent condition is identical.
These admissibility and paired-index facts are inherited, not replaced by
modular inequalities.

## 2. General gap-aware refinement

**Lemma.** Suppose a family of these contracting cylinders has a valid
positive lower bound Delta>=d_min>0. Then

    t_min <= ceil((C+P)/d_min).                  (2)

**Proof.** Write L=ceil((C+P)/d_min). For t>=L,

    x=rho+tM >= (C+P)M/d_min,
    Delta x >= d_min x >= (C+P)M > H.

This is precisely strict descent. The strict carry bound handles equality
at the real-valued threshold; no subtraction of one from the ceiling is
justified in general. The adjacent companion has the same exact cutoff. ∎

A safe modular way to obtain d_min is this: if Delta is a positive integer
congruent to r modulo q, its least possible positive value in that residue
class is d=r for 1<=r<q, and d=q for r=0. One may take the minimum of these
positive representatives across any **proved** set of possible residues.
A residue of zero must never be treated as a positive lower bound zero.
This device does not assume that the minimum is actually attained by a
power difference; a possibly unattainable lower bound is still valid.

## 3. Modulo 8, then modulo 24

Assume N>=3 and h>=1. Modulo 8,

    2^N=0,
    3^h=3 if h is odd, and 1 if h is even.

Thus positive Delta is 5 modulo 8 when h is odd, and 7 modulo 8 when h is
even. In particular Delta>=5 uniformly. The N>=3 hypothesis matters:
2^2-3=1 is an explicit exception below it.

Retaining N's parity and reducing modulo 24 gives the sharper table:

- N even, h odd: Delta=13 modulo 24, so Delta>=13.
- N even, h even: Delta=7 modulo 24, so Delta>=7.
- N odd, h odd: Delta=5 modulo 24, so Delta>=5.
- N odd, h even: Delta=23 modulo 24, so Delta>=23.

Indeed for N>=3, 2^N is 16 modulo 24 for even N and 8 for odd N; for h>=1,
3^h is 3 modulo 24 for odd h and 9 for even h. Subtracting gives exactly
the listed residues. The condition h>=1 is needed for this table: at h=0,
3^h=1, so for example 2^4-1=15 does not occupy the even/even residue 7.
No Diophantine power-gap theorem or unproved separation estimate is used.

## 4. Concrete bank constants and improved heads

For the inherited expanding bank {K,V=K^4F} and reset bank {F,H},

    C=83829391903560839836419813590128166431
          /157489106280555890461927692448373459,
    P=1.

Every eligible composite word has N>=3 and h>=1. Equation (2) therefore
gives the new uniform bound

    t_min<=ceil((C+P)/5)=107.                    (3)

Using the exact exponent parity gives:

- N even, h odd: t_min<=ceil((C+P)/13)=42.
- N even, h even: t_min<=ceil((C+P)/7)=77.
- N odd, h odd: t_min<=ceil((C+P)/5)=107.
- N odd, h even: t_min<=ceil((C+P)/23)=24.

If the expanding word contains k K blocks and v V blocks, and the reset
word contains f F blocks and g H blocks, then

    N=30k+122v+2f+9g,
    h=19k+77v+f+4g.

Hence N has parity g and h has parity k+v+f. This gives a direct class
selection without computing the huge powers in Delta. Orders within the
two banks are arbitrary, though the expansion-then-reset factorization and
full-word source guards are still required.

The earlier optional one-extra-reset theorem already removes all head
exceptions on its enlarged lawful cylinder. The present refinement instead
improves the guaranteed head for the original chosen reset word. It does
not authorize that extra reset on sources outside the enlarged cylinder.

## 5. Scope and exact controls

This sharpens a bound within already constructed sufficient descent
cylinders. It does not prove that additional arbitrary integer sources enter
these cylinders, that reset capacity is available at arbitrary exits, or
that any residue density supplies a deterministic trajectory frequency.
The true optimal uniform head for this concrete bank remains unresolved.
All 396 sampled schedules have actual t_min=0, but this is a finite control,
not a universal zero-head theorem. The Collatz conjecture remains open.

The [read-only checker](../experiments/collatz_modular_bank_head_20261019.py)
uses only Python's standard library and exact integers/rationals. It checks
82,602 positive power differences with 3<=N<=512 and 1<=h<=N, all four
modulo-24 classes, the exact ceilings, 396 bank schedules, each selected
class bound, and the strict endpoint inequality at the certified lift.
There are **167,190 explicit checks**. All four classes occur in these
controls: 91 with minimum 13, 87 with minimum 7, 131 with minimum 5, and
87 with minimum 23. Normal and optimized receipts match byte for byte.

Run from the repository root:

    python3 05-knowledge/experiments/collatz_modular_bank_head_20261019.py > /tmp/collatz-modular-normal.json
    python3 -O 05-knowledge/experiments/collatz_modular_bank_head_20261019.py > /tmp/collatz-modular-optimized.json
    cmp /tmp/collatz-modular-normal.json /tmp/collatz-modular-optimized.json
    cmp /tmp/collatz-modular-normal.json 05-knowledge/results/collatz_modular_bank_head_20261019.json
    cmp /tmp/collatz-modular-optimized.json 05-knowledge/results/collatz_modular_bank_head_20261019.optimized.json

The all-word conclusions follow from the symbolic modular argument and
(2), not extrapolation from finite exponent or bank controls.
