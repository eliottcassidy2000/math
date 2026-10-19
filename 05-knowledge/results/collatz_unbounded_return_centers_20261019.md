# Expanding adjacent returns with unbounded compressed centers

**Status, 2026-10-19: PROVED.** There is a family of genuine oriented
adjacent Collatz returns whose lower sources have no descending prefix,
whose total multipliers lie strictly between 1 and 4/3, and whose expanding
center sizes are unbounded. This is a limitation on treating every such
compressed return as a member of one bounded-center bank. It does not
contradict the fixed-bank reset theorem, prove new source entry, or settle
the Collatz conjecture.

## 1. The actual guarded family

Use the Terras map and the independently guarded adjacent return K from
[the odd-adjacent packet](collatz_odd_adjacent_return_20261019.md):

    K(x)=(Px+C_0)/Q,
    P=3^19=1162261467, Q=2^30=1073741824,
    C_0=4912087879, c=C_0/(P-Q), lambda=P/Q>1.

Its lower and upper parity words are respectively

    110110110110110101011011101010,
    001111111111111101010110000100.

The lower cylinder is 815464443 modulo Q; the upper is one higher.
Every nonempty lower prefix has multiplier greater than one. Let F be the
adjacent contraction (3x+1)/4, with lower word 10, upper word 01 and fixed
point 1. Set mu=3/4 and

    theta=log(lambda)/log(4/3),
    j_L=floor(L theta),
    V_L=K^L F^(j_L), L>=1.                       (1)

Words are read chronologically, so V_L=F^(j_L) composed with K^L. The case
j_L=0 is allowed; every sufficiently large L has j_L>=1. These are
compressed macroreturns, not asserted primitive or first adjacent returns.
In particular the K block boundaries are already earlier adjacent returns.

## 2. Irrationality and the exact multiplier interval

The number theta is irrational. If theta=p/q for positive integers p,q,
then lambda^q=(4/3)^p would give

    3^(19q+p)=2^(30q+2p),

contradicting unique prime factorization. Thus L theta is never an integer.
Writing f_L={L theta}, one has 0<f_L<1 and

    A_L=lambda^L mu^(j_L)=(4/3)^f_L,
    1<A_L<4/3.                                  (2)

An exact integer certificate for any proposed pair (L,j) is therefore

    M<S and 3S<4M,
    S=3^(19L+j), M=2^(30L+2j).                  (3)

These inequalities imply j=floor(L theta), without approximating logarithms.

## 3. Actual joint cylinders and prefix survival

Put

    C_L=C_0(P^L-Q^L)/(P-Q),
    H=3^j C_L+Q^L(4^j-3^j),
    H_upper=3^j(C_L+Q^L-P^L)+2Q^L(4^j-3^j).

These are the integer carries obtained from the two actual concatenated
parity words. Their independently composed upper carry satisfies

    H_upper=H+M-S.                              (4)

The lower full word has its unique lawful cylinder

    x=rho+tM, rho=-H S^(-1) modulo M, t>=0.

The usual parity-word prefix congruences establish lawfulness of every
step. The upper full word has residue rho+1 modulo M. Neither residue is
zero: the lower word and upper word both contain odd steps. Moreover rho
cannot be M-1, since the upper residue would then be zero, which would force
an all-even word. Thus 0<rho<M-1 and the actual upper positive sources are
rho+1+tM. On this joint cylinder the endpoints are V_L(x) and V_L(x)+1.
Every cylinder has arbitrarily large positive sources. This is an explicit
existence statement, not a claim that an arbitrary source follows the word.

Every nonempty lower prefix of V_L has multiplier greater than one. The
prefixes during K^L do, by the inherited K prefix inequality. During the
F segment, after k complete blocks the multiplier is lambda^L mu^k,
which is at least the final A_L>1. The intermediate odd step of each F
block multiplies the preceding boundary multiplier by 3/2. All parity-word
carries are nonnegative. Hence every lower prefix strictly exceeds its
positive original source, and V_L is itself expanding. No such prefix
claim is made for the upper source, which descends immediately.

## 4. Unbounded centers

The normalized affine formula is

    V_L(x)=A_L x+e_L,
    e_L=c A_L+1-(c+1)mu^(j_L).                  (5)

Indeed K^L(x)=lambda^L x+c(lambda^L-1), and
F^j(y)=mu^j y+1-mu^j. The expanding center size is

    c(V_L)=e_L/(A_L-1).                         (6)

There is a sequence L_n tending to infinity for which f_(L_n) tends to
zero from above. Here is an elementary version of the irrational-rotation
fact, so no quantitative equidistribution assertion is needed. By putting
the N+1 fractional parts of 0,theta,...,N theta into N equal subintervals,
some q>0 has distance from q theta to an integer at most 1/N. Irrationality
makes that distance delta positive. If {q theta}=delta, this is already a
positive fractional part at most 1/N. Otherwise {q theta}=1-delta. Put
k=floor(1/delta). Since delta is irrational,

    0<{kq theta}=1-k delta<delta<=1/N.

Thus arbitrarily small positive fractional parts occur at positive indices.
The resulting indices cannot remain bounded, since a finite set of their
strictly positive fractional parts has a positive minimum. Extracting an
increasing subsequence gives the assertion.

Along that subsequence, (2) gives A_L tending to 1 from above. Since
theta>0, j_L tends to infinity, and (5) gives

    e_L -> c+1 >0,
    c(V_L) -> +infinity.                        (7)

This proves the unbounded-center result. It does not depend on a numerical
continued-fraction expansion or a finite list of resonances.

## 5. What this obstructs, and what it does not

There is no finite uniform expanding-center bound covering **all** these
compressed, first-descent-free, adjacent macroreturns. Thus a method which
first compresses arbitrary such returns into single expanding letters and
then assumes one universal C must establish an additional restriction.
The fixed finite bank in [the arbitrary-bank theorem](collatz_bank_reset_20261019.md)
is not refuted: any fixed finite selection has a finite maximum center.

There is also no contradiction with decomposing the present family using
the original small center c of K. In the original K/F alphabet, K^L F^j
has the required expansion-then-reset factorization, but its selected total
multiplier is still greater than one. The descent theorem becomes applicable
when enough actual F blocks are appended to cross below one; indeed j+1
has multiplier below one by (2). Its bounded-carry conclusion retains the
original constant c+1, not the inflated compressed center c(V_L)+1.
Lawfulness of those extra blocks is imposed by their enlarged cylinder;
it is never inferred for every source of the shorter word.

No primitive-return assertion, least-counterexample entry theorem, global
trajectory-frequency statement, or global Collatz conclusion is made.

## 6. Exact resonant controls and replay

The checker verifies the seven proposed pairs

    (L,j)=(1,0),(4,1),(11,3),(40,11),(69,19),
           (3341,920),(6613,1821).

Their center sizes have integer floors

    55, 532, 3949, 12849, 644326, 970938, 1944005,

respectively. For example the final slope gap lies in the rational interval

    29060157/10^12 <= A_6613-1 < 29060158/10^12,

and its center lies in [1944005,1944006). The receipt records rational
enclosures, not floating-point decisions. These witnesses are not claimed
to be consecutive best approximations or a complete resonance census.

The [read-only checker](../experiments/collatz_unbounded_return_centers_20261019.py)
performs **11,581 explicit checks**, including integer certificate (3),
exact formulas (4)-(6), both independent source congruences, source hashes,
and the rational enclosures. For L<=69 it additionally reconstructs both
full parity words and literally replays their trajectories, checking every
lower prefix is above its source and every upper parity guard. For the two
larger controls it verifies exact composite coefficients and congruences,
not a full step-by-step replay. The universal prefix and cylinder assertions
are proved in section 3 independently of those finite controls.

Run from the repository root:

    python3 05-knowledge/experiments/collatz_unbounded_return_centers_20261019.py > /tmp/collatz-centers-normal.json
    python3 -O 05-knowledge/experiments/collatz_unbounded_return_centers_20261019.py > /tmp/collatz-centers-optimized.json
    cmp /tmp/collatz-centers-normal.json /tmp/collatz-centers-optimized.json
    cmp /tmp/collatz-centers-normal.json 05-knowledge/results/collatz_unbounded_return_centers_20261019.json
    cmp /tmp/collatz-centers-optimized.json 05-knowledge/results/collatz_unbounded_return_centers_20261019.optimized.json

Normal and optimized receipts agree byte for byte. The checker is standalone,
uses only Python's standard library, prints to stdout, and makes no network
calls or file writes. The infinite theorem is (7), not finite extrapolation.
