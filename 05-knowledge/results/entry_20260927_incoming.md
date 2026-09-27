# Independent audit of incoming coefficient-descent claims

**PROVED elementary repairs / FINITE-EXACT independent checks / OPEN
global first-descent coverage.** This audits incoming commits4403f0015 and
5a37804ed, especially the original versions of
[THM-4512](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md)
and [the residual note](collatz_precision_residual_20260926.md). Those
documents may subsequently contain the repairs described here.

## 1. Exact words and coarse cylinders are different objects

Let w=(v_1,...,v_j), all v_i>=1, A_t=sum_(i<=t)v_i, A=A_j, and

    S=sum_(t=0)^(j-1) 3^(j-1-t)2^(A_t), A_0=0.

The exact-word cylinder has modulus2^(A+1), with residue

    eta=(2^A-S)*3^(-j) mod2^(A+1).                    (I1)

Its positive members have exactly that valuation word and
U^j(n)=(3^j*n+S)/2^A. Its density among odd integers is2^(-A).

The incoming residue

    rho=-S*3^(-j) mod2^A                             (I2)

instead describes a coarse cylinder: the first j-1 valuations are exact,
and the last valuation is at least v_j. On this cylinder

    Q=(3^j*n+S)/2^A is an integer,
    U^j(n)=oddpart(Q)<=Q.                            (I3)

To prove the prefix claim, reduce the terminal divisibility congruence
modulo2^(A_t+1), for t<j. The contribution at A_t has odd coefficient;
later contributions vanish. Consequently the corresponding intermediate
quotient is odd. The last quotient need not be odd. Requiring that final
quotient to be odd adds exactly the bit in(I1).

**Minimal hostile:** w=(1), A=1, n=1 satisfies(I2), but Q=2 and U(1)=1.
The claimed exact word and affine equality on modulus2^A are false.

**Strongest survivor:** if2^A>3^j and n>S/(2^A-3^j), then Q<n and hence
U^j(n)<n even on the coarse cylinder. Thus the sufficient descent
certificates survive. The affine threshold is an iff only on the exact
word cylinder; the coarse formula supplies an upper bound.

The implemented audit checks6240 exact/coarse instances: all words of
length1..4 with valuations1..5, four positive lifts of each residue in
each interpretation. It independently reproduces both conclusions.

## 2. Finite claims preserved by independent exact computation

The elementary bound

    S<=2^A*((3/2)^j-1)

is valid. The audit uses A=bit_length(3^j), not a floating logarithm, and
checks for every j=1..5000 that

    (3^j-2^j) < (2^A-3^j)*2^j.

The maximal ratio is exactly211/416 at(j,A)=(5,8). Consequently N=S/(2^A-3^j)
is below2^A in this range. Each coarse class has at most its least positive
representative left uncertified by that threshold. Exact classes, with
larger modulus2^(A+1), also have at most their representative uncertified.

For first coefficient descent at j<=14, independent recursion reproduces
606746 labeled coarse cylinders and606746 corresponding exact-word
cylinders, with the sole threshold exception n=1,w=(2),N=1 in either
interpretation. The coarse cylinders can overlap; their number is not a
count of disjoint exact-word classes in the coarse interpretation.

The original script does not run to40 valuations beyond critical as its
output says. Its actual cutoff is total A<=40 for j=1, and A<=41 for
j>=2. This still suffices. Before first coefficient descent,
2^(A_t)<=3^t, so

    S<=j*3^(j-1).

For every j<=14 the audit checks j*3^(j-1)<2^41-3^j. Hence N<1 throughout
the omitted tail A>=41; no positive representative there can fail the
sufficient threshold. This supplies an explicit valid tail proof.

An effective irrationality measure would give a computable eventual
cutoff for the analogous all-j gap inequality, with its constants retained.
It does not automatically cover every j>5000. An explicit estimate and
verification of any intervening finite interval are required before that
claim is promoted beyond the checked range.

The incoming vectorized test reports equality of actual and coefficient
stopping times for odd3<=n<=10^7. This audit does not rerun that universe.
It independently confirms that the floating coefficient comparisons used
through the reported maximum step155 agree with integer thresholds. It
does not independently certify the vectorized int64 orbit arithmetic or
its overflow bounds. Preserve the result as producer-reported finite
computation, separately from the independently checked statements above.

## 3. Residual densities survive; two claimed identifications do not

An integer-only dynamic program, grouping exact valuation words by their
total A and weighting each word by2^(-A) among odd inputs, reproduces

    D(41)=1530343662856563/576460752303423488,
    D(60)=12185976031772265023015553/19807040628566084398385987584.

These are approximately0.002654723078269616 and0.000615234565339227.
The class measure uses(I1), not the twice-as-large coarse cylinder(I2).
This validates the finite density computation after the dictionary repair;
it does not validate an orbit-wise independence claim or a pointwise
termination claim.

The complement of the old finite bank is not the classical no-coefficient-
descent set. The audit independently regenerates the171 bank rows and
verifies that4091 belongs to none of them, whereas its actual orbit is

    4091,6137,4603,6905,5179,7769,5827,8741,1639.

Its first actual and coefficient descents both occur at step8. Thus a
bank-unresolved source can already possess a short classical coefficient
certificate. The bank cutoff, consumed valuation sum, and recreated
precision are distinct coordinates until an explicit map identifies them.

Likewise, ordinary entry into the region of local one-step descents is
unconditional for every odd source>1: the initial valuation1 run has length
v2(n+1)-1 and ends at an integer1 mod4. Such entry does not pay earlier
growth relative to the original source. Equivalence with Collatz applies
to eventual descent below each original source, followed by strong
induction, not to merely visiting a moving local-descent region. See the
[source-preserving entry analysis](entry_20260927_entry.md).

Finally, literal Haar-genericity is not what Collatz asserts for positive
integers. A convergent Syracuse orbit is eventually1, so its valuation
sequence is eventually constant2. This is not a Haar-generic geometric
valuation sequence. A local/global analogy must not identify those
distinct predicates.

## 4. Reproduction and audit boundary

    python 04-computation/experiments/entry_20260927_incoming.py
    python -O 04-computation/experiments/entry_20260927_incoming.py

The [retained output](entry_20260927_incoming.out) is identical in both
modes. No incoming producer module is imported. All validation checks
remain active with optimization enabled. This audit does not claim a
literature-priority result, verify every literature attribution in the
incoming note, or establish its stated Syracuse asymptotic root limit.
Its purpose is to preserve the valid finite results and sufficient
certificates while repairing the demonstrated type and quantifier errors.
