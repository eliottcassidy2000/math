> **Curated publication copy.** Executable proofs and certificate data are in
> `code/`; independent audits are in `audit/`. Run `verify_all.py`.
> Historical exploratory scripts mentioned below are intentionally omitted;
> their finite searches are not dependencies of the certified theorems.

# Ternary tests and the exact word spectrum

Status: PROVED HERE, with exact rational certificates and separate checker.
This is a restriction on two-fixed-ear/transitive-spine accumulators, not a
nonattainment result for Hamiltonian-path counts.

Use the letter convention and36-state transfer in
`TWO_EAR_TRANSFER_AND_COMPILER_OBSTRUCTION.md`.

## 1. One-letter ternary accumulators are impossible

Suppose every actual descendant with H>1 must offer the three actions
H→3H+d, d∈{-2,0,2}, each by appending one spine vertex. Letter3 is a universal
sink and preserves H, so cannot realize any action (H=3H+d impliesH≤1).
Distinct digits require distinct physical letters. Thus letters0,1,2 must
all be available at every descendant, regardless of any phase bookkeeping.

An elementary obstruction is to repeat letter1. In an accepted xy path,
no appended1 can belong to P, since x beats it. In an accepted yx path,
no appended1 can belong to Q, whose last vertex must beat x. If the seed
spine has m0 vertices, this gives

    H(seed followed by1^k) ≤ 2·3^m0·2^k.

But any sequence of the proposed actions satisfies

    H_k ≥ 3^k(H_0−1)+1,

a contradiction forH_0>1. Digit assignments may change by phase or count;
that does not affect either bound.

### Exact finite algebraic certificate, including the boundary H=1

Let c be the output row and
D3=c(M0+M1+M2−9I). The sum of the three digit actions forces D3v=0, and
all descendants force D3M_wv=0 for every word w over{0,1,2}.
Its row-span closure has dimensions1→4→13→20→20. It contains the functional

    ell(v)=sum_{p=0}^2 v[(xy,p,empty,empty)],

which equals1 for EVERY actual state, including the empty spine. Hence no
actual state can satisfy all these sum identities. This certificate is
stronger than a finite seed search and makes no affine-ray assumption.

`ternary_certificates.json` expresses ell as an exact rational linear
combination of18 rows D3M_w, with |w|≤3. The checker independently rebuilds
the transfer by dictionary-based partition-state transitions and verifies
all36 coordinates using Python `Fraction`, without SymPy or producer imports.
It also checks the convention-sensitive control: seedword1 has child counts
[3,1,5,1], so the first sum identity alone holds, but its descendant closure
cannot hold.

## 2. Fixed length-two radix9 banks are also impossible

Here the digits are {-8,-6,-4,-2,0,2,4,6,8}. A fixed physical bank requires
nine distinct length-two words, all available at every descendant; digit
assignments themselves may vary.

If a word has x-bit1 at every position, or y-bit1 at every position, its
repetition gives at most O(2^m) accepted paths, by the same forbidden-run
argument. This contradicts the required radix9 lower bound on repetition.
There are exactly seven such words:

    11,13,31,33,22,23,32.

Consequently the bank is forced to be exactly

    00,01,02,03,10,12,20,21,30.

For T_ab=M_b M_a, put D9=c(sum_bank T_ab−81I). Its closure under the nine
bank matrices has dimensions1→9→20→20 and again contains ell. The second
certificate expresses ell as20 exact defect rows, with suffix length≤4.
It excludes every actual seed, not merely short seeds.

**Scope warning:** excluding a slow word by repetition assumes it stays
available. This does NOT exclude arbitrary phase-dependent changes of the
physical nine-word bank: a slow word may be transiently available and then
unavailable. The certificate does exclude phase-dependent permutations of
digits within the fixed forced bank. Do not upgrade the fixed-bank conclusion
to changing banks.

## 3. A common triangular filtration gives every word's spectrum

For one18-state ear-order block, the letter matrix is the Kronecker sum of
three local append operators, on dimensions3,3,2. Explicit common bases are:

* P: (e_empty−e_good, e_good−e_bad, e_good). Its reset has diagonal(0,0,1).
* Q: (e_good−e_bad, e_good, e_empty). Its diagonal is(0,1,0); the final
  column records whether starting Q is allowed, but does not change it.
* R: (e_started,e_empty). Its diagonal is(1,0).

Each local matrix is upper triangular in its indicated basis, for EVERY
letter. Therefore every full letter matrix is upper triangular in the
corresponding tensor basis, and the diagonal is independent of the letter.
Its per-order multiplicities are encoded by

    (2+t)^2(1+t)=4+8t+5t^2+t^3.

For every nonempty word w of lengthL, the full36-state transfer therefore has
characteristic polynomial

    z^8 (z−1)^16 (z−2^L)^10 (z−3^L)^2.

This statement includes multiplicity but does not claim diagonalizability;
Jordan blocks may contribute polynomial factors to count sequences.
`spectral.py` verifies the explicit common basis, all four letter matrices,
and all sixteen length-two products in exact arithmetic.

### Consequence for stateless repeated-template recursions

Suppose one fixed length-L template, repeatable forever, obeys
H→BH+d with B>1, |d|≤B−1, starting atH>1. Its iterates have form

    H_k = B^k(H_0+d/(B−1)) − d/(B−1),

with strictly positive coefficient of B^k. By Cayley–Hamilton, this scalar
sequence is annihilated by the word matrix's characteristic polynomial.
Substitution forces B to be one of its eigenvalues, hence

    B ∈ {2^L,3^L}.

Thus a STATELESS fixed-radix template bank cannot evade these restrictions
merely by allowing different word lengths: every repeatable template must
have the same length, and the radix must be an integer power of2 or of3.
Combined with the prior theorem, all stateless binary-power accumulators on
this carrier are excluded, even if unequal lengths were initially allowed.
Subsequent exact certificates exclude L=3 and L=4 as well; see
`RADIX27_REPEATED_TEMPLATE_OBSTRUCTION.md`. The untouched longer-macro
frontier in this package is L≥5.

### Necessary cycle law for finite-phase variable-length models

For a phase cycle of k steps, let L be total appended length. Its scalar
composition is H→B^kH+D with |D|≤B^k−1. Repeating the cycle and applying the
same spectral argument gives

    B^k ∈ {2^L,3^L}.

This is a necessary cycle law, NOT a sufficient construction or a complete
obstruction for phase-dependent variable-length models.

## Verification

Run `python -B verify_all.py` from the package root (or invoke its absolute
path). The curated packet includes exact witness data and independent verifiers,
not exploratory certificate generators. The ternary standalone checker reports
72 exact coordinate checks, four rejected mutations, and the seed1 control.
All bundled normal/optimized receipts match; independent audits reconstruct
the row closures and spectral calculations separately.

### Two hostile spectral controls

The exact checker also establishes

* H((01)^k)=(4^k+6k−1)/9 for k≥1, by the state annihilator
  (M_01−4I)(M_01−I)^2 M_01 v_empty=0 and three initial values.
* H(2330^k)=14·3^k−17 for k≥1, using the square-free annihilator
  M0(M0−I)(M0−2I)(M0−3I)=0 and three initial values.

Thus polynomial lower modes really occur; conversely the2^k mode can cancel
on actual states. Neither diagonalizability nor universal noncancellation
is a valid shortcut to the compiler obstruction. The second family iterates
H→3H+34, outside the centered ternary digit range; it is an actual explicit
infinite family, not a coverage theorem.

### Constant-defect strengthening

The coefficients of each conserved-unit certificate sum to ZERO. Therefore,
for either stated physical bank, it is impossible even to have

    c(sum_bank T−B^2 I)v = K

with the SAME constant K at every bank descendant: the certificate would give
1=K·0. This excludes changing the centered digits to any other fixed-sum digit
set while retaining the same physical bank. It does not by itself force that
bank for arbitrary digit bounds or phase-dependent changes of physical words.
The independent checker explicitly verifies the zero coefficient sums.
