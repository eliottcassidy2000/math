# Independent complete support comparison through 2^24-1

**FINITE-EXACT: PASS.** This strengthens the earlier prefix control without
extending the target range. The all-height five-core covering statement
remains OPEN. No classification of every graph embedding follows from a
count-support computation.

The original signed44-signature inverse and an independent positive90-cut
inverse give identical support bytes for every positive odd target through
16,777,215. Each byte retains all five separate core-support bits. Thus the
comparison checks all8,388,608 targets individually, including which core
realizes each one, rather than merely comparing histograms or the two holes.
The complete support SHA256 is

    cfc1065f31b4d76b69ec449df45946404821c3eb8b63590dc6f2ebaf7a68e306

[full_support_audit.json](full_support_audit.json) records both programs'
source hashes, the independently recovered observer hash, and the exact
pair counts:419,808,372 in the old traversal and888,728,504 in the new one.
The original full JSON result is reproduced exactly. Large temporary
support binaries are compared byte for byte and then discarded.

## Three independent implementation choices

1. **Actual observer recovery.** For every one of the90 raw core cuts,
   construct the seven literal children at (A,B)=(0,0),(1,0),(3,0),(2,0),
   (0,1),(0,3),(1,1). Count their Hamilton paths by a first-vertex recurrence,
   without importing the inherited compiler or path-signature extractor.
   The proved seven-coordinate form determines h,d,L,i,o,s,e from these
   seven counts. In particular

       i=(H30-H10)/2, s=H10-H00-i,
       o=(H03-H01)/2, e=H01-H00-o,
       d=(H20-H00-2i)/2, L=H11-H10-H01+H00.

   Independently check recovered d against H(core), all nonnegativity, and
   parity. [positive_raw_signatures.json](positive_raw_signatures.json)
   retains all630 literal counts. These formulas use the already proved
   response theorem; they are not curve fitting in place of that theorem.
2. **Direct digit sum.** The new C++ implementation computes g separately
   for every integer using the binary inversion sum: start at1, scan from
   the most significant end, and at each zero add the preceding binary
   prefix plus the immediately preceding bit. No recursive g-table entry
   is used. The original implementation fills g by its halving recurrence.
3. **Different inverse and no cut quotient.** The original implementation
   loops A first through44 signatures after flank-swap support quotienting.
   The new one retains all90 raw cuts, loops B first, and subtracts the
   nonnegative fixed-B g/parity contribution before bounding A. It evaluates
   only the positive seven-coordinate formula.

The old program gained only an optional binary-mask output argument; its
enumeration and count arithmetic were unchanged.

## Why the independent inverse is complete

Write c=N-h. For a mixed cut, both linear coefficients i,o are positive.
Fix B<=c/o and put

    base(B)=h+d(g(B)-1)+oB+e(B mod2).

The remaining terms are d(g(A)-1)+(LB+i)A+s(A mod2), all nonnegative.
Thus if base(B)<=N, every solution has

    A <= (N-base(B))/(LB+i).

The program evaluates every such A. It never assumes monotonicity of g.

For a uniform cut one slope vanishes. The corresponding free axis has
count d*g(v), because h=d and its parity coefficient also vanishes.
Repeatedly remove trailing ones using g(2q+1)=g(q), obtaining an even
representative2q or zero. Since g(2q)>=q+1, all possibilities are represented
by0<=q<floor(N/d). Off this axis, if o=0 then A>=1 and LAB<=c gives
B<=c/L; the same fixed-B bound on A applies. If i=0, B=0 is precisely
the already handled axis and every B>0 gives positive denominator LB.
This exhausts all cases. No bound on graph order is used.

The construction uses canonical binary address words; arbitrarily many
neutral leading-zero flank vertices do not change the support value.
The source grammar and inherited axis normalization therefore remain exact
for count existence at the finite target bound.

## Reproduction

Run `python -I -B verify_full_support.py` from this directory (C++17/g++
available). It rebuilds both programs in a temporary directory, recovers
all observers from literal graphs, reruns both complete inverses at2^24-1,
compares the complete support arrays and retained result, and prints the
audit receipt. `--bound 4095` supplies a fast smaller control. `--fail`
is an intentional nonzero hostile control.
