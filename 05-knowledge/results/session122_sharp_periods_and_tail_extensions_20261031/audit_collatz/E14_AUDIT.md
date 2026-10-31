# Independent audit of the E14 finite-content closure

**PASS: complete potential replay, independently reconstructed near-tight graph,
exact arithmetic, native-return lift, and unfinished-suffix margin.**
This proves the stated finite-depth content inequality using the inherited
E10/R1 rational certificate. That inherited large certificate was not rerun.
Neither an all-depth bound nor Collatz convergence is claimed.

## 1. Object and potential interpretation

For actual positive valuation words, retain the original clock A, length L,
and affine map mX+t. The content is
min(v5(m-1),v5(t-m+1)); this equals v5(gcd(3^L-2^A,B)).
The E10/R2 state keeps z=t-m+1 modulo5^10 and the multiplier modulo
5^max(2,min(v5(z),10)). Every actual word is represented, including every
needed compatible lift. An artificial path need not be an actual word.

`audit_integer.cpp` is adapted from the pinned earlier independent z-major
checker, not from producer scheduling. It uses repeated modular halving,
checks all letters1,...,64, all compatible multiplier lifts, all630 paid
canonical entries, and both infinite-letter tails. Its optional shortcut
has been checked algebraically: a compatible target residue different from
one fixes v5(m-1) below the common precision; residue one permits the upper
bound10. If the largest possible right side is below the minimum stored
potential, every such edge is proved automatically. The receipt counts those
edges explicitly as bounded, separately from individually evaluated edges.

The full replay passes on245,312,500 states and18,850,000,000 semantic edges:
1,422,905,275 individually evaluated and17,427,094,725 certified by that
inequality. It checks630 paid entries. The decoded range is[-128,-10],
and terminalV(0,1)=-30. All values are read from the eight-byte file and
validated before storage in signed16-bit memory. The unreached sentinel
is rejected. This is a mathematical inequality audit, not trust in producer
convergence or a check of only sampled transitions.

The producer's40-letter coarse omitted-entry bound fails at this range.
It is not used. Explicit checking through64 and the65+ bounds succeed.
The earlier E7/R2 control also passes; raising its terminal value from-42
to0 is rejected on edge(z,m,a,m2)=(0,1,1,39064).

## 2. Completeness of the three deficit layers

Let the actual excess be E=8S10-5A. Potential deficit is V(state)-E.
An edge adds precisely V(target)-V(source)-8charge+5letter, which the
potential audit proves nonnegative. Consequently every actual endpoint
word with H10=30/8,31/8,or32/8 lies in the reverse augmented graph with
terminal deficit0,1,or2. No word-order or translation realizability is
inferred from a clock pair; the graph only supplies necessary pairs.

For the initial all-one phase, write k=bT_ell+s,0<=s<T_ell,
T_ell=2*5^(ell-1). Since floor(s/2)<5^(ell-1), adding bT_ell inserts b
above the low ell-1 digits without carry. The exact identity

    5k-8 sum_(j=0)^(E-1) floor(floor(k/2)/5^j)
      =2s5(floor(k/2) mod5^E)+10floor(floor(k/2)/5^E)+5(k mod2)

therefore gives the canonical-phase loss

    2s5(b mod5^(E-ell+1))+10floor(b/5^(E-ell+1)).

At loss at most2, b is zero or a single power of5 in the retained digit
range. Thus the complete phase bank consists of s and s+T_e,ell<=e<=E.
Noncanonical phases enter only when the canonical deficit is zero.
The same phase multiplier is verified exactly. This is essential: retaining
only canonical entries would omit long clock pairs such as3906592.

Every letter>=65 has edge slack at least

    -128-(-10)-80+325=127>2,

and entry slack at least-128-80+325=117>2. Such letters cannot occur in
these layers. Thus no untested letter or unbounded entry phase remains.

## 3. Independently reconstructed graph

`independent_corridor.cpp` was written separately from the producer-side
corridor implementation. It decodes states from integer indices, reconstructs
inverse affine transitions, enumerates the full phase bank, and propagates
exact(clock,length) sets with Kahn topological processing. The producer uses
coordinate nodes and depth-first processing. The independent implementation
imports no producer module or graph receipt. Resource guards and cycles
produce errors rather than a successful truncated answer.

Both reconstructions agree on the full E10 census:

-3,908,743 augmented states and3,908,754 edges;
-630 canonical entry checks and1,611 exact long-phase controls;
-10 retained entry records; graph acyclic;
-score30:2 paths, pair(342,326);
-score31:no paths;
-score32:38 paths and11 pairs(A,A-16), with
 A=72,352,392,552,592,1592,6592,31592,156592,781592,3906592.

All12 pairs fail3^L=2^A modulo5^10. In fact ten have exact multiplier
height2, and only(352,336),(552,536) have height3, so all already fail
modulo625. This diagnoses the forgotten multiplier digits; it does not
justify changing the graph's retained precision without a new certificate.

The full potential gives H10>=30/8, and the independent graph excludes
all three smaller integer score levels. Hence every actual non-one word
ending at identity modulo5^10 has H10>=33/8.

## 4. Same-word interpolation and the native lift

Inherited E10/R1 supplies lambda12*A-S10>=eta10. Both it and the new
H10 bound hold on the same actual word, without identifying artificial
paths across abstractions. Their convex combination gives

    lambda14=(24/25)*(5/8)+(1/25)*lambda12,
    eta14=(24/25)*(33/8)+eta10/25
         =5017089801/1220703125>4.

For each non-one primitive11 core, its lambda14 fee is therefore at least
eta14-1. Splitting a non-one primitive13 excursion at actual12 and11
returns preserves all native charges. If C>=B>=1 count its non-one11 cores
and12 blocks, its fee is at least(eta14-1)C-B-1, plus nonnegative all-one
fees. Subtracting eta14-3 leaves

    (eta14-2)(C-1)+(C-B)>=0.

Thus every non-one primitive13 excursion pays at least1+alpha,
alpha=134277301/1220703125>0. The all-one11,12,13 fees at lambda14
are exactly31/125,6/25,1/5.

For the fourteenth-level boundary, record the affine tangent(x,y) inF5^2
at actual13 returns and use Phi(x,y)=((-x) mod5)/5. Any tangent increment
has extra-charge plus DeltaPhi at most1. An all-one13 increment is(4,4)
and has this bill at most1/5. The potential starts at zero and finishes
nonnegative, so completed blocks retain surplus at leastalpha whenever
one is non-one. These are checked across all625 tangent transitions.

An unfinished suffix has no13 return. Its S14 equals its native S12,
because its entering map is identity modulo5^13. If non-one, the inherited
global S12 margin61035156/48828125 exceedsalpha and remains valid at
the larger slope. If all-one, its exact floor sum pays its clock. Therefore
no unfinished suffix is dropped, and every non-one word retainsalpha.

The resulting sharp slope is762939453/1220703125. Equality is witnessed
by2*5^13 consecutive ones, with exact S14=(5^14-1)/4. These are finite
floor-sum identities, not a billion-letter enumeration.

## 5. Reproduction and evidence

Given the regenerated eight-byte potential file:

    c++ -O3 -std=c++17 audit_integer.cpp -o audit-potential
    ./audit-potential 10 integer10.bin 2
    c++ -O3 -std=c++17 independent_corridor.cpp -o audit-corridor
    ./audit-corridor 10 integer10.bin
    python3 -B check_closure.py
    python3 -B -O check_closure.py

The arithmetic checker reads the independently generated corridor receipt,
not the producer's. It supplies53,056 phase controls,625 tangent controls,
all12 endpoint exclusions, exact constituent prices, interpolation, margin,
and sharpness. Its normal and optimized receipts are identical. The full potential and corridor replays were each repeated with NDEBUG;
both normal/NDEBUG receipt pairs are byte-identical. The potential SHA256
was independently recomputed and agrees with the commitment below.

The potential is a regenerable1,962,500,000-byte artifact, not a required
checked-in blob. Its expected SHA256 is
9afc878aaaf0ff11b075678493f6ef935c1735afe4f81bdcd9a2e1282c652c31.
A hash alone proves no theorem; the full replays above were performed.
No Lean formalization, external human peer review, literature priority,
all-layer theorem, or convergence conclusion is claimed.
