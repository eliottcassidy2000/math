# Seven physical measurements recover the full coloured-core compiler

**PROVED elementary consequence of inherited theorems; FINITE-EXACT controls.**
The all-height covering statement and Chan–Pak remain OPEN. This note does
not claim that the old alternating identity is new.

## 1. Recover the March mechanism before changing coordinates

The source is [THM-016, claim-b-alternating-sum](../../../../01-canon/theorems/THM-016-claim-b-alternating-sum.md),
proved in March 2026. For any nonempty tournament C and vertex v,

    sum_(S containing v) (-1)^(|S|+1) H(C-S) h_start(C[S],v)
        = h_end(C,v).

Its path-reversed version exchanges start and end. This identity applies
exactly to the residual-H-weighted paths of the
[Session251 compiler](../../session251_creative_completion_20261111/ternary/PROOF.md).
The source/target map groups a directed path by its vertex support. It
preserves the residual Hamiltonian-path weight and distinguished endpoint.
Discarding the endpoint before applying the identity loses the cut data.

Use the inherited notation: O is the set of core vertices beaten by the
new apex, I=C minus O, d=H(C), and k_±,i_±,o_±,j_± are evaluations of the
four weighted path polynomials. Define ordinary endpoint counts

    s = sum_(v in I) h_start(C,v),
    e = sum_(v in O) h_end(C,v).

Grouping by support and applying THM-016 gives, for every coloured core,

    k_-=-d,        i_-=-s,        o_-=-e.                 (1)

The first identity is obtained by summing the endpoint identity over all
start vertices: each directed path is counted once. No parity restriction
on the order of C, regularity, or strongness is required.

By the inherited odd-cycle formula [THM-002, ocf](../../../../01-canon/theorems/THM-002-ocf.md),
or by Session251's compiler at A=B=0, the actual apex extension C+x has

    h := H(C+x) = d+j_++j_-.                            (2)

Here h is a physical path count, not a formal choice of a constant term.

## 2. A nonnegative seven-coordinate normal form

Put L=d+k_+, i=i_+, o=o_+. For arbitrary nonnegative integer flank
addresses A,B and their parities alpha,beta, the exact count is

    Phi(A,B) = h + d(g(A)-1+g(B)-1) + L AB
                   + i A + o B + s alpha + e beta.       (3)

This follows by substituting (1)–(2) in Session251's formula. In particular
the apparent alpha*beta coefficient d+k_- vanishes for EVERY tournament
core, not merely the five chosen cores. All seven coefficients

    (h,d,L,i,o,s,e)

are nonnegative integers, d>=1, h>=d, and L>=2d. The last inequality follows
because the spanning core paths alone contribute d to k_+.

Moreover L is even, i+s is even, and o+e is even. Indeed the positive and
negative evaluations of an integer-coefficient polynomial agree modulo2,
so (1) gives k_+=d, i=s, o=e modulo2. Since g is odd and h is odd, (3)
keeps the output parity visible. On writing A=2a+alpha, B=2b+beta,

    (Phi-1)/2 = (h-1)/2 + d*((g(A)-1)/2+(g(B)-1)/2)
                +(L/2)AB+i a+o b+((i+s)/2)alpha+((o+e)/2)beta.

This is a sum of nonnegative integers. It sharpens the interpretation of
the inverse, not the inherited target bound: the minimum over parity
corrections in that bound already equals h-d by (1). In mixed contact,
Phi<=N implies i A+o B+L AB<=N-h. Uniform zero axes still require the
inherited g-value normalization; they are not silently bounded by i or o.

## 3. Seven small actual contexts determine every flank length

Write p_ab=Phi(a,b). The seven evaluations at

    (0,0), (1,0), (3,0), (2,0), (0,1), (0,3), (1,1)

recover every coefficient exactly:

    h = p_00,
    i = (p_30-p_10)/2,       s = p_10-h-i,
    d = (p_20-h-2i)/2,
    o = (p_03-p_01)/2,       e = p_01-h-o,
    L = p_11-h-i-o-s-e.                                  (4)

The proof uses only g(0)=g(1)=g(3)=1 and g(2)=3. Each evaluation constructs
an actual tournament by Session251's arc constructor; the largest has
|C|+3 vertices. Thus for cores of order at most five, seven contexts of
order at most eight determine every arbitrary-length two-flank response.
This is a complete observer for this response function. It is not a
complete tournament invariant or a state sufficient for arbitrary future
births, substitution, or rerooting.

The seven coordinate functions in (3) are linearly independent on all
nonnegative A,B, since (4) explicitly recovers their coefficients. The
actual coefficient vectors of coloured tournament cores also span a
seven-dimensional rational space: `probe.json` gives seven concrete
linearly independent vectors from the complete core bank through order
five. This rank certificate is independently checkable by rational row
reduction. Hence fewer than seven scalar linear measurements cannot
recover an arbitrary vector in this linear span. This is a linear-observer
minimality statement, not a claim against nonlinear encodings of a finite
bank or a requirement that a coverage proof use seven states.

## 4. Why this is useful for the five-object search

The bank has five underlying cores but its response space has rank seven,
and its integer addresses are unbounded. These are different types of
cardinality. A correspondence to a five-object family must specify how
that family's operations act on (3), or retain a separate coordinate
update. Treating five names as five scalar states erases the required data.

Conversely, an old geometric or combinatorial operation can be tested
cheaply and decisively: construct its proposed image and compare just the
seven contexts (4). Agreement proves equality of the entire fixed-core
two-flank response; a disagreement identifies an exact lost coefficient.
This is a reusable consequence of recovering the March mechanism, without
asserting any new symmetry between the five core types.

## 5. Checks and boundaries

Run `python3 -B 05-knowledge/results/chanpak_five_core_bridges_20261110/algebra/probe.py`.
It enumerates all 20 tournament isomorphism classes of orders 1..5 and all
470 cuts, including uniform cuts. Direct permutation endpoint counts are
compared with (1); actual subset-DP values in seven small contexts recover
(3); the old and new formulas agree for all 25 pairs A,B in0..4, giving
11,750 comparisons. Separate actual-child DP checks use all A,B in0..2.
The rank witness is exact rational arithmetic. Reproduction does not
assume the five-core covering conjecture. No all-height spectrum coverage
follows from interpolation alone.
