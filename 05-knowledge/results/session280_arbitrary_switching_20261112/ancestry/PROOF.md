# Hidden switches are self-delimiting in a thin physical root strip

**Session280, 2026-11-12. PROVED from the inherited active-block transfer.
No Chan–Pak coverage conclusion.**

The fixed-word spectral decoder of Session279 requires a known suffix. Here
we change the question: given a full root packet, can its unknown switching
instructions be recovered? The answer is yes on a physically invariant strip,
for arbitrary active words of arbitrary lengths, even when switching expands
transverse real directions. Together with this session's physical probe chart,
three counts recover an entire hidden switching history, not just a state.

## 1. Objects and inheritance

Use the inherited [Session279 block law](../../session279_interval_tightening_20261112/physical/PROOF.md#8-classifying-a-larger-supply-of-decoding-blocks-and-its-first-obstruction).
An active binary word w starts1, ends0, and has length n>=2. Set

    t=2^n, P=value_2(w), a=P+1,
    Q=sum_j (1-w_j)2^(j-1), b=Q+1, h=g(P).

The recurrence for g is inherited there. The actual rooted tournament
operation sends

    (A,W,d) -> (hA+bW-d, aA+tW, A).                   (1)

Its physical implementation appends the one-apex block and reverses the seam
arc between successive distinguished roots. All claims below concern that
operation, not arbitrary graph switching.

Consider the thin integer strip

    S={ (A,W,d): A positive odd, 0<=d<=A<=W<2A,
                   W=d (mod2) }.                    (2)

The inequalities/parity do NOT characterize physical tournaments. We keep a
specified physical seed. The triangle has packet (3,3,1) and is in S.

Our distant-source inspiration is [recognizable guarded composition](../../support_image_composition_20261026/COMPOSITION.md): a state chart needs an image predicate and initialization, not just rational invertibility. The [doubling-defect arithmetic discussion](../../session116_dimension_owner_towers_20261030/tower/TOWER_DEFECT.md#3-integer-matrices-the-exact-bad-prime-set) suggests keeping real shape and modular algebra separate. No theorem about sedenions is imported here. The connection is explicit: the real ratio locates a dyadic size interval, and an odd modular unit locates the binary word inside it. Neither coordinate alone supplies the predecessor.

## 2. Every active block enters and preserves the thin strip

For even positive P, the recurrence gives

    h=g(P/2)+P/2+(P/2 mod2) >= P/2+1=(a+1)/2.

Because w ends0, Q>=t/2, so b>=t/2+1. Also a>=t/2+1 and a<=t-1.
The broad cone C={0<=d<=A<=W} is preserved by the inherited law.
For any nonzero broad-cone source with A>=1,

    2A'-W'=(2h-a)A+(2b-t)W-2d
           >= A+2W-2d >= A >0.                      (3)

Consequently every active image satisfies W'<2A'. Parity follows from
h,a,b odd and t even: A'=A+(W-d) mod2 is odd, W'=A mod2=d'.
Thus (1) sends broad-cone physical packets into S and preserves S.
No spectral assumption, block-length bound or periodicity appears in (3).

## 3. The length is read from an ordinary real ratio

Let x in S and y=(N,V,D)=M_w x. Then D=A is positive odd and

    V/D = a+t(W/A),  1<=W/A<2,  t/2<a<t.

Therefore

    3t/2 < V/D < 3t.                                (4)

The open intervals (3*2^(n-1),3*2^n), for n>=2, are pairwise disjoint.
Hence y determines at most one possible length n among ALL active words.
The boundaries cannot be physical active outputs. This conclusion relies
on the strip of the INPUT, not only on the fact that the output is thin.
An arbitrary first active block applied outside S is not covered by this
length decoder, though its output enters S by (3).

This is a variable-length code without an externally supplied block delimiter.
The code length comes from an Archimedean scale interval, not from a guessed
number of vertices or from floating logarithms. Integer comparisons in (4)
find it exactly.

## 4. The word is read from an odd modular unit

For the unique candidate t, reduce the second equation in(1) modulo t:

    a = V D^(-1) mod t,  0<=a<t.                     (5)

D is odd, so its inverse modulo2^n exists. The actual a is in this canonical
range. Reject unless a is odd and t/2<a<t. If it passes, P=a-1 has a unique
n-bit representation starting1 and ending0: this is w. Thus the entire word,
not merely its length, is determined. Compute its h,b and put

    A_old=D, W_old=(V-aD)/t,
    d_old=hD+b W_old-N.                              (6)

Reject unless (6) lies in S. Recomposition by (1) is exact.

**Theorem 1 (disjoint active images on S).** The map

    (w,x) -> M_w x,  w ANY active word, x in S,

is injective. Its image has the exact algorithmic predicate (4),(5),(6) and
membership of (6) in S.

Proof: (4) uniquely determines n; (5) uniquely determines a and w; (6)
uniquely determines x. All necessity statements follow from(1). Conversely,
a passing preimage in S recomposes to the supplied y. This proves the image
predicate for FORMAL strip packets. Physicality still requires ancestry to an
actual seed. The proof does not claim that every y in S has a predecessor. ∎

This strengthens the known-schedule inverse: no word, length, depth or common
Perron direction is supplied. It is unrelated to convergence of normalized
ratios. Real expansion under a bad switched product does not affect it.

## 5. Full unknown-itinerary recovery and its physical gate

Each active block satisfies h>=3,b>=3. On C,

    A'=hA+bW-d >=5A.                                (7)

Fix any actual seed packet s in S, together with a construction of its rooted
tournament. Given y, repeatedly apply the unique predecessor until its first
coordinate is at most that of s. Accept iff the final packet equals s.
Reversing the recorded words constructs the unique itinerary from s to y.
The empty itinerary is accepted at y=s. Every nonempty predecessor step
reduces A by at least5, so the algorithm terminates after at most
floor(log_5(A_y/A_s))+1 tests, with exact integer arithmetic.

**Theorem 2.** For a fixed physical seed packet s in S, distinct finite words
over the infinite alphabet of active blocks give distinct output packets.
The deterministic inverse accepts precisely the packet orbit of s and supplies
a literal tournament construction for each acceptance.

Proof: use Theorem1 inductively; (7) excludes a nonempty cycle back to s and
ensures termination. Physicality follows by composing actual blocks with the
stored seed. Different rooted tournaments with the same seed packet are NOT
distinguished; this is packet injectivity, not an isomorphism theorem. ∎

For a bank of seeds, one deterministic ancestry chain suffices. Each encounter
with a bank packet gives a valid owner. Different bank seeds can lie on one
another's orbit, so uniqueness of seed-plus-itinerary is asserted only for a
fixed seed or an ancestry-antichain bank. The decoder must not silently choose
one such owner and claim unique origin.

A seed outside S may first be advanced by one active block. One can store its
actual resulting packet as a seed in S. This is a paid physical initialization,
not a free claim that the original source already satisfied(2).

## 6. Combining with physical probes

[The parallel-probe chart](../wildcard/PROOF.md#2-equal-order-probes-give-a-unimodular-integer-chart)
recovers (A,W,d) from A and the two counts after separately appending100 and110,
without spectral assumptions. Combining that chart with Theorem2 gives:

**Corollary.** From those three exact counts and a fixed physical seed packet
in S, one either recovers the unique entire active switching itinerary or
correctly rejects membership in this orbit. Neither the hidden words nor their
lengths need be given.

Only one source is being probed; the two enlarged tournaments are experiments,
not additional coordinates one may independently prescribe. This is not an
inverse solver from the one desired target N alone. It changes the rigorous
question from “can a switching system preserve a scalar approximation?” to
“which minimal exact observations make its whole construction recognizable?”

## 7. Relevance and limits for count coverage

The [capacity lane](../coverage/PROOF.md) counts presentations. Theorem2 now
shows those presentations are genuinely different ROOT PACKETS for a fixed
thin seed, even when word lengths vary. Thus its growing branch supply is not
merely repeated spelling of the same root interface. It can still collide
under the projection (A,W,d)->A. Controlling these scalar fibers and showing
that each prescribed odd target has a nonempty physical fiber remain open.

No claim about all positive odd values, eventual thresholds, uniqueness of
unrooted tournaments or a complete structural grammar follows. The exact
positive ingredient is an infinite-alphabet, nonperiodic, spectrum-independent
physical ancestry code, with a constant number of physical observation counts.

## Reproduction

`python3 verify.py` and `python3 -O verify.py` from this directory, or supply
an absolute path. Explicit check gates survive optimized mode. The verifier
checks all511 active words through length10 on varied thin packets; literal
small physical examples; all343 three-block schedules over the first7 letters;
120 long unknown schedules; and initial/boundary controls. The proofs above,
not the finite tests, establish all-length claims.
