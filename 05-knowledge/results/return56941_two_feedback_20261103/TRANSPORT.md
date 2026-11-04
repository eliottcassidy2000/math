# A uniform1-ready to3-ready map and its hidden third mode

**PROVED finite-dimensional reduction and exact transport law, with
FINITE-EXACT coefficient certificates and extra-size controls.**
The construction is the actual second-port birth in
[TWO_FEEDBACK F2](TWO_FEEDBACK.md), on the repaired source
H(T_m)=1243*2^m-1061. Keep the first clone r=31 marked. Let U_m=T_m+q
with the literal named second-port cut described there.

## T1. Five samples determine these responses for a proved reason

Put k=m-1: r remains a fixed vertex, and k further uniform clones vary.
For either family, every one of

    H, E_r(1), S_r(1), H(T-r), E_r(2), S_r(2)

is an integer linear combination of0^k,1^k,2^k,3^k,4^k, with0^0=1.
Here S_r(2) counts ordered spanning covers whose first path starts at r.

Proof: delete the two feedback vertices from one or two path bins. There
are at most four increasing base pieces. The variable clones form an
interval in the base order, hence each piece has at most one insertion
gap for them. Fix the words of all other vertices, including the marked r,
the order of the feedback vertices and the bin division. The legality of
each insertion depends only on which candidate gaps are occupied: their
first/last external incidences are uniform. Empty gaps retain their own
adjacency guards; this explicitly includes the orphan-path boundary.
For each permitted occupied subset J, assignments of k distinct clones
onto J number |J|!S(k,|J|), an integer linear combination of j^k for
0<=j<=|J|<=4 by inclusion-exclusion. All internal clone orders are forced.
Summing over the finite external data proves the claim. Ending/starting
constraints and deletion of r restrict that data without increasing four.

The5-by5 sample matrix V_(k,j)=j^k,0<=k,j<=4, is invertible. Therefore
five literal samples determine the entire sequence. This is not unbounded
extrapolation from a fitted pattern: the prior occupancy proof bounds the
function space. The separate m6,m7 computations are additional controls.

## T2. The exact joint response law

Write unprimed fields for T_m and primed fields for U_m. The independent
endpoint counts in `feedback_ports.json`, interpolated in the proved
five-mode space, give for every m>=1:

    H'       = 3H,
    E_r(1)'  = 3 E_r(1),
    S_r(1)'  = 3 S_r(1),
    H(U_m-r) = 3 H(T_m-r),
    E_r(2)'  = 3 E_r(2),
    S_r(2)'  = 3 S_r(2)
                 +135276*(3^(m-1)-2^(m-1)).             (T1)

For the last difference, the exact spectral coefficients at roots0..4
are(0,0,-135276,135276,0). The other five difference vectors are zero.
The saved root profiles give E_r(2)-H=0 at m1 and1 for every m>=2;
therefore E_r(2)'-H'=0 at m1 and3 for every m>=2. Thus this is an actual
**uniform1-ready to3-ready map**, not just a scalar3H family. It does not
return the3-ready state to either ternary-ready input language.

The third-mode term is already nonzero at m2. Its first values are
0,135276,676380,2570244,8792940,... . Four familiar rooted response fields
and H all scale by3 while this fifth response does not. This provides a
specific mechanism for growing alignment memory in the live construction.

## T3. Golden coordinates, F49/F343 and the six-step correction clock

Define the inherited golden vector

    g=(H,X,Y), X=(S_r(2)-E_r(1))/2,
               Y=(S_r(2)+E_r(1))/2.

Then(T1) gives the exact integral transport

    g' = 3g + d_m*(0,1,1),
    d_m=67638*(3^(m-1)-2^(m-1)).                         (T2)

In the inherited F49-plus-charge chart

    (a,b,c)=(X+4(H-Y),Y,H-Y) mod7,

this becomes

    (a',b',c')=3(a,b,c)+d_m*(-3,1,-1) mod7.             (T3)

Because67638=4 mod7 and3/2=5 has order6 in F7*, the correction vanishes
**iff m=1 mod6**. This is an exact growth-clock condition, not a claim
that the physical tournaments or their integer port tensors repeat.

In the current F343 vector codec z=H+Xt+Yt^2, with t^3+t^2+1=0,

    z'=3z+d_m*(t+t^2),
    (z')^(7^i)=3z^(7^i)+d_m*(t^(7^i)+t^(2*7^i)),       (T4)

for i0,1,2. The added coefficient is in F7, so these equations explicitly
retain Frobenius phase. This is the current response codec; the ancestral
interval-clock states remain separate and frozen.

Every one of the five inherited actual Q4 matrices L acts linearly, so
L(g')=3L(g)+d_m L(0,1,1). `feedback_transport.py` checks this against all
five recorded presentations at m1..7, as well as all three F343 frames
and the F49-plus-charge chart. No scalar identification replaces the
underlying tournament or authorizes an unguarded complement.

## T4. Integer sample guards explain the bad characteristics

The same small computation exposes an exact integer boundary. Let S be
the lower-triangular Stirling matrix S_(k,r), and B_(j,r)=binom(j,r).
Both are unimodular, and

    V = S diag(0!,1!,2!,3!,4!) B^transpose.

Thus V has Smith factors1,1,2,6,24 and index288. Integer reconstruction
from arbitrary proposed sample values needs the corresponding divisibility
guards. At characteristic7 the sample map is invertible, so these five
samples recover the modal coefficients without a precision loss. At
characteristics2 and3 it is not invertible; collapsing to field values
discards distinctions that integral observations can retain.

The source and target here are actual matrices, so their samples already
satisfy the realization guards. The Smith calculation concerns the sample
lattice, not a full classification of physically realizable port tensors.
It supplies an explicit link between ordered-run memory, integer precision
and Frobenius transport; it does not identify288 with an unrelated object.

The next operator to seek is a realization-preserving return from this
3-ready family, carrying d_m and the joint endpoint data. The missing
56935+2 lane and the81/27 hostile prevent replacing that problem by a
scalar recurrence or by equality in a finite-field chart.
Incoming [SF3](../second_feedback_18979_20261103/PROOFS.md) additionally
excludes repairing the desired +2 family to one-ready states inside any
fixed uniform-clone exterior; simply adding more feedback while demanding
that same closure is insufficient.
