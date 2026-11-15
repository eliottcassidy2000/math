# Same-host dyadic reader fork

The seeded historical draw is in DRAW_PORTS.txt: Python seed 36006213, Git's ordered
Markdown list under 05-knowledge/results filtered to January–September 2026.
Two operative prompts were reset_20260926_colours.md (same numerical value
with different marked atoms) and glued_xor_20260921_conjugacies.md (a
projection can lose an essential shell coordinate). Here the marked atoms
are actual choices of which vertex the future closing apex contacts. The
quotient coefficient agrees exactly, but the offsets and a higher slot
response do not.

Use the fixed five-vertex tournament with outgoing bit masks
(18,8,3,5,14), indexed 0 through 4. The two readers are O={1,2} and
P={0,1}. Let W count full ordered two-path covers with first end and second
start in the marked set. U counts full ordered two-slot covers with the
first slot optional and second path starting in the marked set. E counts
paths ending in the marked set. Here H=13 and Q=P_B(2)=104; optional covers
permit empty paths. Both readers have U=32 and E=6, but W_O=8 and W_P=6.
Hence r=H-W is 5 or 7. A literal closure appends TT_t ordinally then adds
a single apex whose out-neighborhood is exactly the marked set. For t>=1
its H is r+2^t U, giving consecutive odds 2^(t+5)+5 and 2^(t+5)+7,
including 69 and 71. This is a physical pair not a signed formal state.
The same rows are used for the entire two-element window.

The cyclic-context law is exact: replace B by A⇒B⇒S, keeping marks in B.
H and both-slot Q multiply through ordinal factors. In W, the first path
ends in B, so all A belongs to it, and the second starts in B so all S
belongs to it; W scales by H(A)H(S). In U, the second starts in B, so
all A belongs to first slot and the suffix is distributed over either
slot; U scales by H(A)P_S(2). E scales by H(A) for empty S and vanishes
for nonempty S. Thus the offset gap is exactly 2H(A)H(S), with shared
2^tU=2^t·32H(A)P_S(2), and shared base H,Q,E. A and S may be empty
with both counts defined as one.

## Exact limitation: closed outputs do not share the quotient again

For an arbitrary marked B, write h_D=H(B[D]) (including h_empty=1)
and let (r_D,u_D,e_D) be its subset reader profile. Split an optional
two-slot cover of the closed tournament by which slot contains the
apex, and which B vertices D that slot contains. For that slot's t'
tail vertices its count is r_D+2^t' u_D if t'>0, and is
r_D+u_D-e_D if t'=0. Indeed without a tail the total is H-E for
apex last, and U-W for all other placements. The remaining slot
contributes h_(B-D) because its tail is an ordinal TT. Binomial summation
therefore gives exactly
    P_close(2)=2 sum_D h_(B-D) (2^t r_D+3^t u_D-e_D).
This applies also to the empty D (r=1,u=e=0). For this concrete fork
the three summed differences, O minus P, are -4, 3, and 3. Thus the
closed quotient difference is 2(-4·2^t+3·3^t-3), negative at t=1,
positive at t=2, strictly increasing thereafter, and never zero at
integer t>=1. Same preclosure H,Q,U,E is not a continuation guarantee.
The invariant sufficient to repeat would have to retain subset-weighted
moments (or use a separate repair to match the new Q). The fork is not
a universal selector and proves no cofinite coverage.

The stronger version holds for every common, nonempty ordinal tail S.
If the apex is last, it contributes H(B)H(S). Otherwise its next vertex
must be a marked B vertex. Partition B into ordered paths (first optional,
second starts marked) and S into two optional paths; only the choices with
nonempty first B path ending marked and empty first S path are blocked.
This proves H_close=rH(S)+U P_S(2). Applying the slot split once more gives
    P_close(2)=2 sum_D h_(B-D) [P_S(2)r_D+P_S(3)u_D−H(S)e_D].
In the P_S(3) term the three S slots are the two sides of the apex's path
and the other optional path. Therefore the displayed fork's ΔQ/2 is
−4P_S(2)+3P_S(3)−3H(S)=−2H(S)+5c_2(S)+3c_3(S), where c_j counts
ordered j-path covers with all j paths nonempty. It equals −2 for singleton
S and is strictly positive for |S|≥2 since splitting any Hamilton path
after its first vertex gives an injection proving c_2(S)≥H(S). No common
nonempty S can restore Q; changing the source or matching two different
tails is genuinely needed for this fork.

In fact the individual summed triples (r,u,e) are (76,132,44) and
(80,129,41). Using two different fixed nonempty tails followed by common
TT_t padding cannot make the two quotients coincide identically in t,
either. The constant terms would require 44 H(S_O)=41 H(S_P). The
left is even and the right odd by Rédei parity. This does not preclude
matching at one height, or changing the tails with t.
