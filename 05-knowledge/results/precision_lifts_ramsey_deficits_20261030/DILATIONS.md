# Joint two-tail dilation: only the two-point fibre needs an extra test

**PROVED:** for 0<r<=1/8, coverage for every pair of positive attached speeds
on a fixed core C, together with coverage for every pair of odd attached
speeds on2C, is equivalent to pair coverage on every dilation dC, d>=1.
**FINITE-EXACT:** the required doubled-core banks are certified here.
**OPEN:** coverage for arbitrary primitive core shapes.

This consumes the exact single-tail gcd-fibre result from incoming commit
1da4abafb8,[Session120 LRC.md](../session120_major_frontier_interfaces_20261030/lrc/LRC.md).
Its multiple-tail warning is essential: different tails cannot choose
different fibre labels. We count bad labels in the common fibre instead.
The older [THM-4446, primitive ten-pack descent and dilation rays](../../../01-canon/theorems/THM-4446-lrc14-primitive-ten-pack-descent-and-dilation-rays.md)
already uses common-grid counts for three tails and closes its clock-three
dilation rays. The argument here is an exact two-tail exceptional-fibre
reduction with its own checked cores, not a claim to invent fibre lifting.

## D1. Automatic joint lifting except for two order-two orbits

For a base-safe time s, all t_j=(s+j)/d, 0<=j<d, preserve the core dC.
The phases of speed u form an order-h grid, h=d/gcd(d,u), with each phase
repeated d/h times. An open bad arc of length2r contains at most ceil(2rh)
grid points. Its bad fraction is therefore at most

    ceil(2rh)/h <= ceil(h/4)/h                             (1)

when r<=1/8. This bound is1/2 at h=2 and strictly below1/2 at every h>=3:
check3 directly; for h>3, ceil(h/4)<=h/4+3/4<h/2.

If both attached speeds have nontrivial orbits, not both of size2, their
bad fractions sum to less than1. Some **same** label j works for both.
This holds for every original safe s. It assumes neither independence of
phases nor disjoint bad sets.

If both orbit sizes equal2, write u=(d/2)a and v=(d/2)b with a,b odd.
Dividing the whole speed row by d/2 reduces the question exactly to2C
with attached odd a,b. This is the only nontrivial joint-fibre exception.
It is real: at r=1/14, d=2, u=1, v=7 and s=7/50, neither label works
for both tails although each tail has a successful label. This is a
pointwise fibre obstruction, not a claim that the whole core fails.

## D2. Size-one orbits and the equivalence

If both speeds are divisible by d, divide the row by d and use pair
coverage of C. If just u=d*a is divisible, choose s safe for C union{a}.
Pair coverage of C guarantees this by repeating a as both tails. Its
u-phase is independent of j; the other tail has bad fraction at most1/2,
so a successful common label exists. Otherwise apply D1 or the2C hypothesis.
This proves sufficiency. Necessity specializes the dilation to1 and2.

The theorem is an equivalence of coverage properties. Its witness remains
in the actual time fibre. This proof permits a finite search over j; it
does not claim a logarithmic joint-fibre decoder independent of d.

## D3. The checked banks

[joint_phase.json](joint_phase.json) certifies all positive pairs at1/14 for
four eleven-speed cores. [dilation_pair.json](dilation_pair.json) certifies
all positive pairs on their doubles, stronger than the needed odd-only
case. Its169,029 finite queries pass the same independent one-tail reference
checks and analytic tails. The grid control checks279,070 cases with
2<=d<=44,1<=u<=v<=2d. The separate(2,1,7) hostile is outside that speed
window and is checked explicitly at1/14.

Consequently, for each of the consecutive, Fibonacci-adjacent and prime
cores in[LRC.md](LRC.md), every positive dilation d and every positive pair
u,v satisfies the1/14 clearance condition. The26-times consecutive core
is now redundant: it is in the first core's dilation family. Its doubled
bank completed as a positive control, not as a new primitive family.

The reusable step is precise: retain the two orbit orders, isolate the
small fibre where the union bound loses strictness, and solve that fibre
with its shared label intact. The primitive-core coverage question remains.

## D4. A constructive lcm flag for several tails

There is a stronger pointwise statement for **unequal** orbit orders:
if h_u,h_v>1 and h_u!=h_v, every base phase has a common successful lift
for both tails at every radius0<r<=1/4. It needs only two modular inverses,
without enumerating the d labels. Equal orders are not covered by this rule.

More generally, order k tails so that, with L_0=1 and

    h_i=d/gcd(d,v_i), L_i=lcm(L_(i-1),h_i),
    H_i=L_i/L_(i-1)=h_i/gcd(h_i,L_(i-1)),

every H_i>=2. Then every base-safe phase s has a native lift safe for all
tails at radius1/4. A finer sufficient radius is
min_i (H_i-1)/(2H_i). This is a sufficient arithmetic flag, not a necessary
condition for joint coverage.

**Proof and decoder.** After i-1 stages retain j_0 modulo L_(i-1).
Changing j by L_(i-1)*t fixes every previously chosen tail phase modulo1.
The new tail has an order-H_i grid in this remaining fibre. Its maximum
norm is at least(H_i-1)/(2H_i), by the inherited single-tail formula.
Choose its nearest grid point to1/2 by rational rounding and a modular
inverse. This fixes j modulo L_i while preserving all earlier phases.
Induction gives0<=j<L_k<=d and t=(s+j)/d with the actual original time.
For two distinct orders at least one does not divide the other; put that
order second to ensure H_2>=2. QED.

The first obstruction to iterating this trick is a redundant lcm:
orders2,3,6 admit no ordering with every stage strict. This can reflect an
actual pointwise failure: d=6, tails(3,2,1),s=1/10 have no common safe
label at radius1/4. The core C={3} is safe at this s, so it is a genuine
core fibre. This does not say the complete speed row has no other safe
time. The equal-order hostile in D1 fails even at1/14.

[fibre_flag.py](fibre_flag.py) checks the native decoder directly against
all tail norms in its declared finite universes, tests both hostile fibres,
and constructs a lift with d=2^150*3^100*5^80. The algorithm has no search
loop depending on the numerical value of d: one gcd, rounding and modular
inverse per tail suffice. Input bit lengths still control bit complexity.

## D5. Incoming 78-core bank extends to every dilation

Commit373661611e's
[L4, complete shared-bank phase classification](../clique_palette_joint_ports_20261030/LRC.md)
uses our two-tail conductor and independently certifies every11-subset of
[1,13], checking6,397 finite pairs by closed-tooth and literal-comb paths.
It supersedes fixed-ratio enumeration for these78 cores.

We apply D2 to this incoming result. The additional necessary bank is every
odd pair on each doubled core; [dilation78.py](dilation78.py) proves the
stronger all-positive-pair statement on every doubled core. All21,945
finite queries pass the Euclidean pair compiler and the independently
organized one-tail check after forming C+u. The largest both-large cutoff
is319; every remaining single-small-tail branch has positive width and an
explicit finite cutoff. No extrapolation from the finite queries is used:
the two-tail conductor pays their infinite tails.

Consequently every eleven-subset C of[1,13], every positive integer d, and
every positive pair u,v has a common time at clearance1/14 for

    dC union {u,v}.

Together with the Fibonacci-adjacent and prime cores in D3, this gives80
distinct primitive core shapes with full dilation and two-tail closure.
The base78 assertion is credited to the incoming packet. This result adds
its dilation closure; arbitrary primitive cores remain outside the theorem.
