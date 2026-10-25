# Clock-correct reset bounds and a ternary role ledger

**PROVED elementary inequalities and identities / FINITE-EXACT controls /
OPEN global amortization, 2026-09-26.** No novelty claim. This companion to
[the product inequality](reset_20260926_inequality.md) uses the original
[precision controller](nextforest_20260926_precision.md), including its
extra arithmetic step before resetting.

## 1. What the clock changes

Write U(n)=oddpart(3n+1). A retained state is

    y=q*2^r+u,   q,u positive odd, r>=1,
    K=q*u*2^r=u(y-u).

Set a=v2(3u+1), v=U(u). At a collision a=r, put
b=v2(3q+v)>=1. The actual division exponent is m=r+b>=2,
and x=U(y). The inherited P5 protocol next takes z=U(x), if x>1,
then splits z canonically. Thus this reset costs TWO U steps from y.
If x=1 it stops. Let K_new=0 at terminal1; otherwise any positive
integer split of z has K_new<=floor(z^2/4)=(z^2-1)/4.

Call the old split balanced when y/4<=u<=3y/4. Since a collision has
y=1 mod4, integrality gives

    K >= 3(y-1)(y+3)/16.                            (C1)

Indeed both summands are at least ceil(y/4)=(y+3)/4, and their
product is minimized at an endpoint of that interval.

**P5 reset inequality.** Every such balanced collision satisfies

    K_new < (27/16) K.                              (C2)

If x>1 and the extra step has d=v2(3x+1)>=2, then

    K_new < (27/64) K.                              (C3)

Both constants are sharp over all admissible balanced triples.
Sharpness over the smaller class of triples reached from a specified
canonical initialization is not asserted.

Proof. Since m>=2 and d>=1,

    x <= (3y+1)/4,       z <= (9y+7)/8.

For y>=9, (C1) and the maximum product bound give

    K_new/K <= (81y^2+126y-15)/(48(y-1)(y+3)) <27/16.

The difference after clearing positive denominators is 36y-228>0.
The only smaller possible source is y=5, which collides to1.
For d>=2, use z<=(9y+7)/16 instead. Then

    K_new/K <= (81y^2+126y-207)/(192(y-1)(y+3)) <27/64;

the cleared difference is 36(y-1)>0. Terminal cases have zero product.
These proofs also allow noncanonical fresh splits, so no optimization
assumption is hidden in the reset.

In contrast, an IMMEDIATE split of x has sharp factor below3/4 in
the companion note. The constants differ because an actual extra U
step can increase the endpoint. They are not competing estimates of
the same operation.

## 2. Sharpness and an actual canonical-orbit hostile

For k>=0 set

    y=(2^(6k+8)-31)/9,
    (q,r,u)=(3(y-1)/8, 1, (y+3)/4).

These are positive odd registers with a=r=1, actual collision exponent
m=2, extra exponent d=1, and z=2^(6k+5)-3. A canonical split of z gives

    K_new/K=(81y^2+126y-527)/(48(y-1)(y+3)) ->27/16.

The first member is (9,1,7), y=25, with x=19,z=29 and ratio104/63>1.
It is reached on an actual canonically initialized phase:

    (1,7,49)=177 -> (3,5,37)=133 -> (9,1,7)=25
       -> collision19 -> P5 reset29=(1,4,13).

Therefore even balanced P5 resets can increase K on the prescribed
chronological controller. Replacing27/16 by1 is false.

For k>=1 set z=2^(6k+1)-1, y=(16z-7)/9 and

    (q,r,u)=((3y-11)/8,1,(y+11)/4).

These are balanced admissible triples, with m=2,d=2 and U^2(y)=z.
Their canonical reset products satisfy

    K_new/K=4(z^2-1)/((y+11)(3y-11)) ->27/64.

The two families establish sharpness of (C2),(C3) on their stated
domain. They do not establish the density of either type along one orbit.

## 3. Precision has an exact ledger, not yet a drift bound

For a noncollision write r' for the successor precision. The controller
has r'=|r-a| and m=min(r,a), hence

    2m=r+a-r'.                                      (C4)

At a collision set terminal precision to zero; since a=r and m=r+b,

    2m=r+a+2b.                                      (C5)

For J steps from initial precision r0 through the first collision,
with A the sum of ACTUAL division exponents,

    2A=r0+sum_(j=0)^(J-1) a_j+2b.                  (C6)

For an open prefix, the corresponding identity is
2A=r0+sum a_j-r_J. These telescope through swaps without losing the
source. A subsequent P5 extra step adds its own exponent d; a fresh
precision R is then a new register, not a free positive term in (C6).

Equation (C6) is valuable bookkeeping, but it is not an inequality
controlling sum a_j or regenerated R. Treating it as such would assume
the missing amortization. The product law in the companion note and
(C2),(C3) identify exactly which events need a compensating charge.

## 4. The genuine ternary information: role and run length

After EVERY noncollision, exactly one of the two odd registers q,u is
divisible by3. This follows because U(u) is always a3-adic unit:

| branch | new coefficient | new core | nonzero3-adic valuation |
|---|---|---|---|
| consume |3q|U(u)|v3(q_new)=v3(q)+1|
| swap |U(u)|3q|v3(u_new)=v3(q)+1|

Thus a role bit records which branch just happened. Starting with
canonical q=1, v3(q) equals the length of the current consecutive
consume run. A swap transfers that run length plus one into v3(u_new)
and returns v3(q_new) to zero. This is an exact ternary run clock with
unbounded depth, naturally related to a finite coloured marker PLUS
an integer register. It does not determine the next binary valuation.

The older [ternary sibling clock](ternary_digits_20260925.md), section3,
supplies an exact hostile to replacing that register by finitely many
ternary digits. Put u_k=(4^k-1)/3. For any fixed D>=1, the subsequence
k=4+j*3^D has constant u_k modulo3^D and v3(u_k)=0, but

    (1,1,u_k) -> (1,2k-1,3).

The new precision is unbounded. The congruence follows from
v3(4^t-1)=1+v3(t); v3(u_k)=v3(k)=0. These are admissible noncanonical
states, not a claim that every displayed triple arises from canonical
initialization. Their actual values u_k+2 have the inherited defectE=12
for k>=3. Hence even a small absolute defect and arbitrarily many fixed
ternary digits do not bound a single swap's regenerated precision.

The comparison with [the three-unit colouring](reset_20260926_colours.md)
is consequently precise: an extra marker can encode a missing role;
the3-adic valuation supplies a depth coordinate which a finite marker
alone cannot replace. Neither statement provides a well-founded rank.

## 5. Controls and scope

[The exact script](../../04-computation/experiments/reset_20260926_clock.py)
checks49152 triples (q odd1..63,u odd1..255,r1..12), including4096
collisions and2267 balanced cases. It checks both sharp families,
all2047 canonically initialized phases from odd3..4095 with cap512,
and fixed ternary guards at depths1..6. The finite phases all collide;
the proof does not assume that an arbitrary phase must do so.

Source: exact retained triples. Target: split product plus a precision
and role ledger. Preserved: actual integer, actual step count, and
valuation budget. Lost by the product alone: split imbalance and
ancestral arithmetic. Required sidecars: the exact registers and reset
clock. Hostiles:177->133->25->19->29, the sharp families, and the
fixed-ternary precision recreation family. The remaining task is a
source-preserving inequality paying for complete phases, not another
local assertion of balance or finite-state expressibility.
