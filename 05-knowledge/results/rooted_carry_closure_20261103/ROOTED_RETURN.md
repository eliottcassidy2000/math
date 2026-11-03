# RC: parity-selected H497 entry and exact clocks

**PROVED all k>=1.** The rooted profile z=(h,s,q,F) records H, paths starting
at the root, H after root deletion, and ordered two-bin covers whose first
nonempty path ends at the root; the second bin may be empty.

## RC1. A neutral reset changes the primitive return clock

Use chart(P,R)=(3,4) in
[APEX mechanisms](../apex_wedge_constructor_20261103/mechanisms.py).
There are seven base vertices and physical apexes a=7,b=8. Reverse only
a->b with the inherited two-native-move reset and one helper birth/deletion.
Every old name survives. The full eighteen-polynomial observer gives these
profiles, independently recounted by endpoint subset DP:

|Orientation|Root|h|s|q|F|delta=F-h-q+s|
|---|---:|---:|---:|---:|---:|---:|
|plus|a|497|59|69|549|42|
|plus|b|497|55|61|597|94|
|minus|a|497|58|69|622|114|
|minus|b|497|56|61|528|26|

A root birth adds x with r->x->T-r. R retains r; N selects x. The exact
inherited laws are

    R(h,s,q,F)=(h+F-s,q,q,2F-q),
    N(h,s,q,F)=(h+F-s,h-s,h,h+F).

[Rooted sealing R1/R3](../two_port_return_induction_20261103/ROOTED_RESETS.md)
gives R^k N R=(b,t,t,b). Put a0=h-s,d=F-q. Then

    b=a0+3q+(3*2^k-1)d,    t=a0+q+(2^k-1)d,
    b-3t=2delta,           gcd(b,t)=gcd(delta,2^k d+q).     (RC1.1)

Odd q removes the factor2 in the gcd. Substitution gives:

|Repaired root|b|t|gcd(b,t)|
|---|---|---|---|
|a|1659*2^k+93|553*2^k-45|19 if k=5 mod18, otherwise1|
|b|1401*2^k+157|467*2^k+35|13 if k=8 mod12, otherwise1|

For a, 553*2^k+69 never vanishes mod3. Mod19 its zero condition is
2^(k+1)=7, equivalent to k=5 mod18. For b the mod13 condition is2^k=9,
equivalent to k=8 mod12. In contrast plus a always retains3, and also7
exactly when k=1 mod3. This identifies the mechanism of the repair.

**Select a at even k, b at odd k.** Every selected seal is primitive.
Each fixed port fails infinitely often. Two parity states suffice, and
one fixed choice does not, for this particular selection problem.

## RC2. One C5 always enters the repeatable primitive controller

The actual-source cycles of
[NATIVE_SELECTION](../source_selection_closure_20261103/NATIVE_SELECTION.md)
act on balanced pairs by

    C3(b,t)=(61b-2t,9b),    C5(b,t)=(421b-2t,33b).          (RC2.1)

Their literal construction works on any supplied balanced source. N at
the root creates r; add j children at r, keeping r and the newest z.
Move this physical z to the oldest position with j-1 adjacent module
swaps, each two native moves and one helper birth/deletion. N at z, then
R at its newborn, returns to a balanced root. Costs are j+3 permanent
births,2(j-1) native moves and j-1 helpers. H497 ancestry from C3 is not
assumed. Every path guard and name is recorded in returns.json.

For coprime odd b,t, a C5 common divisor can only be3 or11. Its eleven
guard is3b-2t !=0 mod11. For selected a at even k this is6-2^k mod11,
zero only at k=9 mod10. For selected b at odd k it is2*2^k+5, zero only
at k=8 mod10. Parity excludes both. Mod3 the a seed has(b,t)=(0,1)
and the b seed has(1,0); both give3|t' and3 does not divide b'. Thus C5
is primitive and t'/b'=0 mod11. This enters the inherited eight-state
controller: C3 is always legal, and C5 is legal except at ratio7, resetting
the ratio to0. Hence the source family has repeatable primitive returns.

Including the initial neutral reset, entry costs ten native moves, five
helper births/deletions and k+10 permanent births; final order is k+19.
Independent native code replays the complete source program. Full
eighteen-polynomial transport matches its endpoint by a virtual history
that births z first, retaining it while the other children are born. The
recorded name permutation identifies that graph with the actual native
endpoint; it is an observer computation, not a free physical relabeling.

At k=1 an optional C3 return gives prime H178561 on18 vertices, independently
counted by subset DP. For every even k, the C5 oldest exit is

    E_k=226b-t=7(53483*2^k+3009).                          (RC2.2)

At k=6 this is23981447=7*3425921 on24 vertices. Both claimed primes
are checked by trial division. The larger path count uses the proved
marked block transfer and literal graph replay, not a24-vertex DP census.
No prime infinitude is asserted.

## RC3. Scalar and full golden recurrences

For selected u=b or t, u[k+4]=5u[k+2]-4u[k]. With generating functions
sum_(k>=1)u[k]z^(k-1),

    B(z)=(2959+6729z-3430z^2-7008z^3)/(1-5z^2+4z^4),
    T(z)=(969+2167z-1074z^2-2032z^3)/(1-5z^2+4z^4).

The inherited golden observer is(h,X,Y), with X=(S(2)-E(1))/2 and
Y=(S(2)+E(1))/2, where S/E are ordered-bin starts/ends at the root.
For the selected H497 seal:

|k|h|X|Y-X|
|---|---|---|---:|
|even|1659*2^k+93|12920*3^k+96*2^k-703|14|
|odd|1401*2^k+157|10820*3^k+1832*2^k-345|12|

After C5:

|k|h|X|Y-X|
|---|---|---|---:|
|even|697333*2^k+39243|26944660*3^k+4563888*2^k-1046432|224|
|odd|588887*2^k+66027|22565110*3^k+6123264*2^k-509726|192|

**All-order proof certificate.** Cap the cover states at two unnamed
terminal endpoints, giving46 states. Under indegree-one births the only
consumed incoming endpoint is the named root, so unnamed-end count cannot
decrease. The discarded states can never affect this cap. All golden
observations are linear functions of it. For each initial root and each
suffix(NR, or NR followed by virtual C5), the exact rational operator in
phase_recurrence.py proves that the observation row annihilates

    R(R-I)(R-2I)(R-3I).

Thus the formulas hold for every k>=1; three initial evaluations determine
the coefficients. They are not extrapolations from a fitted prefix. Full
polynomial transport independently checks them at k=1,...,18. The46-state
cap is not claimed minimal. Parity merging gives

    u[k+6]=14u[k+4]-49u[k+2]+36u[k],

with generating-function denominator(1-z^2)(1-4z^2)(1-9z^2).

In the inherited convention theta^2=4theta+2 over F7,

    alpha_G=4(h+2X-Y)+Y theta,    charge=h-Y mod7.

C5 output phases(constant,theta coefficient,charge) for k=1,...,6 are

    (3,6,0), (4,0,1), (2,2,6), (5,2,6), (2,1,1), (0,6,2).

The formulas give exact period6; all six values differ. The receipt retains
all five Q4 presentations individually. q2 multiplies alpha_G by theta
and preserves charge; Frobenius sends theta to4-theta. This is not the
balanced-pair beta coordinate of PC1.

## The remaining attachment obligation

Incoming [READY sections4–6](../two_entry_funnels_20261103/READY.md)
classifies the five R-stable ready families and proves that two successive
root births from a strong source of order>=4 leave no1/2-ready root,
even after arbitrary re-marking. Every final seal/return here ends with
that reset. Thus its endpoints are not ready sources. This does not
exclude a genuine source switch or different attachment.

Incoming E_k has H=3*2^k-1 and a converse choice for both offsets; its
golden phases match iff k=2 mod6, then(alpha_G,charge)=(1,4). Our same
clock position gives(4,1). Equal clock moduli do not identify the phases.

A complete six-command search from each H497 orientation, initially
retaining(7,8), has seven continuation states at H<=1495 and depth at
most1. Attained counts are{497,987,1039} before the reset and
{497,969,1061} afterwards. Neither1493 nor1495 occurs. Strict H growth
makes the cutoff complete in that language. It excludes neither native
source changes nor other port policies or two-polarity births.
