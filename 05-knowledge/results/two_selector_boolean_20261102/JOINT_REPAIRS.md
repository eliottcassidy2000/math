# Three-edge repairs expose a genuinely third-order source law

**Inherited FINITE-EXACT escape theorem, independently replayed two-edge
test; new FINITE-EXACT joint-law witnesses and elementary general identities.**

## What the requested two-edge test establishes

Incoming [prime inverse MOVES](../prime_inverse_ports_20261102/MOVES.md)
already resolved the768 six-vertex strict local traps for target17. We
replayed all768 sources against every pair of edge flips:480 have a two-edge
target;288 do not. The four ordinary isomorphism orbits have representatives
72,140,178,314 and masses240,144,240,144. The two144-element orbits require
three flips, with a sharp score-distance obstruction. These counts and the
three-edge packets are inherited, not a new discovery here.

All four representatives have a common-path-preserving packet whose sequential
counts are15 -> 13 -> 23 -> 17 and whose minimum possible peak squared error
is36. A multi-edge jump must not be called a monotone single-edge descent.

## The entire packet cube is physically realized

Fix one inherited three-edge packet with distinct changed edges e1,e2,e3.
For u in{0,1}^3 let T_u flip exactly the selected edges of T_000. Every T_u is
an actual labelled tournament. The inherited common Hamilton path avoids all
three changed edges, so it remains legal at all eight corners, not just at the
four sequential packet states. We verified every corner independently by
endpoint subset DP and by the full path-subcube count table.

Order the corners by u=0,...,7 with bit0 selecting e1. The H values are

| Source | H at the eight corners |
|---|---|
|72|15,13,27,23,13,9,23,17|
|140|15,13,33,23,13,9,29,17|
|178|15,13,31,23,13,9,29,17|
|314|15,13,33,23,11,9,29,17|

In each cube,111 is the only target17 corner.

## All edge-pair marginals agree, but target support differs

Put uniform probability on the four even-parity corners to obtain mu_even,
and on the four odd-parity corners to obtain mu_odd. Every proper marginal
on the three varying bits is uniform for both laws. All other tournament
edges are fixed. Therefore **every one-edge and every two-edge joint
orientation distribution of the full tournament agrees** between the laws.

Nevertheless,

`Pr_mu_even[H=17]=0,    Pr_mu_odd[H=17]=1/4.`

For E=(H-17)^2 the exact expectations are

| Source | E under even law | E under odd law | Delta_123 E |
|---|---:|---:|---:|
|72|35|33|-8|
|140|62|72|40|
|178|62|57|-20|
|314|62|77|60|

Here Delta_123 E=sum_u (-1)^(3-|u|) E(u), so the difference of the two
expectations is -Delta_123 E/4. These are positive laws on actual common-path
sources, not formal negative weights or fitted moments. The control rules
out a pairwise-edge distributional carrier for the target-support query;
it does not say that the complete deterministic edge table is insufficient.

The source72 cube gives an even sharper control. Both parity laws have mean
H=35/2. Mix each with the SAME atom T_000 of H15, giving weight4/5 to the
parity law and1/5 to that atom. The new laws still have every pair marginal
equal, and **both now have mean exactly17**. In corner order0,...,7 their
weights, divided by5, are

```text
nu0: 2,0,0,1,0,1,1,0,
nu1: 1,1,1,0,1,0,0,1.
```

Yet their target-hit probabilities are0 and1/5 and their squared-error
expectations are144/5 and136/5. Thus even the exact target mean plus all
pairwise edge laws misses whether the support contains a target tournament.
All source weights are nonnegative and all supported orientations are actual.

More generally, the two parity laws on a d-cube agree on all proper
marginals and differ on an observable f by its top mixed difference divided
by2^(d-1), with the appropriate sign. The assertion follows by summing over
one omitted coordinate. This is the classical Boolean finite-difference
identity, applied here to actual tournament attachments.

The older [Liouville–Möbius eleven note](../liouville_mobius_eleven_20261027/ELEVEN.md)
used precisely this proper-face versus full-corner mechanism on squarefree
divisor cubes. A literal indexing map here sends u to2^u1 3^u2 11^u3 and to
T_u. Then mu(d)=lambda(d)=(-1)^|u| on this cube. The alternating divisor sum
is exactly the displayed response difference. No primality, Collatz, or LRC
theorem is transported; the preserved operation is the signed mixed difference.

## Correlated choice removes the expected-energy barrier, at a known cost

For each packet, let one Bernoulli variable Z choose T_000 when Z=0 and
T_111 when Z=1, with Pr(Z=1)=t. This actual two-source law has

`E[(H-17)^2]=4(1-t),    E[H]=15+2t.`

Its expected error decreases to zero. In contrast, independently flipping
each of the three edges with probability t gives these exact polynomials:

```text
source72:   4+120t-116t^2-8t^3,
source140:  4+276t-320t^2+40t^3,
source178:  4+216t-200t^2-20t^3,
source314:  4+296t-360t^2+60t^3.
```

Every derivative at zero is positive. The correlated law uses one common
selector for a **known valid packet**, retaining three-edge dependence. It
does not discover the target, prove a universal inverse, or reduce the peak
error36 when physical execution must flip edges one at a time.

This is the relevant information-saving device: a certified whole-packet
instruction can share one choice, while independent coordinate updates erase
the dependence that makes the packet useful. The joint-source distinction
also explains why the two parity laws above cannot be replaced by matching
marginal statistics.

## Birth, old dead states, and future operations

The inherited [97-state square25-to26 interface](../eleven_wedge_honeycomb_20261102/NATIVE_26.md)
keeps all old sources. Exactly9,184 active G26 orientations arise from old
dead orientations; they contribute58 additional pre-birth classes. The full
interface is97=38+58+1, not the old33-state fixed-order observer.

There are two different losses here. The birth example loses deterministic
edge information needed by a later path; the parity-law example loses joint
probability information although every pair marginal survives. The repair
in each case is operation-specific: refine the actual source before composing.
The numbers26,97,9184 are not a statement about the26th tournament class.
G26 is the square-sum graph on the literal integer interval1,...,26.
