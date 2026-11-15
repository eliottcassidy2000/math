# Recenter a changed reader by dilating both address axes

Two older, tangential sources suggest this move. The [Lonely Runner phase-branch pullback](../../../07-reflections/lrc14-phase-branch-color-pullback-codex-20260629.md) insists on pulling a numerically plausible gate through the correct discrete branch before calling it an obstruction. The [residue/phase/incidence synthesis](../../../07-reflections/residue-phase-incidence-synthesis.md), itself a random historical reading, distinguishes the residue surviving a projection from the phase surviving a decomposition. Here pull back the old Hamiltonian count to two coordinate charts of different aspect ratios, keeping both ordinary and bit-reversed phases. Fixing one tape and repainting only was the too narrow chart.

Use the physical two-tape formula from [the unit-slide core](../chanpak_reservoir_slide_20261115/UNIT_SLIDE.md). Its fixed eight-vertex underlying core is `(4,1,130,7,143,159,63,75)`; its two guards are 1 and 0. At height zero *every* mask m on these eight vertices gives
```
F_m(x,y)=B_m+T_m x+S_m y+5760xy+45[p(x)+p(y)]
p(0)=0, p(2v)=p(v)+v+(v mod 2), p(2v+1)=p(v).
```
The unmarked spine counts, hence 5760 and 45, do not change with the mask, and the guards kill both signed channels and give strongness. The three coefficients are obtained by literal Hamiltonian subset-endpoint DP at (0,0), (1,0), (0,1) since p(1)=0. In particular m=35 gives (2855,4034,3914), and m=195 gives (2793,3914,3944), in (B,T,S) order.

Put M=16384 and b0=11158. From the [row-fold analysis](../chanpak_modular_kernels_20261115/ROW_FOLDS.md), the row x=u has duplicate value G(u)=F_35(u,Mu+b0)=F_35(u,Mu+b0+45). Each of the following is a literal changed-reader return:
```
F_195(4096u+r,4u+2) − G(u) = K_r    for every integer u>=0

r       2840     2842     2845     2847
K_r   -11126   -12118   -13786   -14958
```
Thus a change of mask can preserve absolute height up to a *constant*, despite the earlier linear displacement bound when the addresses themselves were fixed. The product of the new scale factors is exactly M. The derivative and low-bit phase then cancel at the next order as well. All four shifts are even, with gcd 2. This palette alone is finite and all four go the same direction; no independent persistent composition follows, and it does not establish coverage.

## Exact theorem within the common-high-word dilation model

Fix a power of two P and a low word a in [0,P). Evaluate the old pair at x=P v+a, y=M(P v+a)+b0; let v range over all nonnegative integers. Pick dyadic factors d,e>=1 with de=M and any mask other than 35, and try the new pair
```
x'=dP v+dA+r,  y'=eP v+eC+s,
0<=A,C<P, 0<=r<d, 0<=s<e.
```
If the count difference stays bounded as v tends to infinity, then d=4096, e=4, m=195, A=C=a, s=2, and r is one of the four values above. Thus adding arbitrary common high-word depth does not improve this finite palette. The theorem is not about affine addresses with a shifted high prefix, nondyadic scalings, growing mask sets, or a different underlying core.

**Digit proof.** Write rev_q for reversal of exactly log2(q) bits, with rev_1(0)=0. Iterating the recurrence for a fixed low word b gives
```
p(qv+b)=p(v)+(q−1−rev_q(b))v+p(b)+ε_b(v mod 2),
```
where ε_b is 1 iff the leading low bit is 0, and ε_0=0 at q=1. Only the first digit appended can see the parity of v. For the old and new pairs, the 2p(v) terms and quadratic leading product terms cancel. Boundedness is therefore equivalent to zero coefficient of v, with at most a two-valued parity-dependent remainder.

Let R=rev_M(b0), L=A+C−2a, and J=(2rev_P(a)−rev_P(A)−rev_P(C))/P. Splitting the bit reversals at d,e,M and comparing the coefficient of v shows
```
0 = (T_m d+S_m e−T0−S0 M)/45
   + d+e−M−1 + R−rev_d(r)−rev_e(s)
   + 128[ M L+ d s+ e r−b0 ] + J.
```
The unreduced equation forces 45 to divide T_m d+S_m e−T0−S0 M (P is a unit modulo 45). It then forces J to be an integer and hence in {-1,0,1}. Also |L|<=4, using 2880<=S_m,T_m<=5670 from literal DP, d+e<=M+1, and -M<d s+e r−b0<2M. More explicitly, the first term multiplied back by 45 is less than 5760M in absolute value, and the second and third groups outside 128[.] including J contribute less than M to the equivalent equation after division by 128. Thus |L|<4 already suffices; the certificate uses the looser bound.

Now all variables except a,A,C,P have bounded, fixed numerical ranges. The exact integer coefficient certificate in replay_recenter.py checks them, using rev_e(s) in [0,e−1] to restrict candidates s before testing equality. It leaves only the four displayed low motifs, all with L=J=0. The two exact midpoint identities
```
A+C=2a,   rev_P(A)+rev_P(C)=2rev_P(a)
```
force A=C=a: first A,C have the same parity by the ordinary sum. If their parity differs from a's, both reversed words have their leading bit on the opposite side of the half interval from rev_P(a), so their midpoint cannot be rev_P(a). All parities therefore agree, and deleting their common last bit reduces both midpoint identities to width P/2. Induct from P=1. The new addresses now have the form d x+r,e x+2; substituting the digit formula with high prefix x proves the constants in the table, including cancellation of both parity terms. ∎

## Why moving only one axis cannot do this

Keep x, change the mask with S_m != 3914, and allow any new right address z. If z<Mx or z>=M(x+1), the leading-product loss or gain already prevents equality with G(x). For the upper side bound the old surcharge by 45(Mx+b0), the new surcharge below by 0, and use the smallest S_m=2880 at z=M(x+1). For the lower side use the largest S_m=5670 and the largest surcharge 45Mx at z=Mx; the resulting upper bound is below G(x). These inequalities are strict for all x>=0 from the stated coefficient ranges and p(t) in [0,t].

In the remaining 14-bit chart z=Mx+b the digit identity says the difference is exactly a slope times x plus an intercept depending only on x mod 2. The exhaustive coefficient certificate gives minimum absolute slope 14 and solves every affine exact-equality equation with its parity. Its only exact return with S_m changed is x=0, m=233, b=8886: there is none on any positive row. The best-slope chart is m=177, b=9347, where the exact discrepancy is 14x−1184510. It gives +2 at x=84608, but this near-return is isolated rather than repeatable. The two-axis dilation above is what removes the entire linear drift.

## Construction debt

The paired projections show that even exact bounded recentering need not produce a unit: residue reachability, tangency, recentering, and independently iterable response are separate. A wider return can change the whole underlying core or splice a translated high word across the old fold. For the latter, one must retain the digit-dependent jump p(v+c)−p(v), which the present common-prefix certificate deliberately factors out.

Run `python replay_recenter.py` in this directory. It rebuilds all 256 profile triples from literal tournaments, certifies the fixed-row affine reduction, certifies the all-depth low motif table, and checks the exact return at both parities and extreme address sizes.
