# A fourth deletion in the two-pair AND case

**Status: PROVED local structural lemma, conditional on the explicitly stated
legal-state invariants.** This improves the AP2 row of the pinned provisional
elimination proof. It is not independent certification of that proof's global
coverage, affine-disperser construction, or a new circuit lower-bound record.
No literature-priority claim is made.

Source: [provisional elimination manuscript, Sections 1–7](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/4a61f024a4993592d97b49ca7d00979029b7fece/submission/elimination.md).

## 1. Setting, including sharing

The acyclic part of a full-B2 state has one output. No acyclic gate feeds the
auxiliary affine system. Every gate except the output has at least one child,
so following children eventually reaches the output. All gates depend on two
distinct nonconstant parents after the preliminary simplifications.

Let G be an AND-type gate with paired parameter parents p,w. Their distinct
partners are r,h. Each of p,w,r,h has exactly one child, and paired parameters
occur in no affine equation. Thus the unique children H_r of r and H_h of h
are gates, neither is G; they may coincide. The two marked pairs are distinct.
G has a child and is not the output.

Perform the source AP2 restrictions: fix r and h safely, and then fix w to the
value killing G. These restrictions preserve the product-domain semantics,
remove three parameters and two pairs, and decrease safe-core dimension a by
one. Assume the output restricted to the resulting allowed domain is
**nonconstant**. In the proposed disperser application this follows from the
source hypotheses a>=d+2 and a'=a-1>d, not from a finite graph experiment.

A gate with a constant pin is removable as a constant or a literal of its
other pin. Count distinct original gates, including a gate already removed
while processing another restriction. Subsequent simplifications only help.

## 2. Theorem: at least four distinct gates disappear

The source counts three; the following argument forces a fourth.

### Case A: H_r and H_h are distinct

The three distinct seed gates G,H_r,H_h are already removable: G is constant,
and the other two have their r or h pin fixed. Follow a directed child path
from G to the output in the original acyclic graph.

As long as this path stays among the three seed gates, every encountered gate
is constant after the restriction. Indeed G is constant. Each later seed is
H_r or H_h; one of its pins is its fixed paired parameter, and the other pin
on the path is a constant predecessor. Both its pins are therefore fixed.

If the entire path stayed among the three seeds, the output would be constant,
contrary to the hypothesis. Therefore the path first exits the seed set at a
fourth distinct gate K. The preceding gate on the path is constant, so K has
a constant pin and is removable. We have four distinct deletions.

This argument counts a shared descendant only once. It does not require the
children of G to be disjoint from the partner children: those overlaps are
exactly why the constant path argument is needed.

### Case B: H_r=H_h=H

The two pins of H are r,h, so H becomes constant. G also becomes constant.
Neither gate reads the other, since their pins are the specified parameters.
Both have children: if either were the output, its constant value would already
contradict the nonconstant-output hypothesis.

If the union of their distinct child sets contains at least two gates, those
two gates are outside {G,H}, have constant pins, and are removable. Along with
G,H this gives four deletions.

Otherwise that union is a singleton {K}. Both G and H feed K; its two distinct
pins are therefore constant, and K becomes constant. It cannot be the output.
Choose a child L of K. Acyclicity makes L distinct from G,H,K, and its constant
K pin makes it removable. Again there are four deletions.

These cases exhaust all shared-child possibilities. The proof uses only the
binary pin bound, acyclicity, live output ancestry and post-restriction output
nonconstancy. It imposes no formula restriction and does not assume disjoint
subcircuits.

## 3. Revised row and its strictly conditional consequence

The original parameter and usage accounting is unchanged:

    Delta m = -3, Delta q = -2, Delta a = -1,
    Delta I <= 0, Delta s <= -4.

Thus for W=g*s+i*I+p*q, with nonnegative i,p and g>0, AP2 satisfies

    Delta W <= -4g-2p.

At the source weights (g,i,p)=(1,5,3), this yields -10, a surplus of one over
its dimension charge nine. At weights (1,4,2), it yields -8. Combining the
revised row with **all the other source rows and their global invariants**, if
valid, gives the conditional telescope |C|>=4n-8(d+1). This paper proves only
the strengthened AP2 row; the complete revised potential frontier and its
remaining obligations are recorded separately in this session.

## 4. Sharp legal-state witness: four cannot become five unconditionally

Let

    f(p,w,r,z,h,t,a,b) = pw + rz + ht + ab

with binary addition. Use identity affine parameter map and mark (p,r) and
(w,h), each excluding corner (0,0). Set d=4. This quadratic has nondegenerate
polar form on eight variables, pairing p with w, r with z, h with t, and a
with b. Any constant affine flat has a totally isotropic direction of
dimension at most four; a four-dimensional zero flat exists. Hence f is
nonconstant on every affine flat of dimension greater than d=4.

The state has m=8,q=2,a=6=d+2 and no affine equations. Its seven gates are

    G=pw, H_r=rz, H_h=ht, K=G+H_r,
    O=K+H_h, A=ab, F=O+A (output).

Every gate is live and nondegenerate. Every paired parameter is isolated,
and partners have different children. There is no unused paired parameter,
paired high-fanout parameter, paired affine-equation occurrence, affine gate
with affine parents, or other preliminary deletion. G can be selected as the
first minimal AND. (A is another available minimal AND; this witness is about
the stated legal AP2 choice, not every strategy.)

Fix r=h=1 safely and w=0. Exactly the four gates G,H_r,H_h,K disappear. The
remaining circuit is

    O=z+t, A=ab, F=O+A.

All three remaining gates are nondegenerate and live. The output depends on
z,t,a,b, so no B2 circuit for it can use fewer than three gates. The retained
parameter p is unused, as the source invariant permits. The exact changes are

    (Delta s,Delta I,Delta q,Delta a)=(-4,0,-2,-1).

The initial seven-gate circuit is also gate-minimal, even for agreement only on
the marked product domain. Every one of the eight parameters has two legal
same-fiber assignments differing only in that bit and having different output.
For paired-bit flips, choose the partner equal to its safe value one; its
quadratic mate can then be set to one. For z,t choose r,h equal to one, and
for a,b set the other member of that quadratic pair to one. Thus any agreeing
circuit depends on all eight inputs and needs at least seven binary gates.
The residual three-gate circuit is likewise minimal by its four essential
inputs. Sharpness is not an artifact of avoidable gates at either endpoint.

This satisfies the full stated legal-state conditions and actual dispersal
hypothesis, rather than merely a graph-shaped mockup. It does not claim that
this marked state is reached from the original q=0 state by the prescribed
elimination strategy. The example proves sharpness of the local lemma under
its stated invariants, not an obstruction to strategies with extra reachability
information or to a global amortized improvement.

## 5. Finite checks and their precise universe

`check_ap2_surplus.py` checks two independent items:

1. Every prescribed acyclic shape with three through six gates, distinguished
   first G, two external nonfixed signals, isolated paired pins, and either
   distinct or common partner children. It allows arbitrary sharing among
   all earlier gates, and filters to full output ancestry. The earlier distinct
   partner child is named H_r by symmetry. G may be moved first because it
   depends only on parameters. The search uses a deliberately weak simplifier:
   a gate with one fixed pin is treated as an alias of its other pin even when
   its real truth table would make it constant. Equal pins are similarly
   retained as a literal. Thus it is a topology/propagation regression test,
   not an exhaustive enumeration of all B2 truth tables or all auxiliary states.
2. All 256 inputs of the literal seven-gate sharp witness, its restricted
   output, legal same-fiber essentiality witnesses for all eight initial inputs,
   essential residual inputs, and the nondegenerate polar matrix.

Both normal Python and `python -O` runs must match `receipt.json`; explicit
exceptions survive optimization. These finite checks corroborate the all-state
proof, rather than replace the nonconstant-output and safe-core hypotheses.
