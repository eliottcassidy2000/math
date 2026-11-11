# Incoming offset35 closure and an explicit three-level9-grid constructor

**PROVED integration and explicit constructors.** The incoming full35-grid
closure is credited to its source. The present packet supplies a kernel-based
middle level for the9-grid, combined with two elementary product levels.
Neither a full9-grid theorem nor universal seed rearming is asserted.

## 1. Exact incoming source and its contribution

Incoming commit `b42131e609ec7a66a3f41a49e6d87592e7b29069`,
`05-knowledge/results/multiplicand_saturation_chanpak_20261111/chanpak/PROOF.md`,
proves every

    N(35,r,s)=1+2^s+35*2^(r+s),    r,s>=1.

The previously missing s=2 level is realized by the ordinal product
5(1+7*2^(r+2)); levels3,4,5 use literal selected sources with offsets9,17,33.
Levels1 and s>=6 are credited there to the earlier seed-recursion packet.
We read the proof and retain its exact scopes without broadly re-auditing it.
No35-grid closure is claimed as a contribution of the current constructor lane.

The point of contact is physical, not just arithmetic: the general tail law
is D+2^t U, not1+2^t U. The incoming source deliberately retains positive
nonunit D. Our near-rearm family has D=9−4*2^m, eventually negative, and still
produces a valid positive Hamilton count through that same exact graph
operation. A signed offset is a bookkeeping term, not a negative graph.
Neither argument installs an abstract same-q rearm or permits dropping D.

The incoming product lemma says that when2^s+1 divides d, the residual level
factors as

    N(d,r,s)=(2^s+1)[1+(d/(2^s+1))*2^(r+s)],

provided both factor counts have actual tournament constructors. This gives
one of the following two easy levels and clarifies why the middle level
needs our different kernel.

## 2. The9-grid at levels1,2,3

**Theorem.** Every N(9,r,s)=1+2^s+9*2^(r+s), with r>=1 and s=1,2,3,
is the count of an explicitly constructed tournament.

### Level1: ordinal multiplication

    N(9,r,1)=3(1+3*2^(r+1)).

Use C_3 for the factor3. For the inner factor, the selected base C_3 with
two vertices selected has (H,W,U)=(3,2,6). Its tail of length r>=1 has
count1+6*2^r=1+3*2^(r+1). Ordinal join multiplies the counts. This is a
recovered use of the inherited tail/product mechanisms, not a new ray.

### Level2: count-kernel replacement and one actual head

For r>=2 substitute a singleton and TT_(r−1) into ports5,7 of the quotient
from BILINEAR_KERNEL.md:

    Q=(198,196,136,211,167,207,20,64).

Its singleton row gives

    H(Q[1,TT_(r−1)])=72*2^(r−1)+5
                    =5+9*2^(r+2)=N(9,r,2).

At r=1 the target is77. Use the literal seven-vertex tournament

    (32,9,3,69,15,30,55).

This head is inherited from the record `19 2 77 census 7 131216` in
`05-knowledge/results/kps_s21_two_adic_lines_20261110/below_rays.txt`.
Its rows were decoded and its77 Hamilton paths counted directly. It is
credited as an inherited head, not newly discovered numerical coverage.
The constructor has exactly r+6 vertices for every r>=1.

This level is not supplied by the incoming divisibility lemma:5 does not
divide9. The operation really changes from the finite head to our two-module
kernel. There is no assertion that the77 graph itself admits a single port
that generates the entire line.

### Level3: incoming product principle, with explicit factors

    N(9,r,3)=9(1+2^(r+3)).

The factor9 is C_3→C_3. A singleton selected base has (H,W,U)=(1,0,1),
so its tail of length r+3 gives1+2^(r+3). Ordinal join gives the displayed
count. This is the incoming product lemma with d=9 and s=3, whose realization
hypotheses are supplied explicitly.

## 3. Status and reproduction

The twelve literal controls r=1,...,4 and s=1,2,3 are in
`realize9_lower_receipt.json`. The script emits the full graph for arbitrary
valid parameters:

    python3 realize9_lower.py
    python3 realize9_lower.py --r 20 --s 2

The all-r conclusion is proved by the kernel identity and tail/product laws,
not inferred from twelve samples. We make no claim here for s>=4, no claim
that these values were previously missing from the numerical spectrum, and
no global novelty assertion against unsearched historical constructor menus.
The concrete contribution is this compact, uniform physical construction
and its integration with the current offset35 proof.
