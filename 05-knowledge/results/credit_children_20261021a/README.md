# Grounding the scale-credit owners and their paid children

Status: **FINITE-EXACT** for all 380 declared base owners and supplied parent
words; **PROVED application of inherited inverse/debt identities** to every
guarded three-scale ray. Universal entry and the remaining uniform residual
bound are **OPEN**.

This packet connects the source ownership in
[scale_credit_supply_20261021a](../scale_credit_supply_20261021a/README.md) to
literal fully odd Collatz certificates. Its 342 linear-family sources and
38 early-frontier sources have 374 distinct parents. Every parent has a
supplied, replayed route to 1; the longest has 99 odd steps. Production
verification never searches for a route. Bounded discovery is separated in
`discover.py`, with a declared cap of 10,000 steps per parent.

The closest mechanisms are the primitive inverse and scale lifts in
[residual_external_credit_20261020k](../residual_external_credit_20261020k/README.md)
and the literal original-source debt in
[source_gap_bursts_20261020k](../source_gap_bursts_20261020k/README.md).
The corrected near miss is to use a negative interaction as a ROOT
certificate; it is not one. The hostile case is the branch a=1,j=0, which
need not descend from its original source. We require s>u and verify every
valuation. The least-used sidecar is the relation between a primitive edge's
strip clock and its two literal odd Collatz steps.

## The compiler and its payment

Write U(x)=(3x+1)/2^{v_2(3x+1)}, B(x)=4x+1. For a primitive parent u and
an admitted inverse branch a in {1,2}, strip count j>=0, put

\[
s=\frac{2^a(3u+1)4^j-(2^a+3)}9,
\quad z=B^j(u),\quad b=v_2(3u+1).
\]

The source guard is: the fraction is integral, s is primitive, s>u, and its
actual primitive edge is (u,a+2j). Then U(s)=z with valuation a. If j>0,
U(z)=U(u) with valuation b+2j, since 3B^j(u)+1=4^j(3u+1).

Suppose u has supplied ROOT word w=(b,w_2,...). For j>0 the child's word
is (a,b+2j,w_2,...); for j=0 it is (a,w). When u=1 the strict parent word is
empty and the j>0 word stops after (a,2+2j). These are literal identities,
not statistical approximations. The verifier stops exactly at its first 1.

Use inherited debt D'=3D/2^a+max(0,1-4/2^a), initially zero. For a=2,
z=(3s+1)/4<s when s>1 and D after that step is zero: the first child is
strictly paid. For a=1, s>u forces j>=1. The second endpoint is U(u), and

\[
U(u)\le(3u+1)/2<s,
\qquad D_2=1-4/2^{b+2j}<1.
\]

The first inequality follows already from s>= (8u+1)/3 at j>=1. Both s and
U(u) are odd, so their positive gap is at least two. Hence
D_2<s-U(u): payment occurs by step two. The selected child's own suffix is
replayed and independently inverted back from 1.

For the 380 base sources there are **360 first-step payments and 20
second-step payments**. This grounds the children actually used by the new
credit bank. It does not identify them with the terminal charges or move a
charge to a different source.

## Infinite guarded rays, with separate clocks

Replace j by j+9k, k>=0. With F=64^{3k}, the source becomes

\[
s_k=Fs+(2^a+3)(F-1)/9.
\]

The primitive parent is unchanged, its native edge clock increases by 18k,
and its literal second valuation increases by 18k. The same supplied parent
word grounds every member of the ray. Original source/paid child registers
are recomputed rather than inferred from an unmarked residue class. The
source remains in [4N F,256N F), and the three-scale charge phase restores.

There are two different clocks. Native arrival time stops at the fixed
frontier after a declared number of primitive edges; its accumulated
division clock controls the charge sign. The fully odd word goes on to 1
and has its own step count and valuation sum. Both are retained. Grounding a
parent does not show that a growing set of parents arrives by a uniform
native deadline. The rays of a finite parent bank also do not provide an
increasing number of distinct owners within a single annulus.

## Exact scope and reproduction

The map is from `(source, parent, branch, strips, parent ROOT word)` to a
literal paid prefix plus grounded suffix. It preserves the original source,
actual valuations, both modular registers, and native receipt clocks. A
finite word can be compressed to the parent word and two integers without
losing its reconstruction. Dropping the parent word would lose the stopping
certificate; dropping the source would lose the payment comparison.

```text
python -B 05-knowledge/results/credit_children_20261021a/check.py
python -B -O 05-knowledge/results/credit_children_20261021a/check.py
```

The checker pins the upstream credit receipt, verifies the exact 380-source
universe, replays every supplied parent word, reconstructs sources backwards
from 1, verifies strict payment and each child suffix, and checks lifted
instances k=1,2,10 (up to 201 bits). The all-k conclusion is the displayed
identity, not extrapolation from those probes. Boolean branch aliases,
illegal j=0 prefixes, and false valuation words are rejected. Every check
remains active under -O. There are 545,725 body checks; verification adds one
immutable receipt comparison.

Independent API review caught a draft issue before promotion: for j>0 the
compiler replaced the parent's first valuation without authenticating it.
The full receipt had already replayed every parent, so its source claims
were unaffected; the standalone compiler now replays the supplied parent
before slicing. A forged first letter and a ROOT-padded parent are retained
as hostile controls.

These are 380 explicit grounded source families. They do not prove that
every positive integer lies in one, that the new credit supply is timely at
every scale, or that R' satisfies its all-horizon budget.
