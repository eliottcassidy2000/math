# Source-gap payment with interrupted reset bursts

2026-10-20. **PROVED scoped:** a sharp affine emission envelope after forgetting
letter order, a source-specific payment guard, and a common-floor certificate
for a variable block language. **FINITE-EXACT:** supplied strict ROOT histories
for the67 sources2^e-1,14<=e<=80. Universal entry, universal first-exit payment,
and the separate native-annulus153N inequality remain **OPEN**.

## Inheritance and the changed question

The [positive complement packet](../source_gap_selector_20261020j/README.md)
proves K=x+1-D>0 for every actual odd Collatz history, and therefore E=x-1+D<2x.
It pays any half-descent and any descending head followed by three native
valuations at least2. The [stable reset monoid](../grounded_block_children_20261020i/README.md)
already handles arbitrary debt-stable ordered blocks. The present question is
different: can a cheap summary of an interrupted reset burst pay a head whose
entire incoming history is unknown except for K>0 and its original source?

Canonical hostile: source55 has a later descent that does not pay; descent
alone is insufficient. Corrected near miss: a formal valuation permutation
need not be native at the same integer. Least-used sidecar: the total number
of valuation1 injections and the actual valuation cost of the repair.
The live board is source / current child / positive envelope / burst summary /
actual word and clocks / grounded ROOT suffix. This is a companion lane to
the annulus interval-flow calculation, not an identification of its debt.

## 1. A sharp summary of the positive envelope

Use the fully odd map x'=(3x+1)/2^a with its exact native valuation a.
The inherited positive envelope obeys

    E'=(3/2)E+1             when a=1,
    E'=(3/2^a)E             when a>=2.

For a finite word v let l be its length, A its valuation sum, and r its
number of1 letters. Then for every incoming E>=0,

    E_after = mu E+P_v <= mu E+2((3/2)^r-1),
    mu=3^l/2^A.                                             (1)

Proof: exchange adjacent letters1,a with a>=2. Performing1 then a
contributes3/2^a, whereas performing a then1 contributes1; their multiplier
is the same. Moving every1 to the end only increases P. The other letters
have no additive injection, and r final ones inject the geometric sum
1+3/2+...+(3/2)^(r-1)=2((3/2)^r-1). Equality is attained by this formal
reordering. It is sharp for the prescribed multiset of affine maps, not a
claim that those permutations all occur at the original integer source.

The map word -> (l,A,r) retains the multiplier and this worst-order envelope.
It loses the exact injection, valuation legality, original source and child;
those remain in the authenticated word. For example(1,2,2) injects9/16,
while(2,2,1) injects1. A summary may reject an actually paid ordered word.

## 2. The source-specific burst guard

Let a supplied actual history from original odd n>1 reach h<n. Its incoming
envelope satisfies E_head<2h even if earlier descents and resets occurred.
If an actual repair v satisfies

    2h*3^l/2^A + 2((3/2)^r-1) <= n-1,                       (2)

then the concatenation pays the same source: E_after<n-1, equivalently
D_after<n-child. Strictness is supplied by E_head<2h; equality in(2) is
allowed. No imposed valuation, source inflation, or credit reuse occurs.
The repair may contain1 letters and arbitrarily many different reset sizes.
If it reaches ROOT, stop there; padding1 is not an allowed continuation.

This is a sufficient guard. It need not be necessary: actual order can
lower P_v, and the real E_head can be much smaller than2h. It proves no
existence of a future word meeting the guard.

## 3. A common floor permits changing block patterns

Every length-three word with at most one1 satisfies, for E>=0,

    E_after <= (27/32)E+1.                                (3)

The equality case is(2,2,1). Larger non-one valuations only decrease the
result. Therefore any sequence of such blocks, with freely changing
patterns but all actually native, obeys after q blocks

    E_q <= c+(27/32)^q(E_0-c),       c=32/5.                (4)

For original odd n>=9, n-1>c. After any descent to h, a finite q satisfying

    c+(27/32)^q(2h-c) <= n-1                              (5)

therefore pays the original source. It suffices that the actual next3q
letters admit this consecutive three-letter partition; a prior ROOT is
already grounded. This extends the three-consecutive-reset guard to a
language allowing interrupted resets, without prescribing a single repeated
word. It does not establish that an arbitrary orbit enters that language.

There is even a uniform finite bound in this guarded class: **ten blocks,
or30 odd steps, suffice for every original odd n>=9 after any descent**.
Indeed h<=n-2, and

    ((n-1)-c)/(2(n-2)-c)=(5n-37)/(10n-52)>=4/19

for n>=9, while(27/32)^10<4/19. Substitute into(5). The constant ten is
the first permitted by this relaxed common-floor bound at n=9; this is
not a claim of a native orbit requiring thirty steps. Earlier payment or
ROOT terminates the certificate. The unproved condition is still the
actual thirty-letter burst pattern, not its sufficient payment bound.

More generally, a block with exactly r ones and s nonones has
mu<=3^(r+s)/2^(r+2s) and injection at most2((3/2)^r-1).
Whenever this multiplier is below1, the quotient injection/(1-multiplier)
is a common floor for the stated class. The floor must lie strictly below
the original payment threshold. One1 and one2 have multiplier9/8>1: merely
alternating an up-step and a reset is not a contraction certificate.

## 4. Grounded Mersenne controls

[inputs.json](inputs.json) contains every source2^e-1 with14<=e<=80, each
with a complete strict odd ROOT word. A separate bounded [discovery script](discover.py)
used a5000-step cap; none hit the cap, and the largest word has386 letters.
Production [check.py](check.py) searches for no orbit or ROOT suffix: it
replays and independently inverts every supplied word.

At the actual first exit,35 sources satisfy the inherited half-descent
guard. The remaining32 satisfy(2) after an actual repair of length at most24.
Of these,30 lack the previous immediate three-nonone/early-ROOT repair gate.
These are new guard certificates in this declared experiment, not a claim
that the sources were historically unverified or that all Mersenne exponents
have been covered. In particular this finite range has no exponent-lift
argument identifying an unresolved original source with a different one.

Each receipt retains the original source, compatible mod15/mod35 registers,
actual first-exit cut, paid child, that child's literal ROOT suffix, odd-step
counts and valuation-sum clocks. The fully odd clocks here are not the
canonical-P braid/depth clocks in the annulus calculation. Both are kept
in their own carriers. No proof of the153N bound is inferred from this bank.

## 5. Exact controls and remaining obligation

The checker independently verifies(1) on every formal word of lengths0..6
over{1,2,3,4}, with the maximizing reordering, and(3) on all guarded triples
over that alphabet. These are controls for the all-letter proofs above,
not native trajectories. It then verifies the complete67-source input
universe, first ROOT and first exit, strict payment, actual child suffixes,
source registers and both clocks. Checks use explicit exceptions under-O.

    python -B -X utf8 05-knowledge/results/source_gap_bursts_20261020k/check.py
    python -B -O -X utf8 05-knowledge/results/source_gap_bursts_20261020k/check.py --verify

Both runs pass108,061 exact controls. Script and input LF hashes are in
[receipt.json](receipt.json). The outstanding universal statement is actual
entry into a paying language or another grounded descent certificate;
neither a common floor nor a successful finite bank supplies that quantifier.

Independent peer proof and API audit: **PASS**. The peer checked the sharp
histogram envelope, the thirty-step common-floor corollary and its entry
condition, literal source/register types, strict ROOT stopping, and the
retained child suffixes. The optimized replay matches the saved receipt.
