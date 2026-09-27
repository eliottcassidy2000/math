# Three unit colours give an exact marked Zeckendorf lift

**PROVED elementary coding and loss statements / FINITE-EXACT controls /
OPEN global stripe continuation and Collatz reset inequality.** No novelty
claim. The user's35 supplied colours remain authoritative data throughout.
This note specifies one natural interpretation of the extra-unit proposal;
it does not assume that the prose uniquely fixes the indexing convention.

## 1. Inheritance and the construction

The anchor is an exact nonconsecutive representation with three differently
coloured units. The niche is the boundary charge of
[the earlier colour note](crossroads_crossing_20260926_colour.md), C7--C10,
and the wildcard is whether these markers preserve the roles of the
[precision-collision controller](nextforest_20260926_precision.md).
The live board is unit identity, nonadjacency, full charge, row coherence,
normalization, and source-preserving precision history.

The closest proved mechanism is the terminating distinct-word rewrite in
[the Zeckendorf charge note](duck_zeckendorf_20260925.md), Z1--Z9.
The hostile is equating the same visible integer length with the same
lifted colour. The corrected near miss is calling a Fibonacci atom's
ordinary Zeckendorf representation nontrivially splittable. The least-used
sidecar is the distinguished unit's position in the extended alphabet.

Use ordered atoms with weights

    g_0=g_1=g_2=1,  g_i=F_i for i>=3,
    colours(g_0,g_1,g_2)=(K,B,R),                    (ZC1)

where K is black, B blue, R red. The weights are therefore
`1,1,1,2,3,5,8,...`. An extended representation is a finite subset of
indices with no two consecutive indices. The original ordinary
Zeckendorf alphabet is the part indexed2,3,4,... . The two extra units
are exactly g_0 and g_1; the ordinary unit g_2 stays red.

This definition makes nonconsecutive refer to *ordered atom indices*.
It permits red1 and black1 together, since their indices differ by2,
but forbids using black1 and blue1 together. Equal numerical weights
are distinct atoms, as required by the colouring proposal.

## 2. Complete fibres: every positive integer has two or three representations

Let Z(n) be the ordinary Zeckendorf support in indices>=2. Every extended
representation of n>0 is exactly one of

    Z(n),
    {0} union Z(n-1),
    {1} union Z(n-1), provided 2 is absent from Z(n-1).  (ZC2)

**Proof.** If neither extra unit occurs, ordinary uniqueness applies.
If black index0 occurs, index1 is excluded and all remaining indices
are an unrestricted ordinary Zeckendorf support of n-1. If blue index1
occurs, both0 and2 are excluded, so the ordinary support of n-1 must
omit2. These exhaust the possibilities and are all legal.

Consequently the fibre size is exactly2 or3, never unbounded. This is a
positive exact coding theorem for the proposed extra units. It is more
specific than matching a small graph or a number of colours.

In particular every Fibonacci atom g_i, i>=3, has a proper decomposition
into smaller-index nonconsecutive atoms. The preferred one is

    g_i=g_(i-1)+g_(i-3)+...,
    ending at g_0 when i is odd, g_1 when i is even.   (ZC3)

Examples are `2=red1+black1`, `3=2+blue1`,
`5=3+red1+black1`, and `8=5+2+blue1`.
For odd i this is the only proper representation. For even i there is
one other: replace the final blue unit by the black unit. Formula(ZC2)
and the usual alternating decomposition of F_i-1 prove this classification.

## 3. Exactly two sheets preserve the full integer charge

Retain the inherited charge Q(n)=(n-2b(n),b(n)), with
`b(n)=floor((n+1)/phi^2)`. Put

    U=(1,0), V=(-1,1),
    charge(red1)=charge(blue1)=U,
    charge(black1)=V,
    charge(g_i)=Q(F_i) for i>=3.                     (ZC4)

Both U and V have visible length1 under `(A,B)->A+2B`. Red and blue
are distinct typed atoms despite having the same full integer charge.
Their distinction requires the marker, not another interpretation of Q.

For every m>=0,

    Q(m+1)-Q(m)=V if Z(m) contains the ordinary unit,
                U otherwise.                        (ZC5)

There is a direct rewrite proof. If the unit is absent, adjoining it
gives a distinct representation; normalizing adjacent Fibonacci atoms
preserves charge and adds U. If the unit is present, its neighboring
weight2 is absent. Replace the two units by that weight2, whose charge
is U+V, and then normalize the remaining distinct word. Relative to
the old charge this adds V. No estimate of an irrational rotation is
needed for this argument.

Exactly two of the representations(ZC2) therefore have charge Q(n):

    ordinary sheet: Z(n),
    marked sheet: {0}+Z(n-1) if 2 is in Z(n-1),
                  {1}+Z(n-1) otherwise.              (ZC6)

When a third representation exists it is the black version, and its
charge differs by `V-U=(-2,1)`. Thus one extra boundary bit, together
with n, identifies every charge-preserving representation. The permitted
nonordinary marker is determined by n; it is not a freely chosen ternary
digit.

There is an explicit terminating normalizer. On the marked blue sheet,
replace blue1 by red1, which was absent. On the marked black sheet,
combine black1 with the present red1 to make the weight2 atom. Then use
the ordinary distinct-word merges `011->100`. The first operation removes
the marker and the remaining operations reduce the number of atoms.
Charge and visible value are preserved. Keeping the original sheet bit
makes the normalization lossless; dropping it forgets which representation
was supplied. This terminating representation algorithm is not a Collatz
time-descent algorithm: its steps keep n unchanged.

No invariant depending only on the visible integer can assign three
different values to the three units, since all three represent1.
The construction instead distinguishes typed representations and states
exactly which changes preserve Q. It does not identify these colours
with addition in Z/3Z.

## 4. Recursive piece colours, a legal35th blue, and the next-row obstruction

Keep the three coloured units terminal. Recursively expand every larger
Fibonacci atom by(ZC3), listing the pieces in descending index order.
Let W_i be its resulting word. Then

    W_0=K, W_1=B, W_2=R,
    W_3=RK, W_4=RKB,
    W_i=W_(i-1) W_(i-2), i>=5.                       (ZC7)

The last identity follows by comparing the two alternating lists(ZC3).
The words have lengths F_i for i>=2 and are nested from i=2 onward.
Their infinite limit starts

    RKBRKRKBRKBRKRKBRKRKBRKBRKRKBRKBRKR...

The ordinary Zeckendorf expansion of any n, with its Fibonacci pieces
listed in descending order and each replaced by W_i, gives the length-n
prefix of this limit. Indeed if F_j<=n<F_(j+1), the remainder is below
F_(j-1); the concatenation W_(j+1)=W_j W_(j-1) provides exactly the needed
next prefix, and induction treats the remainder. The small unit cases
are immediate. This is a precise prefix-coherent row construction.

The authoritative supplied word is

    RKBRKRKBRKBRKRKBRKRKBRKBRKRKBRKBRKB.            (ZC8)

It matches all first34 symbols of(ZC7), but has B at35 instead of R.
There is a useful distinction missed by simply rejecting a Fibonacci
continuation: **the entire35-symbol supplied word is a valid isolated
row under the extra-unit proposal.** Represent35 as34+blue1 and expand
the34 atom. This row even has the correct full charge Q(35), because
red1 and blue1 both have charge U.

It cannot extend to a prefix-coherent row36 under this same grammar.
The complete representations of36 are

    34+2,
    34+red1+black1.                                  (ZC9)

Both expand to a34-word followed by RK, so position35 must be red.
Replacing it by blue at length35 was a valid boundary choice, but that
boundary choice is incompatible with the next fixed diagonal stripe.

This remains true even if every occurrence of an even-index Fibonacci
atom is allowed to choose its alternate black ending. By induction from
the complete proper-representation classification, the full expansion
language of any g_i consists exactly of W_i with an arbitrary subset of
its blue symbols replaced by black. Those choices never change a red
position. Thus all possible expansions in(ZC9) still force red at35.

**Scope.** This excludes prefix-coherent continuation only in the stated
descending-piece, terminal-unit grammar. It does not exclude a different
diagonal indexing, a piece-order change, a context-dependent unit role,
or a recolouring rule. Such a rule must be stated to choose a continuation;
the authoritative35 symbols alone do not supply it. No supplied symbol
has been silently changed.

## 5. What the couple of marker bits can retain in a Collatz reset

Consider the exact controller's states

    y=q*2^r+u, q,u positive odd, r>=1.

At a *fixed* odd integer y>1, there are exactly `(y-1)/2` such triples.
The bijection chooses an even integer x with0<x<y, writes uniquely
`x=q*2^r` with q odd, and sets u=y-x. Thus the amount of decomposition
history attached to one visible integer is unbounded.

The branch partition itself has a clean fixed-source formula. Put
`t=v2(3y+1)`. The controller consumes precision when r>t, swaps the
roles when r=t, and collides when r<t. To prove this, use

    3u+1=(3y+1)-3q*2^r.

For r!=t the unequal valuations give their minimum; for r=t, the two
odd coefficients subtract to a positive even integer and the valuation
increases. This reproduces exactly the three cases of the precision note.

There is a small decisive comparison at y=41:

| q,r,u | a=v2(3u+1) | controller branch |
|---|---:|---|
|1,1,39|1|collision|
|3,2,29|3|swap|
|1,3,33|2|consume|

But41 has only two extended nonconsecutive representations,

    34+5+2,
    34+5+red1+black1.                                (ZC10)

So these representations cannot even carry all three branch tags
injectively at that fixed integer, much less its20 full triples.
Adding two unrestricted bits gives four states and still cannot retain
all20 triples. More generally an exact encoding of every decomposition
over the same visible y needs at least
`ceil(log2((y-1)/2))` bits of extra information.

The [ternary role/run clock](reset_20260926_clock.md) sharpens this
comparison: after a noncollision exactly one of q,u is divisible by3.
Even within that restricted interface, for any fixed odd y coprime to3
there are exactly `floor(y/3)` legal triples. In the even-x bijection,
they are the x with residue0 or y modulo3. Counting the two classes
gives2k for y=6k+1 and2k+1 for y=6k+5. All three displayed41 witnesses
obey this role invariant, and13 states survive over41. These are legal
states compatible with the noncollision output condition; no claim is
made that every one occurs in one particular canonical-reset orbit.

This does **not** say that n plus a finite marker cannot compute Collatz:
n itself suffices to compute its next iterate. Nor does it exclude a
canonical algorithm choosing just one triple, or an arbitrary nonlinear
potential of the full n. It says that the full precision
decomposition cannot be compressed to those couple of bits while
preserving all its legal states. Which q/u role was swapped can be
marked by a bit; their values and precision still need their exact
unbounded registers or words.

The same distinction explains why the colour normalizer cannot pay the
reset inequality. It terminates at the same integer, whereas the arithmetic
controller changes that integer and can recreate arbitrarily much precision.
The inherited family `(1,1,(4^k-1)/3)->(1,2k-1,3)` is a required hostile.
An inequality that preserves this full decomposition must account for its
unbounded registers; whether a quotient can suffice for descent remains
open. The existence of the finite controller alone gives no sign to it.

## 6. Connection ledger and reproduction

| Source -> target | Map and preserved predicate | Lost data / sidecar | Decisive test |
|---|---|---|---|
| Three-unit independent sets -> n plus marker |(ZC2), all exact representations | marker needed to reconstruct the original word | fibres have2 or3 elements |
| Charge-preserving words -> n plus one bit |(ZC6), full Q and ordinary value | red/blue differ despite equal Q | row35 has two charge-correct endings |
| Fibonacci pieces -> nested coloured rows |(ZC3),(ZC7), lengths and nonadjacency | another row-boundary choice need not extend | legal35-blue row has no row36 continuation |
| Exact precision state -> visible y |q2^r+u, actual integer and its future orbit | chosen decomposition and valuation accounting |20 states over41, only2 row representations |
| Marked normalization -> Collatz rank | no arrow-preserving map established | representation moves keep n fixed | arbitrarily recreated precision in the inherited reset family |

Run [the experiment](../../04-computation/experiments/reset_20260926_colours.py)
normally and with `-O`. The [retained output](reset_20260926_colours.out)
records complete extended representation fibres for n1..6764; preferred
prefix coherence through6764; complete Fibonacci expansion languages through
F9=34; all isolated35/36 row choices; and2,096,128 precision triples across
odd sources3..4095. The minimal41 example is a finite-exact minimality
claim in that searched universe; the stated example and unbounded fibre
formula have direct proofs. All checks remain active under `-O`.
