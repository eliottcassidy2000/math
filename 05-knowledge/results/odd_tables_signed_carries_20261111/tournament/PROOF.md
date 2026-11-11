# A ternary carry actuator and a uniformly beneficial paired root reset

**2026-11-11. PROVED IN NOTE:** an arbitrary-prefix ternary simultaneous
rebase; a target-preserving change of root that improves both witnesses
at every such rebase; a small133/135 reset; and an explicit infinite
beneficial family. **FINITE-EXACT:** the small controls in `receipt.json`.
**OPEN:** completing a prescribed original target rather than a padded
address, and the Chan–Pak odd-spectrum conjecture.

## 1. Inherited mechanisms and the missing coordinate

The inherited [residue-phase note](../../h13_paired_descent_20261111/square_graceful/RESIDUE_PHASE.md)
gives two actual paired families, a multiplier-three carry lift, and an
even doubling ray on which their targets agree. The
[repeated-9 actuator](../../h13_paired_descent_20261111/square_graceful/ADDITIVE_CARRY.md)
gives an arbitrary-prefix completion in a different residue phase, but
its H3 witnesses always have worse future response than the H13 ones.
The [rooted-order theorem](../../h13_paired_descent_20261111/dynamics/PROOF.md)
and [ordered amalgamation](../../h13_paired_descent_20261111/modular_cycles/ADDITIVITY.md)
supply the exact future comparison used below.

The concurrent [odd-multiplier compiler](../../central_modular_chanpak_20261110/compiler/PROOF.md)
already lifts every fixed odd address multiplier; that general lift is
inherited here. Its warning that an address operation is not automatically
a scalar-H Collatz operation remains essential. The new content is a
particular coupled zero-defect language and a physical source/root change.

The closest proved mechanism is a common defect for two literal
tournament witnesses. The hostile is a scalar equality with wrong
endpoint response. The corrected near miss is all-prefix completion
mistaken for all-target completion. The underused coordinate is the
choice of retained root inside an already realized graph. The live board
is multiplication-by3 carry / shared target defect / periodic residue
phase / actual root position / positive future response.

Retain the inherited binary function

    g(0)=g(1)=1,
    g(2n+1)=g(n),
    g(2n)=g(n)+n+(n mod2).

The H13 core Q=(18,8,3,5,14), with apex cuts6 and3 and right address B,
has counts L(B),L(B)+2, where

    L(B)=13g(B)+32B+6(B mod2)+18.

The C3 core=(2,4,1) gives another actual pair:

    F_(C3,cut3)(0,2C+1)=N(C),
    F_(C3,cut1)(0,4C+3)=N(C)+2,
    N(C)=3g(C)+12C+10.                             (1)

These named cuts and addresses specify every tournament and apex arc.
The shared defect for the map C=3B is

    D(B)=L(B)-N(3B)
        =13g(B)-3g(3B)-4B+6(B mod2)+8.             (2)

Every zero of D simultaneously gives actual C3-based witnesses for both
H13 targets. Their residues are1 and0 modulo3. We will construct many
zeros and then improve their future response by changing the retained
root, without changing either realized target.

## 2. A second residue-five actuator, now with multiplication-by3

Define

    F(x)=13g(x)-3g(3x)+x+13(x mod2)+20.

**Lemma1.** For all x>=0,

    D(16x+5)=F(x),
    F(16x+3)=F(x)+4.                              (3)

Proof. The binary recurrences give

    g(16x+5)=g(x)+5x+(x mod2)+2,
    g(48x+15)=g(3x),
    g(16x+3)=g(x)+3x+(x mod2),
    g(48x+9)=g(3x)+18x+4.

Substitute in (2) and F. QED.

Consequently, for m>=0, let

    x_m(q)=q16^m+(16^m-1)/5,
    B_m(q)=16x_m(q)+5.

The integer x_m appends m copies of binary0011 to q; B_m then appends
0101. The exact residual is

    D(B_m(q))=F(q)+4m.                            (4)

This is the operational modular-table connection. The repeating block
0011 represents3/15=1/5, and multiplication by3 sends it to the block
1001 representing9/15=3/5:

    3x_m(q)=3q16^m+3(16^m-1)/5.

There is no carry between the repeated blocks. The two period-four
phases are synchronized by multiplication by3 modulo5, while (4)
retains the integer defect discarded by a residue-only table. The period
is four because 2 has order4 modulo5. This does not assert a Collatz
orbit theorem or identify the address B with the Hamilton count L(B).

## 3. Every prescribed binary prefix has a completed zero

**Lemma2.** For every p>=0, exactly one of the two choices

    P=8p+2,       P=8p+4

satisfies `g(P)+g(3P+2)=0 modulo4`.

Proof. Both g values are odd, so their sum is0 or2 modulo4. Direct
binary expansion gives

    g(8p+2)+g(3(8p+2)+2)
      =g(p)+g(3p+1)+26p+10,
    g(8p+4)+g(3(8p+4)+2)
      =g(p)+g(3p+1)+18p+12.

The difference is8p-2, which is2 modulo4. QED.

Choose this P and set

    q=32P+31,
    m=-F(q)/4.

Five appended1 bits give

    g(q)=g(P),
    g(3q)=g(3P+2)+24(P+1),
    F(q)=13g(P)-3g(3P+2)-40P-8.                  (5)

Thus F(q) is divisible by4. Since P>=2 and the inherited bound is
`g(P)<=P+1`, while `g(3P+2)>=1`,

    F(q)<=-27P+2<0.

So m is a positive integer.

**Theorem1 (all-prefix completion in the missing phase).** For every
finite binary prefix represented by p>=0, the procedure above produces

    B=16[q16^m+(16^m-1)/5]+5

with D(B)=0. For p>0, its binary word is exactly

    prefix(p) · (010 or100) · 11111 · (0011)^m · 0101.     (6)

The empty-prefix case uses p=0. Therefore both targets L(B),L(B)+2 have
the literal C3 witnesses in (1) at C=3B, or right addresses6B+1 and12B+3.
They lie in the residue phase1/0 modulo3, which the inherited
repeated-9 family did not cover.

Proof. Lemma2 supplies the parity repair, (5) supplies sign and
integrality, and (4) gives D(B)=0. Equation(1) gives the actual graphs.
QED. The construction is not claimed to be efficient in input bit length:
m can grow proportionally to the numerical prefix p.

## 4. Re-rooting repairs the future response in every leading phase

Use the standard rooted state `(H,E,U,d)`. Consider any vertex r with
exactly one out-neighbour x.

**Lemma3 (one-out-neighbour equality).**

    U=H+d-E.                                      (7)

Proof. The root-alone distinguished path contributes d. Every other
root-start distinguished path starts r→x. Its complementary path, if
nonempty, ends at a vertex different from x and hence points into r.
Concatenation gives a Hamilton path not ending at r. Conversely split
any such Hamilton path immediately before r. These are inverse maps,
and there are H-E Hamilton paths not ending at r. QED.

For the two C3 graphs in (1), choose r to be the **last right-flank
vertex**, immediately before the old apex in the compiler's vertex
order. It beats only the old apex. This changes no vertices, edges,
or Hamilton counts.

For C>=2 write

    c=bitlength(C),
    C0=C-2^(c-1),
    epsilon=the second leading bit of C.

Deleting r removes the leading1 of the right-flank address. Any remaining
leading zero bits correspond to universal sinks, so they can be stripped
without changing the count. The deletion counts are consequently
N(C0) and N(C0)+2. To count paths ending at r, append r to a deletion
Hamilton path. This succeeds unless that path ends at the old apex.
If epsilon=0, the last remaining flank vertex is a universal sink, so
no path ends at the old apex. If epsilon=1, the inherited endpoint formula
gives exactly3 such paths. The same epsilon applies to both graphs.
Lemma3 now proves the complete states

    lower: (N(C),   N(C0)-3epsilon,   N(C)+3epsilon,   N(C0)),
    upper: (N(C)+2, N(C0)+2-3epsilon, N(C)+2+3epsilon, N(C0)+2).    (8)

The pair difference is exactly `2(1,1,1,1)`. The old apex-root pair had
a pure-H difference. Thus a root change within already realized graphs
transfers that difference to the singleton direction while preserving
both current targets.

**Theorem2 (uniform beneficial simultaneous reset).** At every B>=1
with D(B)=0, replacing each H13 witness by its corresponding C3 graph
rooted as in (8), with C=3B, strictly reduces the last three coordinates
of `q=(H,U,U-d,H+U-E)` and preserves the first coordinate H.

Proof. The old apex-root states have E=d=13. Write U_old=h+delta for
the lower or upper target h. The inherited O/Z recurrence implies

    delta>=8 on the lower branch,
    delta>=6 on the upper branch.                  (9)

For clarity, the invariant before the final digit is U-H-d+E. Its base
values at31/33 are8/6. A Z step preserves it; an O step increases it
by U-E>=0. Every nonempty canonical right word finishes in Z, making
E=d and giving (9).

Write d'=N(C0) on the lower branch and N(C0)+2 on the upper branch.
The exact old-minus-new difference in q is

    (0,
     delta-3epsilon,
     delta+d'-13-3epsilon,
     delta+d'-13-6epsilon).                       (10)

The function N is at least13. Thus the last three entries are at least
5,5,2 on the lower branch and3,5,2 on the upper branch. All are strictly
positive. QED.

By the inherited exact order theorem, every common O/Z/J program and
every fixed right-root amalgamation preserves weak H dominance. A
common word containing Z or J gives strict H improvement. No modular
or spectral approximation is used in this comparison.

This resolves the leading-phase obstruction of the apex-root comparison:
every completed prefix in Theorem1 has a beneficial source/root reset,
regardless of the leading bits of p. The new root deletion count d' can
be large; the claim is an improvement in the exact future-response order,
not a simultaneous decrease of every complexity statistic.

## 5. A small reset and a stable paired continuation

Already D(3)=0. The old H13 states and the new C3 states are

    (133,13,141,13)  →  (133,25,133,25),
    (135,13,141,13)  →  (135,27,135,27).            (11)

The q savings are respectively `(0,8,20,20)` and `(0,6,20,20)`.
The C3 witnesses use right addresses19 and39, with cuts3 and1; their
new roots are vertex7 and vertex8, respectively. The old H13 graphs
have right address3 and retain their apex. All four states in (11)
are checked by literal subset path enumeration.

The earlier211 reset had U saving78 larger than its old deletion13.
Its general converse-switch obstruction does not automatically apply
to (11): these savings are8 and6, both smaller than13. This records
the exact changed hypothesis, not a claim that a compensating switch
has already been constructed.

The difference direction in (8) is fixed by Z:

    Z(1,1,1,1)=(1,1,1,1).

Thus every common all-Z continuation preserves the unit target gap2.
For either branch, if the initial new state is
`(h,d'-3epsilon,h+3epsilon,d')`, then for k>=1 its exact count is

    H_new(k)=d'+2^k(h-d'+3epsilon).                (12)

The two branches have the same coefficient of2^k and intercepts
differing by2. Moreover their old-minus-new count saving is the SAME:

    H_old(k)-H_new(k)
      =2^k(delta+d'-13-3epsilon)-(delta+d'-13).    (13)

Moving from lower to upper subtracts2 from delta and adds2 to d',
leaving (13) unchanged. Thus the construction retains a coupled
paired descent along an infinite continuation language, rather than
two unrelated source replacements.

## 6. A closed infinite family and the alternative apex comparison

For each integer k>=7 set

    q=2^k+7,
    m=15·2^(k-7)-8.

Binary expansion gives

    g(q)=2^(k-3)+1,
    g(3q)=33·2^(k-5)+7,
    F(q)=-4m.

Therefore B=B_m(q) is a zero of D for every k>=7. Theorem2 applies
to this explicit infinite family. Its first member has

    q=135, m=7, B=580679578421,
    targets=21035117728249 and21035117728251.

For comparison, even without the new root change some completed sources
improve the original apex response. Let l=bitlength(q), j=bitlength(3q).
The old H13 U minus the original C3 apex U equals

    (966/5)(3·2^l-2^j)16^m
       +(52·2^l-9·2^j)/5-13rev(q)+3rev(3q).       (14)

This follows by reversing the repeated0011 and1001 words before the
terminal0101/1111 words. For m>=1 it is greater than10 if j=l+1 and
negative if j=l+2: bound each reversed integer by the corresponding
power of2; the leading term dominates the remaining terms. The displayed
closed family has j=l+1. The root change in Section4 is stronger because
it also handles j=l+2, without sacrificing either original target.

## 7. The exact remaining target gate

Theorem1 takes a prefix p and realizes new targets at a **completed**
address B. It does not preserve L(p). Theorem2 preserves the exact two
targets once D(B)=0, and (13) controls what later continuations cost.
These are different quantifiers.

A full target proof still needs, for a prescribed target n, a legal
decomposition or inverse rule placing n in one of these synchronized
charts; alternatively it needs a target-preserving transport from an
arbitrary original address to a zero of D. Merely appending a defect
actuator changes the target and cannot supply that transport. Rerooting
does preserve the target, but does not itself turn a nonzero D into zero.
The exact residual gate is therefore a coupled inverse/repair equation,
not residue reachability or the existence of arbitrary even widths.

The fixed-d25 switch from
[the preceding construction](../../h13_paired_descent_20261111/construction/PROOF.md)
still supplies every even pure-H width. A lawful use of it here must
satisfy the actual baseline and boundary equations of the amalgamation
theorem. Equations(10)-(13) retain those boundaries rather than adding
independent numerical gaps. The small reset(11) changes one previous
obstruction's parameters and supplies a concrete place to seek that
next repair.

## 8. Reproduction and scope

    python3 -B 05-knowledge/results/odd_tables_signed_carries_20261111/tournament/probe.py

`receipt.json` records128 bounded actuator/preparation controls,32
completed arbitrary prefixes with exact q dominance, four literal C3
reroot controls, the tiny H13 controls at B3, and compact certificates
for two explicit completed addresses. The larger graphs are specified
by their named cores, cuts and addresses; their counts use the inherited
exact compiler, not a claimed direct enumeration of all their paths.
The arbitrary quantifiers follow from the proofs above. No finite
coverage range or outside source was broadly re-audited.
