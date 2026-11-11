# A native unit switch from the prime kernel's completion cube

Date: 2026-11-11. **FINITE-EXACT:** the eight literal completions and their
signatures. **PROVED:** arbitrary-flank switch formulas, native even-width
coverage, and the infinite paired C3 rebase. **OPEN:** universal target
coverage or insertion of a switch into an arbitrary given tournament.

This follows the actual three-port, five-path graph recovered in
[PRIME_KERNEL.md](PRIME_KERNEL.md). It is a constructive continuation of
that recovery, rather than an inference from the coincidence of six cycles.

## A small completion cube replaces the parity-only quotient

Use ports a=0,b=1,c=2 and shorten the five paths to

    Q:0→1,  A0:0→2,  A1:2→3→0,
    B0:1→2, B1:1→4→2.

This chooses positive path lengths(1,1,2,1,2) in the same kernel. It is a
new subdivision of the suppressed kernel; it does not claim to preserve the
square labels or lengths of the prime graph. Seven pairs have fixed arcs;
the missing pairs are(0,4),(1,3),(3,4). Orient them with bits b0,b1,b2,
where1 means smaller→larger. There are exactly eight completions:

| Code b0+2b1+4b2 | Adjacency rows | H |
|---:|:---|---:|
|0|(6,20,8,3,13)|13|
|1|(22,20,8,3,12)|11|
|2|(6,28,8,1,13)|9|
|3|(22,28,8,1,12)|9|
|4|(6,20,8,19,5)|15|
|5|(22,20,8,19,4)|9|
|6|(6,28,8,17,5)|13|
|7|(22,28,8,17,4)|9|

Every completion is strong because the fixed directed paths already connect
all vertices in both directions. The three orientation bits are a minimal
record for these labelled completions: they are independently readable arcs.
They are a different eight-state space from the cycle-parity compiler. The
map between the two constructions retains the **five paths and their ports**;
it does not identify subdivision parity with chord orientation.

Let Q− be code7 and Q+ code3. They differ only in the native arc(3,4), and
both have H=9. Their three directed-triangle arc-overlap graphs are P3 and
K3 respectively. Thus the user's distinction between the two same-count
cores is operational: one marked extension separates them while preserving
the continuation response.

## The same root cut produces23 and25

Adjoin a root r beating precisely vertices{0,2,4}, the cut mask21. The two
literal six-vertex tournaments are

    T−=(22,60,8,49,4,21),
    T+=(22,60,8,33,12,21).

Only arc(3,4) changes. Define E as HPs ending at r, d=H(T−r), and U as the
number of ordered two-path covers whose first path starts at r, allowing an
empty second path. Literal counts give

    (H,E,U,d)(T−)=(23,2,50,9),
    (H,E,U,d)(T+)=(25,2,50,9).                       (1)

The finite certificate is independently checked by subset DP and by
permutation/path-cover enumeration. This packet claims the **kernel
derivation and extension law**, not that the values23 and25 or every
isomorphic pair were previously unknown. A targeted search found no earlier
23/25 pair with this stated rooted signature; the exhaustive earlier
small-core compiler could already contain these examples implicitly.

The physical maps inherited from
`05-knowledge/results/h13_paired_descent_20261111/dynamics/PORTABLE_GAPS.md`
are

    O(h,e,u,d)=(h,0,2u,d),
    Z(h,e,u,d)=(h+u−e,d,2u−d,d),
    J(h,e,u,d)=(h+u−e,h−e,h+u,h).

Hence every common O/Z continuation preserves the difference(2,0,0,0).
After a first J, subsequent common continuation W has H-gap2H(WK1), by the
same linearity argument as the inherited portable-gap theorem. A pure-H
difference is sufficient for these operations; no two-sided marked-graph
equivalence is being asserted.

## Complete arbitrary-flank family and native widths divisible by6

Use the fixed-core construction and proved nine-coordinate formula from
`05-knowledge/results/session251_creative_completion_20261111/ternary/PROOF.md`,
§§2–3 (implementation `ternary/compiler.py`): an ordinal transitive left word,
then Q±, then an ordinal transitive right word, plus the root. The left root
out-neighbour word is the binary expansion of A (empty for A=0); the right
word is the bitwise complement of the reversed binary expansion of B.
The cut on Q± stays21. Put α=A mod2, β=B mod2, and

    g(0)=g(1)=1,
    g(2q+1)=g(q),     g(2q)=g(q)+q+(q mod2).

The exact core signatures in the inherited compiler's order
`(d,K+,K−,I+,I−,O+,O−,J+,J−)` are

    Q−: (9,75,−9,22,−4,41,−7,12,2),
    Q+: (9,75,−9,25,−5,41,−7,14,2).

These are sums of actual directed-path weights; `switch_check.py` recomputes
every entry. Substituting them in the fixed-core formula gives, for all
integers A,B≥0,

    H−=9(g(A)+g(B)−1)+84AB+22A+4α+41B+7β+14,
    H+=H−+3A+α+2.                                    (2)

Both constructions have exactly the same continuation sidecars. Write
a=bitlength(A), b=bitlength(B), with empty length0, and reverse exactly b
bits to obtain rev_b(B). Then

    d=9,
    E=2 if B=0, otherwise9,
    U=(84A+50)2^b−9rev_b(B),
    K=84·2^(a+b).                                    (3)

Here K is the ordered two-path-cover count of the root-deleted tournament,
allowing either path empty. To see (3) directly, each singleton adjoined in
ordinal position doubles K. Prepending a left vertex beaten by the root
adds the previous K to U; a vertex beating the root adds0. Reading the left
bits from right to left therefore gives U=50+84A. Right attachments have
the O/Z updates shown above, giving the displayed reversal sum. The last
right bit is0 whenever B>0, which gives E=d. The core d and its E=2 at an
empty right word follow from(1). These arguments also prove that (2) is a
**pure-H** switch for every A,B, not just at the seed.

All these tournaments are strong. A nonempty left word starts with1, so the
root reaches its first vertex and then the entire left chain and core. A
nonempty right word ends with0, so its last vertex returns to the root. The
mixed core cut supplies both directions when either flank is empty.

For A even the width in(2) is2 mod6; for A odd it is0 mod6. In particular
for every k≥1, take A=2k−1,B=0 to obtain the actual strong pair

    h_k=44k+9g(k−1)−4,       h_k+6k,
    common(E,U,d)=(2,168k−34,9).                       (4)

The first pairs are49/55,93/105,155/173. No ordinal triangle multiplication
is required. Together with the inherited H13 width `3A−α+2`, which supplies
widths4 mod6 as well as2 mod6, these native two-core constructions realize
**every positive even width** with matched E,U,d. This is pair-gap coverage,
not a theorem that their absolute counts cover every permitted odd target.

For comparison with the earlier literal6853 source, width860 here has
A=286,B=0 and baseline7683, upper8543. It improves the inherited H13
baseline15449 but still cannot be an induced switch inside a source with
only6853 HPs. This is one exact boundary, not a new census.

## A target-preserving infinite paired rebase onto C3

At A=0 write the lower family as

    N9(B)=9g(B)+41B+7(B mod2)+14.                     (5)

The upper member is always N9(B)+2. For B>0 their common sidecars are
`(E,U,d)=(9,50·2^b−9rev_b(B),9)`. The lower residue satisfies

    N9(B) ≡ floor(B/2)+2 (mod3),

so this one unit-switch family occupies all three residue phases as B varies.

The inherited C3 pair is

    N3(C)=3g(C)+12C+10,       N3(C)+2,                (6)

realized respectively by C3 cut3/right address2C+1 and cut1/right
address4C+3. It has common `(E,U,d)=(3,15·2^c−3rev_c(C),3)` for C>0.
This is the same13/15-seeded family used in the earlier residue-phase work.

**Paired rebase theorem.** Define

    B0=16, C0=58,
    B_(j+1)=4096 B_j+1824,
    C_(j+1)=4096 C_j+3162.                           (7)

Then for every j≥0, `N9(B_j)=N3(C_j)`. Consequently both members of the
unit pair have actual C3 witnesses at the **same target counts**. The first
lower targets are823,3356485,13748120659; add2 for the upper member.

*Proof.* All B_j,C_j are even. For a length-t suffix s, direct iteration of
the g recurrence gives

    g(2^t x+s)=g(x)+(2^t−1−rev_t(s))x
       +[s<2^(t−1)](x mod2)+g(s)−1.                (8)

At t=12, the suffix data are

    s=1824: rev=78,   g(s)−1=1790,
    t=3162: rev=1443, g(t)−1=2050.

Thus for even inputs,

    g(B′)=g(B)+4017B+1790,
    g(C′)=g(C)+2652C+2050.

Keep the two exact defect coordinates

    L=218B−61C+50,       D=N9(B)−N3(C).

Substitution gives the finite recursive law

    L′=4096L,           D′=D+936L.                  (9)

Initially L=0 and D=0 (`g(16)=17`, `g(58)=39`, both targets823). Therefore
both defects stay zero forever. Formula(8) has no finite-length restriction,
so this proves the whole infinite family. QED.

The common lower targets also have the closed form

    N_j=(39149312·4096^j+169513)/47775 + (234/35)j.

Integrality follows either from(7) or(9); the rational expression is only a
compact closed form, not the definition of the witnesses.

This is an actual repeatable two-coordinate exact-target gluing rule. It
does not rebase every arbitrary B: its invariant line L=0 and initial
defect D=0 are essential. Those are concrete sufficient guards, not presumed
necessary conditions for C3 realizability.

## The C3 rebase also reduces every positive continuation cost

For(7), the word lengths are b_j=5+12j and c_j=6+12j. Bit reversal gives

    rev(B_j)=1+2496(4096^j−1)/4095,
    rev(C_j)=23+92352(4096^j−1)/4095.

Therefore, at either matched H target,

    U9−U3=640·4096^j+60+254592(4096^j−1)/4095 ≥700.

Also E9−E3=d9−d3=6. In the positive continuation coordinates
`q=(H,U,U−d,H+U−E)`, the difference is

    (0, ΔU, ΔU−6, ΔU−6),

which is coordinatewise nonnegative. Direct substitution in the physical
maps gives, for q=(h,u,v,t),

    Oq=(h,2u,u+v,h+2u),
    Zq=(t,u+v,2v,t+2v),
    Jq=(t,h+u,u,h+2u).

All three maps have nonnegative coefficients. Thus the C3 replacement weakly
lowers every common future H while
preserving both original target counts and their unit-switch choice. This
is a useful reduction of the marked response as well as a count coincidence.
It supplies no new indegree-one ready source and is not labelled as such.

## Reproduction

    python3 -B switch_check.py --write
    python3 -B switch_check.py
    python3 -O -B switch_check.py

The receipt records all eight completions, both path signatures, direct
rooted DP and independent permutation/path-cover counts, a bounded literal
flank grid, the initial four actual C3/H9 target witnesses at823/825, and
small exact checks of the arbitrary-input defect law(9). The all-j result
is algebraic, not extrapolated from those checks.
