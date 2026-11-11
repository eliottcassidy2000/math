# Native source improvements and the exact periodic-return gate

**2026-11-11. PROVED IN NOTE:** a root wrapper that converts the inherited
H37 arc surgery into a same-H improvement for every future O/Z/J program;
an order-reversing asymmetric birth between the H13 and H25 paired
families; and a necessary-and-sufficient periodic-return criterion.
**FINITE-EXACT:** the explicitly bounded word and seed searches below.
**OPEN:** a compatible repeated return, joint target coverage, and Chan–Pak.
There is no claim of a completed odd-spectrum proof or literature priority.

## 1. Inheritance and live concepts

The closest proved mechanism is the actual H13 rooted switch in
[the preceding history continuation](../../five_core_web_20261111/history/continuation/PROOF.md),
extended to arbitrary left addresses in
[PORTABLE_GAPS.md](../dynamics/PORTABLE_GAPS.md). Its golden P4 triangle
overlap explains where the family arose; its four-coordinate rooted
response is what permits the constructions here. The canonical hostile
is a pair of scalar witnesses whose endpoint data do not agree. The
corrected near miss is a family of portable differences without a way
to install the desired difference at the required absolute baseline.
The least-used inherited sidecar is the full path-colour polynomial Q,
especially its value at2, before reducing to a modular substitution jet.

The live board is: paired endpoint responses; same-count source changes;
positive future-cost order; binary recurrence drift; and exact seed
compatibility. The three constructions below each retain a different part
of this board. None is silently substituted for all five requirements.

Use the state `(H,E,U,d)` from the sibling notes: H is the Hamilton path
count, E counts such paths ending at the marked root, U counts ordered
spanning two-path covers whose distinguished path starts at the root
(the other path may be empty), and d is H after deleting the root.
Write

    q=(H,U,U-d,H+U-E).

The [future-dominance theorem](../dynamics/PROOF.md) proves that componentwise
q order is exactly comparison of H under every common O/Z/J program.
Its physical operators, in q coordinates `(h,u,v,t)`, are

    O(h,u,v,t)=(h,2u,u+v,h+2u),
    Z(h,u,v,t)=(t,u+v,2v,t+2v),
    J(h,u,v,t)=(t,h+u,u,h+2u).

## 2. An inherited native arc surgery becomes a rooted cost reset

Session254, [chan_pak/PROOF.md, section6.1](../../session254_four_frontiers_20261111/chan_pak/PROOF.md),
supplies the literal one-arc reversal

    A =(56,9,51,36,10,18),
    A*=(56,9,51, 4,10,26).

Only the arc3→5 changes. The inherited exact polynomials are

    Q_A(t)  =t^6+14t^4+22t²,
    Q_A*(t) =t^6+12t^4+24t²,
    Q_A*(t)-Q_A(t)=2t²(1-t²).

Here `Q_S(t)=sum_k p_S(k)(t)_k`, where p counts unordered directed-path
covers and `(t)_k` is the falling factorial. In particular H=Q(1), while
Q(2) counts ordered spanning covers by two paths with empty slots allowed.
Thus H(A)=H(A*)=37, and Q_A(2)=376, Q_A*(2)=352.

**Root-wrapper lemma.** For any nonempty tournament S, prepend a new
global-source root r. Its exact rooted state is

    (H,E,U,d)=(H(S),0,Q_S(2),H(S)).                 (1)

Proof. Every full path starts at r and then is a full path of S; no full
path ends at r. Deleting r recovers S. Removing r from its distinguished
path gives an ordered two-path cover of S, with that slot empty if the
original path was just r. Conversely prepend r to the first slot of any
such cover. This is a bijection proving U=Q_S(2). QED.

Apply the same wrapper to A and A*. The actual seven-vertex graphs have
rows `(A,63)` and `(A*,63)`, retaining root6, and states

    (37,0,376,37),             (37,0,352,37).

Their q difference is

    q_old-q_new=(0,24,24,24).                     (2)

**Theorem1.** This literal count-neutral arc reversal preserves the root,
all seven vertices, and every other arc. It preserves H and strictly
reduces H after every common O/Z/J program containing at least one Z or J.
Pure O programs preserve the common H.

Proof. The difference in (2) has nonnegative coordinates. After any number
of O steps its final coordinate remains positive, so the first Z or J
makes its H coordinate positive. Every subsequent operator preserves
positivity of the relevant coordinates and hence the strict H comparison.
Equivalently apply the future-dominance theorem and the displayed maps.
QED.

One Z gives the concrete eight-vertex graphs

    (184,137,179,164,138,146,63,64), H=413;
    (184,137,179,132,138,154,63,64), H=389.

Both are strong: if z is the last vertex, then `r→s→z→r` for every
old vertex s. This also explains why the nonstrong source wrapper is
compatible with later strong constructions.

**Arbitrary ordinal-context extension.** For any tournament R, perform
the arc reversal in `S=A▷R`, where all A vertices beat all R vertices,
and then retain the global-source root. The ordinal identity
`Q_(A▷R)=Q_A Q_R` gives

    H_old=H_new=37H(R),
    U_old-U_new=24Q_R(2),
    q_old-q_new=(0,24Q_R(2),24Q_R(2),24Q_R(2)).      (3)

The same holds with R on the other side or with further ordinal factors.
There is **no restriction on H(R) modulo3**: that restriction belongs to
Session254's curvature normalization, while (3) is an integer identity.
For example R=C3 has H3 and Q_R(2)=12; the rooted states are
`(111,0,4512,111)` and `(111,0,4224,111)`, with improvement288.

**Arbitrary rooted carrier.** A second extension retains a root already
present in an arbitrary rooted tournament X. Append A or A* as an ordinal
suffix to X. The two exact states are

    (37H_X,0,376U_X,37d_X),
    (37H_X,0,352U_X,37d_X),                       (3a)

so the q improvement is `(0,24U_X,24U_X,24U_X)`. To see the U product,
split either ordered spanning two-cover of the ordinal sum at the factor
boundary. The X part is a root-start ordered two-cover counted by U_X;
the A part is an arbitrary ordered two-cover counted by Q_A(2). Empty
slots are allowed, and concatenation in each slot is the inverse map.
This installs the native reset on any rooted carrier after multiplying
its count by37. It does not preserve the unextended count H_X.

The preserved predicate is exact H with an actual retained root. The
changed sidecar is U. This provides a native patch inside a specified
family of sources; it does not assert that an arbitrary source of count
37H(R) contains that patch or can be replaced without changing its other
rooted coordinates. After Z the totals themselves differ, as they should
for a downstream cost improvement.

## 3. The paired birth can reverse a portable difference

The sibling [paired-reset note](../modular_cycles/PROOF.md) discovered
the following physical asymmetric birth. The H13 core and its two root
cuts are

    Q=(18,8,3,5,14),       cuts6 and3.

Insert one vertex y after Q and before the right transitive flank. Let
y lose to the retained root and all left-flank vertices, and beat every
right-flank vertex. On the lower side y beats Q vertices1,2,3; on the
upper side it beats only Q vertex3. The resulting cores and root cuts are

    C-=(50,8,3,5,46,14),     cut38;
    C+=(50,40,35,5,46,8),    cut35.                 (4)

Deleting y gives the named original graph on both sides, at every pair
of flank addresses A,B. Thus this is an actual pair of one-vertex births,
not just equality between compiled formulas.

Let α=A mod2, β=B mod2. The inherited H13 difference is

    H_old,+-H_old,-=3A-α+2.                        (5)

The two signatures in the sibling note give, for (4),

    H_new,- =25(g(A)+g(B)-1)+304AB
              +107A+15α+122B+16β+48,
    H_new,+-H_new,-=-3A-α+2β.                     (6)

For B>0, write b=bitlength(B), r=rev_b(B). The two new rooted sidecars
agree exactly and are

    E=d=25,       U=(304A+147)2^b-25r.             (7)

For completeness, the weighted splitting proof used for the old H13
pair applies verbatim with the new common signature entries
`d=25,k+=279,o+=122`. The root-alone path contributes25; marked starts
in the left flank contribute `(25+279)A2^b`; selected core starts
contribute122·2^b; and right-flank starts contribute
`25(2^b-1-r)`. Their sum is (7). The last right-flank vertex beats the
root, so E=d=25. At B=0 the two E values are9 and7, respectively;
the restriction B>0 must therefore be retained in the pure-H claim.

**Theorem2.** If A is odd and B is positive even, the asymmetric birth
changes a pure-H pair of width `3A+1` into a pure-H pair of the same
width with the named pair order reversed:

    Δ_old=(3A+1,0,0,0),
    Δ_new=(-(3A+1),0,0,0).                         (8)

Proof. Equations(5)-(7), together with the common H13 sidecars proved in
PORTABLE_GAPS.md. QED.

For example A1,B2 gives the actual paired transformation

    (385,13,583,13) → (1097,25,1779,25),
    (389,13,583,13) → (1093,25,1779,25).

This order reversal is impossible under a common fixed-root O/Z word
starting from a positive pure-H difference; those words preserve the
same difference. The birth therefore supplies a real operation outside
that restricted continuation language. It changes the absolute baseline
and is not an unproved operation for subtracting a chosen gap from an
arbitrary target. At A0,B odd it reduces to the sibling unit-pair
construction, including69/71→211/213.

**Corollary: a fixed H25 deletion core realizes every even pure-H width.**
For any positive even w, choose

|w modulo6|A|B|
|---:|---:|---:|
|0|w/3|2|
|2|(w+1)/3|1|
|4|(w-1)/3|2|

Then (6) equals `-w`, and (7) makes the difference pure H. Reversing
the pair order gives width+w. Thus a fixed deletion count25 suffices
for every even width; ordinal multiplication is unnecessary in this
family. Only the flanks grow with the chosen width. This still specifies
its own absolute baseline rather than freely choosing that baseline.

For later composition, the two deletion graphs have equal Q(2) as well:
`Q_C-(2)=Q_C+(2)=d+k+=304`. With flank lengths a,b, deleting the root
leaves `TT_a▷C±▷TT_b`, so its Q(2) is `304·2^(a+b)` on both sides.
The deletion graphs themselves differ. This distinguishes the exact
sidecar retained by the pair from a stronger and unnecessary claim of
identical underlying graphs.

## 4. Exact mixed-word drift and the seed compatibility equations

The next positive question is whether two different source languages can
meet repeatedly, so that a same-count cost reset can be reused. Merely
finding many isolated equalities does not establish that. The following
criterion isolates the three exact obligations.

The inherited binary function satisfies

    g(0)=g(1)=1,
    g(2n+1)=g(n),
    g(2n)=g(n)+n+(n mod2).

Fix a word w of length ℓ and integer value r, including leading zeros.
Set `M=2^ℓ`, `Q=M-1`, β=r mod2, and

    a=Q-rev_ℓ(r),
    b=1 if w starts with0, otherwise0,
    c=g(r)-1.

Induction on the appended bits gives the exact identity

    g(Mn+r)=g(n)+an+b(n mod2)+c.                   (9)

For `n_k=value(prefix p followed by w repeated k times)`, k>=1,

    n_k=(p+r/Q)M^k-r/Q.

Put `D=bβ+c-ar/Q`. Summing (9), with the first prefix parity kept
separately, gives

    g(n_k)=a(p+r/Q)M^k/Q + Dk
             +g(p)-a(p+r/Q)/Q+b((p mod2)-β).       (10)

No asymptotic error or finite truncation occurs. In particular for a
scalar response `F(n)=dg(n)+on+e(n mod2)+h`,

    F(n_k)=C M^k+dDk+E,
    C=(da/Q+o)(p+r/Q),
    C+E=F(p)+(db-e)((p mod2)-β).                  (11)

**Theorem3 (periodic return criterion).** For two responses using words
of the same length ℓ, their outputs agree for every k>=1 if and only if
all three equations hold:

    C_1=C_2;                                    growth match
    d_1D_1=d_2D_2;                              drift match
    F_1(p_1)+(d_1b_1-e_1)(parity(p_1)-β_1)
      =F_2(p_2)+(d_2b_2-e_2)(parity(p_2)-β_2).   adjusted seed match

Proof. Subtract (11). The sequences M^k, k, and1 are linearly independent
for M>1. Necessity follows by successive differences, and sufficiency
by substitution. QED. Actual compiler sources and nonnegative addresses
make every solution a physical repeated equality; formal coefficients
without those sources do not.

## 5. A positive 13/25 drift cancellation, with the remaining gate explicit

The common-count211 source pair from the sibling note has lower responses

    F_old(C)=13g(C)+334C+198,
    F_new(D)=25g(D)+244D+186.                      (12)

Both associated upper responses are exactly2 higher. Thus any repeated
meeting of (12) is simultaneously a meeting of the full paired targets.

The mixed words

    w_new=00000000000101,       r_new=5,
    w_old=00001001110101,       r_old=629,

both have length14 and satisfy

    D_new=18434/16383,
    D_old=35450/16383,
    25D_new=13D_old=460850/16383.                  (13)

The strictly positive mixed drifts therefore cancel between the two
source languages. A positive drift inside one language is not itself
an obstruction to a paired return. This is the surviving useful signal
from the periodic investigation.

The bounded exhaustive word probe `period14.py` enumerates all mixed
words at each length2 through14. There are no such weighted drift pairs
through length13 and1568 ordered pairs at length14 (28 choices on the
new side and56 on the old side). This is a finite first-occurrence
statement, not a general role for the number14.

The [synchronization dataset](../modular_cycles/synchronization.json)
contains all227 meetings of (12) for C,D in1..32767. For example
`C=79,D=102` gives26649/26651. Testing those227 seeds plus `(0,0)`
against all1568 drift-compatible word pairs finds no solution to the
remaining two equations in Theorem3. The parity correction is handled
exactly: since e=0 in (12), equality of the corrections13b_old·ε_old
and25b_new·ε_new forces both to vanish. This is the complete tested
universe, not an assertion that a periodic seed does not exist.

The precise next construction obligation is a physically realized seed
meeting the growth and adjusted-seed equations for a drift-compatible
word pair, or a source-changing operation that compensates one of these
equations. The H37 reset changes future U cost while keeping current H;
the order-reversing paired birth changes the sign of a portable gap;
and the sibling all-even-width theorem supplies every possible gap
magnitude somewhere. Combining these at the required same baseline is
still an open compatibility problem.

## 6. Controls and narrower discarded probes

Run from the repository root:

    python3 -B 05-knowledge/results/h13_paired_descent_20261111/construction/check_positive.py
    python3 -B 05-knowledge/results/h13_paired_descent_20261111/construction/period14.py

`positive_receipt.json` records literal endpoint/subset-DP controls for
the new H37 rooted claims, the nonunit R=C3 extension, six actual paired
births (including deleting the inserted vertex), a nontrivial rooted
carrier, and40 exact loop-formula
controls. The all-address and arbitrary-context statements are proved
above; these small tests are controls, not their justification.

Two narrower exploratory receipts are retained with their scripts:

- `intersections.json`: no positive all-k direct H13/H25 family match in
  `p digit^k suffix`, prefixes0..512, suffix length<=5, k>=1, odd new
  addresses, with empty prefix allowed only for a one loop. The search
  imposes exact coefficient equality and a U decrease.
- `ray_rebases.json`: the lower ray77·2^m-33 had no representation in the
  stated bounded fixed-left five-core template. The upper ray
 77·2^m-31 recovered the sibling history lane's H9/H5 source rebases.
  Those are inherited discoveries, not a new result of this probe.

These negative finite observations motivate the mixed-word and native
source-change operations; they do not constrain the full tournament
spectrum. This note proves new constructions and an exact return gate,
while leaving the final all-odd target-coverage implication open.
