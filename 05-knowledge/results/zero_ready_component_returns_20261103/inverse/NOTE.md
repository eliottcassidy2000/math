# A component reader with an exact target inverse

**PROVED operation laws and terminating inverse for each specified seed;
FINITE-EXACT three-round target acceptance.** This is not an all-target
acceptance theorem or an unbounded ternary induction.

The starting zero-ready source has physical root39, sole in-port31, prefix
50->51, and marked strong component B with (h,F)=(2195,42153). The inherited
[zero-ready decomposition](../../triangle_return_170807_20261103/ZERO_READY.md)
evaluates its count as h+4F=170807. Merely resizing that prefix or repeating
the old fixed-module rewrite cannot give another ternary return.

## 1. A retained-mark binary reader changes the component grammar

Fix a seed B at w. Append a vertex for each bit: bit1 is a global sink;
bit0 beats only w and loses to every other existing vertex. The old mark
and all old relations survive these reader births. Intermediate components
can be weak. If the seed has at least two vertices, w reaches every seed
vertex, and the last bit is0, the completed component is strong: every
non-w vertex reaches the last newborn, which reaches w, and w reaches the
old seed and all appended vertices. Singleton seeds require a separate
strongness test; a singleton followed by0 is transitive2.

Use component fields

    (h,e,t,F,p,q,P2,S3)
      =(H,E_w(1),S_w(2),E_w(2),P_(B-w)(2),H(B-w),P_B(2),S_w(3)).

All path bins are ordered and may be empty, except the designated first bin
in E or S. The exact update laws are

    1: (h,0,2t,F,2p,q,2P2,3S3),
    0: (h+t-e,q,2t-q,F+p,2p,q,2(P2-F+S3),3S3-p).

The [independent interface proof](../interface/NOTE.md#2-two-genuine-finite-scalar-closures-with-distinct-operation-scopes)
derives these identities from the inherited polarized path-bin calculus and
checks them by literal small graphs. The first six fields suffice for the
count inverse below. The last two retain the zero-ready wrapper's golden
phase. A named spare needs additional joint data or the full graph.

## 2. The inverse eliminates the unbounded word-length search

For a nonempty word ending0, let

    r = sum_(i:bit_i=0) 2^i,                    i starts at0.

Then r>=1 and its binary length is exactly the word length. Conversely each
positive r gives a unique word by reading its bits from low to high and
complementing them. Write a(r) for the Hamilton count of the ordinary
one-apex chain word so obtained. The inherited binary reader gives

    a(1)=1,
    a(2u+1)=a(u),
    a(2u)=a(u)+u+(u mod2),                       u>=1.

These are also direct graph identities: a leading0 is a global source;
a leading1 adds precisely the cycles with that first chain vertex. The
one-apex endpoint count is1+r, so

    1 <= a(r) <= r+1.

One may prove the upper bound without a cycle formula: every Hamilton word
injects into a two-bin cover ending at the apex by cutting immediately
after the apex. We use the inherited exact one-apex evaluator from
[Session179](../../session179_niche_methods_20261104/spectrum/REPORT.md),
and independently compare its streaming and low-bit recurrences.

For an arbitrary actual seed, summing the six-field reader yields

    h_out = h + q(a(r)-1) + (t-q)r + (q-e)(r mod2),
    F_out = F + p r.                                      (I1)

The first zero's correction q-e is essential. The formula also follows by
induction using the displayed recurrence of a. For a transitive prefix
clock s>=1 put M=2^s, K=t-q+Mp, d=q-e and D=N-h-MF. The exact target
equation h_out+M F_out=N is therefore

    K r + q(a(r)-1) + d(r mod2) = D.                      (I2)

Here0<=e<=q, t>=q, p>0, q>0. The endpoint bound follows by deleting w
from an ending-at-w Hamilton path; t>=q uses a singleton first bin w.
Consequently K>0, d>=0 and every possible r lies in the finite interval

    max(1, ceil((D-d)/(K+q))) <= r <= floor(D/K).           (I3)

For each integer in this interval, reconstruct its one word and check I2.
This is necessary and sufficient for the declared reader family. Include
the empty word separately when the seed is strong. Since both h_out and
F_out are nondecreasing, only clocks with h+2^s F<=N need examination.
Thus the inverse covers every word length and every clock without a
sampled length limit. Its acceptance for every requested N remains OPEN.

**Witness-size bound (PROVED).** Let the seed have m vertices and let the
nonempty reader word have length L and end in 0. Each of the 2^(m-1)
ordered subset splits of B-w supplies at least one two-bin cover, so
p=P_(B-w)(2)>=2^(m-1). Also r>=2^(L-1). Thus

    N = h_out + 2^s F_out >= 2^s p r >= 2^(s+m+L-2).

The full wrapper has n=m+L+s+1 vertices, giving

    n <= floor(log2 N)+3.

Both first-round witnesses below attain this integer bound: they have
21 vertices and their counts lie between 2^18 and 2^19. Thus the additive
constant cannot be reduced uniformly for this grammar.

This applies to every accepted nonempty word, even if the seed bank is
enlarged. The empty-word case is separate. The bound proves neither
acceptance nor a polylogarithmic running time: the inverse integer interval
can still contain many candidates. Strongness is checked independently.

## 3. The first two named component replacements

The +4 branch uses the ordinary one-apex component word

    10101001010010100,
    (h,F)=(72701,109931),       h+4F=512425.

Its apex keeps physical name31. The +2 branch needs a different seed:
a directed triangle6->8->55->6, with apex31 beating only6. Append the
reader word

    10111010010110,
    (h,F)=(49955,115617),       h+4F=512423.

In both cases the outer root39 and prefix50->51 are retained, and all
seventeen old names survive. Each endpoint has21 vertices; only four new
permanent names55,...,58 are used. The two full named graphs and their
ordered attachment cuts are in `named_targets.json`.

The +4 alignment changes four old component arcs. The +2 alignment changes
26. The inherited [single-helper compiler](../../adaptive_metagraph_reformulation_20261103/METAGRAPH.md#5-one-temporary-vertex-compiles-any-old-arc-reversal)
gives respectively eight and52 legal native complements, with all temporary
helpers deleted and peak order22. The full paths, relations and deletions
are in `native.json`, with a separate independent replay. No minimum native
distance is claimed. H, strongness and port stars need not be preserved at
native intermediate states; all old arcs outside B agree after each
completed helper gadget.

The original 170807 component is unique-start. Its four-vertex triangle
reader replacement has q=3; the one-apex replacement has q=1. Changing this
deletion component is a real operation, not a re-marking of the old graph.

## 4. Finite repeated returns and the next acceptance obligation

The small bank uses all570 labelled strong seeds of orders3,4,5, reduced to
29 distinct eight-field keys, and35 inherited order-six representatives in
`metagraph_primegame_20261102/low_sources.txt`. Together they give219 keys.
The retained representative graph and root accompany every key. Equality
of these keys is used only for this reader and wrapper, not for arbitrary
joint-port queries or native moves.

Exact inversion supplies both branches for three complete rounds:

    170807
      -> 512423, 512425
      -> 1537271, 1537273, 1537277, 1537279
      -> 4611815, 4611817, 4611821, 4611823,
         4611833, 4611835, 4611839, 4611841.

The six-vertex seed is necessary for one branch **within this stated bank**:
the entire order<=5 bank misses1537279. One order-six witness has rows
(6,44,56,49,35,1), mark2, reader10011101110 and prefix clock4. It gives
(h,F)=(66287,91937), so66287+16*91937=1537279.

The named continuation packet retains parent-to-child maps and full graphs.
This establishes finite reuse of the source-changing procedure, not closure
of the219-key bank at arbitrary depth. To prove an induction, one must prove
that the interval inverse accepts each required next target (possibly with
an enlarged constructive seed grammar), and continue to supply the actual
named transport. Field residues alone do not prove this acceptance.
