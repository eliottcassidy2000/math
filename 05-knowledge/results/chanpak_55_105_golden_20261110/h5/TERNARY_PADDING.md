# One active singleton suffices below the fivefold threshold

Status: PROVED IN THIS NOTE, using the incoming Session251 padding theorem
and the fixed-core compiler. FINITE-EXACT controls accompany the proof.
This is a source-normalization and branch-support theorem, not a proof that
either ternary branch exists for every parent count.

## Inheritance and exact scope

The incoming
`05-knowledge/results/session251_creative_completion_20261111/geometry/PADDING_CUT_NORMAL_FORM.md`
at Git commit `5bddb480c6` proves that, below five times the source H,
singleton attachment bits have at most one inversion. Its isolated-triangle
rigidity strengthens this: any inversion gives a separate new triangle
component and exactly3h. Consequently a low-ratio child whose new-vertex
SCC contains the old nontrivial core has a monotone singleton word, and
at most one active flank. That theorem is inherited, not reproved here.

The additional ingredient is `PROOF.md`, Theorem H5.1, together with the
recovered Session197 single-exception identity. Both endpoint-cover sums
are retained. The resulting sharpening bounds the number of active old
singletons by one, and gives an exact common-source preparation rule.

Let C be a fixed nontrivial strong tournament, n=|C|, h=H(C)>=3. A source
is B=TT_a -> C -> TT_b. Add a new vertex x with arbitrary incidences.
On C let O be its out-set and I the complementary in-set. Use the fixed
core coefficients from `PROOF.md`: q=Q_C(2), e is the ordered-cover
end sum on I, s the ordered-cover start sum on O, u the Hamilton-start
sum on I, v the Hamilton-end sum on O, and t=H(C+x_O).

The left incidence integer is A; the complement-reversed right integer
is Bnum (to distinguish it from the source B). The exact count is

    F(A,Bnum)=q A Bnum+e A+s Bnum+h(g(A)+g(Bnum)-2)
              +t+u(A mod2)+v(Bnum mod2).                 (1)

## 1. All-order sharpening

**Theorem TP1.** Suppose H(B+x)<5h and the SCC of x contains C. Delete
only the exterior singleton SCCs of the child. The remaining source is
one of C, TT_1 -> C, or C -> TT_1; the remaining child is strong. In
particular there is at most one active old singleton vertex.

Proof. By the inherited isolated-triangle rigidity, the source's singleton
word has no inversion. The left word is therefore 0^r1^k, and the right
word is 0^l1^s. They cannot have k,l both positive, since this would give
an inversion across C. Initial left zeros and final right ones remain
exterior source/sink SCCs in the child and can be removed without changing
either H count. Thus the normalized integers are

    (A,Bnum)=(2^k-1,0) or (0,2^k-1), k>=0.              (2)

For a left-active word with k>=1, g(2^k-1)=1, so (1) becomes

    F=t+e(2^k-1)+u.                                    (3)

For every vertex i, its individual ordered-cover end response is at least
h: split every Hamilton word of C immediately after i; concatenation
recovers the word and makes the injection explicit. Hence e>=|I|h.
Ordinary insertion monotonicity gives t>=h.

If I has at least two elements and k>=2, (3) gives F>=h+3(2h)=7h,
contradicting the strict threshold. If I={b}, the recovered Session197
single-exception identity is

    t+u=h+e.                                          (4)

For clarity, (4) also has a direct path proof: a path beginning at x
contributes h-u; a path ending at x contributes end_b(C); and a path
with x interior contributes e-end_b(C), because its first piece must end
at b and its second piece cannot contain b. Their sum is t=h-u+e.
Thus (3) is h+2^k e>=5h when k>=2. The right-active argument is the
converse, using s>=|O|h and t+v=h+s when O is a singleton.

It remains to exclude uniform core cuts. If O is empty, the left-active
count is h+(q-h)(2^k-1)+h>=q+h. Cutting each Hamilton word in n+1
positions gives q>=(n+1)h, so this is at least5h because n>=3. Its
right-active count is h and x does not join C. The full cut is the
converse. With no active flank, a uniform cut also leaves x outside C.
Thus the relevant core cut is proper and nonempty, and k is at most1.
Such a cut makes C+x strong; including its one active singleton preserves
strongness by the displayed orientations. QED.

The threshold is sharp. For C=C3, let I be one core vertex, take two
left singleton vertices both beaten by x, and take no right vertices.
Then A=3, e=h=3, and (4) gives F=3+4*3=15=5h. The whole child is strong
and contains both old singleton vertices. This is an actual equality
case, not merely a failure of the proof's estimates.

## 2. Exact ternary branch support and common-source preparation

**Corollary TP2.** Fix C and a,b>=0. For N in {3h+2,3h+4}, an attachment
of N from TT_a -> C -> TT_b exists exactly when some proper nonempty core
cut O satisfies at least one enabled equality:

    N=t;               always enabled,
    N=t+e+u;           enabled iff a>=1,
    N=t+s+v;           enabled iff b>=1.                (5)

Proof. Since h is odd, gcd(h,3h+2)=gcd(h,3h+4)=1. If C lay outside
x's SCC in the child, its count h would divide N by SCC multiplicativity,
which is impossible. Also N<5h for h>=3. Apply TP1 and (1) at the three
integer pairs (0,0),(1,0),(0,1) to obtain necessity. Conversely each of
these three cuts is an actual graph construction. To restore every old
singleton, put zeros at unused left positions and ones at unused right
positions. The active left vertex, if used, is the last left vertex;
the active right vertex is the first right vertex. This preserves the
given old source and supplies sufficiency. QED.

Thus the entire ternary branch support of a fixed core depends on the
padding lengths only through the two predicates a>0 and b>0. At most
3(2^n-2) exact arithmetic candidates settle both branches for that core.
The omitted uniform cuts cannot contribute either target.

**Corollary TP3 (common-source preparation).** The single prepared source

    TT_1 -> C -> TT_1                                 (6)

enables every ternary branch attainable from any singleton padding of C.
If both branches are individually possible using this same C, both can
be constructed from the common source (6), using different attachment
cuts when necessary. More generally, starting from any already named
padding of C, add one new universal source and one new universal sink.
All old names and arcs survive, H remains h, and the full branch support
for C is enabled. No identification of different strong cores is used.

This statement would be false for arbitrary low targets: N=3h can come
from an isolated triangle using two adjacent old singleton vertices and
x, so its existence depends on more than positivity of each flank length.
The coprime ternary targets and core-containing condition are essential.

## 3. Prime and seven-times-prime lanes

Every realized parent count h=p or h=7p, with p an odd prime and h not
7 or21, has exactly one nontrivial strong component. Indeed component
counts are odd positive factors of h, and every nontrivial strong
component has count at least3. For prime h this permits only one such
factor. For h=7p, two such factors would require one to be7 (including
the p=7 case h=49), excluded by inherited THM-1370,
`01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`.
All other components are singletons. Therefore TP2 and TP3 cover exactly
these parent lanes, for every actual strong-core realizer.

The remaining existential problem is concrete: for each desired parent
count, find a strong core C and a proper cut giving the requested equality
in (5). The two branches may need different C. This result eliminates
arbitrarily long neutral padding as a source of new ternary support for
a fixed C; it does not eliminate the need to switch actual cores.

The normalization preserves the scalar parent and child counts. For
later joint-port queries, retain the complete original names, their
positions, and all incidences; formally deleted exterior padding is not
permission to forget a required old port. TP3 preserves those names by
adding vertices, but does not prove closure of the joint17-coordinate
observer or its F49/F343 phase data under repeated ternary returns.

## Controls

`ternary_padding_check.py` checks the finite-cut criterion against the
complete fixed-core inverse on four specified cores: C3, Q4, and the two
displayed H9 order5 cores. It also checks all masks with A,B in0..7,
the named singleton-exception identity, the per-vertex cover lower bounds,
and the strong equality example at15. Graph witnesses for positive
branches are recounted by subset DP; the equality graph is also counted
by literal permutations. This is corroboration of TP1–TP3, not the source
of their all-order quantifiers. No new H9 classification is asserted.

The controls also display why actual core choice remains necessary. The
H9 core (2,4,9,19,7) supports both29 and31 after suitable padding; the
other specified H9 core (2,4,25,3,11) supports29 but not31, even with
arbitrary singleton padding. The negative is a finite proper-cut check
made all-order by TP2, independently agreeing with the complete inverse.
Common-source preparation synchronizes padding for one core; it cannot
repair this difference between the cores.

    python3 05-knowledge/results/chanpak_55_105_golden_20261110/h5/ternary_padding_check.py
    python3 -O 05-knowledge/results/chanpak_55_105_golden_20261110/h5/ternary_padding_check.py
