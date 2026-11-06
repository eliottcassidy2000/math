# The same obstruction for the actual source Hamilton count

Status: **PROVED, independently audited in the stated scope.** This extends the main theorem
from its closing readout F to the actual Hamilton-path count D of the source.
It still concerns only the physical source grammar; it is not a claim about
arbitrary vertex additions or arbitrary tournament surgery.

Use the notation and exact parent sets A_plus and B_minus of PROOF.md:

    A_plus={6*2^k-1,10*2^k-1,6(2^j+1)*2^k-1:j>=1,k>=0},
    B_minus={2^r+1:r>=1} union {3(2^j+1):j>=1}.

The display for A_plus means that j is used only in the third family.

**Corollary.** Let n>=3 be odd. Choose any actual nontrivial source s with
D(s)=n, allowing arbitrary replacement of an initially supplied source by
another such source. Append any physical O/Z/J word. The resulting source
has D=2n+1 exactly when such a procedure exists with n in A_plus; it has
D=2n-1 exactly when such a procedure exists with n in B_minus. Consequently
both children can be produced, allowing separate parent representatives,
exactly at n in {5,9,17}.

Here 'exactly when such a procedure exists' classifies existence over the
choices of s and appended word; it does not say every source at that parent
works.

**Proof.** O preserves D. Both Z and J replace D by

    F=D-E+U >= D-(D-1)/2+D=(3D+1)/2,

using the physical inequalities U>=D and E<=(D-1)/2. After two D-changing
letters, regardless of intervening O letters,

    D_final >= (9n+5)/4 > 2n+1.

Thus exactly one Z or J is possible in a successful word. Any preceding O
letters can be included in the arbitrarily chosen parent representative;
any following O letters preserve the output. It remains to classify an
actual source s with D(s)=n and F(s)=2n±1.

The source s cannot be the seed because n>=3. It cannot end in O: such an
output has E=0 and U>=2D, giving F>=3n>2n+1. Thus in any physical recipe its
last letter is Z or J. Write s=t.c with c in {Z,J}, where this t is an actual
nontrivial source, not necessarily the initial seed. Then

    F(t)=D(s)=n,
    F(s)=F(t.c)=2n±1.

This is exactly the one-letter classification of the main theorem. Conversely,
for every one-letter witness t.c in that theorem, s=t.c has D(s)=F(t)=n;
appending one further Z or J gives D_new=F(s)=2n±1. The main theorem's exact
intersection is {5,9,17}. Its minus parent 3 is also valid here: the triangle
source has D=3 and F=5, while no plus child is present. QED.

The distinction of interfaces does not create an escape: even when the
parent is its actual source-graph path count rather than a prospective
closing count, the same complete binary-parent sets occur. A successful
cofinite construction must use an operation outside this arbitrary-rebase
plus append-word scheme, or a different inductive step from these binary
children.

The companion `verify_d_interface.py` exhausts all relevant physical profiles
for D<=3333 using U<=5000. Indeed U=D+E±1 implies U<=(3D+1)/2<=5000 there.
It compares both complete parent sets, checks the growth inequality, and
counts literal source graphs at the three double-child parents. Finite
controls support the proof; no stabilization is inferred.
