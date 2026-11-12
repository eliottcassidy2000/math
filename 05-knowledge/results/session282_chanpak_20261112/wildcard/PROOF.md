# Rational receiver carries localize the Hamilton-count interaction

**Session282. PROVED, independently reviewed.** All-order statements
below are derived from the inherited physical fusion carrier; finite tests do
not supply their quantifiers. Chan–Pak remains OPEN.

## Inheritance and connection contract

The remote source is [THM-2163, radix relation-carry descent](../../../../01-canon/theorems/THM-2163-radix-relation-carry-descent.md)
and [THM-2231, completion and the exact radix clock](../../../../01-canon/theorems/THM-2231-relation-carry-completion-and-exact-radix-clock.md).
Their carry technique is inherited, not claimed new. The physical source is
Session264's rooted-fusion law, reproduced in
[Session281's finite experiment basis](../../session281_chanpak_20261112/wildcard/PROOF.md).
[THM-2259's signed-interaction lesson](../../../../01-canon/theorems/THM-2259-boolean-continuation-hasse-field-and-signed-interaction-dividends.md)
suggests retaining the exact interaction before forgetting a coordinate.
Here the fusion algebra has at most two strict upper jumps, so a quadratic
interaction is exact, unlike a general Boolean interaction expansion.

Map: integer relation words -> pairs of physical module words. Preserved
predicate: exact ordinary receiver equality. Lost quantity: the central count;
the local action below restores its exact difference. The absolute baseline
remains a separate coordinate. A closed finite carry path is not a universal
scalar-target construction. Concept board: carries, physical realizability,
quadratic interaction, signed parity channel, content, absolute baseline.

## 1. Alphabet and inherited matrix law

Fix a finite alphabet of nontrivial rooted tournaments sharing root-deletion
path count h and ordered two-path-cover count Q, with empty cover slots allowed.
Then h is positive odd and Q>=2h. For each letter put

    a=U-h, c=V-h, b=H-h, s=d-h, t=e-h,
    M=[h a s b; 0 Q 0 c; 0 0 0 t; 0 0 0 h].           (1)

Rooted fusion identifies roots and directs every nonroot vertex of the left
factor to every nonroot vertex of the right. Its matrices multiply. The
singleton is separately the identity, not a nontrivial letter in this alphabet.
Deletion of the root from a weighted starting/ending path injects into an
ordered deletion cover; the singleton-root contribution is h. Therefore
`0<=a,c<=Q-h`. The two empty/nonempty covers give Q>=2h. We do not assert
that every digit in this integer box has a physical realization.

For a word of k letters, matrix multiplication gives

    U_word-h^k=sum_i h^(i-1)a_i Q^(k-i),
    V_word-h^k=sum_i Q^(i-1)c_i h^(k-i),
    d_word=h^(k-1)d_k,  e_word=h^(k-1)e_1.             (2)

Word length means number of modules; equal vertex orders require a common
letter order or an additional order label.

## 2. Exact finite rational carry relation

Set `g=gcd(h,Q), p=h/g, q=Q/g`; gcd(p,q)=1 and q>p. Compare equal-length
words L and R, priming the R coordinates and writing delta for R−L.

**Theorem1.** Equality of V is equivalent to integer carries

    x_0=x_k=0, q*x_i=p*x_(i-1)+delta c_i.              (3)

Equality of U is equivalent to reverse integer carries

    y_0=y_k=0, q*y_(i-1)=p*y_i+delta a_i.              (4)

All these carries satisfy `|x_i|,|y_i|<g`.

Proof. Divide the V difference in(2) by g^(k-1). Substitute
`delta c_i=q*x_i-p*x_(i-1)` in
`sum delta c_i q^(i-1)p^(k-i)`; it telescopes to q^k*x_k.
Conversely reduce this zero sum modulo q, use gcd(p,q)=1, divide, and repeat.
This supplies every integer carry and terminal zero. Reverse the word for U.
Since |delta c_i|<=g(q-p), induction gives
`|x_i|<=g(1-(p/q)^i)<g`; the reversed argument bounds y. ∎

Thus a two-sided carry automaton has at most(2g−1)^2 states. Its left-to-right
transition on the physical letter pair(i,i') is

    q*x'=p*x+delta c_i, p*y'=q*y-delta a_i.           (5)

Start/end at(0,0). Equal first e and equal last d are the remaining endpoint
conditions for full rooted-packet equality except H. An endpoint flag suffices.
When g=1 all carries vanish: U,V equality forces both digits equal letter by
letter. This is not letter isomorphism or H equality. When g>1 actual carries
can trade receiver digits rather than matching them separately.

## 3. Localizing the exact count change

For k>=2 the Hamilton count has the exact expansion

    H_word=h^k+h^(k-1)sum_i b_i
      +sum_(i<j)h^(i-1+k-j)Q^(j-i-1)a_i c_j
      +h^(k-2)sum_(i<k)s_i t_(i+1).                 (6)

A matrix path either jumps directly, visits the Q channel, or visits the
zero-diagonal signed channel at adjacent letters. There is no third strict
upper jump. At k1 the two pair sums are absent and H_word=H_1.

For closed ordinary carries define the integer local action

    J=h sum_i(H'_i-H_i)
       -p sum_(i<k)[a_i*x_i+c'_(i+1)*y_i]
       +sum_(i<k)[s'_i*t'_(i+1)-s_i*t_(i+1)].        (7)

**Theorem2.** `H_R-H_L=h^(k-2)J`. At k1 this reads J=h(H_R−H_L).

Proof. Normalize the cross term in(6) by h^(k-2) and put r=q/p. Expand
`a'_i c'_j-a_i c_j=a_i delta c_j+delta a_i c'_j`.
The closed carries give the exact telescoping identities

    sum_(j>i) r^(j-i-1)delta c_j=-p*x_i,
    sum_(i<j) r^(j-i-1)delta a_i=-p*y_(j-1).

The direct and signed channels supply the other terms of(7). ∎

For k>=2, every receiver-preserving gap is therefore divisible by h^(k-2).
Tournament parity strengthens any nonzero gap to a multiple of2h^(k-2).
For h>1 a fixed-alphabet full-central trade at length>=3 cannot give gap2.
This is an exact length barrier, not a failed search.

**Finite weighted compiler.** Retain the preceding physical letter pair in
the state(5). Each letter contributes site weight h*deltaH and each transition
the adjacent terms of(7). The series counting accepted labelled pairs by
length z and normalized gap J (Laurent variable t) is rational, of form
initial-vector `(I-zB(t))^-1` terminal-vector. This is a normalized-gap
series, not A000568 or an absolute-count selector; program collisions remain.
Finite-state compatibility does not supply all baselines or arbitrary targets.

## 4. Physical full-central +2 at deletion169

Let D have outgoing masks `(18,8,3,5,14)`. Add root0 with outgoing cut c
on D, giving T_c. For A=T_9, B=T_3, C=T_24 the packets are

    (H,h,Q,U,V,d,e):
    A=(37,13,104,46,71,8,5),
    B=(33,13,104,45,63,7,5),
    C=(43,13,104,54,72,8,6).

There is an exact positive physical identity

    M(B)M(C)=M(A)^2+2E_14.                            (8)

The eleven-vertex graphs A*A and B*C have counts2747,2749 and the SAME
remaining full packet

    (h,Q,U,V,d,e)=(169,10816,4030,6955,104,65).        (9)

The reduced base is8; L=AA,R=BC have carries x=(0,−1,0), y=(0,1,0).
The three contributions to(7) are26,−26,2. Discarding the signed path-length
channel incorrectly predicts zero. Every statement in(8)-(9) is independently
checked by literal eleven-vertex counting, not just matrix substitution.

For arbitrary actual common rooted contexts L,R, singleton permitted,

    deltaH=2h(L)h(R),                                (10)

and all other packet coordinates remain equal. Indeed
`M(L)E_14M(R)=h(L)h(R)E_14`. Unit-deletion contexts preserve gap2; nonunit
contexts pay the exact multiplication cost. This does not supply a +2 move
at an arbitrary source. Deletion169 is coprime to the older deletion15
central41/43 trade and therefore provides a different central content.

A common left triangle sends H toH+2V and V to4V−h, e toh. It yields
16657/16659 sharing(h,V,e)=(169,27651,169), a renewable two-bank.

## 5. Arithmetic cost of carry flexibility

For the common(h,Q) alphabet, (6) gives

    g^max(0,k-2) divides H_word.                      (11)

Thus g1 makes receiver digits rigid, while g>1 permits carries but forces
growing content in long words. More generally, for differing letter pairs
(h_i,Q_i), let g_i=gcd(h_i,Q_i). For every prime l and k>=2,

    v_l(H_word)>=sum_i v_l(g_i) - the two largest v_l(g_i).                 (12)

Each matrix path skips at most two diagonal factors. Every surviving term,
and hence their sum, has this divisibility. Zero diagonal terms do not create
an exception. This is a necessary budget, never a converse realization rule.

For repeated A in Section4 and k>=2, the explicit formula is

    H(A^k)/13^(k-2)
       =169+312k+1914*(8^k-7k-1)/49+40(k-1).          (13)

The quotient `(8^k-7k-1)/49=sum_(j=1)^(k-1)(k-j)8^(j-1)` is integral.
Modulo13 the right side is8k−8^k. Since8^4=1mod13, this residue
has period52; checking that finite period gives the following four classes.
Thus the exponentk−2 is attained whenever
k is not in{1,34,44,51} modulo52. At k34 it first gains an extra factor13
after the k>=2 boundary. No exact valuation beyond this first test is claimed.
This strengthens the original finite sharpness control without promoting it
to an invalid all-k equality.

## 6. Exact monomer-dimer fibers and hidden rerooting freedom

Within the THREE-letter alphabet A,B,C, compare a word W with A^k. Its
V-digit differences are0,−8,+1. The carry automaton from0 has exactly

    0 --A--> 0,   0 --B--> −1,   −1 --C--> 0.

All other exits from these reachable states have a nonintegral next carry.
Therefore, for every k>=1,

    V(W)=V(A^k) iff W is tiled by monomers A and dimers BC.                (14)

The identity(8) makes every such word agree with A^k on its complete packet
except H. A word with j dimers has

    H(W)=H(A^k)+2j*13^(k−2), 0<=j<=floor(k/2).         (15)

At k1 there are no dimers and the count is simply37. The empty word is the
singleton identity and is treated separately. At k>=2, expand each dimer
matrix as M(A)^2+2E_14. Two E_14 insertions give zero, and a single insertion
is multiplied by13^(k−2) by the remaining letters. This proves(15).

There are exactly binom(k−j,j) such words with j dimers, with bivariate
program generating function1/(1−z−t*z^2), including the empty word.
D is strong (verified directly). Deleting the distinguished root of any
word gives the ordered SCC chain D->...->D. Any rooted isomorphism must
preserve each SCC position; restricted to that block and the root, it
recovers the letter because A,B,C have distinct H37,33,43. Thus distinct
programs here give distinct ROOTED isomorphism classes. In particular
binom(k−j,j) rooted classes share the SAME FULL rooted packet. This is a
statement within the specified alphabet, not the full packet fiber among
all tournaments or an unrooted class count.

The lost microscopic information becomes usable when rerooting is permitted.
At k3, BCA and ABC both have packet

    (H,h,Q,U,V,d,e)=(264251,2197,1124864,409318,717743,1352,845).

But their vertex-deletion count multisets differ: BCA has minimum nonroot
deletion count30693, while ABC has29047. The separate verifier records and
independently recounts both complete menus. They are therefore nonisomorphic
even without distinguished roots. This is not just a different presentation.

For ABC, reroot at vertex1 in the literal fusion labeling. Its ending mass
is325181 and its deletion count29047. A triangle fused at that NEW root has
count264251+2*325181=914613, which is11mod13. Thus the fixed-language content
obstruction of Section5 is genuinely escapable by a physical root change;
it must not be stated as a theorem about arbitrary continuations. No claim
is made that all odd targets are accessible through this escape.

## 7. Reproduction and scope

`python3 wildcard/verify.py` supplies independent backward subset counts,
permutation/split controls for all32 core cuts, literal eleven-vertex packets,
all1,048,576 two-letter pair comparisons, local-action checks on every accepted
pair, and seeded longer controls in bases8,28/3,352/37. `fiber_extension.py` independently checks the monomer-dimer support through
k8, its exact full-packet formula, and both sixteen-vertex root menus.
Normal and optimized runs must match. The all-order conclusions follow from the displayed proofs.

A preliminary sampler of50 small labelled degree/H/Q signatures found theh13
unittrade and no other baseh; those signatures were only search heuristics,
not an isomorphism classification. Simple transitive blowups of D failed to
preserve(8). These searches are not necessary for the theorem and prove no
global exclusion. The unresolved task is to leave the fixed-base length barrier
while retaining actual carry realizations and controlling the absolute target.
