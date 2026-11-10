# Round 2: exact synchronized-pair continuation and finite-prefix immunity

Status: **PROVED scoped algebraic rigidity and scalar countermodel**.
No claim of all-odd coverage. Physical construction of a complete inverse
invariant remains **OPEN**.

## 1. A chosen missing prime after any finite exact prefix

Let F be any finite set of attained odd values, with 7,21 absent. Choose any
prime p>max(151,max F), and define

    S_p = {positive odd integers} \ {7,21,p}.

Then S_p contains F, is multiplicatively closed, and has every compact gap
pair, long odd interval, infinite gap fibre, and uniform additive-28 anchored
gap property proved in PROOF.md. Indeed S_p contains the previous model;
the same explicit witnesses still work. Multiplicativity is also direct:
a product cannot equal the excluded prime p unless one factor is p, and
the only nontrivial factorization of 21 uses the excluded factor 7.

Thus every finite exact prefix, however long, can be respected while omitting
one specified prime beyond it. There is no inference from such finite data
plus those scalar properties to full coverage. This does not model the
entire physical tournament category or refute a proof using further physical
identities. One may replace max F by any given finite prefix endpoint.

## 2. All common synchronized-face continuations are Z strings

Use the inherited O/Z/J four-coordinate linear maps, and the paired-word
outputs proved in this session's `../composition/PROOF.md`. They have the
same E=d and U=D, with positive difference

    Delta = q(1,0,1,0), q>0.

Apply the **same arbitrary finite word** W in O,Z,J to both rooted states.
If the outputs again satisfy

    Delta E=Delta d=0, Delta U=Delta D,                      (1)

then W is Z^k for some k>=0. Conversely every Z^k satisfies (1), and the
resulting full-count gap is 2^k q.

**Proof, excluding J.** Work with the normalized difference coordinates
(x,y,z,t), initially (1,0,1,0). The cone

    x>y>=0, z>t>=0

is preserved by all three maps. For O this is immediate. For Z the first
strict difference becomes (x-y)+(z-t)>0 and the second becomes 2(z-t)>0.
For J the two strict differences become z>0 and z>0; the new y=x-y>0
and new t=x>0. Once any J has occurred, t is strictly positive forever:
O and Z retain t, and each later J replaces it by a positive x. Hence a word
satisfying Delta d=0 contains no J.

**Proof, excluding O.** For an O/Z word of length k, y=t=0 throughout,
while z=2^k. If positions are numbered i=0,...,k-1, then

    x=1+sum_{i: W_i=Z} 2^i.

The equality z=x holds iff every position is Z. This proves the theorem.
No finite search, scalar fit, or assumption about graph order enters it.

**Consequence.** The new synchronized compiler is stable under exact gap
dilation, but cannot provide a nontrivial same-word synchronized additive
repair, gap halving, or movement at constant gap. In particular, preserving
synchronization on both coordinates while exploiting additional O/J moves
is impossible under a common continuation. This is a restriction of this
alphabet and this synchronized pair, not of asymmetric constructions or
arbitrary tournament operations.

## 3. Asymmetric universal face continuations still have an anchor term

A *universal face-return word* is a word whose linear map carries every
formal state (h,e,h,e) into that same two-dimensional face. Suppose two such
words A,B have equal final endpoint coordinates when applied respectively
to (h,e) and (h+g,e), as an identity in the independent variables h,e,g.
Then both words are Z strings.

Indeed, let the endpoint row for B be (c,d). Equality forces c=0, and also
forces A's endpoint h-coefficient to be zero. For any universal face-return
word, zero endpoint h-coefficient implies zero deletion h-coefficient.
Applying the word to the derivative vector (1,0,1,0), its output therefore
satisfies (1). The preceding theorem forces that word to be Z^k.

Write A=Z^i, B=Z^j. Since

    Z^k(h,e)=(2^k h-(2^k-1)e,e),

the resulting gap is exactly

    2^j g+(2^j-2^i)(h-e).                                  (2)

It is independent of the source anchor h-e iff i=j. Therefore even allowing
different universal face-return words does not provide an anchor-independent
repair while synchronizing endpoints on an arbitrary independent pair.
This universality hypothesis matters: target-dependent words on a restricted
physical subfamily can evade the conclusion and remain a valid research
route. No impossibility is claimed for that broader setting.

## 4. The new face maps impose exact inverse gap congruences

For the actual face maps in `../composition/PROOF.md`,

    F-(h,e)=(20h-3e,6h-e),
    F+(h,e)=(20h-e,6h-e),

suppose two physical face inputs have synchronized output endpoints and
output full-count gap 2N. Input h and e are both odd, so their differences
are even. Solving the exact two-by-two systems gives

    F-: (Delta h,Delta e)=(N,6N),
    F+: (Delta h,Delta e)=(N/7,6N/7).

Consequently same-F- synchronized gaps are divisible by 4, and same-F+
synchronized gaps are divisible by 28. These are necessary congruence
conditions, not sufficient physical realizability criteria. They clarify
why the cross-branch compiler gap 2e with odd e is not reproduced by applying
one common branch to another synchronized pair. Its utility comes from the
branch change, not from a universal integer inverse of determinant 2 or 14.

## 5. A finite algebraic induction certificate to seek

A usable positive replacement can be specified without assuming arbitrary
target realizers. Choose a finite collection of *physical state classes*
C_i, each defined with its actual graph-sidecar interpretation; finite seed
states; and a finite directed list of literal graph operations. An edge
(i,j,a,b) certifies, on **every state in C_i**, both:

    H(output)=a H(input)+b, and output belongs to C_j.

The graph operations may be word maps or proved grafts. For an O/Z/J word,
the scalar identity is a finite matrix equality modulo the defining linear
relations of C_i; closure additionally requires actual physical operations
and proof of the sidecar relations, not just integer inequalities. If each
class is specified by invariant linear equalities, their preservation is
also a finite matrix computation. Any required inequalities must be checked
on the whole class, including strongness or root reachability guards.

To conclude all-odd coverage, separately prove a finite residue-based
inverse selector: each desired target/class pair outside the seed set has
an incoming edge for which (m-b)/a is a positive smaller allowed integer,
and lies in its stated residue class. This last arithmetic induction is
finite-checkable when the selector is a finite congruence partition with
explicit threshold, together with exhaustive base cases.

This certificate format is substantially more structured than asserting a
target graph exists: finite matrices, class equations, physical closure,
residue guards, and a decreasing rank suffice. The five-seed binary tree in
PROOF.md is its simplest arithmetic instance. The new rigidity theorem says
that classes using only a synchronized face pair and common O/Z/J words
cannot supply the needed additive edges; asymmetric restricted fibres or
extra graft coordinates are necessary. **No completed certificate of this
kind is supplied here.** Proving one remains a substantive construction
problem, not a consequence of the existing every-gap theorem.

## 6. Verification and a caught sign correction

Run `python3 -I -B verify_round2.py` and compare with `ROUND2_CHECK.json`;
repeat with `-O`. The checker exhausts all words through length 10, verifies
the invariant difference cone and synchronized-word classification, checks
all resulting universal face maps for the zero endpoint coefficient rule,
and checks the inverse congruences for N=1,...,10000.

The direct inverse check caught a sign reversal in the initially supplied
composition inverse formula. The correct inverse of F- is

    h=(h'-3e')/2, e=3h'-10e'.

For example (9,3) maps to (171,51), which this inverse returns to (9,3).
The correction was sent to the parent for the composition packet. It does
not change the forward maps or the new pair compiler. This packet's final
inverse differences have the positive signs shown in §4.

## 7. The tempting branch-cancellation repair acquires content 6

For the newly compiled pair, let its lower face state be (h,e), its upper
state (h+g,e), and write the original source coordinates as (D,E,U,d). Then

    g=2E, e=3D-E+3U.

Try applying F+ to the lower member and F- to the upper member, the branch
choice whose subtraction can reduce the gap. The full difference becomes

    (20g-2e,6g,20g-2e,6g)
      = (6(7E-D-U),12E,6(7E-D-U),12E).

Every coordinate is divisible by 6. Because every O/Z/J transfer matrix is
integral, *every subsequent common continuation* retains that divisibility.
Its total-count gap can therefore never be 2 or -2. This is an exact content
obstruction to that specified branch-cancellation repair, not to other
asymmetric operations or source constructions. The forward-branch choice
and more general noncommon continuations are not covered by this claim.
