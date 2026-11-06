# Exact binary-child classification after arbitrary physical source rebasing

Status: **PROVED, independently audited in the stated scope.** The theorem classifies every
physical source parent that admits either binary child by an append-only word,
even after replacing the parent source by any other actual source with the same
path readout. Both children are possible precisely at parent counts 5, 9 and 17.
This obstructs that particular induction strategy. It is not a nonexistence
theorem for child tournaments and does not prove or refute eventual odd coverage.

## 1. Physical domain and the proposed induction

Use actual nontrivial rooted start-one sources, with profile
s=(D,E,U,d), the complete chronological O/Z/J grammar from
`dependencies/PHYSICAL_GRAMMAR.md`, and seed t=(1,0,2,1):

    O(s)=(D,0,2U,d),
    Z(s)=(D-E+U,d,2U-d,d),
    J(s)=(D-E+U,D-E,D+U,D).

A final J closes a source into a strong C1 tournament whose Hamilton-path
count is F(s)=D-E+U. The nonnegative physical coordinates are

    y=(a,f,e,p)=(D-E-d,U-D+E-d,E,d-E),
    D=a+2e+p, E=e, U=a+f+2e+2p, d=e+p.

In particular d>=1 and D is odd. The inherited balanced inverse in
`dependencies/BALANCED_SOURCE.md` says that every nontrivial source with
D=U is canonically OJ of an actual source W, allowing W singleton, with

    (D,E,U,d)=(D_W+2U_W,D_W,D_W+2U_W,D_W).

Only this physical inverse, not inversion of a free integer vector, is used.
The source-domain completeness is supplied by the physical generator in
`dependencies/COMPLETE_SOURCES.md`, Section 1, and its equivalent one-vertex
grammar in `dependencies/BALANCED_SOURCE.md`, Section 1. Alternatively, every
statement below can be read just on the explicitly defined physical orbit;
no all-tournament classification is needed for the orbit theorem.

The proposed induction is deliberately permissive: given a parent count n,
choose ANY actual source s with F(s)=n, independently for the two signs, then
append ANY nonempty physical word w and require F(s.w)=2n+1 or 2n-1.
The chosen s need not be related to an originally supplied source, so this
includes arbitrary target-preserving source rebasing before the extension.
It does not include replacing a terminal segment while changing the parent
readout, more general tournament surgeries, or a change to a C2/C3 interface.

Closest mechanism: the complete physical grammar and balanced inverse.
Hostile parent: 25, which is represented but has neither binary child by this
procedure. Corrected near miss: a signed offset identity does not supply the
source sidecar needed at a prescribed parent. Least-used sidecar: the exact
integer h=f-d, which Z preserves. Live concepts are physical inverse, binary
scalar induction, source rebasing, endpoint saturation and append length.

## 2. Two or more letters are already too large

Direct substitution gives

    F=(2,1,3,3)y,
    F(O s)=F(J s)=(3,2,6,5)y,
    F(Z s)=(4,3,5,5)y.

Every output is at least (3/2)F(s), coefficient by coefficient on the
nonnegative physical cone. Consequently a word of length at least two has

    F(s.w)>=(9/4)n>2n+1, when n>=5.

Thus for every odd parent n>=5, an append-only binary child exists if and
only if a ONE-letter child exists after a permissible rebase. No bounded
word-length assumption is being imposed; the inequality excludes every
longer suffix.

Put h=U-D+E-2d=f-d. The exact one-letter readouts are

    F(Z s)=2F(s)+h,
    F(O s)=F(J s)=2F(s)-(D-2E).

Since D-2E=a+p is a positive odd integer, O and J can never give the plus
child. Z gives the plus/minus child exactly at h=+1/-1; O and J give the
minus child exactly at a+p=1.

The coordinate moves needed to classify those slices are

    O(a,f,e,p)=(a+e,a+2f+e+2p,0,e+p),
    Z(a,f,e,p)=(2a+f+e+p,f,e+p,0),
    J(a,f,e,p)=(f+p,0,a+e+p,e).

## 3. Complete physical classification of h=-1 and h=+1

Write s_j=t.Z^j for j>=1. Its actual profile is

    s_j=(2^j+1,1,2^j+1,1).

Sources with D=1 are precisely the transitive sources t.O^r, r>=0,
whose cover count is 2^(r+1), plus the separately allowed singleton with
cover count 1. Indeed O preserves D, whereas Z/J on a nontrivial source
create D=F>=3. This also follows from the source graph interpretation.

### Negative slice

**Lemma 1.** Every physical source with h=-1 has the numerical profile of
exactly one of t.Z^m, m>=0. Its readouts are exactly

    A_minus={2^r+1:r>=1}.

Proof. Z preserves h, so remove all final Z letters from a physical recipe.
An O output has h=a+2f+p>=0. A J output has h=-D_old; hence it can have h=-1
only from a physical source with D_old=1. Such a source is t.O^r, and
J(t.O^r) has the same profile as s_(r+1). The alternative terminal object
is the seed t itself. Restoring the final Z letters proves exhaustion.
Every displayed profile is physical, and direct iteration gives the stated
readouts, including F(t)=3. QED.

### Positive slice

**Lemma 2.** Every physical source with h=1 has the numerical profile of
one of the following actual words, and all listed words have h=1:

    t.O.Z^k,                  k>=0;
    t.Z.O.Z^k,                k>=0;
    t.Z^j.J.O.Z^k,            j>=1, k>=0.

Their respective readouts are

    6*2^k-1,  10*2^k-1,  6(2^j+1)*2^k-1.                 (1)

Proof. Strip final Z letters. Neither t nor a J output can have h=1, so
the last remaining letter is O. Its predecessor obeys

    a+2f+p=1,

hence f=0 and (a,p) is either (1,0) or (0,1).

In the first case the predecessor is (D,E,U,d)=(2e+1,e,2e+1,e), e>=1.
It is balanced. The canonical physical inverse has D_W=e and
U_W=(e+1)/2. Since U_W>=D_W, e<=1, so e=1 and the predecessor is s_1.

In the second case the predecessor is (2e+1,e,2e+2,e+1). If e=0 it is t.
If e>0, a last O is impossible because E>0, and a last Z is impossible
because p=1. It therefore ends in J. Its actual predecessor has
D=U=e+1 and E=1. The canonical balanced inverse has D_W=1, U_W=e/2.
The D=1 classification forces e=2^j for some j>=1. Thus the profile before
O is exactly J(s_j). Conversely each of these profiles is supplied by the
displayed physical words.

Their post-O readouts are 5, 9 and 6*2^j+5. Since Z preserves h=1 and sends
F to 2F+1, k more Z letters give (1). QED.

These are classifications of numerical profiles; no uniqueness of source
graphs or source recipes is asserted.

## 4. Complete physical classification of O/J minus parents

**Lemma 3.** The readouts of actual sources with D-2E=a+p=1 are exactly

    B_minus={2^r+1:r>=1} union {3(2^j+1):j>=1}.           (2)

Proof. If E=0, then D=1 and the sources are t.O^r, giving the first set.
Assume E>0. A last O is impossible. For a last Z, its output has
p_new=0 and a_new=2a+f+d. Equality a_new=1 forces a=f=0,d=1. The predecessor
is either t (e=0,p=1) or (D,E,U,d)=(2,1,2,1) (e=1,p=0); the latter has even
D and is impossible. The output is therefore the triangle profile s_1.

For a last J, a_new+p_new=f+p+e=f+d=1. Thus f=0,d=1. If e=0 its predecessor
is t, and the output is again s_1. If e=1,p=0 its predecessor is balanced
with E=1. As in Lemma 2 the balanced inverse has D_W=1, forcing predecessor
s_j for some j>=1. The resulting J(s_j) has F=3(2^j+1).
All listed values are supplied by these literal words. QED.

In the J argument with e=0, the conditions f=0,d=1 imply U=D+1, E=0.
An E=0 source cannot end in Z or J, since either gives positive E.
A last O would give f_new=a+2f+e+2p>=e+p>=1, a contradiction.
Hence the only such source is the seed t. This explicitly handles the weak boundary rather
than assuming E=0 implies D=1 without the extra condition f=0.

## 5. Exact theorem, intersection and hostile 25

**Theorem.** For every odd n>=5, after arbitrary target-preserving rebasing
to an actual nontrivial source and any nonempty O/Z/J extension:

- The plus child 2n+1 is available iff n belongs to the union (1).
- The minus child 2n-1 is available iff n belongs to (2).
- Both children are available, allowing different rebases, iff
  n is one of 5, 9, 17.

The first two assertions combine Sections 2--4; the negative Z slice is
already contained in (2). For the intersection, first take n=2^r+1.
For r=1, n=3 is not a plus parent. For r>=2, n+1=2(2^(r-1)+1) has exactly
one factor of two, forcing k=0 in (1). The equations give n=5 or 9 in the
first two cases. The third requires

    2^(r-1)=3*2^j+2.

If j=1 this gives r=4 and n=17. If j>=2 the right side has 2-adic valuation
one, forcing r=2, incompatible with its size. Thus only 5,9,17 occur.

Next take n=3(2^j+1), j>=1. The two families in (1) with coefficient 6
are 2 modulo 3 and cannot occur. Equality with 10*2^k-1 gives

    3*2^j=10*2^k-4.

At k=0 this forces j=1 and n=9; k=1 gives no power-of-two solution. For
k>=2 the right side has valuation two, forcing j=2, but then its size is
already too large. This completes the intersection proof.

For positive controls, the minus/plus source words at parent counts 5,9,17
may respectively be chosen as

    5: Z / O;
    9: ZZ / ZO;
    17: ZZZ / ZJO.

Appending Z gives both desired children. Source words are from t, and the
final tournament compiler appends one further J.

The source word OJOO gives profile (5,0,20,1) and F=25, so parent 25 is an
actual C1 count. It belongs to neither set (1) nor (2). Thus NO physical
source representation of 25, and NO nonempty append word from it, can give
49 or 51. Those child counts themselves are realized: words JZJJZ and OOZZ
have F=49 and F=51 respectively. This distinction between a forbidden move
and an absent target is essential.

The two exact parent sets have O((log N)^2) values through N, and hence
density zero. This follows directly by bounding the nonnegative exponents
in (1) and (2). The theorem rules out an entire prescribed-parent binary
induction on this interface, not just a fixed finite menu or a synchronous
selector. Other target-dependent surgeries, different interfaces and
unrestricted tournament constructions remain outside scope. No threshold
for eventual odd coverage follows.

## 6. Verification

`python3 verify.py` and `python3 -O verify.py` use explicit exception gates,
independent positive-coordinate and raw-coordinate transfer calculations,
complete physical-source enumeration through U<=5000, exact generated slice
profiles, literal subset Hamilton-path/cycle counts for small positive and
hostile controls, and the arithmetic intersection through a declared range.
The all-order statements follow from the proofs, not finite stabilization.
The independent audit is supplied separately.
