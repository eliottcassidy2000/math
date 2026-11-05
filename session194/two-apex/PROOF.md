# Unavoidable periodic counters and the short-block radix boundary

Session 194. **PROVED with finite exact invariant-space certificates.** The
Hamilton-path spectrum conjecture with exceptions `{7,21}` remains **OPEN**.
The bounds of this theorem concern the length of a repeated attachment block,
not the size of the initial tournament or the number of repetitions.

## 1. Physical definitions and inherited mechanisms

Fix named apexes a,b with a->b and an increasing chain x_1->...->x_m. The
letter c_i in {0,1,2,3} records a->x_i in bit 0 and b->x_i in bit 1. Let T(w)
be this tournament and H(w) its directed Hamilton-path count. A **physical
seed** is any finite word u, including the empty word. No strongness,
minimality, or restriction on the seed length is imposed.

A **periodic block** d is a nonempty word on this alphabet; its repetitions
mean literal suffix concatenation u,ud,udd,... . Old apex names, their mutual
arc, the old chain order, and all old incidences are unchanged. A block is
**mixed** when it uses at least two distinct letters. This is not a theorem
about arbitrary tournament extensions, internal edits, changing apex marks,
or state-dependent choices of blocks.

The closest proved tools are S191's two-apex transfer (36 nonnegative
transfer states, reachable rational dimension 28) and S192's common tensor
flag for all attachment matrices, routed through the published
[Session 192 guide](https://github.com/eliottcassidy2000/math/blob/8d04f8557ea50c985521f039d5f58aa335617c9a/00-navigation/S192.md). The
[Session 193 guide](https://github.com/eliottcassidy2000/math/blob/169792fb5eeb483e351112b5a09706e6d7bd89f0/00-navigation/S193.md)
routes S193's full-response decoder and its fixed binary/ternary suffix
grammar obstruction, which explicitly leaves larger radix blocks
open. Here the retained sidecar is the actual polynomial part of the mixed
period response, with an exact modular certificate that it cannot vanish.

The hostile is still H(120)=13 and the all-order failure of an end-growth
ancestry 13 -> 41 -> {125,127}; this result neither repairs that ancestry nor
claims that the target values are absent. The corrected near miss is fitting
an exponential scalar recurrence on a few seeds while discarding a small
polynomial counter. The concept board is periodic words, invariant reachable
space, nilpotent counters, modular normalization, and exact-target branching.

## 2. The all-prefix periodic-counter theorem

For a block d of length L, write B_d for its 36-state transfer matrix. The
common tensor flag, proved below in the inherited S192 coordinates, gives
for every prefix u and every k>=3,

    H(u d^k) = A(u,d) + J(u,d) k + K(d) k^2
             + (D(u,d) + E(u,d) k) (2^L)^k
             + G(u,d) (3^L)^k.                     (1)

Initially the general theorem allows K to depend on u. The certificate below
proves its independence for all blocks used here. The six functions in (1)
are linearly independent over Q: their annihilating roots 1,2^L,3^L are
pairwise distinct, with respective multiplicities 3,2,1. Thus the six
coefficients are uniquely determined by k=3,...,8. The all-k formula comes
from the common tensor flag, not from the six samples.

**Theorem 1 (all physical prefixes, bounded mixed period).** If d is mixed
and 2<=|d|<=6, then for every finite physical seed u the polynomial
A+Jk+Kk^2 in (1) is nonconstant. The same conclusion holds for each block

    d = c^L e^L,  c!=e in {0,1,2,3}, 1<=L<=6.       (2)

In particular the periods in (2) may have length as large as twelve.

### Self-contained common-flag and Jordan justification

For one apex order the transfer space is Q^3 tensor Q^3 tensor Q^2.
In a bounded-run factor choose E,B,D=G-B (empty, bad-last, and the
signed difference of good-last and bad-last). Give these coordinates levels
0,1,2; give the final E,O factor levels 0,1. A single-run reset is

    E -> alpha (B + delta D),
    B -> B + delta D,
    D -> 0,

where alpha is the permission to start that run and delta indicates a good
last edge. In the final factor E->alpha O and O->O. The full letter matrix
is the sum of the three local resets tensored with identities on the other
factors. Thus it only increases coordinates, and its diagonal at a tensor
state t is j(t), the number of level-one factors, independent of the letter.
The verifier checks this change-of-basis identity on all 144 letter/basis
column combinations across the two apex orders.

A length-L product therefore has diagonal j(t)^L. Among states with j>=1,
the total coordinate level is j+2a, where a is the number of D factors.
It lies between j and 6-j. Hence a strictly increasing coordinate path visits
at most 4-j states with diagonal j^L. For diagonal zero it visits at most
three states, since the final factor must be empty and each bounded factor
can only move from E to D.

To turn these counts into a Jordan bound, write the block matrix as diagonal
D_0 plus strictly triangular N and expand its resolvent as the finite series

    (zI-B_d)^(-1) = sum_a ((zI-D_0)^(-1) N)^a (zI-D_0)^(-1).

Each term is a sum over increasing coordinate paths. Its pole at j^L has
order at most the number of visited states with that diagonal. Consequently

    B_d^3 (B_d-I)^3 (B_d-2^L I)^2 (B_d-3^L I) = 0.

The zero component vanishes by k=3; the three positive components give
exactly (1), with degree bounds 2,1,0. Distinct rational eigenvalues give
rational coefficients. This proves the all-period representation within
this packet, without extrapolating any scalar samples or needing an external
linear-algebra certificate.

The complete certificate counts are:

| Period length | All mixed blocks | Nonzero quadratic | Forced nonzero linear |
|---:|---:|---:|---:|
| 2 | 12 | 0 | 12 |
| 3 | 60 | 36 | 24 |
| 4 | 252 | 192 | 60 |
| 5 | 1020 | 900 | 120 |
| 6 | 4092 | 3816 | 276 |

There are 72 separately verified composition presentations in (2), including
some duplicates of the first census. No isomorphism or necklace filter is
used. Every physical word in the stated universes is checked.

### Exact proof certificate, including the unbounded-prefix quantifier

`verify.py` constructs the actual run-assignment transfer, not fitted
matrices. In each apex order a path is L_0,p,L_1,q,L_2. A bounded run retains
empty/good-last/bad-last status; the final run retains empty/nonempty. First
boundary edges are checked on birth; last boundaries are checked at readout.
The two orders give 36 states. This is S191's path bijection.

Starting from the actual empty vector v_empty, the verifier builds an
independent rational spanning family of reachable physical vectors. Each
letter image of each added vector is tested. All 28 accepted vectors are
therefore an invariant reachable basis, so identities checked on them hold
for every finite prefix, not merely for prefixes up to a sampled length.

For each d, the six exact rows f B_d^k, k=3,...,8, are multiplied by the
inverse of the six-by-six interpolation matrix from (1). Integer numerator
rows and a common denominator are retained; all arithmetic is exact. The
quadratic row is checked to have one constant value on all 28 physical basis
vectors. Every physical vector has normalization

    v_0 + v_6 + v_12 = 1,

so this proves prefix independence on the entire reachable space.

If that value is nonzero, the block is certified. If it is zero, let r/d be
the linear coefficient row, reduced so that r is a primitive integer row
and the nonzero rational scale is retained. The verifier finds a prime p
and a residue t!=0 such that

    r v(u) = t (mod p) for every physical prefix u.  (3)

This is not inferred from testing many prefixes. For the chosen p, it builds
an invariant reachable spanning family over F_p, again testing all four
letter images of every added vector until closure. It checks (3) on that
family, using the same normalization coordinate. Linear algebra then proves
(3) on every reachable vector. Thus r v(u) cannot be zero as an integer, and
J(u,d) cannot vanish. This completes the theorem for all prefixes.

The deterministic receipt lists every block, whether its obstruction is
quadratic or linear, and respectively its exact quadratic coefficient or
linear scale/prime/residue. A SHA-256 digest pins all primitive linear rows.
The script regenerates the rows rather than loading them from the receipt.

### Small transparent examples

For d=01,

    J(u,01) = (2/3) [v_0+v_6+v_12-v_32-v_34].        (4)

The bracket is always odd; the modular reachable-space certificate uses
p=2. Its value can have either sign, so positivity would be the wrong proof.
For example prefixes empty,120,2330 give J=2/3,2/3,-10/3.

For d=03 and d=12, the primitive numerator row has value 1 modulo 3; its
scale is 1/9. At the empty seed J=4/9. These blocks have K=0, but their
linear counters cannot be removed by a clever prefix.

For d=003, K=-8/49 for every prefix. For d=012, K=13/49 for every prefix.
The counter itself need not be positive; the full Hamilton count remains
positive because the exponential terms dominate.

The exceptional zero-quadratic cases at length six include genuinely
nonconstant patterns such as 010323. Its linear coefficient is a primitive
integer row times 1/351; that row is 6 modulo 13 on every physical prefix.
Thus checking only the quadratic coefficient would miss a necessary case.

## 3. Consequence for exact affine iteration

**Corollary 2.** For any block covered by Theorem 1, any physical seed u,
any real q>1 and any real b, the sequence h_k=H(u d^k) cannot satisfy

    h_(k+1) = q h_k + b  for every k>=0.             (5)

It cannot satisfy this recurrence eventually either.

Indeed (5) gives h_k=C q^k + constant. If C=0 the sequence is constant. If
C!=0, the matrix recurrence forces q to be 2^L or 3^L: the other possible
eigenvalues are 0 and 1. In either case uniqueness of exponential-polynomial
expansions forces J=K=0, contradicting Theorem 1. For an eventual recurrence,
discard finitely many terms; a nonzero polynomial counter remains nonzero.

Constant-letter blocks are intentionally excluded. They really can be
exactly affine. For all k>=0,

    H(2 0^k) = 2^(k+1)+1,
    H(2330 0^k) = 42*3^k-17.

Thus the theorem does not confuse periodic growth with a universal no-go.

## 4. Consequence for short fixed-block digit grammars

**Theorem 3.** Let F be any nonempty family of physical words. There do not
exist two fixed nonempty blocks d_0,d_1, each of length at most six, a common
real q>1, and distinct offsets b_0,b_1, such that F is closed under suffix
concatenation by either block and, for every u in F,

    H(u d_i) = q H(u) + b_i,  i=0,1.                (6)

There is no bound on the seeds, family size, word lengths, or hidden state.
The two macro lengths may initially be different.

**Proof.** Each branch has some seed whose scalar value is not its affine
fixed point: otherwise all values in F would equal that fixed point, and
the other branch would force equal offsets. Repeating the branch and using
the common tensor flag forces q in {2^|d_i|,3^|d_i|}. Equality of these
powers forces both lengths equal, say L, because no positive power of two
is a positive power of three.

If either block is mixed, its repeated sequence contradicts Corollary 2.
Thus d_i=c_i^L. Distinct offsets force c_0!=c_1. Repeating the composition
d_0d_1 gives an affine scalar recurrence with multiplier q^2>1 and offset
q*b_0+b_1. Its period is c_0^L c_1^L, covered by (2), another contradiction.
QED.

This closes the larger-radix fixed-macro loophole through six attachment
letters per branch. It does not eliminate arbitrary longer macros or
finite-state rules that choose different macros in different states.
State-dependent rules, internal surgery, re-marking apexes, and rebuilding
a target tournament remain outside its scope.

## 5. What this does and does not advance

The full-response decoder in S193 is too rigid to permit a scalar correction
that survives every suffix. This theorem uses much less future data, just
one repeated block, yet shows that short mixed periodic blocks still carry
an unavoidable lower-order counter. The modular normalization is what makes
the statement uniform in the physical seed.

It does not imply a missing-density theorem for unrestricted two-apex words.
It does not imply any target is absent, any fixed-order interval is full, or
that a larger finite prefix proves a spectrum theorem. The inherited positive
missing density for one-apex constructions still blocks replacing exact
integer landing by modular coverage. The H13 ancestry obstruction remains
applicable to pure end-growth; the present theorem is complementary and
requires neither that target nor a strong-order bound.

A constructive route must control, rather than discard, these polynomial
coordinates, or change the physical operation. The certificate is a fast
hostile test for future proposed fixed-block digit constructions.

## 6. Reproduction and independent controls

Run from the repository-independent packet directory:

    python3 verify.py > normal.stdout
    cp receipt.json receipt.normal.json
    python3 -O verify.py > optimized.stdout
    cmp normal.stdout optimized.stdout
    cmp receipt.normal.json receipt.json

The verifier uses only the Python standard library. It checks every mixed
word of length at most six, all 72 constant-branch compositions, invariant
reachable closure over Q and the stated finite fields, 341 complete words
of chain length at most four by independent vertex-subset Hamilton-path DP,
and five selected larger physical words. It also verifies the genuine
constant-letter affine controls and the full-vector constant-letter
annihilator. All checks use explicit exceptions, so optimization cannot
remove them.

The all-period exponential-polynomial representation is proved directly in
section 2 using the inherited S192 common tensor coordinates; the verifier
checks the local change of basis. No external literature bound or unsupported
interpolation assumption is needed. This packet makes no literature-priority claim.
