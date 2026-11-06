# An intrinsic physical section of the 6/7 signed-response fiber

**PROVED, independently audited in the stated scope.** A two-equation slice of the
actual source orbit has an exact physical inverse. On this slice the signed
response `6U-7D` is automatically positive and is an actual unsigned source
path count. Every slice point either already has a canonical five-letter
physical presentation or admits a two-letter, target-preserving physical
rebase to one. The inverse contracts cover height by at least a factor six.

**OPEN:** proving that a prescribed signed fiber intersects this physical
slice. Signed-response surjectivity alone does not provide that intersection.
No universal C1 coverage or full `{7,21}` theorem is claimed.

## 1. Physical domain, inheritance, and the formal boundary

Let S be the complete nontrivial physical rooted start-one source orbit of
`t=(1,0,2,1)` under chronological operations

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

Also allow the singleton `e=(1,1,1,1)` and write S*={e} union S. Its first
O gives t. The response coordinates and graph constructions are those of
the frozen physical-adjoint and balanced-normal-form dependencies. Write

    F(s)=D-E+U, G(s)=D+2U,
    b=(3,1,3,1), z=(-1,-7,-1,-7),
    Q=[[0,0,0,1],[0,1,0,-1],[0,0,1,-1],[1,-1,-1,1]],
    K(v,w)=v^T Qw.

The adjoint reverses words and swaps O with Z while fixing J. It gives
`K(z,s)=6U_s-7D_s`, `K(e,s)=D_s`, `K(t,s)=F(s)`, and `K(b,s)=G(s)`.
The vector z is formal and nonphysical. No graph is assigned to it.

The distinction between a signed target and a positive target is visible
in a fixed rank-two physical transfer.
Direct exact iteration gives

    z.JJJOJ=e, z.JJJOJO=t, z.JJJOJOZ=b.

Consequently, on every formal input s,

    6U(s.JZJJJ)-7D(s.JZJJJ)=D_s,
    6U(s.ZJZJJJ)-7D(s.ZJZJJJ)=F(s),
    6U(s.OZJZJJJ)-7D(s.OZJZJJJ)=G(s).            (1)

All words applied to an actual source in (1) are actual source extensions.
Only the boundary calculation uses a formal vector. The global chronological
identity `JZ=OJ` rewrites the first suffix as

    W=OJ^4.                                     (2)

The new issue is the converse: can the image of this fixed physical suffix
be recognized inside an actual signed-response fiber without assuming that
its source word already ends in W?

## 2. The exact rank-two transfer

Put `n=D_s` and `m=G(s)=D_s+2U_s`. Then

    W(s)=(6m-n, 2m, 7m-n, 3m).                   (3)

Equivalently its matrix, in raw order D,E,U,d, is

    [[5,0,12,0],
     [2,0,4,0],
     [6,0,14,0],
     [3,0,6,0]].

Thus the old E,d are forgotten, but the old D,U are retained injectively.
The signed target of the output is exactly n. Since `U_s>=D_s` on S*,
`m>=3n>0`. Equality means the precursor is balanced, including the singleton.

Define the **intrinsic slice** T to consist of actual q in S satisfying

    E_q=2(U_q-D_q), d_q=3(U_q-D_q).               (4)

Membership in T requires actual source membership; equations (4) alone are
not a physicality criterion. For q in T set

    m=U_q-D_q, n=6U_q-7D_q.

Then q has the form on the right of (3). The physical parity invariants give
m odd and n odd. Positive coordinates give m>=1 and n<=m, but by themselves
allow negative n and do not give a physical inverse.

## 3. Intrinsic inverse and the exceptional branch

**Theorem 1 (physical section theorem).** Every q in T satisfies exactly one
of these alternatives:

1. `m>=3n>0`. There exists an actual precursor s in S* with

       D_s=n, U_s=(m-n)/2, q=W(s).

   The two displayed precursor coordinates are uniquely determined by q.
   Its forgotten E,d need not be unique at the numerical-profile level.

2. `m=n>0`. There exists a balanced physical source v in S*, with
   `D_v=U_v=n`, and

       q=v.OJJ=(5n,2n,6n,3n).

In particular, the positive-cone chamber `n<m<3n` and every nonpositive n
are empty on the actual slice.

**Proof.** Use the physical positive coordinates

    (a,f,e,p)=(D-E-d,U-D+E-d,E,d-E)>=0.

For q in (4) these are `(m-n,0,2m,m)`. Thus m>0 and n<=m. A source with
positive E but E!=d cannot end in O or Z, so every presentation of q ends
in J. Its physical predecessor q_1 has

    D_1=3m, E_1=m, U_1=4m-n.                    (5)

Its d_1 has not yet been inferred; J forgets the old deletion coordinate.

If n=m, (5) is balanced. Balance and nonnegativity force d_1=E_1=m.
The canonical balanced inverse supplies an actual balanced v in S* with
D_v=U_v=m and `q_1=v.OJ`. This proves alternative2.

Suppose n<m. The predecessor q_1 has E_1>0 and so ends in Z or J. If it
ended in Z, its old deletion coordinate would be d_1=E_1=m. For every
physical Z output one has

    2D-U=2(D_old-E_old)+d_old>=3d_old=3E.

Applied to (5), this would give `2m+n>=3m`, contrary to n<m. Hence q_1
ends in J. A J output has f=0, so its previously unknown deletion count is
now forced to be

    d_1=U_1-D_1+E_1=2m-n.

Its physical predecessor q_2 therefore has

    D_2=2m-n, E_2=m-n, U_2=2m.                  (6)

Again E_2>0. If q_2 ended in Z, the same Z inequality would require
`2(m-n)>=3(m-n)`, impossible. Thus q_2 also ends in J, and f=0 forces
`d_2=m`. Its physical predecessor q_3 has

    D_3=m, E_3=n, U_3=m.

Nonnegative E_3 and odd n give n>0. Balance forces d_3=n. The canonical
balanced inverse supplies s in S* with D_s=n, U_s=(m-n)/2 and
`q_3=s.OJ`. Since q=q_3.J^3, it is W(s). The inequality U_s>=D_s gives
m>=3n. This proves alternative1 and all exclusions. QED.

The proof never arbitrarily fills in a deletion coordinate forgotten by J.
It forces the predecessor's terminal letter before using that letter's
zero coordinate. This is the physical step missing from a mere matrix inverse.

## 4. Target-preserving rebase and quantitative height reduction

**Corollary 2.** Every q in T can be rebased to a canonical W output without
changing the signed target n. In alternative1, leave q unchanged. In
alternative2, append JJ:

    (5n,2n,6n,3n).JJ=(17n,6n,20n,9n)=W(v),       (7)

where v is the balanced precursor in Theorem1. This is a literal physical
extension and preserves `6U-7D=n` exactly on that exceptional branch.
No universal target-preservation identity is asserted off the branch.

The inverse produces an actual source with unsigned D=n. In alternative1,

    14 U_s < U_q,

because `U_q=6D_s+14U_s` and D_s>0. In alternative2 it produces balanced
v with `6U_v=U_q`. Therefore the entire intrinsic inverse has cover height
at most U_q/6, with equality on the exceptional branch. This is a contraction
for an already supplied actual slice witness, not a bound proving that such
a witness exists at every prescribed target.

The physical graph inverse is constructive: peel the forced terminal J
births in the proof and apply the canonical balanced graph inverse. A
profile-only implementation can instead enumerate actual source profiles
through the proved smaller cover bound, retaining a witness with the required
D,U. Neither implementation replaces source membership by cone membership.

## 5. Exact positive spectrum and the remaining occupancy lemma

The complete source grammar immediately gives

    {D_s:s in S*}={1} union {F(s):s in S}.        (8)

Indeed O preserves D, while each Z or J birth replaces D by F of its physical
parent. Conversely `D(Zs)=F(s)`; the singleton supplies1. The strong C1
spectrum is F(S) by the inherited physical classification.

Combining (3), Theorem1 and (8) gives the exact scalar equality

    {6U_q-7D_q:q in T}={1} union F(S).           (9)

Thus the positive C1 problem can equivalently be posed as intersection of
the globally surjective signed-response fibers with the intrinsic physical
slice (4). The theorem proves an actual section and inverse, not just an
existential encoding of physicality.

For any proposed target set A of positive odd integers, the following is a
precise sufficient, and in this C1 setting necessary, arithmetic lemma:

> For every n in A, there is an actual source profile q satisfying
> E_q=2(U_q-D_q), d_q=3(U_q-D_q), and 6U_q-7D_q=n.

If this occupancy lemma holds, (9) realizes every n in A\{1} by a strong
C1 tournament. All source and positivity conditions must be retained.
The full signed-response theorem does not establish this simultaneous
intersection. For example the actual unit source `t.OZ=(5,1,7,1)` already
has signed response7, but it does not satisfy (4).

There is also a complete all-height hostile at7 for a mistaken intersection
claim. Any source with D=7 can have its O tail stripped. Its positive-endpoint
core would have `7<=U<=11`: a last Z gives `2D-U>=3d>=3`, and a last J gives
`2D-U>=3` as well (use D_old>=3, or separately D_old=1). Every core with
U<=11 is covered by the complete source bank at that height. That bank has
no D=7 core. Hence D=7 is absent from S*, and (9) rules out any slice witness
with response7, despite the explicit ordinary signed witness above.
This small exact inventory is an all-height D=7 proof after the core bound;
it is not a claim that an arbitrary finite response search proves absence.

## 6. Verification

The verifier checks all rank-two identities on a full formal basis, both
formal-boundary transfers, the complete physical slice through a declared
cover cap, its decomposition into the two constructive families, the unique
D,U inverse, both sharp height statements, the forbidden intermediate chamber,
and the complete D=7 core control. It separately counts all four source
response coordinates and the root-start count on literal small graphs.
Normal and optimized runs must agree. Independent proof review and separate checks are preserved in audit/.

The new progress is an intrinsic physical rebase with a bounded inverse.
The unresolved existence problem is visible as the two simultaneous linear
constraints (4) on each signed fiber; it has not been removed.
