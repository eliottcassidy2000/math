# Canonical sink cores and an exact deficit-fiber compiler

Session 200. **PROVED in the stated scope; independently audited.**

Every physical joint fiber admits a finite-core/dyadic-tail normal form at a
fixed deficit `r=m-h`. For a supplied physical source profile, the sink-tail
length and the core's three retained coordinates are canonical; the joint
pair alone need not determine them. The new compiler produces actual source words,
not formal cone solutions. It gives explicit infinite positive target families
and a quantitative obstruction to covering targets with uniformly small
physical deficit. Universal `{7,21}` coverage remains **OPEN**.

## 1. Inheritance and the arithmetic coordinate

Use the complete physical source orbit S from `t=(1,0,2,1)` under

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

Words are chronological. The physical cone supplies

    D-E>=d>=1, D<=U, U>=2d, 0<=E<=d,

and D,d are odd. Put

    F=D-E+U, G=D+2U, H=2U-d, r=G-H=D+d.

The inherited joint compiler proves `3<=H<G<2H` and

    F(S)={3} union {G(s)+2(2^k-1)H(s):s in S, k>=0}.       (1)

At k=0 this is the balanced readout G. The source word `s.O.Z^k`, followed
by a closing J, realizes the displayed target. Conversely terminal-run
analysis proves (1), including when a last J is rebased to a balanced word.
See `CARMICHAEL-DEPENDENCY.md`, Sections 3 and 7, and the independently
audited common-word inverse in `JOINT-DEPENDENCY.md`.

Set `M_k=2^(k+1)-1`. The target equation becomes

    n=M_k h+r, h=H(s), r=D_s+d_s.                       (2)

The retained deficit is not an independent integer. For example the
unsynchronized pair `(m,h)=(15,13)` passes the open joint cone but is not
physical. The present construction gives a direct arithmetic reason.

Concept board: common-word synchronization; source-fiber occupancy; physical
boundary faces; dyadic tail; finite core; target-dependent deficit. The
Carmichael transfer is to remove a universally available fiber operation and
identify its exact residual chamber. The Diophantine transfer is to keep a
literal physical core certificate behind the integer parameters.

## 2. A sharp physical boundary and a canonical tail

**Lemma 1 (endpoint/height dichotomy).** Every source s in S satisfies

    E=0 iff U>=2D.

If E>0, then

    D<=U<=2D-3.                                      (3)

The seed has E=0 and U=2D. In particular there are no physical sources with
U equal to `2D-2` or `2D-1`.

**Proof.** A last O gives E=0 and `U'=2U>=2D`; conversely O is the only
letter with zero output endpoint, and the seed has the stated equality.
A last Z gives

    2D'-U'=2(D-E)+d>=3d>=3.

For a last J,

    2D'-U'=D+U-2E>=D,

because U>=2d>=2E. If the old D>=3 this proves the bound. If old D=1,
then `D-E>=d>=1` forces E=0,d=1 and U>=2, so the same difference is at
least 3. The inherited `D'<=U'` supplies the lower bound. QED.

A *core triple* is `(D,u,d)` with a physical realization `(D,E,u,d)` having
E>0, together with the exceptional seed triple `(1,2,1)`. Store one literal
physical recipe for every retained triple. Call the set of these triples C.
No claim is made that an arbitrary integer triple satisfying (3) is in C.
Nor is uniqueness of the forgotten endpoint coordinate asserted.

**Theorem 2 (canonical sink-core projection).** Every source s has a physical
presentation

    s=q.O^a, a>=0, (D_q,U_q,d_q) in C.                (4)

For D_s>1, every presentation of s has exactly

    a=floor(log_2(U_s/D_s))                          (5)

terminal O letters. The core projection is uniquely

    (D_s, U_s/2^a, d_s).                             (6)

For D_s=1 the source is transitive: `s=t.O^a`, `U_s=2^(a+1)`, d_s=1.
The triple and a in this exceptional case are again unique.

**Proof.** Strip the maximal O suffix from any physical word. A nonempty
remaining word ends in Z or J, so its endpoint is positive; an empty
remaining word is the seed. O leaves D,d unchanged and doubles U. If
D>1, the core cannot be the seed, and (3) puts its cover count in `[D,2D)`.
There is exactly one dyadic rescaling of U in that half-open interval,
proving (5)-(6), independently of the chosen word. If D=1, no Z or J can
occur: each raises D by `U-E>0`, and no letter decreases D. Thus only the
seed and O letters occur. QED.

The theorem is canonical on the *three-coordinate projection*. O forgets E;
it would be an unjustified strengthening to claim a canonical complete
profile or a unique graph from this argument.

## 3. Complete finite cores at a fixed deficit

Let `C_r={(D,u,d) in C:D+d=r}` and `C_<=R=union_(r<=R) C_r`.
Only positive even r occur. For a nonseed core,

    u<=2D-3<=2r-5.                                  (7)

**Lemma 3 (complete core inventory).** For integer R>=2, C_<=R consists
exactly of the seed triple and the projections of

    Z(p), J(p), where p in S has U_p<=R-2,
    and the output has D+d<=R.                      (8)

Duplicates of triples can be merged while keeping a physical recipe.

**Proof.** Every generated output is physical and has positive E. Conversely
write any nonseed core as a word whose last letter is Z or J. In both cases
its predecessor obeys

    D_q=D_p-E_p+U_p>=d_p+U_p,

so `U_p<=D_q-1<=R-2`. This includes every core. QED.

Both cutoffs are sharp along `q=t.O^j.Z`, j>=0: with
`r=2^(j+1)+2`, its predecessor has `U=r-2` and the core has `u=2r-5`.
This is an all-order sharpness family, not a finite maximization claim.

A coarse uniform cardinality estimate useful below is

    |C_<=R| <= sum_(D=1)^(R-1) D(R-D)
             = (R^3-R)/6.                          (9)

For D>1 there are at most R-D possible positive d values and at most D
possible integer u values in (3). These bounds ignore parity and other
physical constraints, so are upper bounds. At D=1 the single exceptional
seed triple is covered by the term R-1. The endpoint coordinate is not
counted because it has been legitimately removed from this compiler.

## 4. An exact dyadic-ray description of the physical joint image

**Theorem 4 (fixed-deficit joint fibers).** For integers m,h, put r=m-h.
Then `(m,h)=(G(s),H(s))` for an actual source s if and only if there exist
`(D,u,d) in C_r` and a>=0 such that

    h=2^(a+1)u-d, m=h+r.                           (10)

In the constructive direction, choose a stored physical core q and use
`s=q.O^a`. Conversely Theorem 2 gives exactly these equations.

This is an existential description of each joint pair, not a unique
factorization of that pair. For example `t.ZZZJJ=(27,9,35,17)` and
`t.ZZJJJ=(29,10,34,15)` both have `(m,h)=(97,53)` and a=0, but their core
projections `(27,35,17)` and `(29,34,15)` are different. The canonical
projection in Theorem 2 requires the supplied source's D,U,d, not only G,H.

Thus at each fixed r the entire joint image is a finite union of explicit
dyadic rays. Core physicality is checked at height `R-2`, not at a height
proportional to h. This can be smaller than the inherited square-root inverse
when r is small. It is not uniformly smaller for arbitrary pairs.

The same statement preserves minimum J cost for the *joint pair*: store the
minimum J cost over all physical realizations of each core triple. Stripping
or adding O costs zero. Merging cores with the same triple loses no G/H data,
and the subsequent k>=1 target recipe begins with O, which also removes E.
This does not license merging by a scalar or using independently chosen
core coordinates.

The first complete slices are

    r=2: C_2={(1,2,1)},       h=2^(a+2)-1;
    r=4: C_4={(3,3,1)},       h=6*2^a-1;
    r=6: C_6={(5,5,1),(5,7,1)}, h=10*2^a-1 or 14*2^a-1. (11)

Their recipes are respectively empty, Z, ZZ, OZ. Completeness for r=2 uses
the transitive clause of Theorem 2. For r=4,6, positive E and `D-E>=d`
exclude D=d, leaving d=1 and D=3 or 5. Then E=1, u is odd, and (3) gives
exactly the displayed u choices.

The pair (15,13) has r=2, but 13 is not `2^(a+2)-1`; it is rejected without
an unsynchronized boundary search. The parity-correct positive-cone tuple
`(D,E,U,d)=(3,0,4,1)` is also nonphysical: it violates Lemma 1. Its putative
pair `(11,7)` has r=4 and fails the second ray in (11). Cone feasibility is
therefore strictly weaker even with the usual parity constraints.

## 5. Constructive target normal form and its exact limitation

**Theorem 5 (core/double-dyadic target compiler).** The complete C1 scalar
image is

    F(S)={3} union {
      (2^(k+1)-1)(2^(a+1)u-d)+D+d:
      (D,u,d) in C, a,k>=0 }.                      (12)

For every choice, a stored core word w compiles the target to the literal
source `t.w.O^(a+1).Z^k` and a closing J. Conversely combine (1) and
Theorem 4. This is a necessary-and-sufficient arithmetic normal form with
physical certificates, not a claim of solvability for arbitrary n.

In particular (11) supplies explicit infinite target families, for a,k>=0:

    n=M_k(2^(a+2)-1)+2,
    n=M_k(6*2^a-1)+4,
    n=M_k(10*2^a-1)+6,
    n=M_k(14*2^a-1)+6.                             (13)

The first is the transitive/Mersenne family. The next three retain distinct
physical cores rather than pretending that arbitrary h values with the same
small deficit are available.

For an actual source profile x whose scalar F(x) lies outside the balanced
scalar image `B={3} union G(S)`, the inherited canonical Carmichael descent
gives

    k=v_2(U_x-d_x), h=(U_x-d_x)/2^k,
    m=F(x)-2(2^k-1)h.

The parent projection is then uniquely reconstructed as

    d=d_x, U=(h+d)/2, D=m-h-d.                     (14)

Theorem 2 supplies its unique sink-tail length and core projection. Thus the
two dyadic exponents and the physical core projection are canonical for each
actual source profile x satisfying `F(x) not in B`. They need not be unique
for the scalar F(x), whose fiber can contain different profiles. For example
443 has the two inherited canonical outer exponents 5 and 2, not 1.

The hypothesis is scalar nonmembership in B, not the profile condition
`D_x!=U_x`. The actual source `x=t.ZJ=(5,2,6,3)` has unequal D and U but
`F(x)=9=G(t.O)`, hence lies in the balanced scalar image. Here
`v_2(U_x-d_x)=v_2(3)=0`; it is not an instance of the positive-outer-exponent
descent. Theorem 2's sink-core normal form itself still applies to every
physical source without this scalar restriction.

**Theorem 6 (small-deficit target obstruction).** Let A_R be {3} together
with the values of (12) from cores with D+d<=R. For N>=3 and integer R>=2,

    |A_R intersect [1,N]|
       <=1+(R^3-R)/6 * (1+floor(log_2 N))^2.        (15)

For each fixed R, A_R has natural density zero. More generally, if
`R=R(N)=o((N/(log N)^2)^(1/3))`, then A_R contains o(N) targets through N.
Consequently uniformly bounded physical deficit, or that subcritical growing
budget, cannot cover a positive fraction of odd integers. Any all-target
proof using (12) must allow the core to grow beyond that budget for almost
all targets, whether it proves C1 coverage or only a cofinite subfamily.

**Proof.** For any term in (12), `h=2^(a+1)u-d>=2^(a+1)` because
u>=2d and d>=1. Also `M_k>=2^k` and h>=3. Thus n<=N forces both a and k
at most floor(log_2 N). Count choices using (9); collisions only lower the
number of targets. The asymptotic conclusions follow. QED.

This does not obstruct the unrestricted source spectrum, and says nothing
by itself about C2/C3 repairs. It pinpoints a lost quantifier: a finite bank
of genuine cores supplies infinitely many physical targets, but not
cofinite target coverage. Arithmetic definability or synchronized modular
solutions cannot fill this deficit automatically.

## 6. Verification and remaining frontier

`verify.py` uses explicit exceptions, checks every word through a declared
length, compares complete forward joint and F inventories with an independently
constructed bounded core bank, retains exact minimum J costs, checks (11),
replays both nonphysical hostiles, and counts paths/cycles on literal graphs.
Run normal and optimized Python serially; their JSON receipts must agree.
Finite checks are implementation evidence, not the proof of the all-order
lemmas. The independent audit in `audit/REPORT.md` passed after the two scope
clarifications recorded there.

The remaining constructive problem is now concrete: given an odd n, choose a
*physical* core of target-dependent deficit and exponents a,k solving (12).
The synchronization requirement has been moved into a finite core with an
honest height bound, not removed. Universal coverage and the all-order
singleton-fiber conjecture both remain open.
