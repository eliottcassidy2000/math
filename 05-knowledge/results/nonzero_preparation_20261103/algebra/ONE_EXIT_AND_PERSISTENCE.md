# One-exit joint preparation, persistent old names, and one physical clone return

**PROVED — [independent audit](AUDIT.md) passed.** The results
are the exact one-exit joint identity in §2 and an all-clock, all-depth application
of an inherited permanent-readiness certificate in §§3–4. Section 5 is a new
application of the inherited clone formula, not a new clone mechanism. It gives
one actual further `+4` return and proves that its particular enlargement cannot
iterate. Neither nonzero two-branch closure nor global spectrum coverage is
claimed.

## 1. Inheritance, objects, and the missing coordinate

The closest proved operator is
[construction/RESULT.md §§1–5](../../multiport_defect_grammar_20261103/construction/RESULT.md):
the square-zero vertex path algebra, the two-birth C2 operator, its transported
zero channel, and its transitive-clone formula. C2 is the converse of the older
ready reset. Its arbitrary-re-marking obstruction already excludes immediate
reuse; we do not reprove that as a new result.

The strongest relevant sidecar is the actual named certificate in
[interfaces/NOTE.md §§I1–I2](../../balanced_ternary_returns_20261103/interfaces/NOTE.md).
It proves conditional endpoint monotonicity and sufficient certificates for
nonreadiness under **every induced extension**, not only under a chosen birth
menu. We apply its converse, with every witness and triangle named below.

The actual starting family is
[PREPARED_FOUR.md](../../multiport_grammar_20261103/component/PREPARED_FOUR.md),
including its fixed retained root. The canonical hostile is the mismatch between
Hamiltonian offset and retained-root response: `sigma>0` alone does not pay
`tau=0`. The corrected near miss is enlarging a triangle to get one more `+4`:
the clone exponent then fails the next required factor. The least-used sidecar
is an ordered cover by two paths with **two named final endpoints**.

The live concepts are: (i) the actual old adjacency relation; (ii) two named
marks; (iii) the joint ending count; (iv) endpoint witnesses surviving growth;
(v) a transitive clone clock; (vi) fresh versus already-saved auxiliaries.
The map in §2 sends an actual one-exit tournament to a deletion graph with two
named endpoints, retaining the exact C2 guard. Forgetting the second endpoint
loses the joint obstruction. The cheapest decisive test is literal path
enumeration; §7 exhausts that test through order five.

Use out-row encoding and compose paths in their directed order. For a tournament
T, let H(T) count directed Hamilton paths. For a retained root w and auxiliary a,
use the square-zero vertex ring with P including the empty path, A_a counting
nonempty paths starting a, W_w counting nonempty paths ending w, and Q_aw counting
paths starting a and ending w. Full-vertex coefficient extraction gives

    F_w=[W_w P], S_a=[P A_a],
    G_(a,w)=[W_w P A_a]+[P² Q_aw],
    sigma_a=S_a−H, tau_(a,w)=G_(a,w)−F_w.

The proved C2 formula is

    (H,F_w) -> (3H+2 sigma_a, 3F_w+2 tau_(a,w)).       (1)

All slots are ordered; P has constant coefficient one. These conventions matter
even for the smallest deletion graphs.

## 2. Exact one-exit preparation identity

**Theorem.** Let a,w be distinct vertices of any finite tournament T.
Suppose a has exactly one out-neighbour b, and put U=T−a. Then

    sigma_a=e_b(U),
    tau_(a,w)=J^end_(w,b)(U):=[W_w E_b]_U.             (2)

Here E_b is the path polynomial ending b in U. Thus J counts ordered spanning
covers by two **nonempty** disjoint paths ending w and b. In particular

    b=w  => tau=0;
    b!=w => tau>=sigma.                               (3)

Consequently positive preparation `sigma>0,tau=0` at an auxiliary distinct from
the retained root forces its unique exit to be the retained root. Conversely,
if b=w, the exact preparation value is the deletion ending count e_w(U).
There is no strongness assumption. This is an iff for this one-exit interface,
not a classification of arbitrary auxiliaries; the case a=w is not asserted.

Proof. All the following polynomials live on U. Write A_b for paths starting b
and Q_bw for paths starting b and ending w. The actual one-exit extension has

    P_T=P+z_a(P−E_b+P A_b),
    W_(w,T)=W_w+z_a P Q_bw,
    A_(a,T)=z_a(1+A_b), Q_(aw,T)=z_a Q_bw.              (4)

A path containing a either ends there, or continues immediately to b; this
proves (4), including the singleton a. Extracting gives

    H(T)=H(U)−e_b(U)+[P A_b], S_a=H(U)+[P A_b],
    F_w(T)=[W_w P]−[W_w E_b]+[W_w P A_b]+[P²Q_bw],
    G_(a,w)(T)=[W_w P]+[W_w P A_b]+[P²Q_bw].           (5)

Subtract to obtain (2). If b=w, repeated use of vertex w makes W_w E_b=0.
If b!=w, take each HP of U ending b and cut immediately after w. Both parts
are nonempty, ending w and b; restoring the removed join recovers the path.
This injects the e_b(U) paths into J^end and proves (3).

The equality tau=0 does not follow from outdegree one alone. Nor does the
transported zero auxiliary, which exits to the previous auxiliary rather than
the retained root, become positively prepared by changing its label.

## 3. The inherited permanent certificate, in the required orientation

For induced K⊂V containing a, the converse of NOTE §I1 says

    e_a(V)=0 or e_a(V)>=e_a(K).                         (6)

The zero alternative is essential. No strongness of V or K is assumed. The proof
is the endpoint-preserving insertion injection in that source, read backwards.

By the converse of NOTE §I2, any of these finite named certificates for an old
vertex v in T forbids sigma_v∈{1,2} in every induced tournament extension V⊃T:

1. v has at least three old out-neighbours.
2. It has exactly two old out-neighbours a,b, and T−{v,a,b} retains a directed
   triangle.
3. It has one old out-neighbour a, a triangle survives in T−{v,a}, and an actual
   induced K⊂T−v containing a has e_a(K)>=3.

For completeness, if outdegree reaches three, the inherited pair-defect bound
is at least six. With two exits and the triangle, the pair-deletion lower bound
is again at least six. With exactly one exit, (2)'s sigma identity and (6) give
zero or at least three. If a new second exit appears, deleting it leaves the
old certified triangle. These exhaust every future outdegree. The certificate
is sufficient, not necessary, and constrains sigma independently of the root
observer and its tau.

## 4. Every saved vertex of the actual prepared zero channel has a certificate

Let k=t−2>=1. Use the fixed actual six-vertex tournament

    D=(2,4,57,3,11,27), port w_D=4,
    (H_D,P2_D,F_D,G_D,J_D,e_D,s_D)=(17,224,28,44,0,2,4).

Let B=B_t be the actual bulk44 source of PREPARED_FOUR: the ordered suffix is
TT_k→D→u; the retained root r beats TT_k and D−{w_D}, and receives the arcs
w_D→r and u→r. In particular u has only the exit r. The inherited formulas are

    h=H(B)=63*2^t−29, f=F_r(B)=63*2^t+17,
    H(B−u)=28*2^k+17, e_r(B−u)=2,
    e_u(B)=H(B−u)−e_r(B−u)=28*2^k+15.                 (7)

The endpoint subtraction in (7) is literal: any HP of B−u can append u except
those ending r. The directed triangle D_0→D_1→D_2→D_0 will be retained throughout.

Append x with only exit u, then y with all B→y→x. The actual output T_0 has
terminal triangle u→y→x→u and

    H_0=189*2^t−83, F_0=189*2^t+51.                   (8)

Every old vertex of B other than u gains both exits x,y. Since B is strong,
each already had an exit, so certificate 1 applies. Vertex u has exits r,y;
deleting u,r,y preserves the D triangle, giving certificate 2.

For x use certificate 3 with exit u and K=B, whose ending count in (7) is at
least 71. For y use exit x and K=B+x. Every HP of K ending x is a path of B
whose final vertex is not u, so

    e_x(B+x)=h−e_u(B)=224*2^k−44>=404.                (9)

Both residual graphs preserve the same D triangle. Thus **every physical
vertex of T_0 is permanently nonready under arbitrary induced growth**.

Now perform j further C2 steps, always transporting the auxiliary to the new x
and retaining r. The inherited exact operator gives actual strong T_j with

    (H_j,F_j)=3^j(H_0,F_0), |T_j|=t+8+2j.            (10)

Old certificates persist. For each new pair x',y', x' exits only the previous
auxiliary a_j, with witness T_j. Its ending count is

    e_(a_0)(T_0)=h,
    e_(a_j)(T_j)=H_(j−1) for j>=1.                    (11)

Indeed a Hamilton path ending at the newly born x' in the C2 output must place
the other newborn y' immediately before x', preceded by an arbitrary HP of the
old T_j. This proves (11), also for the first preparation from B.
For y', use witness T_j+x', with ending count H_j−e_(a_j)(T_j). It is 2h+4>=954
at j=0, and 2H_(j−1) thereafter. The fixed D triangle survives all deletions.
Certificate 3 therefore applies to both new vertices.

Induction proves the claimed **all-clock, all-depth permanence**. Given any
finite saved T_j, arbitrary later vertex insertions and arbitrary incidence
choices cannot make any vertex already in T_j have sigma 1 or2. By (1), no such
auxiliary can supply a next C2 digit 2 or4, regardless of root re-marking.

This does not rule out a fresh auxiliary outside T_j, deletion of an old vertex,
or reversal of an old arc: each can evade the retained named certificate. It
also does not exclude every positive even offset. It is an operation/interface
obstruction, not an obstruction to a Hamiltonian count itself.

## 5. A literal non-C2 further +4, and its exact failure boundary

Keep B,x and replace the singleton y by the nonempty transitive block Y_a with

    B→Y_a→x, x→u, B−{u}→x, a>=1.

The old y remains the first block member; append each next member with only
exit x. Every old arc and name survives. Each birth is strong: the old strong
graph reaches the newborn, which reaches x and then the old graph. This is an
actual graph construction. No additional native macro/guard compilation is
claimed here.

The inherited all-size clone identity, applied to S_u(B)=h+2 and G_(u,r)(B)=f,
gives

    H_a=h+2^a(h+2), F_a=(1+2^a)f.                    (12)

In particular T_1=T_0 in the notation of §4, and two further clones give

    (H_3,F_3)=(3H_1+4,3F_1).                         (13)

At t=3 this is the literal step (1429,1563)→(4291,4689), following the earlier
157→475→1429 chain. It is not a C2 step at a rearmed old auxiliary. The terminal
module u→TT_a→x→u has H=1+2^a; its size is a+2 (a=3 gives H=9).

Within this specified enlargement family, simultaneous tripling of F requires

    1+2^(a+ell)=3(1+2^a),
    2^a(2^ell−3)=2.

For integers a,ell>=1 the unique solution is a=1,ell=2. Thus this physical
nonzero step cannot repeat by enlarging the same block. This proves only its
own failure boundary; it does not prove that a different larger operation
cannot continue.

All vertices in every T_a also have permanent certificates. Vertices in B−{u}
have at least three exits. Vertex u has a+1 exits (use certificate 2 at a=1,
certificate 1 thereafter). Vertex x uses witness B. A Y vertex with at least
three exits uses certificate 1; the penultimate Y vertex uses certificate 2
with exits lastY,x; the last Y vertex uses certificate 3 with witness B+x and
count (9). The D triangle is disjoint from all these deleted names.

## 6. What has and has not been repaired

The joint identity identifies exactly why a one-exit auxiliary aimed at an
older auxiliary cannot be a positive prepared port with the old root held
fixed: its two-endpoint cover debt cannot be ignored. The permanence theorem
strengthens immediate reset failure to arbitrary induced growth of the saved
family, using an already proved certificate rather than an extrapolation of
small cases. The clone step achieves one more actual digit while avoiding the
old-auxiliary interface, but consumes its one exact exponent coincidence.

The next unclosed operation must retain more than the old scalar response:
it may change an old named relation, or provide a genuinely fresh auxiliary
with both deletion-ending and two-ending-cover data paid. The unchanged-D
tail resource problem is being handled separately; no implication from its
scalar feasibility to physical construction is made here.

## 7. Exact controls, reproduction, and status

The standalone source is
`04-computation/nonzero_preparation_20261103/algebra_one_exit.py`.
It imports no producer and reads no old receipt. It enumerates all labelled
tournaments of orders 2–5 and every one-exit a and retained w!=a: the per-order
marked counts are 2,24,288,5120 (5,434 total). Every simple directed path is
enumerated literally, and square-zero convolution independently checks both
identities and the cut bound. This universe includes weak tournaments, empty
slots, the repeated-endpoint zero, strong positive sigma1/sigma2 controls and
a strong wrong-exit positive joint hostile.

A separate forward subset endpoint DP checks actual t=3,4,5, zero levels
j=0,1,2, all named certificates and their endpoint witnesses; clone block
lengths a=1,...,5; strongness at each birth; and the exact further +4. Its largest
graph has 17 vertices; stored path counts are at most n!<2^64, and convolutions
use unbounded Python integers. It also enumerates all 2^11 one-birth cuts from
the first T_0, checking the actual retained witness graphs and the appropriate
degree/triangle/witness certificate at every old vertex; this is not
misrepresented as 2^11 independent Hamilton-response counts.

From the repository root run both modes and compare their JSON bytes:

```text
python3 -B 04-computation/nonzero_preparation_20261103/algebra_one_exit.py --output /tmp/algebra_one_exit.normal.json
python3 -B -O 04-computation/nonzero_preparation_20261103/algebra_one_exit.py --output /tmp/algebra_one_exit.optimized.json
```

Both normal and optimized replays pass **44,726 always-active gates** and
produce identical 114,078-byte receipts. Frozen source SHA256:
`a2fa146f99ed10102d4a41cc699967575d8d2d2bdf726d820b72df4a253218c0`.
Frozen receipt SHA256:
`9e6a4dceca0bdaaf72697966313440331cc27764f68146ac49fc26cbb70c6570`.

Finite controls do not replace the all-order proofs or inherited quantified
certificates. The independent audit passed; see [AUDIT.md](AUDIT.md).
