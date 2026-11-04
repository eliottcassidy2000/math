# Incoming Session181/182: exact interactions with the C2 channel

**PROVED scoped bridge boundary; FINITE-EXACT permutation certificates.**
Reviewed incoming commit `b486271d24`. This addendum leaves the previously
pinned `RESULT.md` and `algebra.py` unchanged.

## 1. Provenance corrections and compatible scopes

The complete prototype two-birth menus in our `RESULT.md` section 6 are already
proved and independently audited in incoming
`session181_multiport_defect_20261104/grammar/PROOF.md`, section 3, including
all moving output marks. They should be presented as an independent replay of
that now-shared result, not as a new discovery. The raw weak-suffix H59 hits
in the companion transport packet do not contradict either strong-filtered
menu: both contain a global sink.

Session182's general clone-run formula in `defect/PROOF.md`, section 3,
contains our clone-block H formula as its two-run specialization. For the
six-vertex C2 prototype and its y module, `A1=17,A2=6,A_r=0` for r>=3 gives
`H_m=5+6*2^m`. Our additional retained-mark formula is `F_m=7+7*2^m`.
Its isolated m=1 to m=3 coincidence is compatible with the incoming fixed-module
obstructions; it is not a returning +2 rule.

Incoming `closure/PROOF.md` and `closure_generation/PROOF.md` continue the
actual S181 order-17 multiport H512425 state with defect3690. They do not
identify it with our earlier order-18 zero-ready H512425 state, or with any
small prototype. Their physical grandchildren are genuine additional source
lineages, with separately retained root defects and native guards. No correction
to the C2 operator is needed because of these new finite realizations.

## 2. The two H9-to-H31 constructions are not the same graphs

Incoming `session182_physical_returns_20261104/induction/PROOF.md` gives

    T_t=M(TT_t; TT3, {b0,b1,b2}),       H(T_t)=7*2^t-5,
    U_t=M(TT_t; TT5, {b0,b2,b4}),       H(U_t)=21*2^t-11.

Thus `H(U_t)=3H(T_t)+4`. The suffixes here are transitive; terminal port b2
or b4 makes the full attachment strong. This is outside a grammar that insists
the suffix itself be strong. Its seven-coordinate signature is also not our
four-coordinate C2 signature. In this section T_t is treated as a *whole core*
when computing the ordinary single-vertex responses `E_w(T_t,2)`.

At t=1, use source index order `(b0,b1,b2,a0,r)` and target index order
`(b0,b1,b2,b3,b4,a0,r)`. The incoming rows are

    T1=[22,20,16,7,8],
    U1=[94,28,88,16,64,31,42].

Our C2 +4 preparation and output are

    B=[20,1,2,7,14],       C2(B,2)=[116,97,66,103,110,4,32].

Equal H values do not give isomorphisms. The multisets of the actual one-port
two-slot responses are invariant under any unrooted tournament isomorphism:

| Graph | H | Multiset of E_w(2) over every vertex |
|---|---:|---|
| incoming T1 | 9 | {9,10,15,16,25} |
| our B | 9 | {9,11,11,21,23} |
| incoming U1 | 31 | {31,39,43,71,89,157,183} |
| our C2 output | 31 | {31,37,39,63,87,145,163} |

Both the sources and the targets are therefore nonisomorphic, even before
requiring named-port preservation. These values are independently recomputed
by complete permutation enumeration in `incoming_bridge.py`.

## 3. No member of the incoming source family has a C2 +4 preparation

For a physical auxiliary a in T_t, put

    sigma_a=S_a(T_t,2)-H(T_t).

For a separately retained physical mark w, let `tau_(w,a)=G_(w,a)-F_w`, with
the exact conventions of `RESULT.md`. A C2 all-clock +4 preparation requires
`sigma_a=2,tau_(w,a)=0`.

**Claim.** For every t>=1 and every pair of physical vertices w,a of T_t,
this conjunction fails.

Write the prefix `a0 -> ... -> a_(t-1)`, suffix `b0 -> b1 -> b2`, and root r.
The remaining orientations are prefix→suffix→r→prefix.
There are five auxiliary types:

1. Every prefix vertex has at least three out-neighbors, and so does b0.
   The dual rooted-defect bound in
   `two_entry_funnels_20261103/READY.md`, section 1, gives
   `sigma_a >= d_out(a)(d_out(a)-1) >= 6`.
2. Vertex b1 has out-neighbors b2,r. Deleting b1 leaves the strong cyclic
   blowup with blocks TT_t, TT2 and {r}. In the dual defect formula
   `sigma_b1=delta+q-e_b1`, delta>=2. Also q−e_b1>=1: a Hamilton path of the
   deletion ending at b2 exists by strongness, and cannot be extended by b1
   because b1→b2. Thus `sigma_b1>=3`.
3. Vertex b2 has sole out-neighbor r. The one-out-neighbor defect identity
   gives `sigma_b2=e_r(T_t-b2,1)=1`. Indeed a Hamilton path ending at r is
   forced to traverse the transitive tournament TT_t→TT2, then r.
4. Root r has out-neighbor set equal to the whole prefix. Deleting r leaves
   the transitive order prefix→suffix. In its bad-join decomposition, the
   suffix must lie in the path starting at a prefix vertex. The prefix splits
   between two nonempty slots in `2^t-2` ways. Therefore

       sigma_r=2^t-2.

   It equals 2 only at t=2.

This exhausts every auxiliary, independently of w. At the sole remaining
clock t=2, the source rows, in order `(b0,b1,b2,a0,a1,r)`, are

    [38,36,32,23,7,24],       H=23.

With a=r, the six possible retained marks have

    tau_(w,r)=[4,8,19,1,1,3],

in that order. All are positive. This finite last case is checked by complete
permutation counts on the source and its literal C2 child; the latter has H73
and old-mark E2 values `[128,196,299,74,101,93]`, whereas the source has
`[40,60,87,24,33,29]`. Subtraction of three times the source vector gives
twice the displayed tau vector. This proves the claim at every t.

Even arbitrary output re-marking does not repair this particular C2 branch.
At t=2 its complete endpoint-response list is
`[128,196,299,74,101,93,563,601]`, including the two newborns. It is disjoint
from three times the source list. Thus no choice of input and output marks
restores the missing response-pair identity.

The proof does not exclude a +4 construction from T_t. The incoming operation
already supplies one by changing an old root incidence and the suffix port
set. It proves only that this incoming family cannot be reinterpreted as a
C2 +4 preparation with an unchanged old mark.

### A positive bridge: the incoming source family is a C2 +2 preparation

The conclusion is specific to +4. At the retained physical root w=r and
auxiliary a=b2, every T_t has `sigma=1,tau=0`. Applying C2 replaces the old
terminal suffix vertex b2 by the triangle `b2 -> y -> x -> b2`, with b2 the
only port inside that triangle. The actual suffix becomes `TT2 -> C3`, with
port set `{b0,b1,b2}`. Its seven-coordinate signature is

    (h,p,f,g,j,e,s)=(3,48,21,39,16,1,3).

In particular e=1: every full suffix path ends inside the terminal triangle,
whose sole marked vertex has one Hamilton path ending there.
The whole-state responses are

    H(C2(T_t,b2))=-13+21*2^t=3H(T_t)+2,
    E_r(C2(T_t,b2),2)=3+21*2^t=3E_r(T_t,2).

Together with the incoming U_t rewrite this gives both digits from every
member of this incoming source family, by two different literal programs.
The +2 program is precisely the converse ready reset, not an alternative to
C2. Its outer root defect changes 6 to 16. The new auxiliary x has start
defect zero, a different physical gauge. The physical strong schedule and
shared-clock proof are recorded in `../transport/TWO_DIGIT_ENTRY.md`.
Return of either output to this source family is still unproved.

## 4. The same H73 source/target comparison exposes the missing coordinate

At t=2 the source T2 is literally the same in the two operations. C2 at its
root auxiliary gives H73, as does the incoming port-changing operation.
But at the retained physical root r,

    E_r(T2,2)=29,
    E_r(C2(T2,r),2)=93,
    E_r(U2,2)=85.

The incoming target has defect12; the C2 target has defect20. More generally,
the incoming family's root response satisfies

    E_r(T_t,2)=1+7*2^t,
    E_r(U_t,2)=1+21*2^t=3E_r(T_t,2)-2.

This negative endpoint correction is unavailable to C2 at a fixed old mark,
whose endpoint update is `F'=3F+2tau` with tau>=0. The incoming rewrite is
therefore a concrete example of the different operation needed to escape
that positivity restriction. It is not itself a proof that both output
branches return to a common recursively closed class.

## 5. Different tripling mechanisms have different obstructions

Incoming `gadgets/PROOF.md` triples a TT2 prefix into a homogeneous C3 module
with one birth. Its free order-three automorphism forces divisibility by 3
until the module is destroyed. Our closed C2 channel instead transports its
auxiliary to a fresh vertex and retains a joint response guard.

In a C2 channel state, the old auxiliary a has a sole out-neighbor b. At the
next step `{a,x,y}` is a directed triangle, but it is **not homogeneous**:
`a->b`, while `b->x` and `b->y`. Thus the incoming homogeneous-module theorem
cannot simply be assigned to this channel. Its relevant obstruction to
nonzero C2 digits is the inherited READY re-marking theorem. Neither mechanism
supplies nonzero-digit closure merely by supplying multiplication.

## 6. Zero-ready plus unit start is a different, rigid boundary

Let R be a **nontrivial strong** tournament with physical mark w satisfying

    E_w(R,2)=H(R),       s_w(R,1)=1.

The inherited zero-ready classification gives w a sole in-neighbor a, with
`D=R-w` weak and a in its terminal SCC. No Hamilton path of D starts at a,
since all such paths start in the first SCC. Deleting the initial w therefore
gives `s_w(R,1)=H(D)`. Thus H(D)=1. A tournament has H=1 iff it is transitive:
a directed triangle otherwise contributes three different paths, which can be
extended injectively by deterministic insertion of every other vertex.

Write the transitive order of D. Its terminal vertex must be a; otherwise that
vertex is a global sink in R. Conversely, a transitive D with a terminal, and
w beating every other vertex, has exactly the required strong, zero-ready,
unit-start interface. Hence

    R is C3[TT_m,1,1],       H(R)=2^m+1,       m>=1.        (B1)

The scalar count follows by splitting the transitive block before and after
the other two cycle vertices, as in the incoming triangle-blowup proof.
Consequently an exact tripling within this boundary requires

    2^(m')=3*2^m+2.

For m>=2 the right side has 2-adic valuation one and is greater than2, which
is impossible for a power of two. For m=1 it is8, giving only `(m,m')=(1,3)`.
Thus **3 to 9 is the only tripling step inside this specified boundary**, and
there is no R of H27 satisfying it. This is not a general ternary-closure
obstruction.

The distinction is physical: C2 applied to C3 rows `[2,4,1]` with w=0,a=2
gives `(H,F_w,s_w)=(9,9,3)`. The incoming T1 at physical prefix vertex a0
instead has `(9,9,1)`. The latter belongs to (B1), whereas the former already
left its unit-start boundary. The singleton is deliberately excluded from
the classification; it satisfies the two scalar conditions but not m>=1.

Reproduction:

```text
python 05-knowledge/results/multiport_defect_grammar_20261103/construction/incoming_bridge.py
python -O 05-knowledge/results/multiport_defect_grammar_20261103/construction/incoming_bridge.py
```

The recorded four clock controls verify the formulas and matrix conventions.
The all-clock exclusion follows the auxiliary-type proof, not extrapolation
from those controls. No native complement guards are claimed for the incoming
TT3-to-TT5 family beyond its own stated adjacency-level schedule.
