# Exact continuation fibres and a physical one-step-hidden repair

2026-11-12. **PROVED:** the R/N observation identities, positive-cone
separation, and exact equality language below. **FINITE-EXACT:** the small
physical-source probe. The positive infinite repair theorem is recorded
separately in [EXACT_REPAIR.md](../EXACT_REPAIR.md). Chan–Pak remains open.

## 1. Inherited response, new relative question

For a nonempty rooted tournament use z=(h,s,q,F), with h the full Hamilton
count, s the root-start count, q the root-deletion count, and F the number
of ordered spanning covers (P,Q) with P nonempty ending at the root and
Q possibly empty. The actual root-child commands are inherited from
[ROOTED_RESETS](../../two_port_return_induction_20261103/ROOTED_RESETS.md):

    R(h,s,q,F)=(h+F−s,q,q,2F−q),
    N(h,s,q,F)=(h+F−s,h−s,h,h+F).

Both add x with root→x→all other old vertices; R retains the root and N
marks x. Words are executed left to right. Four-coordinate sufficiency
for the inherited interface is not claimed anew. The present question is
which physical source changes preserve specified *future* full counts.

The older exact-fibre model is July's **THM-2230**,
`01-canon/theorems/THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient.md`:
under its unit-mate hypothesis, a fixed planar Jacobian response has exactly
the target-shear fibre Q_R+k[P]. Here the response is a different, elementary
linear map, and physical realizability of a point in its fibre is an
additional obligation.

## 2. Four literal count observations recover the rooted profile

Write H_w(z) for the full count after word w, including H_empty=h. Direct
substitution gives

    H_R=H_N=h−s+F,
    H_RR=H_RN=h−s−2q+3F,
    H_NR=H_NN=h+2F.

Consequently the four actual observations h,H_R,H_RR,H_NR recover

    F=(H_NR−h)/2,
    s=h+F−H_R,
    q=F+(H_R−H_RR)/2.                              (1)

Thus two actual rooted sources agree on every R/N continuation count iff
their four rooted coordinates agree. Agreement on all words implies
agreement on these four observations; the converse follows by the state
maps. This classifies the observer, not rooted graph isomorphism.

In particular a difference preserving both the current count and its first
birth target has exactly the form

    Delta=(0,v,w,v).                                (2)

The two next branch responses are

    Delta H_RR=2(v−w),      Delta H_NR=2v.           (3)

Both vanish only when v=w=0. A nonzero preparation in this relative fibre
necessarily becomes visible in one of the two next branch counts. The
formula does not forbid preparing for one chosen future target; it forbids
treating a nonzero hidden preparation as neutral for both branches.

The formal two-dimensional fibre is not automatically physically filled.
The small two-arc witness below realizes its pure deletion direction w.
The independent [marked-block lane](../ports/PROOF.md) realizes a different
direction v at its own baseline. Combining their numerical directions does
not by itself combine their graph guards or furnish an arbitrary selector.

## 3. A positive coordinate cone separates every nonempty continuation

Put a=h−s and d=F−q. Each coordinate of (h,a,q,d) is nonnegative on an
actual source: s≤h, and the cover with first path equal to the singleton
root injects root-deletion paths into the covers counted by F. The maps are

    R(h,a,q,d)=(a+q+d, a+d, q, 2d),
    N(h,a,q,d)=(a+q+d, q+d, h, q+d).                (4)

These are nonnegative integer maps. Suppose a pair of rooted sources has
difference (Delta h,Delta a,Delta q,Delta d)=(0,A,0,D), with A,D>0.
After either first command, Delta h>0 and Delta a,Delta d>0. Equation(4)
preserves those strict a,d inequalities and all coordinate inequalities.
Therefore **every nonempty R/N word strictly increases its full count**.

This applies to moving a fixed root forward inside W(TT_n): h and q stay
fixed, s strictly decreases, and F strictly increases, by the inherited
marked block formula and the [constructor proof](../constructors/PROOF.md).
In particular, the count-neutral H17 preparation cannot be commuted past
a selective birth while retaining that same next exact target. This is an
all-word statement for the declared alphabet, not a prohibition on other
native repairs or on mixed-sign differences such as(2).

## 4. A two-arc physical point in the hidden fibre

Use root0 and outgoing rows

    K =(16,9,3,5,14),       z =(15,3,3,21),
    K'=(18,24,3,5,12),      z'=(15,3,5,21).

Reverse the two old arcs 1→0 and 4→1. Both graphs are strong, and their
difference is exactly Delta=(0,0,2,0), realizing(2) with v=0,w=2.
The sources are the inherited cyclic triangle wrapper and regular
five-tournament of [RESET_JUNCTION J2](../../coverage_machine_induction_20261103/RESET_JUNCTION.md),
with the wrapper rooted at its sink-before-the-new-root. The new interface
uses a two-arc change and retains the same physical root and first child:

    Rz =(33,3,3,39),
    Rz'=(33,5,5,37).                                (5)

The new child has identical incidence in both graphs. Their exact named-pair
continuation divisors are3 and1, respectively, by the inherited root-child
gcd theorem. Thus(5) is an actual repair at the unchanged next count33.
The root's old indegree is not preserved.

For every k≥1 the repaired R^k N R returns the inherited balanced family

    b=48*2^k+11,       t=16*2^k+1,       gcd(b,t)=1.

The unrepaired source has delta=6,d=18,q=3, so its corresponding gcd is
gcd(6,18*2^k+3)=3 at every k≥1. The all-length counts themselves were
already recorded in RESET_JUNCTION J4. The gain here is their preparation
inside an exact one-birth fibre, which the parent theorem lifts to all
admissible nonempty left factors, including the infinite target family
21*2^a+15.

## 5. Exact future equality language for a pure deletion repair

For any two rooted profiles differing by w e_q=(0,0,w,0), w>0, the words
whose full counts remain equal are **exactly**

    {empty,R} union N{R,N}*.                         (6)

Proof: N e_q=0, so a first N erases the entire difference, and every later
state agrees. A single R sends e_q to (0,1,1,−1), preserving the count but
retaining a different rooted response. For j≥1, direct iteration gives

    R^j e_q=(2−2^j,1,1,1−2^j)

in (h,s,q,F). Its h coordinate is zero only at j=1. In the coordinates
of(4) this is

    (2−2^j,1−2^j,1,−2^j).

If the first N follows j≥1 copies of R, the difference becomes

    (2−2^(j+1),1−2^j,2−2^j,1−2^j),

which is coordinatewise nonpositive, with h,a,d strictly negative.
Equation(4) preserves strict negativity of the ensuing h,a,d coordinates.
Thus every word beginning R with length at least2 strictly lowers the
repaired full count. Multiplication by w proves(6).

This gives three distinct operational cases: N forgets the repair; R keeps
it at unchanged target; a further R/N command reveals its signed response.
An observation kernel and a continuation-stable kernel are different
objects. This is the useful transition-sidecar lesson of the old Jordan
and principal-parts work; no identification with an algebraic Hodge group
or a Prüfer module is asserted.

## 6. Finite probe and scope

`hidden_probe.py` examines every labelled tournament of orders1 through5
with root0. Relabeling covers any other choice of root. It checks every
single-arc and simultaneous two-arc change, retaining pairs whose H and
first birth H are equal but whose rooted profiles differ. Literal vertex
permutations and cuts compute the profiles; an independent subset recurrence
recounts the selected witness. The receipt records the precise universe,
profile-transition classes and witnesses. No statement about larger orders
or complete target coverage is inferred from the finite probe.

The formal fibre equations and positivity theorem hold at every order.
The two-arc witness supplies one physical generator, not arbitrary motion
inside every formal fibre. Recursive exact repair still needs a physical
source and compatible operation at the next chosen absolute target.
