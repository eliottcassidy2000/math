# A one-arc neutral preparation with accepted return clocks at every block size

**PROVED in this note; FINITE-EXACT literal controls.** For every n≥3, a
physical one-arc change preserves the exact count `2^n+1` and prepares an
infinite explicit family of primitive balanced returns. Every certified
return enters the inherited repeatable controller. This generalizes the
H17 internal-root preparation without assuming that arbitrary target counts
admit these sources. Chan–Pak and target-directed closure remain **OPEN**.

## 1. Inheritance, objects, and the actual neutral move

The closest mechanism is the
[H17 six-state preparation](../../eight_nine_twentyseven_20261112/constructors/SIX_RETURN.md).
Its hostile is an initially primitive rooted source with zero birth-return
defect: the subsequent `R^k N R` becomes imprimitive. The corrected near
miss is identifying a primitive final root-child pair with a primitive
returned root, or with exact target coverage. The least-used sidecar here
is the physical root's position in a homogeneous transitive module.

We inherit the marked block law
[NATIVE_SELECTION (5.3)](../../source_selection_closure_20261103/NATIVE_SELECTION.md),
the exact sealing law
[RESET_JUNCTION J3](../../coverage_machine_induction_20261103/RESET_JUNCTION.md),
and the native adjacent-swap compiler
[METAGRAPH §5](../../adaptive_metagraph_reformulation_20261103/METAGRAPH.md).
Their general formulas are credited; the new combination below is the
all-n preparation, factor-free accepted clock, and controller entry.

Incoming [Session298](../../session298_exact_target_geometry_20261113/README.md)
independently proves the same all-n W_n position formula and one-position
neutral preparation, with a generic multiplicative-order clock. Credit is
shared for that preparation. The explicit factor-free clock `k≡0 mod(n−3)`,
all-n C3 entry, native implementation, and distinction between the existing
pair and the two available roots are the additional refinements here. The
separate [retained-chain renewal](CHAIN_REUSE.md) below changes the ongoing
construction language and supplies an indefinitely reusable exact state.

Let W_n be the actual tournament

    X → Y → TT_n → X,

with module names `0,...,n−1`, exterior names X=n, Y=n+1, and retained
physical root0. Initially the module order is `0,1,...,n−1`. Reverse only
the arc `0→1`. Its new order is `1,0,2,...,n−1`. The entire named graph
changes, but remains isomorphic to W_n; therefore H, the whole unmarked
path-cover polynomial, and every uniform exterior response are unchanged.
In particular its exact count remains `h=2^n+1`.

The source can be supplied with a literal two-port history. Start with the
singleton X; append Y as a sink, then add n successive children at Y in
the order `n−1,...,0`. Each child is beaten only by Y. Retain Y and the
latest child0. The preparation fixes both names, so the first subsequent
birth can choose the available port0. No forgotten vertex is silently
reintroduced as a port.

This old-arc change has the inherited native implementation with one
temporary vertex τ beating only0. Its first Hamilton word is
`(0,1,...,n−1,X,Y,τ)`. Apply native complement off that spine, then use
`(0,τ,1,...,n−1,X,Y)` for a second native move and delete τ. Only the old
arc01 has changed. The script checks both legal words on their successive
actual inputs. The cost is one temporary birth, two native moves, one
deletion, no permanent vertices, and peak order n+3. Count-neutrality is
claimed at the macro endpoints, not at the temporary stages.

Both endpoints are strong. The root's indegree changes **1→2**; it does
not retain its former unique incoming port. The readiness defect `q−s`
below is distinct from the birth-return defect delta.

## 2. Exact response at every physical module position

Use `(h,s,q,F)`, where s counts Hamilton paths starting at the root,
q=H(W_n−root), and F counts ordered spanning covers `(P,Q)` with P
nonempty ending at the root and Q possibly empty. At zero-based module
position j, for `0≤j<n`, the exact rooted profile is

    h = 2^n+1,
    s = 2^(n−1−j),
    q = 2^(n−1)+1,
    F = 2^(n−j)3^j + 2^(j+1) − 1_(j=n−1).          (1)

This is precisely NATIVE_SELECTION (5.3), specialized to the actual
balanced singleton `(b,t)=(1,1)` before its N birth and n-child block.
That formula does not require a strong initial source; strongness is
needed only for the separate obstruction in that source's following
paragraph. The construction in §1 identifies its resulting graph with
W_n, so (1) carries the actual labels and port. The endpoint correction
`−1_(j=n−1)` is retained, rather than extrapolating an interior formula.

Write L=`2^(n−2)`. At the original first position and the prepared second
position, respectively, (1) gives

    z0 = (4L+1, 2L, 2L+1, 4L+2),
    z1 = (4L+1,  L, 2L+1, 6L+4).                   (2)

For n≥3 these positions are interior to the formula's boundary correction.
The birth-return defect `delta=F−h−q+s` changes from0 to `L+2`, while
the readiness defect `q−s` changes from1 to `L+1`. Physical root0 and
the exact full count are preserved; neither defect is asserted neutral.

At n=4, the one-arc operation is

    H17: (17,8,9,18) → (17,4,9,28),    delta:0→6.

This differs from the earlier two-phase move to `(17,2,9,44)`. It chooses
a different actual root position and has its own accepted return family.

## 3. Exact return and a factor-free infinite clock

A birth at root r adds x with `r→x→all other old vertices`. R retains r;
N marks x. The inherited actual state maps are

    R(h,s,q,F) = (h+F−s,q,q,2F−q),
    N(h,s,q,F) = (h+F−s,h−s,h,h+F).

For k≥1, apply `R^k N R` to z1. The final balanced profile is `(b,t,t,b)`
with

    t = (4L+3)2^k + L−1,
    b = (12L+9)2^k + 5L+1,
    b−3t = 2(L+2).                                  (3)

Indeed a=h−s=3L+1, d=F−q=4L+3, and J3 gives
`t=a+q+(2^k−1)d`, `b=a+3q+(3·2^k−1)d`. Every output is an actual
strong tournament of order n+k+4. The root stays physically0 during
preparation and the k R moves; N then deliberately changes it to the
newborn. The last R retains that newborn and seals the final physical pair.

**Accepted-clock theorem.** For n=3, every k≥1 makes `(b,t)` coprime.
For every n≥4, every positive multiple k of n−3 makes `(b,t)` coprime.
No factorization or search for a multiplicative order is needed.

Proof: since t is odd, J3 gives

    gcd(b,t) = gcd(L+2, 2^k(4L+3)+(2L+1)).           (4)

At n=3, L+2=4 and the second argument is odd. For n≥4 set m=n−3,
so `L=2^(m+1)` and `L+2=2(2^m+1)`. Modulo the odd number `2^m+1`,
the second argument of (4) is `−5·2^k−3`. If k is a multiple of m,
then `2^k≡±1`, so that argument is either−8 or2. Either is coprime to
the odd modulus. The factor2 in (4) is harmless because its other argument
is odd. This proves all asserted clocks.

The smallest declared clock k=n−3 for n≥4 has a compact exact form. Put
x=`2^(n−3)`; then

    b = 24x²+19x+1,
    t =  8x²+ 5x−1,    gcd(b,t)=1,                  (5)

on 2n+1 vertices. These are exact totals, not merely congruence classes.
For n=4,5,6,7 the totals are135,461,1689,6449. At n=3, k=1 instead
gives `(b,t)=(77,23)` on eight vertices.

The clock condition is sufficient, not necessary. It cannot simply be
dropped: at n=7,k=6, the same literal construction returns
`(b,t)=(25313,8415)` with gcd17. The formula, not a large graph census,
certifies this hostile. The boundary k=0 is also excluded: `N R` returns
`(16L+9,4L+1,4L+1,16L+9)`, not (3) evaluated at zero.

## 4. Existing-pair content, future sealing, and repeatable entry

The existing final root-child pair is primitive for **every n≥3 and k≥1**,
even for the off-clock hostile. Immediately before the final R, its parent
profile has coordinate gcd

    gcd(t,q,2^k d),    q=2L+1, d=4L+3=2q+1.

This gcd is1 because q is odd and coprime to d. The exact root-child
divisor theorem [ROOTED_RESETS R3](../../two_port_return_induction_20261103/ROOTED_RESETS.md)
therefore gives content1 for the pair just created. On the certified
clocks, (4) separately gives a primitive returned root, so **every fresh
root-child seal at that root** is primitive as well. At the hostile, the
existing pair still has content1, while a fresh seal has content17.

There is also an available primitive root at the **same final exact H for
every clock**, without the restriction on k: choose the other retained
endpoint, the final child. Its actual rooted profile is

    (b, t+L+2, 2t−L+1, 4t+3).                       (6)

To see its gcd is1, subtract twice the third coordinate from the fourth;
the result is q=`2L+1`. The second coordinate is `2^k d+q`, so a common
divisor must also divide `2^k d`. Since q is odd and d=2q+1, that common
divisor is1. Formula (6) follows by applying N, instead of R, to the same
pre-final-birth parent: it describes a different available port of the
same graph, with no extra vertices or changed full count. Choosing it for
the next birth does not retrieve a forgotten vertex. Its defect is zero,
however, so its primitiveness alone does not protect a subsequent singular
`R^ell N R`. The accepted clock ensures the stronger **balanced-root**
interface needed by the next paragraph.

Every certified clock enters the known repeatable primitive controller
after one C3 cycle. In the notation of
[NATIVE_SELECTION §§6–7](../../source_selection_closure_20261103/NATIVE_SELECTION.md),

    C3(b,t)=(61b−2t,9b).

For a coprime odd pair, a common prime of these coordinates could only be3.
But from (3)

    b−2t = (4L+3)2^k +3L+3 ≡ L·2^k ≠0 (mod3).

Thus C3 preserves primitivity and reaches `3|t'`, `3∤b'`. It can then
repeat indefinitely. To enter the inherited eight-state C3/C5 controller,
apply C5 if `r=t'/b'` in P1(F11) is not7. If r=7, one C3 sends it to3,
then C5 is permitted and sends it to0. These cycles have the already
proved native selectors and literal marked reentry, so the entry is an
actual source program rather than a scalar orbit assertion.

For example H17 prepares `(17,4,9,28)`, then `R N R` gives
`(135,41,41,135)`. One C3 cycle gives `(8153,1215,1215,8153)`, in the
primitive three-divisibility guard. The first three intermediate counts
after the H17 preparation are41,79,135.

Without preparation, z0 has a=d=q=2L+1 and zero defect. Every k≥1 then
returns `t=q(2^k+1), b=3t`, and the existing final pair has content q>1.
This is the repaired capability. Both initial root profiles already have
gcd1; the operation does not create initial primitivity from nothing.

The original homogeneous-module swap is asserted neutral only before the
selective births. After a child singles out root0, the old module is no
longer uniformly seen. Repeatability comes from the accepted balanced
controller, not from silently reapplying the initial neutral swap.

## 5. Reproduction and remaining exact-target gate

Run from the repository root:

    python3 -B 05-knowledge/results/h17_exact_repairs_20261112/constructors/wrapper_repair.py

`wrapper_repair_receipt.json` records all25 module-position profiles for
n=3,...,7; all five actual temporary cuts and native word pairs; eight
literal returned graphs through order15; their pre-seal profiles; and the
off-clock hostile with its precise formula-only scope. Endpoint subset DP
checks the literal rooted states. Independent permutation/cut enumeration
checks all positions at n=3,4. A small arithmetic control checks the proved
clock for n≤40 and three multiples each; it is not an unbounded computation
or a source inventory. Ordinary and optimized Python both retain explicit
checks via `require`, rather than removable assertions.

The source/target map is: actual cyclic wrapper with two retained ports →
one physical neutral arc change at the same exact H → chosen root defect →
explicit accepted return length → primitive balanced interface with a
repeatable marked continuation. It destroys the old root's indegree-one
property and the initial module's uniformity after selective attachment;
it preserves physical labels at the neutral step and stores the actual
native cut and birth sequence.

The remaining gate is target-directed: the exact count of the prepared
source is restricted to `2^n+1`, and the declared return totals are (3).
Nothing here supplies every prescribed odd total, preserves an arbitrarily
preselected next target, or proves density of the controller's image.
