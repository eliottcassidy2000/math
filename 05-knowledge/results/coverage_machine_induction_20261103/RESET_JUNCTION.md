# Cyclic gluing and a native reset before singular sealing

**PROVED scoped construction and divisor laws. FINITE-EXACT named controls.
Full prime and 7p coverage remains OPEN.** The new contribution is the
composable junction and its exact guard; elementary five-vertex tournaments,
the inherited birth maps, and classical automorphism divisibility are not
claimed as new discoveries.

## J1 The two cyclic interfaces

For a nonempty tournament T let h=H(T), and let f=F_2(T) count ordered pairs
of spanning disjoint paths, allowing empty slots. Write G(T)=T→a→b→T.
The inherited [S159 cyclic profile operator](../session159_constructions_20261103/profile_operator/NOTE.md)
gives H(G(T))=h+f. At the new root b its full rooted profile is

    (h_G,s_b,q_b,F_b) = (h+f,h,h,h+f).                 (J1)

Here rooted F counts two-slot spanning covers with the first nonempty path
ending at the specified root, including an empty second path. It is not f.
The two actual named old ports can instead be retained by the exact
[polynomial gluing transfer](gluing/NOTE.md).

To prove (J1), a path starting at b must proceed through a Hamilton path of
T and end at a, so s_b=h. Deleting b leaves T→a, giving q_b=h. In a cover
whose first path ends at b, either that path is the singleton b (h choices
for the remaining T→a), or a immediately precedes b. In the latter case
the two T pieces are an arbitrary ordered two-slot cover, including empty
pieces, giving f. Thus F_b=h+f.

Consequently the new-root interface is primitive exactly when gcd(h,f)=1.
Ordinal padding multiplies f by a power of two and leaves h unchanged, so
it cannot repair this gcd. By contrast, choosing an old root changes the
pointed response and must be evaluated separately.

The counter-machine seed [W(0,0)](axiom/COUNTER_CARRIER.md)=G(C3) illustrates
the difference. Its new-root profile is (15,3,3,15), whereas each original
triangle root has profile (15,3,5,19), of gcd one. Merely reassigning the
root to a forgotten vertex is not an allowed two-port operation.

## J2 A literal native junction on five vertices

Start with C3 on 0→1→2→0. Add a=3 as a sink, then b=4 losing only to a.
This is an actual cyclic wrapper; retain the physical ordered ports (a,b).
Rows list outgoing-neighbor bitsets:

    T = (10,12,9,16,7),       H(T)=15.
    P = (4,0,1,2,3).

P is a directed Hamilton path. Complement every arc except its four spine
arcs, obtaining

    U = (6,20,24,3,9),       H(U)=15.                 (J2)

U is the regular cyclic tournament on five vertices. Both sources are
strong, and the native move preserves all physical names. The exact
seventeen-query continuation contents at the *same* pair (3,4) are

    D(T;3,4)=3,             D(U;3,4)=1.

Every root of U has profile (15,3,5,21). Thus the native move supplies a
primitive pointed interface without requiring access to a forgotten old
triangle vertex. This differs from the inherited six-vertex H15 reset in
[NATIVE_RESET_INDUCTION N5](../two_port_return_induction_20261103/NATIVE_RESET_INDUCTION.md).
The common observer applies universally; no C3-adaptive ancestry of U is
asserted. Its entry certificate is exactly the gluing and native program.

This reset improves the *named* continuation divisor while worsening the
uniform-gluing divisor from [G3](gluing/NOTE.md). The complete unmarked slot
polynomials are

    F_T(r)=r^5+8r^3+6r,       d_slot(T)=3,
    F_U(r)=r^5+10r^3+4r,      d_slot(U)=15.

For five vertices the odd-cycle expansion gives the displayed form:
there are four versus five directed triangles, H is15 in both, and no
two disjoint nontrivial odd cycles fit. Their values at r=1,...,5 give
the two stated gcds by binomial inversion. The receipt independently
enumerates the slot covers. Uniform wrappers after U therefore retain15;
the pointed singular steps below are essential. Primitivity must always
name the allowed continuation operations.

## J3 The defect controls a singular sealing clock

Use the inherited [rooted birth maps](../two_port_return_induction_20261103/ROOTED_RESETS.md):
R adds a child at the root and retains that root; N adds the same child but
switches the root to the newborn. For z=(h,s,q,F),

    Rz=(h+F-s,q,q,2F-q),
    Nz=(h+F-s,h-s,h,h+F).

Both linear maps have rank three. Their discarded information is different.
Put

    a=h-s, d=F-q, delta=F-h-q+s=d-a.

Then delta(Rz)=delta(z), while delta(Nz)=0. The actual-source cone has
delta≥0. Each R step also forces s=q. It is therefore material when the
singular N choice is made, even if all current scalar counts are known.

**Exact sealing theorem.** For any actual rooted source and integer k≥1,
the literal word R^k N R ends at the balanced rooted profile (b,t,t,b), with

    t=a+q+(2^k-1)d,
    b=a+3q+(3*2^k-1)d,
    b-3t=2 delta,
    gcd(b,t)=gcd(delta,2^k d+q).                     (J3)

In particular, any positive power-of-two defect makes every such return
primitive. If the source is strong and delta=0, every return has b=3t and
t>1, so this particular sealing schedule is never primitive.

Proof: after k R steps the profile is (t,q,q,2^k d+q). Apply N and R,
using the general identity R(N(z))=(h+2F,h,h,h+2F). Operations in a
program word are read left to right. Subtraction gives
b-3t=2delta. Since t is odd, gcd(b,t)=gcd(delta,t), and t=2^k d+q-delta.
For a strong source, a=h-s>0; if delta=0 then d=a and t=2^k a+q>1.
No strongness assumption is needed for the formula itself.

For delta>0 the exact primitive-clock decision is finite: remove the power
of two from delta, call the remaining odd integer g, and follow 2^k modulo
g. The displayed gcd is periodic, with period dividing ord_g(2), using
period one when g=1. A prime p|g blocks a clock precisely when
2^k d+q=0 mod p. This retains nonunit d and q; no illegal division is made.

## J4 An infinite composable family in both target lanes

Apply (J3) to U at physical root b=4. Its defect is four. For every k≥1,

    (b_k,t_k)=(48*2^k+11,16*2^k+1),
    b_k-3t_k=8,       gcd(b_k,t_k)=1.                (J4)

The whole program is C3, two cyclic-gluing births, one native move, then
k+2 adaptive births. It uses k+7 vertices, no temporary vertex, and every
completed wrapper/native/return source is strong. The intermediate sink
addition is weak. At each adaptive step the root and newest child occupy
the two slots; R replaces the other slot, and N makes the newborn the root.
Thus the final pair is the last root and its last child with an explicit
joint entry program, not just an asserted profile.

At k=1 this gives H=107 on eight vertices. At k=2 it gives
H=203=7*29 on nine vertices. Both are independently recounted by permutation
enumeration; 107 and 29 are prime by bounded trial division. The k=3 control
gives H=395 on ten vertices. No infinitude of prime values is asserted.

These sources can enter the previously proved
[native balanced renewal family](../source_selection_closure_20261103/NATIVE_SELECTION.md).
For example its C3 return sends (b,t) to (61b-2t,9b). On (J4), b=2 mod3,
while t is 0 or2 mod3, so the new b is nonzero mod3. Primitivity is preserved:
any common odd prime in 61b-2t and 9b would divide 3, since gcd(b,t)=1;
the stated mod3 check excludes it. This gives the controller's arithmetic
guard, but not yet membership in its eight projective states. For example
k=2 gives (203,65)→(12253,1827), whose ratio t/b modulo11 is10, outside
that state set.

There is a short exact entry bridge. After the first C3 return, apply C5
if t/b is not7 in the projective line over F11. If it is7, first apply
one more C3, which changes the ratio to3, and then C5. The inherited C5
guard is now satisfied and its output ratio is0. This enters the old
eight-state controller while preserving gcd one, 3|t and 3 not dividing b,
with all names and entry programs. The ratio condition includes b=0 mod11
as infinity; then the C5 guard is automatically satisfied. The receipt
checks all ten k-classes governing the starting mod3/mod11 residues.

This is a closed repeatable construction, not complete target supply. Even
allowing every k here and all later block lengths of that balanced renewal
architecture still gives O(X^(2/3)) declared outputs below X: apply the
inherited bound 30(X/b_k)^(2/3) and sum over k, since b_k grows geometrically.
Fresh cyclic wrappers or other native moves at later stages are outside
that bound and remain legitimate avenues for enlargement.

## J5 Symmetry explains the obstruction but does not determine the reset

Let A be a subgroup of Aut(T) fixing every active port. Births uniform on
the forgotten vertices, with either polarity and any port mask, extend A
by fixing the newborn. Cyclic wrappers likewise extend A by fixing both
new vertices. Forgetting a port imposes no new restriction. As long as the
old moved vertices persist, A still acts faithfully. It acts freely on
Hamilton paths, so |A| divides every future count in this language.
The free action is classical and inherited from
[THM-4747](../../../01-canon/theorems/THM-4747-hamiltonian-cycles-under-automorphisms-fermat-wilson-congruences-for-circulant-and-paley-tournaments.md)
and its cited LEM-003; the extension assertion follows directly from arcs.

**Native intersection lemma.** For any actual Hamilton path P,

    Aut(T) intersect Aut(K_P(T)) = {identity}.       (J5)

This is the inherited simultaneous-stabilizer assertion in
[METAGRAPH section 3](../adaptive_metagraph_reformulation_20261103/METAGRAPH.md),
restated here to identify its role in the junction.
A common automorphism preserves the arcs on which the two tournaments
agree, exactly the directed Hamilton spine P. Its initial vertex and every
successor are fixed, proving the assertion. This kills the *old* common
symmetry; it does not say that the target is rigid. In (J2), automorphism
orders change from three to five, yet the named pair has trivial stabilizer
in the target. The numerical H remains fifteen.

A second control shows why native source changes and gluing need joint
data. The H31 source B=(42,4,1,22,7,94,31) and its neutral native partner
along (4,1,2,0,5,6,3) have F2=644 and596. Their cyclic wrappers therefore
have counts675 and627. At the inherited order-six H31 native pair, both
F2 values are328, so that simpler unmarked gluing does not distinguish them.
Neither observation classifies arbitrary native components.

## J6 Verification and remaining induction obligation

Run `python3 -I -B junction.py` and `python3 -I -B -O junction.py`.
The frozen JSON records actual graphs, paths, named programs and observer
values. Small controls use all labelled tournaments of orders two and
three, the stated H5 source, and T,U; test every root with k=1,2,3; and
check (J5) for every tournament/native-spine pair through order four.
The complete source universe and counts are in the receipt. No larger
integer-spectrum census is used.

The missing global step is a target-directed choice of source, native
path, and singular clock that realizes the requested N and leaves a legal
primitive continuation. The defect register identifies exactly one part
of that obligation. Scalar coverage assumed as an axiom supplies witness
existence, but cannot silently discharge this specified local selector.
