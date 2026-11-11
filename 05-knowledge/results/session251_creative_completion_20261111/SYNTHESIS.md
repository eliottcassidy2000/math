# Five small cores and a rigid low-ratio attachment law

**Session251, 11 November 2026. COMPLETE research deliverable.** New
proved-in-note mechanisms and finite exact construction certificates.
No independent canon audit or literature-priority claim. **Chan–Pak
remains OPEN.** The session builds new constructions rather than repeating
verification of inherited large censuses.

## 1. The main constructive change: replace a growing parent by five fixed cores

The inherited sufficient route asks for one-vertex pairs `(h,3h+2)` and
`(h,3h+4)` for every realizable odd h≥95. This session opens a different
route: use one of just five explicit cores, surround it by transitive
flanks, and add one apex. The cores have Hamiltonian-path counts

    1, 3, 5, 9, 9.

Their outgoing row masks are `(0)`, `(2,4,1)`, `(2,4,9,3)`,
`(2,4,9,19,7)` and `(2,4,25,3,11)`. Flank lengths are unbounded, and every
core cut is allowed. This is an actual graph grammar, not an operation
that depends only on scalar H.

For a fixed core C and cut O, retain nine numbers obtained by evaluating
four residual-H-weighted directed-path polynomials at +1 and −1, together
with d=H(C). Let A encode the left flank cut, and B the complemented,
reversed right flank cut. With α=A mod2, β=B mod2, the exact child count is

    d(g(A)+g(B)-1)+(d+k_+)AB+(d+k_-)αβ
      +i_+A-i_-α+o_+B-o_-β+j_++j_-,

where the inherited positive one-apex function satisfies

    g(0)=g(1)=1,
    g(2q+1)=g(q),
    g(2q)=g(q)+q+(q mod2), q≥1.

The new formula works for every finite core, including uniform apex contact.
Its proof counts the unique odd cycle using the apex and retains the
Hamilton count of the unused core. It does not incorrectly collapse a
module that a Hamiltonian path can revisit. A positive lower bound gives
a terminating inverse for every target, with the uniform axes handled
separately. Every successful inverse outputs all arcs of a witness.

**Finite result:** the five-core bank realizes every allowed odd through
262,143. This is a compact-grammar result, not a new unrestricted coverage
record; the older, richer one-apex-join work already reaches much farther.
The sixth diagnostic source C3▷C3 adds no target in the tested range.

**New open sufficient statement:** the same five cores cover every positive
odd except 7 and 21 at all heights. Proving this would establish Chan–Pak
without proving the prescribed-parent ternary branches. It is a concrete
finite union of binary arithmetic images, with a literal graph constructor
for every image. Finite coverage does not prove it, and it is not claimed
equivalent to the original conjecture.

[Full formula, proof, inverse, masks and exact experiments](ternary/PROOF.md).

## 2. A sharp new small-parent obstruction, and why both H9 charts matter

The first three cores give the **complete all-order one-birth image of all
sources with H≤5**. This follows from inherited SCC multiplication and
the strong minimum: such sources are transitive, a padded C3, or a padded
unique strong four-vertex tournament.

Their first extra missing value is **173**. Consequently every vertex
deletion of every H173 tournament has at least nine Hamiltonian paths.
The bound is sharp. The child with rows

    (126,68,72,82,230,78,128,47)

has H=173 and deletion counts

    (77,71,51,27,51,71,27,9).

This advances the inherited H55 deletion threshold from five to a sharp
nine at a new target. It is a statement about all realizers of 173, not
just about the particular displayed graph. The negative part uses the
proved terminating inverse, not an arbitrary order cutoff.

Through 262,143 there are exactly 38 additional holes in the H1/H3/H5
parent language. The two H9 core types fill them, and neither type can
simply be discarded from this bank:

- Without `(2,4,25,3,11)`, the target 21,901 is missing from the remaining
  four-core grammar at every flank length.
- Without `(2,4,9,19,7)`, the target 86,641 is missing in the same sense.

Both H9 cores have the same uncoloured path totals `(k_+,k_-)=(75,-9)`.
Their coloured endpoint-path coordinates distinguish them. This is a
concrete reason to retain two charts instead of replacing them by an
abstract state labelled H9.

The full hole list, both repair searches and 74 directly recounted repair
graphs are retained. Simple hole rules h→2h±1 or 2h±3 fail in the tested
range. The exact binary coefficient pullbacks are proved, but their slope
grows; no finite carry automaton is asserted. The next useful step is a
proved normalization retaining the chart and target, rather than a larger
finite scan.

[Complete hole and repair certificates](ternary/hole-repairs.json).

## 3. A new rigidity theorem for the original ternary route

There is also structural progress on the inherited h→3h+2/4 route.
Write a source B as its ordered strongly connected components, and let
h=H(B). A singleton SCC is a neutral vertex: removing it does not change h.
List these singleton vertices in source order, putting bit1 when the new
apex beats the vertex and bit0 otherwise.

**Theorem:** if the child count is below 5h and this singleton word has a
forward 1-before-0 inversion, the entire child is

    L ▷ C3 ▷ R,

where the new apex and two adjacent singleton vertices form that C3.
In particular its count is exactly **3h**.

Thus both desired branches 3h+2 and 3h+4, for h≥3, force the singleton
word to be globally monotone:

    0^a 1^b.

If there are s singleton vertices and c vertices in nontrivial SCCs, the
attachment search shrinks from `2^(s+c)` to `(s+1)2^c` for these branches.
This leaves the core choices unbounded; it is not a proof of branch existence.

The mechanism is a positive capacity budget. The inherited odd-cycle
expansion gives

    H(B+apex)=h+2 sum_(i beaten, j not beaten) c_B(i,j),

where c_B(i,j) counts even-vertex directed paths from i to j weighted by
H of their complements. Every inverted singleton pair costs at least h.
For two singleton endpoints around an intervening tournament M, the new
exact formula is

    c(u,w)=H(L)H(R) sum_(S⊆M, |S| even) H(M[S])H(M−S)
           ≥ (floor(|M|/2)+1)h.

This first forces a unique adjacent inversion. Session250's path-cut
injection then sharpens it: any nonuniform apex contact before or after
that pair costs another h. Below5h the triangle must therefore be isolated
as an ordinal factor. This is a genuine synthesis of the earlier rooted
boundary inequality with a different, positive attachment carrier.

**One-core consequence:** if the parent is TT_a▷C▷TT_b with C nontrivial
strong, a ternary-branch witness can be stripped to a strong child having
padding on at most one side of C. Both active flanks would create an
inversion spanning the core. This is an exact simplification of that
source grammar, not a claim that every source has only one strong component.

[Proof, capacity identities and 33,866 labelled attachment controls](geometry/PADDING_CUT_NORMAL_FORM.md).

## 4. The niche bridge: safe inflation creates exact Mersenne arithmetic

An older intercyclic classification immediately rules out the tempting
strategy of obtaining every odd count from pairwise-intersecting cycles:
that subclass also misses H29, H35, H41, H43, and other values. Rather than
repeat that failed search, the session identifies exactly when its clean
arithmetic is legal.

Replace each vertex i of a template R by a transitive block of size m_i.
The inflated tournament has no two disjoint directed cycles **if and only
if every pair of template directed triangles, including a repeated pair,
shares a vertex whose block has size one**.

Under this condition, a cycle cannot reenter a block: two visits would
split into two disjoint cycles. Therefore

    H = 1 + sum_C ( product_(i∈C)(2^m_i−1) − (−1)^|C| ),

summing all simple directed template cycles, including even ones. Even
core cycles matter because their inflated lifts can have odd length.
The result is a positive-coefficient polynomial in Mersenne coordinates
apart from its constant term, with an exact admissibility test.

For all non-one-apex five-vertex templates, this yields the all-size
families `9·2^m−3` and `1+3·2^(a+b)`; the remaining regular core cannot be
safely inflated. The sequence values are not claimed newly discovered.
The new value is the common criterion and the no-reentry proof. The hostile
C3[TT2,TT2,TT2] has H45, whereas applying the formula illegally gives29.

The next productive departure is controlled packing number two, retaining
`H=1+2c1+4c2`, rather than pretending the intercyclic restriction is universal.
No construction controlling that second packing term is claimed here.

[Proof, all-size five-vertex analysis and 2,331 blowup checks](cycles/README.md).

## 5. A recovered modular engine explains what the arithmetic proof still needs

Session197 already supplies joint modular universality for the positive
one-apex compiler. That strong tool was recovered rather than re-proved.
For any physical rooted prefix `(h,e,R)`, its one-arc-join counts are

    F(L)=hG(L)+RL−e(L mod2),  G(L)=g(L−1).

The new complete local classification says that modulo even M its image
is exactly the two classes h and h+R−e modulo `gcd(2h,2R,M)`. Modulo odd m,
it is the two classes0 and−e modulo `gcd(h,R,m)`.

Thus gcd(h,R)=1 removes every modular obstruction. It does not control the
integer value: the singleton prefix has every permitted congruence class
but permanently misses exact13. A second proved observation is that a
constructor hitting all M²/2 joint residue pairs must sometimes use an
address at least M²/2−1. Demanding linear address height for all joint
pairs is therefore impossible.

The actionable change is to optimize the address **over the scalar target
fibre**, leaving the auxiliary cover residue free. It is not useful to
solve another residue census or prescribe both coordinates unnecessarily.
This supports the five-core carry problem while explicitly preserving the
missing height requirement. The new shortest-word experiment covers all
even moduli through128; no all-modulus short-word theorem is claimed.

[Local theorem, constructive residue witnesses and height boundary](gems/README.md).

## 6. What was achieved, and what remains

Proved in this packet:

1. An arbitrary fixed-core, two-flank attachment formula and complete target
   inverse, including uniform-contact cases.
2. The all-order H173 minimum-deletion threshold nine, with an actual witness.
3. Low-ratio singleton-inversion rigidity and monotone ternary padding.
4. Safe transitive-inflation criterion and its full directed-cycle polynomial.
5. Exact local prefix images and the joint-address counting barrier.

Finite exact outcomes include five-core coverage through262,143, the38-hole
small-parent complement, and the two nonredundant H9 repairs. Inherited
large-range coverage, the strong minimum, parity, odd-cycle expansion,
intercyclic omissions and joint modular universality remain credited inputs.

The preferred next attack is the **five-core covering conjecture**:
classify the complement of the first three explicit binary images and
prove that the two H9 charts cover it. The main alternative is the original
ternary induction with the newly proved monotone singleton language.
Neither universal existence statement is established in this session.

## 7. Reproduction and concurrency

Run `python3 -B run_all.py` in this directory. It checks the six new engines
against their saved outputs; `--optimized` repeats with Python's `-O`.
This is validation of new results, not a renewed audit of inherited claims.
Hashes cover the deliverable; navigation files provide a short current route.
No theorem ID was reserved or promoted without independent audit.

The incoming [odd-line seed packet](../odd_c_seed_symmetry_20261110/README.md)
at `2a692c035fa` supplies 18 additional unit-offset ray seeds, including
line245, and proves that offset one forces a trivial selected-set
stabilizer. These are inherited concurrent results, not new claims of this
session. They complement the strong deletion-core search and reinforce
why a coloured core, rather than its bare H value, is the right carrier.
They do not establish the five-core covering conjecture or all-line coverage.

The incoming Fortunate/Landau packet at `a06eb56edac` was preserved during
publication. Its finite-congruence warning is compatible with Section5,
but no primality theorem is transferred to tournament counts. Its objects
and unproved infinitude questions remain separate.
