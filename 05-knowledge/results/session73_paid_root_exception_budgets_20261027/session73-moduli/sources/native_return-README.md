# Native return maps, rooted necklaces and exact reset memory

2026-10-26. **PROVED elementary results N1–N6 below; CITED classical inputs;
FINITE-EXACT controls; OPEN applications.** Local labels belong to this file.
Internal proof review and separate computational formulations, not an
independent agent audit or a claim of historical novelty.

The main construction realizes all four binary return maps by actual marked
tournament operations. A reset appends one output vertex sharing the old
anchor; it preserves the named source edge. The next result identifies an
exact uniform-bound criterion: for cyclic affine controllers containing a
reset, the norm of the real response operator is controlled precisely by
the largest spacing between resets. This does not establish such spacings,
or a paid substitute for them, for the physical Collatz source.

## Inheritance, incoming work and the live board

Anchor: native composition of `I,F,C0,C1` with shared ports. Niche: LCA
grammars and legal path changes. Wildcard: normal-basis field coordinates,
reset-memory spectra and the field-with-one-element interpretation.

Closest proved mechanism: [marked composition M1–M5](../tournament_composition_classification_20261026/metagraph/RESULT.md)
defines the ordered join `J`, its complement conjugate `star`, and their
different relocation laws. The canonical hostile is its pair of order-four
blue half-edges, masks 2 and 5, whose singleton extensions have Hamilton-path
counts 11 and 13. The corrected near miss is the purely set-theoretic
[four-state dictionary](../orientation_monodromy_rule30_20261026/README.md):
it did not supply native composition. The least-used sidecar is a named
source edge retained after the active output port moves.

Recovered stronger prior results, not new claims here:

- **THM-4047**, [Rule 30 left-front affine monodromy](../../../01-canon/theorems/THM-4047-rule30-left-front-affine-monodromy-clock.md),
  already proves the physical fixed-column reset/flip/phase theorem. Its
  moving-center uniformity boundary remains.
- [Prime-circulant necklace bijection, §§4–5](../circulant_necklace_bijection_20261020/CIRCULANT_NECKLACE_BIJECTION.md)
  already proves the cyclic derivative's zero/one/two-predecessor law and
  the bounded Fitting cut rule. Our trace interpretation recovers that
  mechanism; N5 adds variable coordinate selectors and exact response norms.
- **THM-4515**, [fair cycle-necklace splits](../../../01-canon/theorems/THM-4515-fair-consecutive-splits-of-cycle-necklaces-are-circulant.md),
  retains the affine carry and integer clock in a different circulant
  equation. Equality of cyclic indexing does not identify its arithmetic
  integrality test with our binary trace test.
- Incoming commits `d761e4cf0f` and `c9a8823436` were read: the former gives
  the sharp `4p` stabilizer boundary; the latter's
  [carousel classification, C–G](../cyclictriple_composition_20261026/README.md)
  gives complete decorated necklace fibres, physical flip multiplicities,
  spectra, and marked-gap growth. Its two growth children of `011 ~ 101`
  are the analogue of our cut-dependent returns, not the same move.
- The later pull through `8970493290` supplied **THM-4804**,
  [module shattering and normal-basis dictionary](../../../01-canon/theorems/THM-4804-complement-lines-shatter-every-module-and-the-carousel-fibre-as-frobenius-orbits-of-a-normal-basis.md).
  Its actual module-set proof and Frobenius dictionary were checked and
  integrated; N4's basic dictionary is now inherited concurrent work.
  Its source corollary overreached from labelled sets to unlabelled tree
  shapes. The exact order-four witness below repairs that clause without
  changing the module-set theorem. The incoming pre-audit checkpoint is
  not imported as an independently audited package. The new
  [norm-seven taxonomy integration](../norm7_taxonomy_zero_20261026/INTEGRATION.md)
  likewise distinguishes the normal-basis map from whole-surface F1 claims.

| Live concept | Exact map or operation | Preserved predicate; required sidecar |
|---|---|---|
| Source and active port | Append a vertex while retaining anchor 0 | Binary controller readout; original edge and flip parity recover the source |
| LCA composition | Ordered binary tree to free arc bits | A marked tournament; path and leaf order cannot be discarded |
| Necklace / normal basis | Binary word to a finite-field element | Rotation/Frobenius orbit; cut and primitive period retain physical phase |
| Reset memory | Cyclic difference minus a coordinate selector | Exact response and norm; reset spacings control uniformity |
| Collatz address | Exact valuation cylinder to an odd affine address map | Source, endpoint and integer clock; no resets appear on a fixed address cylinder |
| F1 skeleton | Pointed basis set and partial injections | Permutations/selectors; addition and field multiplication require extra structure |

The computation changed the board twice: the missing LCA presentations are
exactly the incoming half-edge hostile, and invertibility of every reset
controller still permits linearly growing response norms. Thus both the
path sidecar and the scale sidecar are essential. No new global meta-pattern
is promoted on the strength of this single packet.

## N1. A native compiler retaining its source

Let `X_n` consist of tournaments on named positions `0,...,n-1` with the
Hamilton path `0->1->...->n-1`. Write a free bit `t_ij=1` for `i->j`,
`i+2<=j`. The actual complement-line involution `tau` flips every free bit
and preserves the path edges. It is not full tournament converse.

For consecutive blocks, `J(T,U)` directs all cross arcs forward. The
inherited operation

    T star U = tau J(tau T,tau U)

preserves both interiors, directs all cross arcs backward, and retains the
one forward path hinge. Both operations are associative.

Start with any `T in X_n`, `n>=3`, named anchor `u=0`, cursor `v=n-1`, and
readout `b=[u->v]`. Let `e` be a fresh singleton. Compile:

| Return map | Actual operation | New active ordered port | Primitive cost |
|---|---|---|---:|
| `I(b)=b` | Leave `T` unchanged | `(u,v)` | 0 |
| `F(b)=1+b` | `tau T` | `(u,v)` | 1 |
| `C1(b)=1` | `J(T,e)` | `(u,new)` | 1 |
| `C0(b)=0` | `T star e = tau J(tau T,e)` | `(u,new)` | 3 |

Here a primitive is one whole `tau` or one singleton ordered join, not one
arc flip. Clock cost in a different native move alphabet has not been
proved. The idle `I` still consumes one controller tick.

**N1 (PROVED).** For every finite controller word, this gives its exact
binary composition, with one shared anchor and a fresh cursor at each reset.
Starting from the two `X_3` sources distinguished by `(0,2)`, if the word has
`r` resets and `f` flips, the final tournament has `3+r` vertices and

    initial bit = final t_(0,2) + f mod 2.                 (1)

Its controller clock is the word length, and its primitive clock is

    #F + #C1 + 3 #C0.                                    (2)

*Proof.* A join gives the anchor-to-new edge value 1. In a star append the
anchor is not the old last vertex, so that cross edge is reversed and has
value 0. Each reset preserves every old internal arc. A flip toggles every
nonconsecutive named pair, including the source and current output pair.
Induction proves the four actions, (1), and the two resource counts. ∎

The source map is injective even when its readout is constant. Any weights
assigned to the two sources remain attached to distinct full states; only
projection to the output bit merges them. This preserves weights without
proving any new sign or positivity bound for them. The archive is not a
record of the entire word: idle steps and pairs of flips make that impossible
without the explicit word/clock sidecar.

The compiler's `X_3` outputs have an even simpler representation: for each
`j>=2`, all free edges ending at `j` have the same bit `z_j`. Thus there are
`2^(n-2)` such presentations at order `n`. Append writes one bit and `tau`
complements the whole finite list. LCA data are unnecessary in this smaller
carrier; they become useful for general block composition in N3.

Converse is retained as well. The inherited `sigma` is reversal of the path
positions together with full converse. It sends a right append to a left
prepend, and reverses the order of factors in both `J` and `star`. Read the
mirrored bit from the left cursor to the right anchor; it equals the original
bit. Source names travel with this reversal, rather than being discarded.

**Boundary.** N1 is a growing, marked-edge realization of four functions.
It does not make the four unmarked order-four isomorphism classes into this
semigroup, and it does not realize resets on an unchanged distinguishing
edge. The actual reset changes the cursor while retaining the source.

## N2. Why a fixed-port synchronous word cannot reset

Use literal vertex labels and retain a Hamilton path `P` as part of the
state. A legal relocation changes `P`, not the tournament's labelled arcs.
`tau_P` toggles exactly the pairs not adjacent in `P`.

**N2 (PROVED).** Apply the same named sequence of relocations and complements
to two states, starting with the same `P`, and assume each relocation is
legal for both. Their labelled arc-XOR difference is invariant. In
particular a named pair that distinguishes the inputs still distinguishes
the outputs.

*Proof.* The path lists remain the same in the two runs. A relocation
changes neither labelled arc relation; a complement applies the same XOR
mask to both. Their difference is therefore unchanged at every step. ∎

The result does not cover input-dependent choices of paths or different
output pairs. N1 takes precisely the latter permitted route. Erasing the
labels after relocation would erase the premise needed for this obstruction.

## N3. LCA grammar, separable permutations, and its first escape

Let `G` be generated from singletons using `J` and `star`. In any ordered
binary expression tree, an arc between nonconsecutive leaves is forward
exactly when their LCA is a `J` node. Every consecutive pair is forward
regardless of the LCA label. This follows by induction at the cut between
the two children; it is an exact pairwise observable with no ties.

Define recursively separable permutations to be those generated from a
singleton by direct sum and skew sum. No external forbidden-pattern theorem
is needed for the following map. For a permutation `pi`, put

    rho(pi): i->j iff j=i+1 or pi_i<pi_j, for i<j.        (3)

**N3 (PROVED).** `rho` takes direct sum to `J`, skew sum to `star`, and
value complementation to `tau`. Consequently `G_n` is exactly the image
of the recursively separable permutations of length `n`, and is closed
under the actual complement involution on the marked path carrier.

*Proof.* Within each factor, relative comparisons do not change. All cross
comparisons in a direct sum are forward, and all in a skew sum backward;
the single adjacent cut is forced forward by (3). Complementation reverses
all nonadjacent comparisons. Induction proves both image inclusions. ∎

This is a new grammar choice that repairs complement closure, not a claim
that the incoming carousel or arbitrary modular-composition families are
complement-closed. Expression trees are not unique, and (3) loses comparisons
on adjacent pairs. They must not be treated as canonical modular trees.

For `n=4` use free-bit order `(0,2),(0,3),(1,3)`, least significant first.
The only absent presentations are masks 2 and 5. Indeed mask 2 forces
`pi_2 < pi_0 < pi_3 < pi_1`, hence the unique permutation `2413`;
mask 5 forces its complement `3142`. Neither admits a first direct/skew cut.
The other six are obtained by the recursive operations. These are precisely
the two incoming strong blue-loop presentations, now identified as the
minimal obstruction to this particular grammar.

**Hostile (FINITE-EXACT, literal four-vertex witness).** Mask 1 is in `G_4`.
Moving position 1 to the end gives the legal path `(0,2,3,1)` and mask 2,
which is outside `G_4`. Thus complement closure does not give relocation
closure. Keeping both an expression and its current path is indispensable.
This witness also shows why no all-order metagraph coverage follows from N3.

**N3 interface corollary (PROVED).** Every expression-subtree interval
`M=[a,b]` has a uniform baseline relation to each outside vertex. The actual
relation differs from that baseline only possibly at the two pairs
`(a,a-1)` and `(b,b+1)` that exist. Hence its cross-interface error matrix
has at most two nonzero entries and rank at most two, over any field.

*Proof.* For any outside leaf `x`, `LCA(x,y)` is the same ancestor for every
leaf `y` of the subtree. Thus its J/star label gives one uniform baseline.
Only a forced consecutive path edge can override this rule. The only such
edges crossing a contiguous subtree's boundary are its two endpoint pairs. ∎

This is a precise rerouting around module shattering: retain an interval,
two endpoint ports, and the baseline rather than demand zero interface
error. Under `tau` the J/star labels switch and the hinge exceptions are
recomputed, so the same expression topology continues to describe the
presentation. These are expression intervals, not canonical modules. The
bound does not give closure under relocation; the escape above still applies.

**Correction to incoming THM-4804's scope.** Its module-shattering theorem
does not force a change of unlabelled tree shape. Masks 2 and 5 have unique
proper non-singleton modules `{0,3}` and `{1,2}` respectively, yet they are
isomorphic via the vertex map `(1,3,0,2)`. Thus the old module disappears
while the same unlabelled tournament class returns. This is the same
incoming blue-loop pair, checked by literal arcs and all subset modules.
The source note, theorem summary and mistakes ledger have been narrowed
to the valid fixed-labelled-module statement. A common tree shape is not
the sidecar that the proof needs.

The final incoming [Session67 synthesis](../session67_metagraph_enumeration_information_20261026/README.md),
commit `69cffefe07`, already provides complete inherited all-order structural
recurrences and ordinary-class recognition by relocation components. Its
small scalar/family grammar remains open. Our grammar is a smaller native
constructor family, not a replacement for those all-order algorithms. Its
future-extension hostile is also relevant: equal present core and full
current automorphism cycle index need not give the same child counting
mass. Thus the actual attachment alignment remains necessary even after
retaining much more information than an automorphism order.

The exact counts for orders 1 through 10 are

    1, 1, 2, 6, 20, 70, 254, 948, 3618, 14058.

These count path presentations in this restricted grammar, not A000568
isomorphism classes. Direct permutation images independently agree through
order 7. No conjectured all-tournament recurrence is inferred from them.

## N4. Rooted necklaces and their finite-field coordinates

Fix a normal basis `beta,beta^2,...,beta^(2^(m-1))` of `F_(2^m)` over `F_2`.
The standard normal-basis and finite-field facts used here are **CITED** from
[Huczynska–Neunhoffer, Finite Fields, §§10–12](https://www.math.rwth-aachen.de/homes/neunhoef/Teaching/ff2013/ff2013.pdf).
The following deductions also have direct coordinate proofs.

For a bit word `w`, set `alpha=sum_i w_i beta^(2^i)`. Frobenius rotates the
coordinates one place to the right. Its orbit length is the word's least
period `d`, and is also `[F_2(alpha):F_2]`. The rotation stabilizer has order
`m/d`. The trace of the normal generator is 1: a zero trace would be a
nontrivial linear dependence of its conjugates. Therefore

    Tr(alpha) = sum_i w_i mod 2.                         (4)

For the incoming carousel family specifically, its proven group formula
becomes `|Aut(T_w)|=3^m m/deg(alpha)` for odd `m>=3`. The restriction to that
family matters: a general prime quotient may have additional automorphisms,
as the incoming order-21 Paley-seven control demonstrates.

Now encode a four-letter word by two bit words `(m_i,c_i)` with update

    z_(i+1)=m_i z_i+c_i over F_2.                        (5)

Their two field elements are subject to **simultaneous** Frobenius: separate
necklaces would lose the alignment of multipliers and translations. A
marked cut retains which port is the physical input. If every `m_i=1`,
the return is translation by `Tr(gamma)`, where `gamma` encodes `c`.
For a word of primitive length `d`, repeating it `m/d` times multiplies
this primitive trace by `m/d mod 2`. Thus a double traversal can conceal a
flip; retaining the primitive clock prevents that error.

If a reset occurs, the return is constant, but its value depends on the cut:
the chronological words `C0,F` and `F,C0` are rotations of one necklace and
return respectively `C1` and `C0`. The unrooted word retains a cyclic
compatibility problem, not a named-input return value.

There are three distinct operations here. Frobenius is a rotation of a
fixed basis. The compiler's `tau` complements all its column bits and thus
is `alpha -> alpha+1` in a normal basis of that fixed length. Appending a
column changes the length. It is **not** a field embedding
`F_(2^m)->F_(2^(m+1))` for `m>1`: such an embedding requires `m | m+1`.
An ordered-basis vector-space injection survives, while compatible field
multiplication generally does not. Nor is rotation of compiler column bits
automatically an automorphism or a legal relocation of its tournament.

## N5. A reset gives invertibility; reset spacing gives uniform bounds

Let `S` send the coordinate vector `z` to `(z_(i+1))_i`, indices modulo `m`,
and let `D=diag(m_i)`. The cyclic equations (5) are

    L z = c,   L=S-D.                                   (6)

In normal-basis coordinates `S` is inverse Frobenius. `D` is coordinatewise
multiplication by the bit mask, **not multiplication in the field**. A
mask supported at one position has rank one; multiplication by its nonzero
field element has rank `m`. This rank hostile is decisive for `m>1`.

If all multipliers are one, (6) is the Artin–Schreier difference equation
`sigma^(-1) Z+Z=gamma`. It has two solutions exactly when `Tr(gamma)=0`,
and none otherwise. This is the prior cyclic-derivative mechanism in field
coordinates: the kernel consists of the two constant words and the image
is the trace-zero hyperplane. It also reproduces the Rule 30 phase/flip
criterion; it is not a stronger physical Rule 30 theorem.

If at least one multiplier is zero, let `ell_1,...,ell_r` be the positive
cyclic spacings between its reset positions, so their sum is `m`.

**N5 (PROVED).** Over `F_2`, (6) has exactly one solution for every `c`.
Moreover, over the reals with the same 0/1 matrices, `L` is invertible and,
up to separate row and column permutations, its inverse is the block sum
of the lower triangular all-ones matrices `H_(ell_j)`. Writing
`ell=max ell_j`, one has exactly

    ||L^(-1)||_infinity = ell,
    ||L^(-1)||_F^2 = sum_j ell_j(ell_j+1)/2,
    ||L^(-1)||_2 = 1 / (2 sin(pi/(4ell+2))).              (7)

In particular a family of reset controllers has uniformly bounded real
response operators if and only if its reset spacings are uniformly bounded.

*Proof.* At a reset, `z_(i+1)=c_i`. Until the next reset, successive states
are successive partial sums of the intervening `c` values. Ordering the
inputs and outputs by these reset segments gives exactly the asserted
blocks. This proves the binary uniqueness and the real block formula,
including signs, without treating XOR as real addition.

Equivalently, `N=S^(-1)D` is nilpotent, `N^ell=0`, and
`L^(-1)=(I+N+...+N^(ell-1))S^(-1)`. The row sums and entry squares give
the first two norms. For a block, `H_s^(-1)` has diagonal 1 and subdiagonal
−1. Its Gram matrix has diagonal `(2,...,2,1)` and off-diagonal −1. The
eigenvectors `sin(j theta)` satisfy the left condition `v_0=0` and right
condition `v_(s+1)=v_s`; hence

    theta_k=(2k-1)pi/(2s+1),  k=1,...,s,
    eigenvalue_k=4 sin^2(theta_k/2).

Taking reciprocal square roots, and then the largest block, proves (7). ∎

The two alphabets have the same integer support matrices but different
arithmetic. Equation (7) concerns the real lift of the controller; it is
not a norm on a finite field and is not yet a ROOT-energy estimate.

**Uniformity hostile.** One reset followed by `m-1` unit multipliers is
invertible at every length, but has infinity norm `m` and 2-norm asymptotic
to `2m/pi`. Thus unique completion, or complete finite coverage, does not
by itself supply a depth-uniform bound. This makes the remaining hypothesis
quantitative rather than merely existential.

### A controlled route around long gaps, with its exact price

**Clock-weighted corollary (PROVED).** On each reset segment label the inputs
and outputs by their actual increasing integer times `t_1<...<t_s`. For
`0<rho<1`, use the same weight `w_i=rho^(2t_i)` on both sides of the real
partial-sum equation. The weighted operator has matrix

    K_(j,i)=rho^(t_j-t_i) 1_(i<=j),
    ||K||_2 <= 1/(1-rho),                               (10)

uniformly in segment length and in the integer time increments. Indeed
`t_j-t_i>=j-i`, so every row and column sum is at most the geometric sum
`1/(1-rho)`; the Schur bound proves (10). Unit-spaced long segments make
the bound sharp in the limit (use long almost-constant test vectors).
Separate row and column permutations retain all named ports and times.

This changes the norm, not the source, the map, or the clock. It gives a
well-controlled comparison model even when reset gaps grow. Returning to
undiscounted physical weights is the explicit remaining invoice: comparing
weights across a time span `Delta` costs up to `rho^(-Delta)` in norm.
There is no depth-uniform comparison without another argument. Applying
(10) to ROOT therefore requires the actual source/terminal weights or a
proved comparison with them; choosing convenient weights is not a proof of
the requested physical energy bound.

### Golden and seventh-root controls from the same boundary problem

For a reset spacing of two, `H_2^T H_2` has polynomial `x^2-3x+1`, so
`||H_2||_2=phi`. For spacing three, the polynomial is

    x^3-6x^2+5x-1,

and `||H_3||_2=1/(2 sin(pi/14))=1+2 cos(2pi/7)`. Its square belongs to the
real cubic subfield of `Q(zeta_7)` and is a unit of norm 1. These are exact
spectral appearances of golden and seventh-root arithmetic in a common
reset family, not identities inferred from matching numerical counts.

For comparison, a periodic seven-cycle Laplacian has the eigenvalue
`2-zeta_7-zeta_7^(-1)`, whose norm in the real cubic field is 7: multiply
the six factors `(1-zeta_7^j)` in conjugate pairs and use `Phi_7(1)=7`.
The open reset block has determinant 1. Thus the same cyclotomic field can
carry either a unit or a norm-seven element depending on the boundary
condition. This is a concrete boundary comparison with the
[norm-seven arithmetic work](../normseven_arithmetic_f1_20261026/README.md);
no map to its K3 surfaces, twists or conductor ideals is asserted.

## N6. Collatz retains an invertible address and its integer clock

Let `T(x)=(3x+1)/2^v2(3x+1)` on positive odd integers. Fix a nonempty
valuation word `k_1,...,k_h`, each positive. Set `s_0=B_0=0` and

    s_j=s_(j-1)+k_j,   B_j=3 B_(j-1)+2^s_(j-1).

Put `s=s_h`, `B=B_h`. Choose the unique residue `c` in `[0,2^(s+1))` solving

    3^h c+B = 2^s mod 2^(s+1),

and put `d=(3^h c+B)/2^s`, an odd positive integer.

**N6 (PROVED).** The exact positive source cylinder for this valuation word is

    x=c+2^(s+1)q,  q>=0,
    T^h(x)=d+2*3^h q.                                  (8)

Consequently its normalized endpoint address is

    (T^h(x)-1)/2 = (d-1)/2+3^h q.                       (9)

This is a bijection modulo every `2^r`. Modulo two it is `I` or `F`, never
a reset. The exact standard Collatz clock (one `3x+1` operation and then
individual halvings) is `h+s`.

*Proof.* Successive affine substitution gives the numerator recurrence.
Oddness of the final quotient specifies precisely the displayed residue
modulo `2^(s+1)`. To verify all intermediate valuations, write the numerator
at any earlier step `j` by unwinding the later recurrence. Its difference
from an odd multiple of `3^j x+B_j` has first term `2^s_j` times an odd
power of 3; all later terms have larger 2-adic valuation. The final
congruence therefore forces `v2(3^j x+B_j)=s_j` for each `j<h`, as well as
for `j=h`. Dividing consecutively gives exactly the specified `k_j`.
Conversely those valuations imply the congruence. The source residue is
odd by reduction modulo two. Substitution gives (8)–(9); `3^h` is odd,
so it is invertible modulo every power of two. Clock counting is literal. ∎

This is a source-preserving local arithmetic bridge, not a new convergence
theorem. Its sidecar is `(c,d,s,h,k_1,...,k_h,q)`, not only the final bit.
Each affine address bit can be compiled by N1, but that readout does not
transport the entire integer trajectory or prove a descent inequality.
The arithmetic clock `h+s`, controller clock and tournament primitive
clock (2) remain distinct integers.

For Rule 30 the inherited fixed-left-front update is

    ell_r(t+1) = ell_(r-2)(t)+ell_(r-1)(t)
                 +(1+ell_(r-1)(t))ell_r(t) over F_2.

It supplies actual multiplier/translation bits for N1 and N5 at every
integer time. In contrast, a fixed Collatz valuation cylinder supplies
only unit multipliers on its normalized source address. Any proposed
reset mechanism for Collatz must therefore justify a changed readout,
an independently certified terminal, or a paid forgetting of source
information. Calling a quotient constant does not pay that obligation.

## F1 interpretation and precise next targets

Use the specific **CITED** convention of
[Szczesny, Representations of quivers over F1, §1](https://arxiv.org/pdf/1006.0912):
finite pointed sets with maps injective away from the inverse image of the
basepoint. In this model rotations are permutations of the nonbasepoint
basis, and coordinate selectors are partial identities. Linearizing gives
`S` and `D`. Forming their difference in (6) uses the addition supplied by
the coefficient field; the F1 skeleton alone does not supply it. This is
one explicit F1 model, not a claim that a one-element field has been built.

If the two bit states are two nonbasepoint elements, neither constant return
to a bit is such a partial injection. N1 instead injects the two sources
into two different archived full states and projects their readouts later.
If one declares bit 0 the basepoint, `C0` is legal but `F` fails to preserve
that point; the whole four-map action still does not fit. The missing memory
coordinate explains the obstruction.

The most useful remaining questions are now explicit:

1. **OPEN — native closure with ports.** Can one add a finite family of
   port-response types to the LCA grammar so that legal relocation admits
   a bounded-cost update? The order-four escape is the first required
   hostile; an unmarked class or decorated tree alone fails.
2. **OPEN — quantitative paid resets.** For the actual weighted transport,
   exhibit a source-preserving renewal decomposition with a bounded weighted
   analogue of (7), retaining the original integer clock. Collatz's fixed
   address cylinders cannot provide resets by themselves. A random or
   unweighted average gap bound would not establish this fixed-source claim.
3. **OPEN — compatible changes of field chart.** Keep the ordered basis,
   cut and transition tensor when the controller length changes. Nested
   finite fields require divisibility of degrees; arbitrary one-vertex
   growth must use vector-space charts or a larger explicit common field.

These questions preserve the actual remaining obligations. Universal
Collatz entry, the Rule 30 moving-center questions, all-order black-core
connectivity, LRC(14), and six-colour forcing remain open here.

## Reproduction and scope of verification

Run from the repository root:

    python3 05-knowledge/results/native_return_necklaces_20261026/verify.py --output /tmp/native-normal.json
    python3 -O 05-knowledge/results/native_return_necklaces_20261026/verify.py --output /tmp/native-optimized.json
    cmp /tmp/native-normal.json /tmp/native-optimized.json

The script uses the standard library and explicit exceptions, not removable
assertions. Its exact universe is all controller words through length 7
on both initial bits (21,845 words, 291,272 prefix states); all pairs of
marked presentations through order 5 for common legal relocation controls;
the `J/star` grammar through order 10, independently checked against
recursive permutation images through order 7; all 28,813 subtree interfaces
of 2,321 recursively separable permutations through order 7; every binary word in chosen
normal bases in degrees 2–9; all four-letter cyclic words of lengths 1–6;
all 1,013 non-all-one multiplier masks of lengths 1–9; 780 valuation words
of lengths 1–4 with exponents 1–5, with 16 actual integer sources each;
1,093 increasing-clock patterns of segment lengths 1–7 and time increments
1–3 for the exact `rho=1/2` Schur bounds; and 64 actual Rule 30 updates on
each of 17 fixed left-front columns.

The finite-field implementation checks irreducibility by trial polynomial
division and finds normal generators by binary rank. Cyclic solutions are
enumerated independently of affine monodromy. Reset inverses are computed
by literal impulse propagation and checked by exact integer matrix
multiplication, then compared with segment norms. The golden and cubic
Gram polynomials have independent determinant checks. The permutation and
tournament paths share no membership recursion. These are separate
computational formulations within one session, not an external proof audit.

Positive controls are all-unit identity/flip words, all-reset words, constant
decorations, and exact arithmetic cylinders. Hostiles are the relocation
escape, the cut-rotation constant swap, rank-one selector versus nonzero
field multiplication, and the single-reset family with unbounded memory.
The full proof gives all-order quantifiers; the finite receipts do not.
