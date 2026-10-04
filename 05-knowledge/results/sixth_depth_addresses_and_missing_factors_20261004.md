# Six, sixty-three, and the coordinates needed for recursion

2026-10-04. **PROVED** for the constructions and obstructions linked below;
**FINITE-EXACT** for the explicitly bounded audits. **OPEN:** universal
Collatz return coverage and LRC(14). This session develops the incoming
sixth-clock results rather than claiming their calculations again.

The strongest outcome is a lossless recursive representation with a visible
failure boundary. The 63-state finite-field clock splits into two ternary
digits and a seven-state phase. Cubic extension introduces an additional
19-state phase that its original generator cannot reach. A general formula
identifies every later omitted factor and supplies a separate register for
it. Alongside that arithmetic construction, an executable Collatz grammar
stores completed inverse routes without expanding their integers.

## 1. Inheritance and research portfolio

The clean worktree integrated incoming commit `676c0c1fe` before deriving
new statements. Its [sixth-clock note](sixth_clock_branches_20261004.md)
already establishes the user's ordered Collatz exponents, the sextic
arising from lambda, the norm19 comparison, and the degree-six residue
clock. Its earlier sources already prove the divisor classification and
inverse-ray odometer. The present work preserves that lineage.

The closest proved mechanisms were the guarded inverse-ray compiler,
prime-power return orders, the tetrahedral pair map, and the sextic's
modulo-two field. The canonical hostiles were the root self-return,
the difference between field and multiplicative generators, and the
failure of a set bijection to preserve divisor order. The corrected near
miss was calling every source in a residue row a power-of-two predecessor.
The least-used sidecars were the second ternary digit, the prime ideal
branch, the divisor complement, and an unbounded block-height quotient.

| Portfolio | Object | What this session resolved |
|---|---|---|
| Anchor | Collatz inverse channels | An injective checked route grammar, digit lift and symbolic evaluator |
| Niche | Divisor layers and tetrahedral pairs | A symmetry obstruction and a replacement rank/charge geometry |
| Wildcard | The 63-state sextic residue clock | Exact norm loss, tensor coordinates, carry, and Fourier observation |
| Wildcard promoted | Cubic cyclotomic depth | The first missing19, all-height missing factors, and full phase decoder |
| Cross-check | Golden prime ideals | Pairwise coprime norm recurrence with two distinct ideal clocks |

The working board was revisited after each decisive pull: prime support
versus depth, inverse address versus height, norm versus branch, rank versus
label, field generator versus full clock, and completed route versus coverage.
Each change produced a concrete map or a hostile example rather than an
identification based only on matching numbers.

## 2. What the sixth term remembers

With `x_0=0` and `x_(n+1)=2*x_n+1`, one has `x_n=2^n-1`. Thus

\[
 x_6=63=3^2\cdot7,\qquad
 \operatorname{ord}_3(2)=2,\quad\operatorname{ord}_9(2)=6.
\]

There is no new prime at index6, but there is a new power of3. The radical
21 forgets that distinction. In the inherited field

\[
 K=\mathbb F_2[t]/(t^6+t^3+1),\qquad b=1+t,
\]

the 63 nonzero states are powers of b. Their two proper-subfield norms
`x -> (x^21,x^9)` have only21 possible values, with three states per fibre.
They recover the cube exactly:

\[
 x^3=x^{21}(x^9)^{-2}.
\]

They also preserve whether x is primitive, but they do not recover x.
An exact code consists of `k mod9=i+3j` and `k mod7=ell`. Multiplication
by b increments i, carries into j when i=2, and increments ell. After21
steps the norms repeat while the missing j digit has advanced; after63
the full state repeats. A multiplicative choice of one representative in
each norm fibre is impossible: the target has an order-three element all
of whose lifts have order nine.

The [complete field note](sextic_subfield_clock_20261004.md) also gives
`F4 tensor_F2 F8 = F64`, with a binary 2-by-3 coefficient matrix. Its
nonzero rank-one states are exactly the21 cubes; its42 rank-two states
include six elements of order9, so matrix rank alone does not imply a
63-step orbit. Six successive trace bits recover the entire field state;
all63 nonzero six-bit words occur as the cyclic windows of the trace clock.
This is a faithful Fourier/linear-algebra interface, with its observation
basis specified.

## 3. Repeating the construction exposes and repairs a missing phase

Take a primitive `3^a`-th root zeta in the field of dimension
`2*3^(a-1)` over F2 and set `b_a=1+zeta`. Write

\[
 q=2^{3^{a-1}},\quad c_a=b_a\zeta^{(3^a-1)/2}.
\]

Then c_a lies in F_q and

\[
 \operatorname{ord}(b_a)=3^a\operatorname{ord}(c_a),\qquad
 D_a=\frac{q+1}{3^a}
   =\prod_{j=2}^{a-1}\frac{\Phi_{2\cdot3^j}(2)}3
\]

is an unavoidable factor of its index in the full multiplicative group.
At dimensions2 and6 the factor is1. At dimension18 it is19:

\[
 |\mathbb F_{2^{18}}^*|=262143,
 \qquad \operatorname{ord}(1+\zeta)=13797=262143/19.
\]

The equality in that last order is **FINITE-EXACT**; the compulsory missing
factor and failure at every larger dimension are **PROVED**. The extra
factor73 in the dimension18 group is retained. Every prime other than3
in `Phi_(2*3^j)(2)` has base-two order exactly `2*3^j`. Thus the omitted
phases are specified by successive primitive-prime layers, not by a vague
recurrence of the numeral6.

There is a complete multiplicative repair:

\[
 \mathbb F_{q^2}^*
   \simeq\mathbb F_q^*\times\mu_{3^a}\times\mu_{D_a}.
\]

The [explicit decoder](cyclotomic_depth_towers_20261004.md) uses the
relative norm, its unique subfield square root, and two CRT projectors.
All three registers together recover every state. The b_a orbit has
third register1. Allowing only that register to vary still leaves the
first restricted to `<c_a>`; full coverage also requires the entire
subfield register. Primitivity of c_a is checked through a=6 but is
not claimed at every height. Both b_a and c_a remain compatible under
relative norms throughout the cubic tower.

The user's `63=3*19+6` fits this precise calculation because
`Phi_18(2)=57=3*19`. Its19 is exactly the first compulsory phase loss
above the63-state construction.

## 4. The golden ideal equality really does recur, with a different base

For `d=3^(k-1)`, let

\[
 A_k=\phi^{2d}+\phi^d+1,\quad
 B_k=\phi^{2d}-\phi^d+1,\quad N_k=L_d^2+3.
\]

The norms of both A_k and B_k equal N_k, and

\[
 N_1=4,\qquad N_{k+1}=N_k^3-3N_k^2+3.
\]

For l>k one has `N_l=3 modN_k`, while every N_k is1 modulo3. Hence
the norms are pairwise coprime. The first four are

\[
 4,\quad19,\quad5779,\quad
 192900153619=3079\cdot62650261.
\]

For every k>=2, the two principal ideals are coprime and

\[
 (A_k)(B_k)=(N_k),\qquad
 \mathbb Z[\phi]/(N_k)\simeq
       (\mathbb Z/N_k)\times(\mathbb Z/N_k).
\]

The two copies carry phi with exact orders `3^k` and `2*3^k`.
This extends the earlier norm19 ideal equality to every level k>=2,
including composite norms. Recording which ideal is being used matters:
the scalar norm alone does not distinguish these two clocks.

The agreement with the binary tower stops at the next test. At k=3,
the binary factor is87211 while the golden norm is5779; base2 has
order5778 modulo5779, while the golden branches have orders27 and54.
The preserved operation is cyclotomic evaluation with a named base,
not an equality of all resulting clocks. The proofs and exact controls
are in [the tower note](cyclotomic_depth_towers_20261004.md).

## 5. A practical Collatz carrier for the three rows

Use `H(n)=(3n+1)/2` and the fully accelerated odd map U. At the root1,
the selected sources in the owner's order are

| Source row mod6 | First nonroot source | H-image |
|---|---:|---|
|3|21|`2^(5+6j)`, j>=0|
|5|5|`2^(3+6j)`, j>=0|
|1|85|`2^(7+6j)`, j>=0|

The last row originally starts with1 and exponent1. Removing that
self-return explains the offset to exponent7. Each row contains many
other integers; the table parameterizes the direct power-of-two sources.

For a certified positive odd hub u not divisible by3 and a source row
r modulo3, there is a unique `kappa in1..6` with
`2^kappa*u=1+3r mod9`. Every inverse source in that row is

\[
 n_b=\frac{2^{\kappa+6b}u-1}{3},\qquad b\ge0.
\]

The stored node is just `(parent certificate,r,b)`. It has exact forward
image u, so its route is inherited. Row0 is a completed leaf, while rows1
and2 can be extended. The sole root self-return is rejected. A finite
ternary address is selected one digit at a time using

\[
 n_{b+d3^{a-1}}\equiv n_b+d3^a\pmod{3^{a+1}},
 \qquad d\in\{0,1,2\}.
\]

Binary and ternary residues can be evaluated recursively without forming
n_b. The tested21-step object represents an integer with more than
`10^102` binary digits, while retaining exact first-hit ranks and legal
extensions. This is an actual finitely described integer and a proved
route, not an expanded enumeration of its digits.

A concrete completion of the previous representation's miss7 is

\[
 \frac{22\cdot4^j-1}{3}\ \longrightarrow\ 11\longrightarrow17
 \longrightarrow13\longrightarrow5\longrightarrow1.
\]

Every j>=0 gives exactly five odd steps and `16+2j` ordinary steps.
This entire fibre is outside the twice-repeated favorable-chart atlas.
It is not claimed as new density beyond the older frozen bank: that bank
already has a first-descent cylinder containing7. The
[codec proof and implementation](inverse_ray_ternary_addresses_20261004.md)
make both distinctions explicit. The grammar is bijective onto the odd
basin of1; finding its code for every positive odd input is still the
global coverage problem.

## 6. The p-squared-q-r analogy has a sharp geometric boundary

For N=p^2*q*r the inherited proper-divisor counts are `F=10,S=7,U=3`.
These are counts, not disjoint literal sets: primes are included among
squarefree divisors. Doubling63 to126 supplies the missing third prime
and realizes this family. Keeping the three primes but extending the
p-height gives the actual rank layers

\[
 p^kqr:\qquad3,\underbrace{4,\ldots,4}_{k-1},3.
\]

The coordinates `(Omega(d),v_q(d),v_r(d))` retain the divisor and turn
complement into a central reflection. The two Boolean exponents are an
arithmetic four-state charge: multiplication adds their parities. Every
middle four-point layer has zero sum for all nontrivial Walsh characters;
the surviving modes come from the two boundary layers.

By contrast, the ten-point tetrahedral pair alphabet has a different
four-state charge. It retains its labelled set and the7+3 split, but cannot
turn divisor complement into an affine tetrahedral symmetry under *any*
bijection. Complement has no fixed divisors here; every tetrahedral
involution fixes at least two degree-two points. This rules out that
stronger identification and explains which coordinate geometry repairs it.
See [the divisor proof](divisor_layers_pair_geometry_20261004.md).

## 7. What is now usable, and what remains worth testing

The reusable design is to specify a quotient, compute its actual fibre,
and retain a coordinate that restores exactly the desired operation.
The fibre is a ternary carry for the two sextic norms, an independent
cyclotomic phase in the larger fields, an ideal branch for the golden
norm, and a block-height quotient plus certified parent for a Collatz
address. These are separate exact constructions with a shared design
principle; no one of their dynamics has been silently transferred to another.

The immediate next research targets are precise:

1. Determine whether c_a is primitive beyond the verified six levels,
   or find its first missing subfield factor. The all-height storage
   decomposition remains valid either way.
2. Use the route codec as the common backend for known return families,
   preserving shared parent pointers and comparing actual integer coverage.
   Ternary support alone is not an arrival theorem.
3. Keep the arithmetic exponent charge when composing divisor letters;
   use the geometric pair charge only for operations proved to respect it.

Four proof notes, four standard-library experiments and their outputs
accompany this synthesis. Independent cross-reviews checked the field
and divisor maps, the golden branch orders and the guarded route grammar.
The final reproduction record is
[sixth_depth_session_audit_20261004.out](sixth_depth_session_audit_20261004.out).
