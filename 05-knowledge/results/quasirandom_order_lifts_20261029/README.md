# Sylvester symmetry, Wieferich delay, and the correlations that make quasirandomness

2026-10-29. **PROVED (elementary, self-audited):** P1–P6 below.
**FINITE-EXACT:** the declared controls in `verify.py` / `exact.json`.
**CITED:** the external group and large-sieve comparisons in §7.
No priority claim for the general Fourier/probabilistic methods. No universal
Collatz, LRC, or A000568 enumeration consequence is asserted. This is a result
packet, not a new audited canon theorem or a promotion of an open hypothesis.

The central connection is an exact correlation identity, not the recurrence
of a numerical label. The same orbit constraints can permit quasirandom
tournaments and strongly biased tournaments. A prime-adic order delay tells
us precisely how long one useful correlation identity survives refinement.

**Concurrent credit.** During the push, commits966ab75d90 (Session107) and
6ce7838d11 arrived. The latter independently proves exactly the19-adic
period law P5, including the same return-defect calculation; that result
is shared concurrent work, not an exclusive contribution of this packet.
Both incoming proofs, their audits, the character-fibre theorem, the Euclid
syndromes and the group boundaries were read. Section10 integrates their
stronger comparisons. The orbit-ensemble theorem P1–P4 and deterministic
Sylvester-indexed construction P3b concern different carriers from their
tournament on the Sylvester terms themselves.

## 1. Inheritance and research board

**Anchor:** explain quasirandomness while keeping the actual arithmetic action.
**Niche:** lift the Sylvester four-cycle at 19 through every prime power.
**Wildcard:** recover the stronger divisor-height repair of the Monster chain.

The six live concepts are: odd orientation symbols; Fourier correlations;
orbit-pair choices; return defects and derivatives; source-residue clocks;
typed divisor-height dynamics.

Closest proved mechanisms:

* [LEM-004, odd-function dictionary](../../../01-canon/theorems/LEM-004-tournaments-are-odd-functions.md): a circulant tournament is an odd sign function. Its higher correlations need spectral information; oddness alone is insufficient.
* [THM-4644, Fermat/Wieferich circulants](../../../01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md): odd-order multiplier orbits, their sign choices, and the 6564/3512 choice-bit counts at the two prime squares. The present packet supplies correlation and completion bounds for those existing families.
* [The complete order tower and source compiler](../fermat_wieferich_order_20261014e.md), §2–4, already prove the all-depth order law and source-residue lifting. We reuse them, rather than claim them anew.
* [THM-4631, compressed twisted transfer](../../../01-canon/theorems/THM-4631-level-factorisation-of-the-compressed-twisted-transfer-operator-of-the-px-plus-1-maps-fast-modes-are-the-doubling-cycles-on-the-units-mod-p-i.md): the graded fast modes retain actual doubling cycles and their multiplicities. Its closure to individual integers remains open.
* Concurrent [tournament proof routes](../tournament_open_problem_routes_20261029/PROOFS.md), P1, use the quadratic large sieve to prove prime-Legendre quasirandomness, with numerical height retained. Concurrent [Session106](../session106_collatz_transfer_20261029/README.md) proves that the native killed Collatz operator has unit norm at every finite depth, even when testing an everywhere-terminating model. Neither result gives a uniform all-source contraction.

Canonical hostiles: the rotation circulant from LEM-004; 17 (non-Wieferich
but not an Artin prime for base 2); a locally admissible prime whose actual
Sylvester orbit misses both roots; and the ramified multiplier modulo 49.
Corrected near miss: the historical June30 reflection inferred attraction
from a shared fixed point. Its corrected real basin is in §8.
Least-used sidecars recovered here: the first nonzero return defect,
the derivative of a polynomial return, and composite-base divisor height.

## 2. P1: an exact isotropy theorem for constrained orientations

Let N>1 be odd, H a subgroup of (Z/NZ)^× of odd order r, and assume

    gcd(h−1,N)=1 for every h∈H other than 1.                 (G)

This is a condition on every prime divisor of N, not just on the order of H.
Reduction of H modulo each such prime is injective. Its image has odd order,
so it cannot contain −1. Therefore h+1 is also a unit for every h∈H.

The action of H on nonzero residues is free, and its orbits occur in distinct
pairs O,−O. Choose one representative of each pair, independently put
σ=ε_O on O and σ=−ε_O on −O, with ε_O uniformly ±1, and put σ(0)=0.
There are m=(N−1)/(2r) choice bits. Define

    K(x,y)=σ(y−x),       x→y iff K(x,y)=1.

This declares vertices, orientation gauge and pairwise observable explicitly.
There are no off-diagonal ties. Complement is global sign reversal; both
translations and H are automorphisms. We count labelled connection sets,
not tournament isomorphism classes or Hamilton path presentations.

**Claim.** Exactly, without a limit or a primality assumption on N,

    E(K Kᵀ) = N I − J.                                    (1)

**Proof by incidence.** For nonzero x,y,
E[σ(x)σ(y)] is 1 if y∈Hx, −1 if y∈−Hx, and zero otherwise.
For a nonzero displacement d, sum this expression at y=x+d. Each h≠1
contributes one positive solution of (h−1)x=d; each h contributes one
negative solution of (h+1)x=−d. All these solutions have x,y≠0; within
either sign there is no duplicate h because the action is free. The sum
is (r−1)−r=−1. At d=0 it is N−1. These are precisely the cyclic row
correlations in (1). ∎

Thus these ensembles have the same *mean* Gram matrix as a Paley tournament.
An individual member need not satisfy that identity. The proof works also
for any finite abelian group of odd order and an odd automorphism group H
for which h−1 is an automorphism whenever h≠1: the same two equations have
unique solutions. That is the actual common incidence structure.

## 3. P2: quasirandom completions, including prescribed choice bits

Circulant matrices are unitarily diagonalized by the additive characters.
Their eigenvalues here are purely imaginary. Write

    λ_k = i Σ_O ε_O c_O(k),     c_O(k) real,    k≠0.

Equation (1) says Σ_O c_O(k)^2=N for every nontrivial frequency, including
imprimitive characters. Independence and cosh(t)≤exp(t²/2) imply

    P(|λ_k|≥t) ≤ 2 exp(−t²/(2N)).

Taking a union bound over N−1 frequencies proves that at least half the
2^m connection sets satisfy

    ||K||₂ ≤ sqrt(2N log(4(N−1))).                         (2)

The displayed logarithm is natural. The zero eigenvalue is exactly zero.
This is an existence/count statement; it does not claim that a fixed
deterministic choice rule succeeds or that enumerating all choices is cheap.

For vertex subsets A,B,

    |Σ_(x∈A,y∈B) K(x,y)| ≤ ||K||₂ sqrt(|A||B|).            (3)

Hence along N→∞ the good completions have cut discrepancy o(N²).
They also have the random density 2^(−binom(v,2)) of every fixed labelled
tournament. To see this directly, use W=(J+K)/2 (diagonal 1/2), telescope
one edge at a time, and fix all other vertex variables. The factors depending
on the two exposed vertices separate into two bounded vectors. The normalized
bilinear error is at most ||K||₂/(2N). Sum over the binom(v,2) edges;
discarding repeated vertices changes the density by O(v²/N). Thus the
operator estimate itself proves the required notion of quasirandomness.
Unlabelled classes have weights proportional to v!/|Aut(F)|, not equal weights.

**Pinned version.** Fix any subset of the orbit-pair signs, and let K₀ be
their partial signed matrix, with zeros on the unassigned packets. For the
remaining independent choices, the frequency variance is at most N. Thus
at least half the completions obey

    ||K||₂ ≤ ||K₀||₂ + sqrt(2N log(4(N−1))).               (4)

In particular a pinning of spectral cost o(N) can be extended quasirandomly.
This is a sufficient condition, not a necessary one: the unfixed part can
sometimes cancel a large pinned mode. Pins are whole H-orbit pairs; an
arbitrary list of edge requests must first be checked for conflicts under
translation, H and converse. A monochromatic rectangle of positive relative
area is an immediate hostile. Ports cannot be replaced by mere counts.

## 4. P3: Sylvester gives both good completions and coherent hostiles

Put f(a)=a²−a+1=Φ₆(a), N=f(a), and suppose a≥3 and 3∤N.
The latter condition is a≢2 mod3. In Z/NZ set g=−a. Then

    g²+g+1=0,   g³=1,
    gcd(g−1,N)=gcd(a+1,3)=1,   gcd(g+1,N)=1.

Consequently g has order exactly 3 and H={1,g,g²} satisfies (G).
We get 2^((N−1)/6) legal connection sets and the completion bound (2).
For the Sylvester sequence s₀=2, s_(j+1)=f(s_j), every N=s_j≥7 qualifies.
The first bit counts are 1 at N=7, 7 at N=43, and 301 at N=1807.
This is an infinite family of quasirandom tournament completions carrying
the prescribed Sylvester order-three symmetry, even when N is composite.

**The lattice map behind the formula.** Let ω²+ω+1=0. Then

    Norm(a+ω)=a²−a+1,
    Z[ω]/(a+ω) ≅ Z/NZ,        u+vω ↦ u−av modN.

The isomorphism follows by substituting ω=−a into its defining polynomial;
it preserves addition and multiplication, not merely cardinality. The
Eisenstein unit ω acts as g, and its six signed unit rotations give the
packets ±{x,gx,g²x}. This makes the hexagonal-lattice analogy an actual map.
The norm identity was already present in the June29
[Eisenstein reflection](../../../07-reflections/the-covering-min-is-an-eisenstein-norm-phi6-and-the-honest-lrc-barrier-opus-20260629.md).
Only that algebraic identity is inherited: its covering-min/LRC discussion
is historical and subject to the subsequent
[covering correction](../../../07-reflections/CORRECTION-the-covering-min-is-1-over-n-the-even-block-not-the-construction-opus-20260630.md).
The new application here is the exact orbit covariance and its completion
bound, not a revived geometric covering claim.

**A hostile inside exactly the same family.** Fix a nonzero frequency k.
Each packet has six entries, so |c_O(k)|≤6. By (1),

    Σ_O |c_O(k)| ≥ Σ_O c_O(k)²/6 = N/6.

Choose ε_O to agree with the sign of c_O(k), with either value at zero.
Then |λ_k|≥N/6. Thus the same Sylvester family contains sequences that
are not quasirandom. For completeness, a bounded circulant symbol has a
large Fourier coefficient only if some cut has a comparable N² discrepancy:
write the character vectors as real/imaginary parts, split positive/negative
parts, and express [0,1]-vectors as integrals of level-set indicators. This
loses only an absolute constant. Symmetry, regularity and exact mean isotropy
therefore do not force the behavior of an individual orientation.

**Prime-factor information.** If p≠3 is prime and p|f(a), then a has exact
order 6 modulo p: a³=−1, and order 2 would force a=−1 and p=3.
Thus p≡1 mod6. This explains the eligibility of both factors

    1807 = 13·139.

Their roles are not identical: 13≡1 mod4 has an even quadratic character,
whereas 139≡3 mod4 supports a Paley tournament. The order-three construction
above works for both. A Fermat prime p>3 is 5 mod6, so cannot divide a
Sylvester term from 7 onward. Separately, its multiplier group is a 2-group,
which has no nontrivial odd subgroup usable in (G). This is the exact
arithmetic separation behind the Fermat/circulant rigidity in THM-4644.
The trivial multiplier group H={1} still satisfies P1; Fermat rigidity is
therefore not an obstruction to quasirandom circulant completions.

### P3b: a deterministic Sylvester family, with an explicit decay rate

Random choices are not necessary to get convergence to quasirandomness here.
For j≥2, N=s_j is3 mod4. Choose, deterministically, the least prime divisor
p of N which is3 mod4; at least one exists. Let M=N/p and g=−s_(j−1).
For d not divisible by p put σ(d)=(d/p), the quadratic character modulo p.
On the kernel pZ/NZ, choose any odd g-invariant sign function, for example
orient each previously unoriented orbit by its least residue. P3 ensures
that these kernel orbits and their negatives are distinct. The action of g
modulo p has order3, so its quadratic character is1; the entire symbol is
g-invariant.

Order vertices as r+pt, 0≤r<p and 0≤t<M. The matrix is exactly

    K = K_Paley,p ⊗ J_M + I_p ⊗ K_inner.

Both matrices on the right have zero row sums on their tournament factor.
On block-constant vectors the first term has norm M sqrt(p); on their
orthogonal complement only the second term acts, with norm at most M−1.
Consequently

    ||K||₂ = M sqrt(p) = N/sqrt(p).                       (3b)

The Paley identity used here is elementary: for d≠0,
Σ_x χ(x)χ(x+d)=−1. After scaling it reduces to Σ_t χ(t(t+1))=−1;
counting y²=t(t+1), or equivalently (2t+1−2y)(2t+1+2y)=1, gives p−1
solutions and proves the identity. Oddness of χ gives the tournament.

There is an explicit index bound. If p|s_j, the residues s_0,...,s_j are
distinct and avoid1: an earlier repeat would prevent a first arrival at0,
and1 is fixed. The j+1 distinct residues s_1,...,s_j,1 all lie in the image
of the quadratic f over F_p. That image has (p+1)/2 elements. Hence

    p≥2j+1,       ||K||₂/N ≤ 1/sqrt(2j+1) → 0.           (3c)

This is a deterministic quasirandom sequence with the prescribed cubic
symmetry; it does not assume that the Sylvester terms are prime or squarefree.
Factoring N to select p is a finite definition, not a claim of efficient
factorization. At1807 the chosen p is139 and M=13: the outer Paley structure
and the inner cubic orientation give ||K||₂²=23491. This assigns the two
factors concrete, different roles. It does not erase the source labels or
identify the13-vertex inner tournament with a Paley tournament.

The asymptotic statement is much weaker than an O(sqrt(N)) norm estimate.
It also says nothing about the random density of individual integer orbits.

## 5. P4: Wieferich delay is an isotropy window; fixed clocks eventually obstruct mixing

Recall the already-proved order law. For an odd prime p, p∤b, put
d=ord_p(b), e=v_p(b^d−1), and let r be the odd part of d. For
g=b^(d/r),

    ord_(p^k)(g) = r p^max(0,k−e).                        (5)

For every k≤e, H=<g modp^k> has order r and injective reduction modulo p.
It satisfies (G), so (1)–(4) hold throughout this window. At the two
specified base-2 Wieferich primes e=2:

| p | d | r | g | choice bits at N=p² |
|---:|---:|---:|---:|---:|
|1093|364|91|16|6564|
|3511|1755|1755|2|3512|

The bit counts are inherited from THM-4644. The new consequence is that at
least 2^6563 and 2^3511 of these connection sets, respectively, satisfy (2).
At 1093 raw doubling is not an admissible odd multiplier: 2^182=−1 mod1093.
Replacing it by 16 must be stated; it is not the same clock.

A smaller exact control connects 3 and 11:

    3^5−1=2·11²,       ord_11(3)=5.

Thus base 3 has e=2 at 11, while base 4 has the same initial order 5 but
e=1. The former has isotropic mean Gram matrix at 121; the latter does not.
Being Wieferich is base-dependent. This does not identify the role of 11
in golden arithmetic, elliptic curves or moonshine.

**All-depth fixed-clock obstruction.** Keep p,b,g,e fixed and let k≥e.
The group <g modp^k> contains the entire principal-unit subgroup
1+p^e Z/p^k Z: its p-part has order p^(k−e), and this is the unique subgroup
of that order in the cyclic p-primary unit group. Hence on unit differences
an invariant symbol is constant on every fibre modulo p^e. Choose distinct
residue classes A,B modulo p^e whose difference is a unit. All |A||B|
cross-edges have one direction, with |A|=|B|=p^(k−e). Therefore every such
invariant tournament satisfies

    ||K||₂ ≥ p^(k−e)=N/p^e.                              (6)

This rules out quasirandomness as k→∞ with p and the native generator fixed,
regardless of how the remaining signs are selected. It explains a genuine
obstruction rather than interpreting a growing orbit count as random mixing.

The mean Gram identity also fails immediately when k>e: at every unit
displacement d the correlation is exactly −p^(k−e). To count it, only unit
x and x+d can share a signed H-orbit. There are (r−1)p^(k−e) elements h
with h−1 a unit, giving that many positive solutions, and r p^(k−e)
elements with h+1 a unit, giving the negative solutions. Their difference
is the claimed value. Thus e is the exact isotropy cutoff for the native
odd multiplier, not merely a sufficient window.

**What changing the generator repairs and loses.** Replace g at level k by
g_k=g^(p^max(0,k−e)). It has order r at every level, so (G) and the good
completion theorem return. This retains the prime-level odd torsion action
but discards the principal-unit motion of the original clock. It is not an
allowed replacement when that actual clock is the observable. Independently,
forcing an old orientation to pull back on every fibre can reintroduce a
monochromatic rectangle even for these repaired groups. Group compatibility
and orientation compatibility are distinct obligations.

Exact hostile controls: at (N,g)=(49,2), the mean correlation at every unit
shift is −7 rather than −1; at (121,4) it is −11. At (1331,3), the base-3
delay has ended and 1220 of 1330 nonzero shifts fail (1). The script also
checks (21,16), where 3|Φ₆(5) breaks the Sylvester guard.

## 6. P5: the complete 19-adic law of the displayed Sylvester cycle

For f(x)=x²−x+1, the orbit of 2 modulo19 is exactly

    2→3→7→5→2.

Let F=f⁴ and H=f¹². The following small, exact calculations suffice:

    F(2)=1807=2+5·19²,
    F'(2)=3·5·13·85=16575 ≡7 (mod19),
    ord_19(7)=3,
    H(2) ≡ 2+17·19³ (mod19⁴),      H'(2)≡1 (mod19).

**Claim.** For every k≥1, the exact period of 2 modulo19^k is

    4                         if k=1,2;
    12·19^(k−3)               if k≥3.                     (7)

**Proof.** The four distinct residues establish the period modulo19.
F fixes 2 modulo19². On x=2+19²t, the next digit modulo19³ obeys
t↦7t+5, whose orbit at 0 is 0→5→2→0. Thus the period at level3 is12.
Now apply the following elementary growth lemma to H at 2.

If P∈Z[x], p odd, P'(x)≡1 modp and v_p(P(x)−x)=a≥2, then

    v_p(P^t(x)−x)=a+v_p(t),          t≥1.                 (8)

To verify it, write δ=P(x)−x, A=P'(x). Taylor expansion along the orbit
gives P^t(x)−x=δ(1+A+...+A^(t−1)) modulo p^(2a) for each fixed t.
For t prime to p this has valuation a. For t=p the geometric sum has
valuation1, while 2a≥a+2, so the valuation is a+1. The return polynomial
P^p still has derivative1 modp. Iterate this step, then the prime-to-p step,
to prove (8) for every t. Applying a=3, p=19 gives exactly the required
additional powers of19. A return of f at level3 must already be a multiple
of12, so no other shorter periods have been omitted. ∎

The direct first-return loop checks k=1,...,6, obtaining
4,4,12,228,4332,82308. These checks support the proof; (8), not enumeration,
covers arbitrary k.

This is the structural analogy with Wieferich lifting: keep a return time,
its defect valuation, and its derivative. For multiplication the derivative
of a return is automatically1 modp. Here it is7, requiring an extra
threefold return before the growing phase. The numerical equality
3·5·13·9=1755=(3511−1)/2 uses canonical representatives of the mod19 cycle;
the actual integer derivative uses85, not9. Only its residue is intrinsic.
No connection to the base-2 Wieferich congruence at3511 follows from that
representative-dependent equality.

There is also a sharp local-entry hostile: f has roots 8 and12 modulo19,
but the displayed orbit misses both. Therefore 19 divides no Sylvester
term even though 19≡1 mod6. Existence of legal roots and actual entry into
their basin are different statements.

## 7. P6 and cited comparisons: Collatz clocks and group representations

For the usual **3n+1** map (interpreting the prompt's `3n=1` in this sense),
fix an odd positive endpoint u with 3∤u. An inverse odd edge has

    n=(u2^a−1)/3,   a≥1,   u2^a≡1 mod3.

The order law at3 gives ord_(3^(k+1))(2)=2·3^k. Among exponent classes
modulo2·3^k, exactly 3^k obey the integrality guard, and

    a ↦ (u2^a−1)/3 mod3^k                              (9)

is a bijection onto all source residues modulo3^k. Indeed 2 generates the
units and the admissible targets are exactly 1+3c. Representatives can be
made positive by adding the period; the resulting source is odd.
Within one fixed a mod6, incrementing a by6 has exact source period3^(k−1).
This recovers the three six-step charts while retaining the division by3
and its extra digit. It is the precise Fermat/order connection, not a
derivation of the multiplier3 or a proof of forward convergence. A specified
integer source must still equal the inverse expression; congruence coverage
and quasirandom averages do not settle its actual trajectory.

The [prior mod19 resonance theorem](../mod19_resonance_depth_20261004.md)
already proves that ord_19(2)=18 creates one shared ternary exponent digit
when six-step addresses are coupled by CRT. That order18, the polynomial
period4 in (7), and its multiplier order3 are different clocks at the same
prime. The common operation is retention of the return data; there is no
asserted conjugacy between their orbits.

**CITED quasirandom-group comparison.** Gowers, *Quasirandom groups*,
[Theorem4.5](https://arxiv.org/html/0710.3877), relates mixing of group
convolution to absence of low-dimensional nontrivial complex representations.
P1–P2 instead make particular symbols on a cyclic group have small Fourier
coefficients. Cyclic groups themselves still have one-dimensional characters
and are not highly quasirandom groups in that sense. The common framework is
to identify which nonconstant representations can carry a large correlation,
then prove they are suppressed. The mechanisms and quantifiers differ.

**Recovered Monster/Baby Monster cycle.** The actual orders, checked against
the [Monster ATLAS](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/M/) and
[Baby Monster ATLAS](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/B/), give
the small-prime core

    C=2^5 3^7 5^3 7^4 11 13²

of |M|/|B|. Its prime-exponent chain is13→2→5→3→7→4. The stronger existing
[Session43 repair](../session43_observers_phase_guards_20261023/cycle-geometry/PROOF.md),
§2, defines h_C(b)=max{j:b^j|C} for composite bases as well. It proves the
literal cycle4→2→5→3→7→4, with13→2 as an additional incoming edge. This
preserves divisibility height, but does not identify prime13 with composite4
as an arithmetic label. Replacing |B| by the actual centralizer order|2.B|
changes C to C/2 and yields2→4→2; normalization matters. Every prescribed
prime-state permutation can be encoded by suitable exponents, so cycle
shape alone is not a selective fingerprint of these groups.

The [prior level filtration](../monster_levels_square_sum_7_21_20261022.md)
and [Monster–Paley action](../session32_square_monster_20261023/paley/PALEY_MONSTER_BRIDGE.md)
are stronger recovered bridges: the latter identifies the 7A normalizer
quotient action with AGL(1,7), retaining orientation and converse. We do not
replace them by a new numerical analogy. Group orders and modular matrices
alone do not determine the least complex representation degree needed for
Gowers' quantitative bound. The subsequent incoming primary character-table
verification supplies the stronger representation comparison in §10.

Concurrent prime-Legendre quasirandomness uses a third mechanism, a quadratic
large-sieve bound. Its retained prime/squarefree labels and height cutoff are
essential. Neither its averaged estimate nor (2) supplies source-specific
ROOT payment. The remaining Collatz obligation stays an independently
positive, source-weighted estimate with the integer clock intact.

## 8. Historical correction, connection ledger, and next questions

The June30 [recursive-mindset reflection](../../../07-reflections/the-recursive-mindset-phi6-iteration-is-sylvester-the-observers-one-is-the-apex-egyptian-fraction-opus-20260630.md)
correctly identifies Φ₆ iteration and its Egyptian-fraction telescope, but
incorrectly calls1 a universal attractor of all these recursions. In fact
f(x)−x=(x−1)². For0<x<1, x<f(x)<1 so the orbit increases to1; 0 maps to1.
For x>1 it grows to infinity, and x<0 maps above1. Thus the exact real basin
of1 is[0,1]. Iterating x+1 also plainly diverges. The telescope sums
reciprocals to1 while its integer denominators grow without bound; a limit
of an observable is not attraction of the original states. A maintained
warning and a mistake entry preserve this correction lineage. No old LRC
covering claim in that reflection is used.

| source → target / map | preserved predicate | lost information / sidecar | cheapest decisive test |
|---|---|---|---|
| Φ₆(a) → multiplier H=<−a> | legal orientation, cubic symmetry | specific sign choices; retain all packet bits | gcd(a+1,N), then exact row covariance |
| orbit-pair ensemble → Fourier spectrum | expectation and tail bound | a specific realization can align a frequency | choose signs of one Fourier row |
| Wieferich order tower → isotropy window | actual fixed multiplier through depth e | after e, unit fibres acquire common orientation | N=49,g=2; N=121,g=3 versus4 |
| changing generator → restored free action | odd torsion at the prime | principal-unit clock is removed | compare orders at p^(e+1) |
| polynomial return → lifted periods | marked source, multiplier and defect | a bare cycle diagram forgets both | f^12(2) mod19⁴ |
| group-order exponents → divisor heights | divisibility, literal forward cycle | incoming prime label and normalization |13 versus4; C versus C/2 |
| inverse exponent → source residue | exact division and clock | residue class is not a fixed integer orbit | test (9), then require integer equality |

Next questions, deliberately narrower than universal Collatz:

1. **OPEN:** improve the deterministic Sylvester rule in P3b to O(sqrt(N))
   norm, or obtain a substantially better unconditional lower bound for its
   selected prime factor. P3b proves only normalized decay O(j^(−1/2)); the
   probabilistic proof gives O(sqrt(N logN)) absolute norm. The aligned-row
   hostile must be excluded by an actual cancellation rule.
2. **OPEN:** characterize prescribed shared-port relations whose partial
   symbol has spectral cost o(N), without checking an entire large matrix.
   Equation(4) supplies the completion step once that cost is paid.
3. **OPEN:** for native arithmetic transport, can a conductor-resolved bound
   retain the exceptional fibres explicitly and charge their contribution
   to the marked source? P4 shows why deleting the principal-unit clock or
   averaging the source away would answer a different question.

## 9. Reproduction and audit scope

`python3 -B 05-knowledge/results/quasirandom_order_lifts_20261029/verify.py`

The verifier declares every finite universe. It checks the six Sylvester
seeds and factorizations; exact 19-adic first returns through k=6; both
methods of covariance at N=7,13,19,43 (packet incidence versus exhaustive
actual orientations); packet incidence also at121 and1807; four ramified
hostiles; the four order towers through k=4; and all admissible source
residues for u∈{1,5,7,11}, k≤5. It verifies the inherited height-cycle repair
and flags the representative-dependent derivative coincidence. It also
checks the first-arrival prime bound through s_5 and the deterministic
construction at N=7,43,1807, including exact Paley correlations and the
outer/inner orientation constraints. A further complete census checks all
admissible multiplier subgroups for every odd3≤N≤101. Such a subgroup
injects into F_p^× for any p|N and is therefore cyclic, so enumerating the
distinct cyclic subgroups really is complete for the stated guard. Normal and
optimized Python must agree; checks do not depend on `assert`.

The proofs were self-audited for: labelled versus unlabelled counts; all
nonzero rather than just unit differences; imprimitive characters; fixed
versus changing generators; exact periods rather than return multiples;
source division by3; expectation versus individual matrices; and unchanged
status of universal coverage. No independent-agent audit is claimed.

## 10. Integration of the two incoming packets

The [concurrent arithmetic-structure proofs](../quasirandom_arithmetic_structure_20261029/PROOFS.md),
§3, give an especially important opposite reading of Sylvester. If the
vertices are the odd terms themselves and the pairwise observable is their
Jacobi symbol, then for i<j, (s_j/s_i)=1 and (s_i/s_j)=−1. The tournament
is transitive at every size, including composite terms. This is an exact
hostile, not a contradiction to P3b: our vertices are residues modulo one
term, with a different pairwise observable. Thus the SAME sequence indexes
a deterministic quasirandom family while its terms form a completely
ordered family under reciprocity. Object, observer and asymptotic scale
must all be specified before saying that arithmetic is quasirandom.

Their §4 proves P5 independently and also identifies the unique compatible
19-adic four-periodic point above2. Its representative modulo19³ is5417,
whereas the fixed integer source2 has period12 there. Our proof concerns
the latter source. The affine digit map t↦7t+5 fixes15, giving
5417=2+15·19², and the Hensel derivative7−1 is invertible. An infinite
compatible periodic lift is therefore not the marked integer we started
with. This is an exact analogue of the source-retention obligation.

[Session107's audited character theorem](../session107_quasirandom_arithmetic_20261029/characters/CHARACTER_STRUCTURE.md)
adds a second useful way to preserve structure: fixed square-class tests
produce an affine binary code, while moving prime edges are quasirandom
inside its attainable colour fibres. Its unequal-fibre bound is
(B_X/2) Σ_(i<j)(|S_i||S_j|)^(−1/2), plus an explicit collision term.
This retains named ports but requires sufficiently large fibres. The
companion arithmetic packet identifies the earlier seven-prime Hamming
bank's character code with its simplex dual: pairwise independent bits
still obey exact Fano triple relations. Neither result permits a growing
bank without new uniform estimates. P2's pinned-packet spectral budget
is a complementary completion criterion, not a replacement for those
arithmetic fibre-density hypotheses.

The [incoming group note](../session107_quasirandom_arithmetic_20261029/groups/GROUPS_AND_RECOVERED_PRIMES.md)
correctly separates group quasirandomness from a single Paley symbol.
The local groups7:3 and7:6 have cyclic quotients and minimum nontrivial
complex degree1. Hence the actual Monster-local-to-Paley map does not
transport the Monster's quasirandomness parameter. The whole Monster and
Baby Monster cannot act nontrivially as finite-tournament automorphisms:
an involution would reverse an arc, and simplicity would make a nontrivial
action faithful. Ordinary graphs or bipartite incidence operators remain
valid carriers.

A primary-source check, prompted by that incoming note, verifies the
representation restriction

    196883 ↓_(2.B) = 1 + 4371 + 96255 + 96256

as a decomposition into the specified irreducible constituents, not just
an equality of dimensions. See Breuer–Magaard–Wilson,
[computational supplement, §2, printed p.3](https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/doc/ctblm.pdf).
The96256 constituent sees the central involution; the other three do not.
The same low-degree list certifies the Baby minimum4371. This is a concrete
reason restriction can create fixed vectors and smaller correlation modes.
The [verified Monster table](https://arxiv.org/html/2412.12182v1) and the
incoming audit supply the ordinary-character context; no table was
recomputed in this session. This representation map is much more informative
than identifying the valuation-cycle labels with group dimensions.
