# Joint spectrum and the exact two-tail self-converse order image

**PROVED elementary consequences of the inherited ordered-SCC mechanism;
FINITE-EXACT independent controls through order six; CITED Busch cutoff.
The unrestricted all-order Hamiltonian-count classification remains OPEN.**
No novelty claim is made for the inherited series, multiplicities, or gaps.

## Inheritance and correction boundary

The closest proved mechanism is ordered strong-component multiplication:
[THM-1862, order-join reduction](../../../01-canon/theorems/THM-1862-order-join-reduction-principle.md).
The canonical hostile is the transitive triple versus the directed triangle:
both ordinary classes are self-converse, but their H values are 1 and 3.
The corrected near miss is “7 and 21 are the only omitted values”: the
all-order exclusions hold, but completeness is OPEN; see
[THM-1370, omissions](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md),
MISTAKE-559, and the later explicit corrections in `01-canon/MISTAKES.md`.
The least-used sidecar here is **the parity of the central SCC's order**.

The exact fixed-H ordinary/self-converse/converse-merged rational series and
eventual polynomial/quasipolynomial behavior were already derived in
[Session24, multiplicities §§1–2](../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md).
The present square-central formula is a recovery, not a new theorem. The
additional extraction below identifies the **exact order support** of each
self-converse H stratum as two independently padded tails, including empty
tails. Session136 already gives the ordinary minimum-order min-plus rule in
[its spectrum structure §5](../session136_full_tournament_spectrum_20261102/spectrum_structure/THEORY.md).

Anchor: exact joint class/H image. Niche: converse-fixed order support.
Wildcard: whether an H-preserving symmetry can explain A000568 parity.
The live board is **SCC words / strong count inventories / square-central
transport / order-parity thresholds / presentation fibers**. The square
map doubles the half's vertex cost and squares its H; it does not retain
source path positions, attachment ports, or arbitrary module-substitution
responses. No claim about full attachment-carrier minimality is used.

## 1. The requested joint enumerator has a typed SCC recurrence

Let H(T) count directed Hamiltonian vertex orders, with H(empty)=1. Define

```
F_n(z) = sum_[T], |T|=n z^H(T),
S_n(z) = sum_[C], |C|=n, C strong z^H(C),
```

where singleton tournaments are strong and S_0=0. These count **ordinary**
isomorphism classes; converse classes remain separate. On polynomials with
positive-integer exponent labels define the bilinear multiplication

```
z^a star z^b = z^(ab).
```

Its multiplicative identity is z, not the ordinary polynomial constant 1.
Then, for n>=1,

```
F_0(z)=z,
F_n(z)=sum_(k=1..n) S_k(z) star F_(n-k)(z).             (1)
```

**Proof.** Every nonempty tournament has a unique first strong component and
a unique remaining tournament in the dominance order. Every path concatenates
one path from each component. Isomorphisms preserve this ordered decomposition,
so each pair of ordinary classes contributes once. QED.

Equation (1) couples support and multiplicity exactly; at z=1 it becomes the
ordinary strong-component count convolution. It is not a recurrence from
previous scalar a_j alone: the unknown strong H inventories remain inputs.
It is not legal to use ordinary polynomial multiplication in (1), which
would add H instead of multiplying it.

Fresh all-class controls give

```
F_3=z+z^3,
F_4=z+2z^3+z^5,
F_5=z+3z^3+2z^5+2z^9+z^11+z^13+2z^15.
```

Their evaluations at one are 2,4,12; the order-six polynomial and all
converse-fixed refinements are emitted by the checker.

## 2. The complete fixed-H order image under converse symmetry

Use x for the vertex-order variable to avoid confusing it with F_n's H marker.
Write

```
G_h(x)=sum_n #{ordinary classes T: |T|=n,H(T)=h} x^n,
Q_h(x)=sum_n #{self-converse classes T: |T|=n,H(T)=h} x^n.
```

For h>1 let R_h(x) enumerate **nontrivial strong self-converse** classes
with H=h, by order. Set R_1(x)=1+x, encoding an empty or singleton center.
The inherited square-central identity is

```
Q_h(x) = sum_(u^2 divides h) G_u(x^2) R_(h/u^2)(x).    (2)
```

Here a zero G_u or R_v is permitted. Each R_h is a finite polynomial by the
strong-order lower bound; no strong-spectrum completeness is assumed.

**Exact map.** A converse-fixed SCC word has the unique form

```
U -> C -> U^op,
```

where U is any ordinary tournament (possibly empty), and C is either empty,
a singleton, or a nontrivial strong self-converse tournament. Its right-hand
SCC word is forced by its left-hand word; there is no quotient by swapping
halves. The map preserves ordinary isomorphism class and gives

```
|T|=2|U|+|C|,       H(T)=H(U)^2 H(C).
```

Uniqueness comes from taking the first half of the unique ordered SCC word,
not by arbitrarily splitting the vertex set. A self-converse tournament with
an even number of SCCs has empty center; with an odd number it has one central
strong SCC. This proves (2) and its injectivity/surjectivity.

### Minimum-order and exact-tail consequence

Let mu(u) be the least order of an ordinary witness of H=u, with mu(1)=0;
use infinity if no witness exists. For epsilon in {0,1}, let theta_epsilon(v)
be the least order of a permitted center with H=v and order parity epsilon.
Thus theta_0(1)=0 and theta_1(1)=1. For v>1 this minimum runs only over actual
strong self-converse classes. Define nu_epsilon(h) similarly over all
self-converse tournaments of H=h and order parity epsilon. Then

```
nu_epsilon(h)
   = min_(u^2 divides h) [2 mu(u)+theta_epsilon(h/u^2)]. (3)
```

The infinity conventions are important: formally admissible arithmetic
factors are not substitutes for actual witnesses.

**Proof.** Decomposing any witness via the exact map bounds its size below by
the appropriate term of (3). Conversely, minimum witnesses U and C for any
finite term construct a self-converse witness of that size. There are finitely
many square divisors, so the minimum is attained when finite. QED.

More strongly, for every n>=0,

```
there is a self-converse n-tournament with H=h
  iff n>=nu_(n mod 2)(h).                              (4)
```

The forward implication is the definition of the minimum. For the reverse,
start with a minimum witness and add a singleton source and singleton sink
in converse positions. This changes n by 2 and leaves H unchanged, and can
be repeated. Thus each parity's support is exactly one tail or the empty
set. The two thresholds need not agree or both be finite.

This is an exact image theorem **relative to the actual strong-center
inventories**. It is an effective criterion for each fixed h because the
classical strong floor bounds the finite inventories needed. It does not
classify all h for which the center input exists, and it does not solve the
unrestricted spectrum conjecture.

## 3. Complete self-converse order support for odd H<=21

The fresh checker enumerates every ordinary class through order six and
independently counts Hamilton paths by subset DP and by literal vertex
permutations. Busch gives f(7)=25 and increasing f(n), so the resulting strong
inventory for H<=21 is complete at **every order**. This is a CITED cutoff
plus FINITE-EXACT base, not an extrapolation of a census.

| H | least even order nu_0 | least odd order nu_1 | all subsequent allowed orders |
|---:|---:|---:|---|
|1|0|1|every order (empty convention at zero)|
|3|infinity|3|odd n>=3|
|5|4|infinity|even n>=4|
|7|infinity|infinity|none|
|9|6|5|every n>=5|
|11|infinity|5|odd n>=5|
|13|infinity|5|odd n>=5|
|15|infinity|5|odd n>=5|
|17|6|infinity|even n>=6|
|19|infinity|infinity|none|
|21|infinity|infinity|none|

For H=9, the odd minimum is a strong order-five center. The even minimum
is U=C3, C=empty, producing C3->C3. These mechanisms have the same H but
different SCC words and different order parities.

For H=19, ordinary witnesses exist at order six, but **no self-converse
witness exists at any order**. Since 19 has no nontrivial square divisor,
(2) requires a self-converse strong H19 center; the complete strong inventory
has two order-six classes forming a converse pair, and no other orders.
This obstruction was already recorded in Session24 and is independently
replayed here. H15 provides a finer hostile: its two strong order-five
classes are self-converse; its two strong order-six classes are not.

## 4. A universal obstruction to an H-preserving free pairing

**PROVED.** At every n>=0, the coefficient of z in F_n(z) is exactly one.
Consequently no fixed-point-free involution on ordinary order-n classes can
preserve H, even at orders where their total cardinality is even.

**Elementary proof of the boundary.** A transitive tournament has its unique
increasing Hamiltonian order, so H=1. Conversely, a nontransitive tournament
contains a directed triangle. The triangle has three Hamiltonian orders.
Each can be extended, one new vertex at a time, by inserting the vertex at
a feasible position in the existing path. Existence follows from the first
out-neighbor/last in-neighbor binary transition argument. The extension is
injective because deleting the inserted vertex recovers the old order.
Thus H>=3. There is one transitive class at each n. An H-preserving involution
must fix its singleton stratum. QED.

At n=3 both classes are also fixed by converse. Hence the parity explanation
must either change H, work on an enlarged presentation carrier, or retain
fixed points. Reversing a path under converse preserves H but does not remove
the fixed-class locus. Even cardinality alone supplies no particular diagram.

## 5. Exact remaining interfaces and reproduction

The inherited ordinary minimum-order recurrence is

```
mu(h)=min(sigma(h), min_(1<d<h,d|h)[mu(d)+mu(h/d)]),
```

where sigma(h) is the least order of a strong witness. Its support at order n
is exactly {h:mu(h)<=n}, by neutral singleton padding. The new parity-specific
formula (3) is the converse-fixed counterpart of that cost-sensitive map.
Both retain an actual realization guard in their strong inputs.

The unrestricted omissions 7 and 21 are all-order theorems. Completeness
reduces to the two infinite strong lanes p (odd prime p!=7) and 7p (odd prime
p!=3), with the inherited carry witnesses 49,63,343; see
[THM-4094, two-prime-lane completeness](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md).
Session142 further proves all prime-H witnesses above 3 asymmetric, and
restricts Aut for H=7p; see
[its minimal-gap proof](../session142_tournament_odd_spectrum_20261102/minimal_gap/PROOF.md).
These sharpen structural obligations without supplying universal witnesses.

The current finite unrestricted coverage through 999999 is inherited from
[Session124 via the square-sum synthesis](../square_sum_spectrum_overnight_20261031b/README.md);
that large bank was not rerun here. It must not be weakened to the older
80557 or 80265 bounds, nor strengthened into all-order coverage. The fixed-H
strong inventory through 623 in THM-4744 is a different assertion: it counts
all strong witnesses, rather than merely furnishing one witness per value.

Run from this packet directory:

```sh
python3 -I -B spectrum_probe.py > /tmp/spectrum-control.json
python3 -I -B -O spectrum_probe.py > /tmp/spectrum-control-optimized.json
cmp /tmp/spectrum-control.json /tmp/spectrum-control-optimized.json
```

No imported catalogue, repository computation module, or external dependency
is used. Universes: all ordinary classes at orders 0..6, obtained by adjoining
all 2^n incident-arc orientations to every order-n class, then canonicalizing
under all degree-compatible vertex orders. This is complete because every
larger tournament has a vertex deletion in the preceding universe and any
isomorphism of that deletion transports all its incident arcs. No graph
symmetry, score sequence, strongness, H, or converse filter is imposed before
generation. Score-compatible canonicalization is exact because isomorphisms
preserve vertex degrees. Controls compare Hamilton DP/permutations, class
orbit mass against all 2^binom(n,2) labelled tournaments, SCC convolution,
converse-fixed square-central convolution, and the positive/hostile examples
above. Correctness checks are explicit exceptions and survive `python -O`.

Primary sources re-opened 2026-11-02:
[Busch, EJC 13 (2006), N3](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/)
and [Chan–Pak, Conjecture 6.2](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf).
The latter states the complete-spectrum conjecture; no claim that it has
been solved follows from a finite certificate or from the present exact
image formulas.
