# Finite classification of the ternary additive actuators

**PROVED in this note:** an exact bounded coefficient graph contains every
additive return of the stated defect. **FINITE COMPLETE:** its reachable
graphs and elementary cycles. The integer translation is retained without
a bound. This classifies a recursive operation; it is not a classification
of all zero addresses or all tournament witnesses.

## 1. Recovering the retained coordinate

The binary recursion and the ternary carry channels are inherited from
[RESIDUE_PHASE.md](../../h13_paired_descent_20261111/square_graceful/RESIDUE_PHASE.md).
The new [paired reset](../tournament/PROOF.md) uses the defect

    D(n)=13g(n)−3g(3n)−4n+6(n mod2)+8.

The older [guarded affine work](../prime_layers/PROOF.md) separates a finite
multiplicative phase from an integer translation. Here the analogous
separation is exact and particularly small. Put

    E(j,u,v,w;n)=13g(n)−3g(mn+j)+un+v(n mod2)+w,

where m is positive and odd, and 0≤j<m. Substituting n=2q+a, let
J,b be the quotient and remainder of ma+j on division by two. Then

    j' = J,
    u' = 2u+13[a=0]−3m[b=0],
    v' = 13[a=0]−3(−1)^J[b=0],
    w' = w+(u+v)a−3(J+(J mod2))[b=0].                 (1)

These formulas follow immediately by applying the defining recursion of g
to n and mn+j. In particular the state is not reduced modulo any number;
w is an exact unbounded integer.

## 2. Every return loop lives in a finite band

**Theorem.** Every nonempty return of the coefficient state (j,u,v)
under (1) lies inside

    0≤j<m,  −13≤u≤3m,  v∈{−3,0,3,10,13,16}.          (2)

**Proof.** The additive part of u' lies between −3m and13. If u>3m,
then u'≥2u−3m>u; it cannot return. If u<−13, then u'≤2u+13<u;
it also cannot return. These inequalities hold at every following step,
so an excursion outside the band cannot reenter it. The six possible v
values follow directly from (1). The j range is invariant. Thus every
return lies in a finite set of at most6m(3m+14) states. ∎

The bound is a consequence of the exact affine recurrence, not a search
cutoff. It gives a finite, exhaustive synthesis algorithm for every fixed
odd multiplier. Each edge retains its integer w increment. Each directed
cycle is therefore a literal repeatable additive actuator, and its drift
is the sum of those increments. Distinct local cycles can only be composed
when their coefficient states agree or an actual connecting path is kept.

## 3. The ternary channel forces the denominator five

For m=3, v after the first step belongs to {−3,0,10,13,16}; hence the
return universe has at most345 states. The bounded graph reachable from
(j,u,v)=(0,−4,6) has exactly20 vertices and22 edges. Its elementary cycles
are precisely:

| State (j,u,v) | Digits read from low to high | Address map | Drift |
|---|---|---|---|
| (0,−4,10) | 0 | n→2n | 0 |
| (0,1,13) | 1100 | n→16n+3 | 4 |

The proof certificate is the complete edge list in [receipt.json](receipt.json).
[compile.py](compile.py) generates this list directly from (1), explores
until closure, and exhaustively lists every elementary cycle. Any closed
walk decomposes into elementary cycles, so the table is a complete finite
certificate, independent of address length. The cycles are in separate
strong components, so no omitted mixed cycle changes their drifts.

The translation at the second state is exactly

    F(n)=13g(n)−3g(3n)+n+13(n mod2)+20,
    D(16n+5)=F(n),
    F(16n+3)=F(n)+4.

The periodic address of its nonzero cycle is

    0.00110011…₂=1/5,
    3/5=0.10011001…₂.

Thus the multiplication-by-three carry preserves a common period-four
clock and changes its phase from1/5 to3/5. The denominator five is forced
by the reachable nonconstant cycle of this defect, rather than chosen
because a different mathematical family has five members.

The other initial phase defects

    D_j(n)=13g(n)−3g(3n+j)−4n+6(n mod2)+8−12j

for j=1,2 have only2 and6 reachable states in the return band and have no
cycles. This does not exclude isolated zeros or nonlinear address changes.
It says that no terminal suffix starting from these states reaches an
additive loop with the same retained coefficient data.

The phase coordinate is essential even as a function, not merely a formal
label. For any j≠k in {0,1,2}, g(3n+j)−g(3n+k) is not affine in
n,n mod2,1: its values at n=0,2,4,6 have a nonzero second difference,
as recorded in the certificate. Once j is fixed, independence of n,
n mod2,1 is immediate. Hence a functional additive return in this family
must return the displayed coefficient state.

## 4. A whole family of denominator five actuators

The cycle has a useful extension beyond the specific13/3 weights. Let
a,b be arbitrary integers with a≡b modulo5, and put

    F_(a,b)(n)=a g(n)−b g(3n)+(6b−a)n/5+a(n mod2)+w.

For every n≥0 and every integer w,

    F_(a,b)(16n+3)=F_(a,b)(n)+2(a−b)/5.              (3)

Indeed g(16n+3)=g(n)+3n+(n mod2), whereas
g(48n+9)=g(3n)+18n+4. Substitution proves (3).
The divisibility condition is exactly what makes the linear coefficient
integral. The drift is controlled by the difference of the two weights.
In particular, weights25/5 give drift8; this is an arithmetic actuator
candidate for the H25/H5 charts, requiring a separate actual-chart match.
It does not assert that any two sources with those deletion counts can
be exchanged.

More generally let R=2^k and choose an odd m and a block r∈[0,R−1].
Fix a carry phase j such that s=mr−(R−1)j also lies in[0,R−1]. The two
addresses n and mn+j then follow synchronized repeated blocks r and s:

    m(Rn+r)+j=R(mn+j)+s.

Write c_z=R−1−rev_k(z), d_z=1−top_k(z), e_z=g(z)−1, and η=j mod2.
The inherited concatenation law gives

    g(Rn+z)=g(n)+c_z n+d_z(n mod2)+e_z.

For F(n)=a g(n)−b g(mn+j)+un+v(n mod2)+w, the difference
F(Rn+r)−F(n) has coefficients

    n:       a c_r−bm c_s+u(R−1),
    n mod2:  a d_r−b d_s(−1)^j−v,
    constant:a e_r−b(c_s j+d_s η+e_s)+ur+v(r mod2). (4)

Vanishing of the first two expressions gives an exact additive actuator;
the third is its drift. This is a constructive way to lift the modular
period and its affine carry to the tournament defect. It keeps the integer
phase j that a mere period-order calculation would discard. For m=7,
r=9 and j=4, the phase changes3/5 to1/5; the sibling investigation uses
this channel for another actual paired rebase.

## 5. What the finite object permits next

The all-prefix reset uses this sole nonconstant loop to discharge an exact
integer debt. Its repetition count can be unbounded even though its
controller has only four cyclic states. This is a finite description of
an infinite construction with the unbounded arithmetic retained.

To obtain different drift or a different target relationship, change an
actual ingredient: the selected C3 phase, a source chart, the leading
binary split, or the odd multiplier. Merely searching for a longer loop
in this same ternary channel cannot produce a new primitive actuator.
The [other multiplier investigation](../prime_layers/PROOF.md) makes those
changes at the level of literal count formulas.

The original target is still changed by terminal-suffix preparation.
The finite return theorem therefore supplies a complete actuator search,
not a proof that every prescribed H13 target has a rebase.

## 6. Reproduction

Run `python3 -B compile.py` here. The complete graph is the proof
certificate; 7,168 direct recurrence comparisons provide a separate
arithmetic control. They include both outgoing digits even when an edge
leaves the return band. The phase second-difference controls guard against
silently equating different carry channels.
