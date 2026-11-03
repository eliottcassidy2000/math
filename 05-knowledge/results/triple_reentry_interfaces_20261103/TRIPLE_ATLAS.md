# A typed atlas of triples: addresses, residues, valuations, and returns

This is a targeted recovery of September–November mechanisms, not a claim
to have reread every session containing the number three. The recent
[S163 survey ledger](../session163_ternary_triples_20261103/survey/META_LEDGER.md)
already surveyed the principal bridges. The present pass extends it with
older native phase clocks, path-cover substitution, valuation-balanced
Gaussian triples, reflected incidence signatures, and the precise
six-root return target. Those inherited results are credited, not re-proved.

A useful triple must have a type. Three **branch choices**, three **coordinates**,
three **residue classes**, and three **exceptional families** are different
objects. The recurring question is whether a map retains what the *next*
operation will inspect. The following contracts record source, target, map,
preserved predicate, lost data, needed sidecar, and a decisive test.

## 1. Tournament ternary arithmetic and actual source reentry

**Source → target:** a marked tournament in a proved source class → another
marked tournament in that class. **Desired scalar map:** `H -> 3H+2j`.
**Retain:** raw H, source predicate, observer port, guard, and any shared pins.
**Loss under the scalar projection:** future path-start and path-cover counts.
**Sidecar:** the marked observer triple and start-conditioned guard charge.
**Decisive test:** two consecutive nonzero returns from the same actual source,
with all marks transported. See [moving guard](../session163_ternary_triples_20261103/strong_returns/MOVING_GUARD.md)
and [first-step obstruction](../session163_ternary_triples_20261103/strong_returns/FIRST_STEP.md).

The new [six-root forest](RETURN_FOREST.md) supplies the arithmetic skeleton.
It does not supply the lift. This is the anchor against which the other
triples are compared.

## 2. Collatz residues 1,3,5 modulo six

For the positive odd accelerated map, an inverse predecessor is
`x=(2^k y-1)/3`. Integrality says `2^k y=1 mod3`: target `y=3 mod6` has
no odd predecessor; `y=5 mod6` requires odd k; `y=1 mod6` requires even k.
These are three **target residue types**, not three independently available
inverse children. Even after the parity constraint, the actual exponent
is unbounded and must be retained as an integer clock.

**Map:** an actual inverse edge → its target residue and exponent parity.
**Preserves:** legality of division by three. **Loses:** the predecessor,
valuation size, and descent/payment. **Sidecar:** full k and source y.
**Test:** enumerate actual predecessors, not just formal residue arrows.
The exact bridge and its limits are already in
[S163 strands](../session163_ternary_triples_20261103/strands/NOTE.md).

For tournaments, the analogous repair is to keep the source-dependent
return charge, not to equate the three residues with the three digits.

## 3. The two legs and hypotenuse; the complete Berggren return

The [September Berggren packet](../ternary_berggren_20260925.md) uses odd
coprime roots `s>t` and

    (A,B,C)=(st,(s^2-t^2)/2,(s^2+t^2)/2).

Its three parameter branches are `(s+2t,t)`, `(2s+t,s)`, `(2s-t,s)`.
They have intrinsic inverse chambers `s>3t`, `2t<s<3t`, and `t<s<2t`.
Except for root `(3,1)`, exactly one chamber applies and the first root
coordinate decreases under inversion. This is the successful model:
branch alphabet **plus** an invariant carrier **plus** a total parent rule.

**Preserves:** ordered primitive triples and root data. **Loses under an
address-only view:** the ordinary arithmetic/Collatz clock. **Sidecar:** the
actual roots and edge count. **Hostile:** legal Collatz edge `7->11` does
not become the ordinary Berggren parent edge; allowing arbitrary ancestry
distance still does not identify all plus-Collatz edges with ancestry.

The triple of numbers `3h,3h+2,3h+4` cannot be three sides of a right
triangle for odd h: all three are odd, already impossible modulo four.
The arithmetic-progression triangle occurs only at h=2, giving 6,8,10.
What transfers is the chamber/parent architecture, not the numerical sides.

## 4. Rational Collatz entry into the integer boundary

The [older rational bridge](../ternary_bridge_20260925.md) sends reduced odd
`s/t` to the signed triangle above. A plus update reduces the denominator
`t=3^r q` to `t/3` while 3 divides t; the prime-to-three q stays fixed.
At those steps the hypotenuse contracts by a factor below 1/3. If q=1,
after exactly r such denominator reductions the orbit enters the integer
stratum. The contraction stops being guaranteed there.

**Source → target:** rational odd charts → integer or prime-to-three-denominator
boundary. **Retain:** reduced roots, denominator valuation, sign of B.
**Loss:** forgetting q would falsely assert universal integer entry.
**Hostile:** `13/3 -> 7 -> 11`, with hypotenuses `89 -> 25 -> 61`.
**Use here:** a proved entry mechanism can terminate exactly where the hard
return obligation begins. A decreasing auxiliary coordinate is not enough
unless every boundary has its own legal continuation rule.

## 5. Three-color Zeckendorf recursion: a two-dimensional carrier

The [zero-slot Fibonacci study](../session41_owned_attachments_fibonacci_20261023/COLORED_GRID_AND_PHASE_QUOTIENT.md)
fixes the substitutions `C->CA`, `A->B`, `B->CA`. Its three-color incidence
matrix factors through a two-register Fibonacci carrier: `q(c,a,b)=(c+b,a)`
and `j(x,y)=(x,x,y)`, with `qj` the Fibonacci matrix. Modulo two the induced
map cycles the three nonzero vectors. The third phase is a quotient feature,
not an arbitrary permutation of the three original color types.

**Preserves:** canonical expansion and the macro-tile recursion.
**Loses:** the C-minus-B contrast in the quotient. **Sidecar:** expansion
phase and the literal zero slot. **Hostile:** deleting the zero slot and
joining its former neighbors falsely forbids the legal `11` choice.

The more recent [Zeckendorf inverse grammar](../inverse65_zeckendorf_20261103/ZECKENDORF.md)
uses two registers `(X,Y)->(Y+d,X+Y+2d)` and a previous-digit guard.
Its append-only prefix obstruction and zero-cycle length restrictions are
already stronger than a scalar Fibonacci analogy. A finite-state tournament
normalizer must retain an analogous carry/guard, and then be tested against
these existing obstructions before another small automaton search.

## 6. The static classification {p,p^3,p^2qr}

The [divisor-balance classification](../collatz_mod6_20260917_divisor_balance_family.md)
proves, for N>=2, `F=U+S` exactly for primes p, prime cubes p^3, and p^2qr
with three distinct primes. F counts proper nontrivial divisors, S those
that are squarefree, U those that are prime. The prime case has U=0,
not the prime-support size one.

**Map:** exponent box → its squarefree Boolean subcube and prime atoms.
**Preserves:** the divisor-count equation when full exponents are retained.
**Loses under radical:** which prime is squared and the entire balance defect.
**New exact application:** the three promotions of1335 in README have the
same unit-interface relation but all have balance `(10,7,3)`, unlike1335.
**Test:** use p^2 or p^2q as negative balance controls; do not infer equality
from the three-bit support alone.

## 7. Three prime choices, four Gaussian charts

The inherited [triple codec](../session163_ternary_triples_20261103/triple_codec/NOTE.md)
uses N=p^2qr with distinct primes all congruent to one modulo four. Gaussian
conjugate choices at three primes, quotiented by simultaneous conjugation,
give `F_2^3/<111>`: four charts, with representatives empty or a singleton.
The same addresses can label the four primitive Pythagorean classes at N.

**Preserves:** prime-choice addresses with a chosen ordering and Gaussian
factors. **Loses:** conjugation sign; the chosen chart must restore the
leg convention. **Hostile:** Gaussian normalization followed by rational
Collatz normalization need not preserve divisor balance; the inherited
36725-to-145 example changes from p^2qr to a two-prime number.
**Use:** shared addresses can be exact even when the arithmetic operations
do not commute. Require a checked square of maps before importing a return.

## 8. A directed triangle needs full path-cover memory

The [October ternary-memory packet](../session77_ternary_color_20261027/TERNARY_MEMORY.md)
uses THM-3121's path-cover substitution kernel. For a directed triangle
of factor ports, the walk-content series is

    (x+y+z+xy+yz+zx+3xyz)/(1-xyz).

In particular `H(C3[T,1,1])=3H(T)+2 pc_T(2)`, where pc counts unordered
spanning path covers. Two strong five-vertex sources both have H=15 but
two-path-cover counts39 and45; the same attachment gives123 and135.

**Preserves:** composition after retaining the complete run/path-cover
profile. **Loses under H:** the way global paths split into factor runs.
**Sidecar:** actual path ownership before dividing by factorial symmetry.
**Test:** attach the same named singleton ports to both H=15 sources.
**Use:** this is a concrete candidate for a larger reentry carrier. It is
closed under substitution, but its dimension grows with source order and
it is not automatically closed under re-marking or nonuniform cuts.

## 9. Native three-phase clocks and noninvertible division

The [native fraction packet](../session44_guarded_native_moves_20261023/fractions-ternary/README.md)
counts fixed-endpoint binary difference words modulo three using `I+Q`.
Its Jacobsthal continuation ratios have an ordinary native clock r but
exponentially larger Berggren address depth. The polynomial lift `I+zQ`
retains exact displacement; reducing z-information to a phase loses it.

**Map:** integer transfer → phase residue. **Preserves:** endpoint phase,
not displacement or execution cost. **Sidecar:** z degree and integer clock.
**Test:** to divide a numerator by three and recover the result modulo m,
keep the numerator modulo3m. At `(r,delta)=(2,0)` and `(6,0)`, numerators3
and66 agree modulo9 while quotients1 and22 do not. This is a particularly
useful rule for nonuniform interfaces over arbitrary composite moduli.
It differs from the unit predicate, which discards prime exponents safely.

## 10. Three reflected incidence colors: joint data matters

The [order-eight incidence study](../session38_incidence_ports_electrical_20261023/incidence/TERNARY_SIGNATURES.md)
colors `u->v` by `(c_v+c_uv) mod3`. Across all6880 classes, this three-color
signature has6864 labels, while the authenticated signature/converse pair
separates classes. Reflection-constrained forest elimination repairs four
middle-edge ambiguities which a singleton-only rule cannot repair.

**Preserves:** the audited port table under its reflection constraints.
**Loses:** native source paths if only the unlabelled signature is kept.
**Sidecar:** authentic converse and actual class edges, including parallel
edges. **Test:** removing a shared owner or merging parallel edges must fail.
**Use:** reconstructing a return from separate marginals needs compatibility,
just as primewise unit witnesses need one shared cut. The three-color
minimum is only within the stated quotient/forest-certification model.

## 11. Valuation-balanced Gaussian triples: a local condition with a gap

The [Gaussian gap-triple study](../session58_hourglass_complement_geometry_20261025/tropical-gap-triples/RESULT.md)
indexes gaps by signed exponent vectors. A signed three-gap relation forces
the largest absolute exponent at each prime to occur at least twice.
For squarefree centers the local row alphabet has21 elements (weights0,2,3
in F_3^3), and the balanced-candidate counts are

    (21^k-6*7^k-3*5^k+14*3^k-6)/48.

**Preserves:** a necessary valuation condition. **Loses:** the signed
numerical gap equation. **Sidecar:** actual Gaussian indices, signs and
normalization. **Hostile:** a balanced candidate is not thereby an hourglass.
**Use:** this is the older analogue of our composite-modulus warning:
locally compatible projections need not share a global realization.
The numbers7 and21 here count local words, not forbidden tournament counts.

## 12. Three norm phases separating 21 into seven-point geometries

The [integer-return arithmetic bridge](../ternary_integer_returns_20261029/RANKS_AND_CONNECTIONS.md)
constructs `PG(2,4)=F_64^*/F_4^*` and the well-defined norm phase `x->x^21`.
Its three seven-point fibres are Fano subplanes. Conditioning the uniform
21-point model on one norm phase recovers a seven-point model.

**Preserves:** an exact finite incidence geometry and conditional distribution.
**Loses:** the number of discarded zero digit blocks and the scalar lift.
**Sidecar:** phase, scalar and original integer clock. **Test:** enumerate
all63 nonzero field elements and retain their actual quotient representatives.
**Use:** it illustrates a genuine three-sheeted refinement. It neither
proves pointwise Collatz coverage nor identifies these geometries with raw
Hamiltonian counts7 and21. Equal cardinalities alone cannot transfer a ban.

## Synthesis for the next experiment

The successful mechanisms have a common specification, not one common
mathematical object: **a legal carrier, a decoding rule, and a boundary rule
that remains valid after the next operation**. Berggren provides all three;
Zeckendorf adds a carry; path-cover substitution adds run ownership; the
unit-interface compiler adds joint prime supports and named cuts.

The most useful proposed experiment is therefore not another matching of
three numerical labels. Give a tournament return source an expanded state
consisting of its marked path-cover profile plus the nonuniform interface
cost profile. Test whether nonzero returns can change observer chart while
preserving shared ports and then return a second time. The full profile is
an intentionally generous candidate, to be compressed only after hostile
sources cease to be distinguishable by the allowed future operations.
