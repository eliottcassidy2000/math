# The 101 clock tower: exact power-free limits and closure under joint ports

Session `clock101_joint_ports_20261029`; starting tree `7c707a4ca50f`.

**Status:** PROVED elementary statements P1–P5 below; FINITE-EXACT arithmetic
receipts in [results.json](results.json), produced by [verify.py](verify.py).
The full-rank and three-deck mechanisms are classical linear/Fourier algebra;
no literature-first claim is made. This packet applies them to the recovered
Collatz clocks and derives the stated all-exponent formulas. These are
same-agent proof and independent-algorithm checks, not an independent referee
report. General Collatz convergence and arbitrary-word exclusion remain OPEN.

The main new consequence is an exact distinction between a small survival
formula and a carrier closed under future joint tests. The combined native
phase/residue clock has `3*101^e` states. Its unlabelled linear alarm loses two
modes, but **two timed alarm ports recover a point indicator**. Closing under
those joint tests therefore restores the full function algebra, of dimension
`3*101^e`. This is a fixed-word finite-residue theorem, not a minimum carrier
for unrestricted integer Collatz or arbitrary appended words.

## Inheritance, scope, and portfolio

The closest mechanism is the independently audited
[Session 101 squarefree calculation](../session101_denominator_connections_20261029/squarefree/README.md):
same forbidden residues modulo `101^2`, different visiting orders, sharp
limits 4,985 and 8,600. Its older source is
[September 21's squarefree guards](../collatz_guards_20260921_squarefree.md).
Neither limit is newly discovered here. The all-prime sieve in that proof
is load-bearing: a local surviving residue alone does not prove an integer
squarefree realization. We extend that proof to every power-free exponent.

Canonical hostiles: omitting the terminal node changes a sharp cutoff; the
two denominator-101 words give denominators 35 and 175 after the same append;
and the 81/27 tournament pair shares 67 labelled cyclic triples. The corrected
near miss is that a denominator prime need not give one full prime-power
cycle: Session 101's `(1,13)` word has carry divisible by 5. Our hypotheses
explicitly require a unit carry. The underused sidecar is the **individual
phase of each forbidden zero**, recovered in the joint clock below.

Anchor: prime-power clock lifting and exact power-free windows. Niche:
minimal observable algebras under joint ports. Wildcard: the older
`7x+169` cycle containing 79 and `19^2`. The live concept board consists of
affine carry, native phase, zero-alarm support, Fourier rank, and tournament
orientation. [RECOVERY.md](RECOVERY.md) records several waves of older sources
on each of 101,19,79,169, including negative and merely numerical connections.
Reserved THM-4848/4849 are not used as proved dependencies.

## P1. Lift the clock instead of enumerating its orbit

For the three positive valuation words in this packet the block return is
`F(x)=(A*x+C)/D`, with `D` a power of 2:

| Native odd map and word | A | C | D | Rational cycle |
|---|---:|---:|---:|---|
| `3x+1`, `(1,1,5)` | 27 | 19 | 128 | `19/101,79/101,169/101` |
| `3x+1`, `(1,3,3)` | 27 | 31 | 128 | `31/101,97/101,49/101` |
| `7x+1`, `(1,3,5)` | 343 | 79 | 512 | `79/169,361/169,337/169` |

The relevant prime is 101 in the first two rows and 13 in the last. In each,
`alpha=A/D=1 mod p` and `beta=C/D` is a p-adic unit. The inherited lifting
identity gives `v_p(F^j(x)-x)=v_p(j)`. Thus F is one full cycle modulo every
`p^e`, and is conjugate to the adding machine on the p-adic integers.

For a phase map `(a_j*x+c_j)/d_j`, its forbidden source residue is
`z_j=-c_j/a_j mod p^e`. Let `i_e(z)` be the unique clock index satisfying
`F^(i_e)(0)=z mod p^e`, with `0<=i_e<p^e`.

**Digit recurrence (PROVED).** Start `i_0=0`. Having `i_e`, calculate F and z
modulo `p^(e+1)` and set

`u_e = ((z-F^(i_e)(0))/p^e)*beta^(-1) mod p`,

`i_(e+1) = i_e+p^e*u_e`.

*Proof.* The inherited orbit identity, or the geometric sum, gives
`F^(p^e)(x)=x+p^e*beta mod p^(e+1)` for every x. Iterating this step u times
gives the displayed unique lift. The quotient in the formula is integral
because the lower e digits already agree. This holds also when
`v_p(D-A)>1`, as in `512-343=13^2`. QED.

One can also write the clock as a p-adic logarithmic ratio
`log(1-(D-A)*x/C)/log(A/D)`, but no numerical logarithms or discrete-log
search are needed. The checker uses exact digit lifting, binary affine
composition, and independently a geometric numerator computed at extra
precision before division by `D-A`.

For both 3-words the three named roots are
`(0,-1/3,-5/9)`; modulo 10201 these are `(0,3400,9067)`. Their named clock
indices are `(0,3795,5214)` for 115 and `(0,8601,8780)` for 133. Scaling
`x -> (31/19)*x` conjugates the two **unmarked** return maps. It does not
preserve their shared forbidden set. This identifies exactly what that
conjugacy loses.

### The same number of coordinates does not preserve the clock

The recovered unit-distance work uses the **field** `F_169`; our 7-word uses
the **ring** `Z/169`. Replacing the ring by the field destroys the return
period. The integer coefficients become `alpha=1,beta=8` over `F_169`, so
the return is translation by 8, with 13 cycles of length 13. Over `Z/169`
the return is translation by 151 and has one cycle of length 169. No set
bijection can intertwine these two permutations. Generally the corresponding
field `F_(p^e)` has `p^(e-1)` return cycles of length p, whereas the ring
has one of length `p^e`.

The missing coordinate operation is explicit. Modulo `p^2`, write
`alpha=1+p*c`, `beta=beta_0+p*beta_1`, and `x=u+p*v`, using digits in
`0,...,p-1`. Then

`u'=(u+beta_0) mod p`,

`v'=v+beta_1+c*u+floor((u+beta_0)/p) mod p`.

Even when c=0, the integer carry in the final term advances the upper digit
after a lower-digit circuit. Componentwise field addition omits it. Thus a
finite-field or Frobenius-inspired compression needs a carry operation as
well as matching cardinality. The checker verifies this digit law for every
state modulo `101^2` and `13^2` in all three examples.

## P2. Re-marking requires the native three-phase clock

Put `N=p^e`. The local state is `(j,x)`, with `j in Z/3` and `x in Z/N`.
Its transition is

`(j,x) -> (j+1, (a*x+1)/2^(v_j) mod N)`.

Each of the three phase maps is bijective modulo N. Since the three-step
return is one N-cycle, this combined transition is one cycle of length
`M=3N`. Starting at `(0,0)`, the three zero events occur at the named times

`T_j=3*i_e(z_j)+j`, `j=0,1,2`.

This is a genuine finite affine clock. It does not assert that its modular
cycle lifts to a periodic positive integer orbit. Finite integer valuation
guards are retained separately by the usual exact cylinder modulo `2^(H*r+1)`.

Sort the T's as `0=t_0<t_1<t_2<M` and set `t_3=M`. For a start phase j,
define three capacities

`K_(j,l)=floor((t_(l+1)-1-j)/3)-floor((t_l-j)/3)`, `l=0,1,2`.

**Endpoint-aware recurrence (PROVED).** The number of source residues whose
`r` full blocks have all `3r+1` displayed nodes nonzero modulo `p^e` is

`h_(e,j)(r)=sum_l (K_(j,l)-r)_+`.

The last allowable block count is `R_(e,j)=max_l K_(j,l)-1`. Equivalently,
`h(r+1)=h(r)-#{l:K_(j,l)>r}`; its second differences recover the capacity
histogram. This recurrence handles overlaps and the terminal node exactly.

*Proof.* In a gap `(t_l,t_(l+1))`, the start time t must satisfy
`t_l+1<=t<=t_(l+1)-3r-1` and `t=j mod3`. Counting those integers gives one
positive-part term. Summing the three disjoint gaps proves the formula.
For different j each capacity changes by at most one, so the three R values
differ by at most one. QED.

Here are the exact phase-zero limits; “e-free” means no prime e-th power
divides any displayed node. Higher e is a **weaker** arithmetic restriction.

| e | 115: last e-free repetition | 133: last e-free repetition | Larger limit |
|---:|---:|---:|---|
| 2 | 4,985 | 8,600 | 133 |
| 3 | 576,469 | 979,116 | 133 |
| 4 | 64,332,491 | 45,324,462 | 115 |
| 5 | 7,110,702,653 | 5,365,816,387 | 115 |

Thus the comparison of the two words reverses by e=4. A squarefree ranking
does not determine the ranking higher in the prime-power tower.

Native re-marking first matters in the displayed data at e=3:

`115 phases 0,1,2: (576469,576470,576469)`,

`133 phases 0,1,2: (979116,979115,979115)`.

Phase j means the rotated word beginning with its j-th valuation. Each
existence statement permits its own source in the named arithmetic row;
it does not claim that one original integer remains unchanged under re-marking.

### Why these are global power-free existence limits

**Theorem (PROVED).** For either 3-word, any `e>=2`, phase j, and fixed
integer `r>=1`, there are positive odd sources `n=5 mod9` realizing r copies
of the rotated word with all `3r+1` nodes e-free **iff** `r<=R_(e,j)`.
When they exist, they have positive relative density in the exact
word-and-row progression. The analogous statement for the 7-word holds
with source row `n=1 mod7` and p=13.

The local obstruction proves necessity. For sufficiency, at every other
prime use the rational cycle itself as an avoiding residue. The only
square exceptions are `13^2` in 115, `7^2` in 133, and `19^2` in the 7-word.
At e=2 replace the fixed point f by `f+q mod q^2` for that exceptional prime q:
its formerly square-divisible phase now has valuation exactly one, and all
other phases remain units modulo q. The return slope is a q-unit, so this
repair persists for arbitrarily many blocks. At e>=3 even that repair is
unnecessary: every numerator in the three displayed rational cycles has
prime valuation at most two. The source rows exclude the multiplier prime;
exact odd valuation words exclude 2. CRT makes all these local conditions
compatible with the exact word cylinder for every fixed r.

For completeness, write the displayed nodes on that progression as fixed
positive linear forms `L_i(s)=A_i*s+B_i`, `0<=i<=3r`. Away from 2 and the
multiplier prime, the slopes are units. If `nu(q)` is the number of forbidden
residues modulo `q^e`, then `nu(q)<q^e` locally and `nu(q)<=3r+1`. The product
`product_q(1-nu(q)/q^e)` is positive. For `0<=s<X`, primes q>z dividing
some `L_i(s)` to power e remove at most
`O_r(X*z^(1-e)+X^(1/e))` parameters. Finite CRT, followed by X tending to
infinity and then z tending to infinity, proves the density. All limits keep
e and r fixed. No uniform density, infinite all-power-free orbit, or general
convergence conclusion is asserted.

## P3. A compact alarm formula can require every state

Fix one block phase and let `S subset Z/N` be its three forbidden clock
indices. Observing the bit `f(t)=1_S(t)` after any number of **read-only**
block advances gives N distinct future behaviors.

More generally, for `N=p^e` and an integer-valued mask with total mass m not
divisible by p, its cyclic convolution matrix has full rank N over Q and
over `F_p`. Let P be the mask polynomial and J the cyclic shift. Modulo p,
`(J-I)^N=0` and

`P(J)=m*I+(J-I)*Q(J)`.

The second term divided by m is nilpotent, so the finite geometric inverse
proves invertibility. For our binary three-point masks, p is 101 or 13.
Equivalently the Fourier transform has no zeros: any cyclotomic factor
`Phi_(p^a)` would force p to divide `P(1)=3`.

Consequences (PROVED): every state is distinguished by some future alarm;
the Hankel rank is N; and the minimal homogeneous linear recurrence for the
periodic alarm has degree N. The compact Hensel formula is an arithmetic
representation with `ceil(log2 N)` state bits, not a small-dimensional exact
linear realization. The three-piece formula for **aggregate survival** does
not contradict this result. A process that irreversibly stops at its first
alarm has a different observation language and can use its distance to that
alarm instead; the read-only hypothesis matters.

Prime-power scope is essential for the rank claim. In `Z/105`, the binary mask
`{0,35,70,1,22,43,64,85}` has mass 8 coprime to 105, and 105 distinct future
states, but rational convolution rank **49**. Its polynomial has precisely
`Phi_15*Phi_105` in common with `X^105-1`. The two disjoint regular polygons
in that support create the missing Fourier modes. Coprime mass still rules
out a translation stabilizer, but does not give full linear rank at a
composite clock.

## P4. The smallest algebra closed under timed two-port tests

Return to the combined native clock `Z/M`, `M=3p^e`, with scalar alarm
`g(t)=1_{ {T_0,T_1,T_2} }(t)`. Its translates span a space of dimension
**M-2**, exactly. The missing modes are the two nontrivial pure phase
characters of order 3: each phase contains precisely one marked zero.

*Proof.* The support polynomial is `G(X)=1+X^(T_1)+X^(T_2)`, where named
times are used even if their order switches. Its exponents are 0,1,2 mod3,
so `Phi_3` divides it. A sum of three complex numbers of absolute value one
vanishes only when they form an equilateral triple. Thus a further root of
unity of order `3*p^a` would require `p^a` to divide both T_1 and T_2. Direct
reduction at e=1 rules that out for all three words, and compatible lifting
preserves that exclusion. No p-power root can be a zero. Therefore
`gcd(G,X^M-1)=Phi_3`, proving rank M-2. QED.

The minimal alarm recurrence has denominator

`(1-z^M)/(1+z+z^2)=(1-z)*(1+z^3+...+z^(3(N-1)))`.

This is also visible as equality of the three phase sums, each equal to one.
But this linear compression is **not closed under joint attachments**.

In each of our three clocks, all six ordered differences between distinct
zero times are different. It suffices to check this modulo `3p` (the checker
also checks every computed higher lift). With `d=T_1`, this gives the exact
two-port identity

`g(t)*g(t+d) = 1_{ {0} }(t)`.

Shifting this point indicator gives all M coordinate indicators. Hence:

**Theorem (PROVED).** The unital Q-algebra generated by the read-only cursor,
the scalar zero alarm, and pointwise joint-port products is exactly `Q^(Z/M)`.
Its dimension is M; all M deterministic states are necessary; degree-two
products suffice. Degree one cannot suffice, since its span has dimension
M-2. This algebra is closed under output complement and native re-marking.

The marked clock already has trivial translation stabilizer: a stabilizer
would preserve the unique ordered pair at offset d and therefore fix its
first point. Nevertheless the linear alarm still misses two modes. Thus
breaking all clock symmetry is not enough to make a linear observer faithful;
the timed product supplies strictly more than the stabilizer calculation.

At p=101,e=2, the numbers are `M=30603` and single-alarm rank `30601`.
One valid offset is 11386 for 115 and 25804 for 133. The isolated state is the
named source `(phase0,residue0)`. The same construction works at every e,
including the separate 13-adic family. Named zero-port observations also
restore full rank immediately. This is a concrete closure theorem, not merely
a statement that more information might be needed.

The state algebra is local to one fixed periodic word and prime-power modulus.
Appending a valuation changes the closure problem; the inherited
`19/101 -> 37/35` versus `31/101 -> 221/175` test remains a separate hostile
against claiming this is an all-word carrier.

## P5. Averaging over sources changes what joint ports can reconstruct

For a block mask f define the two- and three-port averages

`A(u)=sum_t f(t)f(t+u)`,

`B(u,v)=sum_t f(t)f(t+u)f(t+v)`.

A always agrees for S and its reflection -S. In contrast, if all Fourier
coefficients of f are nonzero, B determines f up to translation. This is a
classical three-deck/phase-retrieval result; see
[Jaming–Kolountzakis (2003)](https://nyjm.albany.edu/j/2003/9-12.pdf).
A short proof here: the Fourier transform of B is
`fhat(k)*fhat(l)*fhat(-k-l)`. Two masks with the same B have the same mass
and Fourier magnitudes. Their nonzero Fourier ratios therefore satisfy
`r(k+l)=r(k)r(l)`, so they form a character, exactly a cyclic translation.
P3 supplies the nonvanishing hypothesis on our prime-power block clocks.
For three-point sets the full three-deck also reconstructs them directly.

At N=10201, a **single additional third-port bit** separates each mask from
its reflection while all pair averages agree:

| Mask | Third-port offsets | B for mask | B for reflection |
|---|---|---:|---:|
| 115 indices | `(1419,6406)` | 1 | 0 |
| 133 indices | `(179,1600)` | 1 | 0 |

For these examples, enumerating candidates from the six pair differences
finds exactly two translation orbits: the mask and its mirror. Thus that bit
is a sufficient and minimal nonconstant selector for the remaining binary
ambiguity. This does not claim a one-bit selector for arbitrary tournaments.
The masks of 115 and 133 themselves are already separated by their pair
averages; the mirror is the necessary hostile for the third-port claim.

All such averaged correlations forget absolute translation, even at higher
orders. The named source, phase and terminal port must still be retained.
This is the same boundary emphasized by THM-2312
[`sparse-root-bispectrum-positive-word-current`](../../../01-canon/theorems/THM-2312-sparse-root-bispectrum-positive-word-current.md):
phase-balanced cubic data deliberately lose the absolute root address.
THM-2319
[`crt-unit-bispectrum-needle-and-mixed-polarization-no-go`](../../../01-canon/theorems/THM-2319-crt-unit-bispectrum-needle-and-mixed-polarization-no-go.md)
also forbids turning a positive aggregated current into a same-source
incidence claim. P4 uses **unaveraged** products at one shared source; P5
uses translation averages. Confusing those operations would invalidate the
closure conclusion.

## A legitimate tournament view and its limitation

On the odd joint clock orient `t -> u` when the positive native time from t
to u lies in `1,...,(M-1)/2`. Vertices are phase/residue states, the observable
is directed waiting time, and the gauge is forward native clock motion.
There are no ties. This is the rotational tournament, not a Paley tournament
(101 is 1 mod4). A set of three zero events is a cyclic triple exactly when
all three cyclic gaps are at most `(M-1)/2`.

At e=2, the 115 zero events form a cyclic triple and the 133 events form a
transitive triple. At e=4 the roles reverse. This supplies a coarse geometric
test behind the changing long gap; the exact gap lengths, source, and endpoint
remain necessary to recover h(r). Passing to an unlabelled three-vertex
tournament discards those lengths and cannot determine an exact threshold.
Clock reversal gives tournament converse; it is not an unguarded inverse
map on positive integer Collatz trajectories.

The inherited 81/27 control is recomputed: automorphism orders 81 and 27,
67 common labelled cyclic triples, and two-run response ranks 24 and 22
modulo 3. None of the present three-deck claims says that ordinary cyclic
triples reconstruct arbitrary tournaments. Translation-averaged clock decks
and a tournament's labelled cyclic-triple set are different objects.

## Verification and next decisive test

The checker computes all named clocks through e=12, verifies each lift by two
exact arithmetic paths, and literally traverses every phase/residue state
for e=1,2 in all three families. It compares the entire survival function at
all lengths in that finite universe, checks ranks by polynomial gcd modulo
2, tests reflection/translation losses, exhausts the stated small prime-power
mask families, checks the composite-clock hostile, exceptional-prime repairs,
and the 35/175 and 81/27 controls. Formula proofs, not a finite table, give
the all-e statements. The commands are:

Normal and optimized runs agree byte-for-byte on **157,783** explicit gates,
plus the invoked inherited tournament checker. Each phase capacity is also
checked by independently recompiling the rotated word from its own zero
source, including the phase-order changes at e=3. The result SHA256 is
`e8f578f779e9997abed93e0ca97543de186f82d6b7f689bb9edb256b9c66c897`.
The historical recovery pins 28 source documents.

```text
python -X utf8 -B 05-knowledge/results/clock101_joint_ports_20261029/verify.py --output 05-knowledge/results/clock101_joint_ports_20261029/results.json
python -X utf8 -B -O 05-knowledge/results/clock101_joint_ports_20261029/verify.py --output <temporary-file.json>
```

Next decisive boundary: allow transitions between different word blocks while
preserving the integer source and native valuation guards. P4 closes the
algebra for a supplied periodic word; it does not prove that an arbitrary
orbit enters such a block or that a finite bank of blocks covers every source.
The useful candidate is a source-preserving product of these local clock
algebras with explicit change-of-word seams. Its first test must be the
35/175 append pair, not an averaged clock signature.
