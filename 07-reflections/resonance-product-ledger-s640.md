# Resonance Product-Ledger Atlas S640

HYP-2215 names the parent move: jackknife retained side channels instead of
trusting a raw scalar.  S640 asks how to make that move computable across
domains.  The answer is a local product ledger.

The current cluster is no longer just "these hard problems rhyme."  The rhyme
has a useful engineering form:

```text
scalar quotient + retained side channel + residual obstruction.
```

`H=21` is the cleanest warning.  The scalar `21` is not mysterious by itself.
It becomes impossible as a tournament Hamiltonian-path count only after the
strong-component and odd-cycle/OCF side conditions are retained.  Unit distance
at `n=21` is the mirror warning: the same integer is perfectly legal as a
vertex count, and the exact row has `57` unit edges, because the geometry keeps
a `20`-edge spine and a `37=C_hex(3)` bulk.  Collapse those side channels and
you hallucinate a forbidden tournament value where none exists.

LRC `n=14` is the live version.  The scalar clock is

```text
C = 2n - 1 = 27 = 3^3.
```

But THM-407 says the real object is the folded shell packet: gcd strata
`1,3,9`, unit/nonunit shells, lift/carry data, and endpoint owners.  The
least-positive `Res_27` quotient is now classified, but full LRC still asks
whether lift/CRT branch data can introduce a hidden floor row.  Again: scalar
done, side channel alive.

A000568 is the counting-sequence version.  The raw next term is too thin.
S633 made the adjacent shadows visible: self-converse fixed layer, merged
layer, nonfixed layer, q-deformation, bisection, transporter quotient.  This
looks surprisingly like a sieve.  Burnside cycle types are local data; fixed
and merged layers are residue survivors; q-deformation is pressure.

The pi/e-Schanuel thread is the field-of-definition version.  `S=e+pi` and
`P=e*pi` are trace and norm; `D=e-pi` is the discriminant sheet.  HYP-2212
proves at least two of `S,P,D` are transcendental because any two reconstruct
the hidden source.  S638 goes further: if one algebraic shadow existed, the
transverse shadows would become lonely.  This is the same story with algebraic
independence replacing geometric embedding or residue survival.

Twin primes and Goldbach are the analytic-number-theory version.  The seductive
scalars are `gap=2` and `N=p+q`.  The real side channel is local admissibility,
singular series, sieve weights, major/minor arcs, and the parity barrier.  The
current outside status fits the atlas: twin primes remain open, bounded gaps
are known; weak/ternary Goldbach is proved, strong/binary Goldbach remains
open.  In both cases, local residue survival is not decoration.  It is the
object the scalar depends on.

So the new application is to make a repo-native singular-series object.

Not an analytic theorem at first.  A data type:

```text
domain
local prime/modulus
forbidden residues
surviving residues
local weight
lost side channel
```

For twin primes, this is classical admissibility.  For Goldbach, it is parity
and local congruence survival.  For LRC, it is shell survival under `C=2n-1`.
For unit distance, it is direction-support admissibility before embedding.  For
A000568, it is Burnside cycle-type survival under fixed/merged/transporter
actions.

That would let us compare hard residuals by the same kind of product:

```text
global viability ~= product of local side-channel survivals,
```

with the warning that the product is a ledger, not automatically a proof.

This is where `n=14` LRC and twin primes feel genuinely close.  Both are
blocked by a parity/prime-local phenomenon.  Twin primes run into the sieve
parity problem: ordinary sieves struggle to distinguish primes from products
of two primes.  LRC `n=14` runs into a lift/carry side channel after the
least-positive shell quotient is classified.  In both cases, the local product
gets you very far, then a transverse branch carries the real obstruction.

Goldbach adds the positive version.  The circle method succeeds when the major
arc product dominates the minor arcs.  That suggests a question for LRC: can
the safe-box measure be split into "major shell" and "minor shell" pieces, with
the prime-power residual isolated as a minor-arc problem?  HYP-2164 already
does this in finite form: most `Res_27` rows have strict pair-sum certificates;
the floor rows are a tiny structured family.

Unit distance adds a geometric version of the same product.  A candidate ear
on a `21`-core should first pass local direction/support masks, endpoint
compatibility, and gain budget.  That is an admissible tuple test before the
expensive embedding search.  If we build that ledger, `n=22` might become less
of a raw beam and more of a local-to-global problem.

The most reusable invention from this session is therefore:

```text
local_obstruction_product_ledger
```

It should accept rows from totally different domains and return the same
questions:

1. Which local residues/channels are forbidden?
2. Which survive?
3. What scalar would be unsafe if these channels were forgotten?
4. What transverse branch remains after the product says "probably viable"?

That last question is the important one.  The atlas does not say all these
problems are the same.  It says their false simplifications are the same: each
has a scalar that pretends to be the object.  The cure is to keep the side
channel until the predicate has genuinely been preserved.
