# Carmichael inheritance: exact completion, sparse carry barriers, and character loss

2026-10-19. Repository inheritance checked at `75473bbbee6` after the incoming
Apollonius session was synchronized. No Git or checkout changes were made by
this worker. All repository paths below are relative to the repository root.

**Status.** The elementary propositions in §§2–6 are **PROVED here**, with no
historical-novelty claim. Their computational controls are **FINITE-EXACT**.
The paper's exponent-one conclusion is **conditional in this note on its
analytic inputs (S),(A)**; its author-reported formal verification is not an
independent build. The incoming carry heuristic is neither proved nor silently
adopted. This note identifies two rigorous barriers to its proposed method and
repairs its counting interface. Collatz convergence and general integer-cycle
exclusion remain open.

## 1. Inheritance ledger and the incoming delta

### Closest proved mechanism

`05-knowledge/results/carmichael_successor_blocks_20261014b.md` already contains:

* the exact interface `(N,lambda,e)`, `lambda=lcm_{p|N}(p-1)`, `N-1=e lambda`;
* the complete one-prime criterion `q=1+d lambda`, `d|e`, with a distinct-prime
  guard, and the update `(Nq,d lambda,N+e/d)`;
* a complete two-prime classifier, not merely a search under a prime cutoff;
* the finite-bank product-one compiler with its CRT conductor guard;
* marked polygonal and digit-sum versions of every local Korselt condition;
* the exact distinction between constructing a new product and covering every
  prescribed source in another dynamical problem.

The two-prime enumeration gives exactly 1, 13, 21 extensions of 561, 1105,
1729, respectively. For 561, the sole pair is `(401,3041)`:

`561*401*3041=684106401`, with `lambda=15200`, `e=45007`.

Neither intermediate product is Carmichael. The final `e` is prime and
`lambda+1=15201` is composite, so the one-prime obstruction repeats at the
child. This is a successful simultaneous repair, not an infinite extension
theorem.

### Canonical hostiles

1. `91=7*13`: product one modulo 6 and modulo 2 does not imply product one
   modulo 12. In the prior compiler, dropping `gcd(k,L)=1` at `k=2,L=6`
   makes the assertion false.
2. `561` is Carmichael but not Giuga; `30` is Giuga but not Carmichael.
3. `7 -> 11 -> 17` returns modulo 5 and grows. Even the full modular affine
   transformation being the identity does not establish integer descent.
4. A repeated word has fewer rotations than its length. At shape `(4,2)`,
   the two zero-carry words `0101,1010` form one necklace, not `2/4` necklaces.

### Corrected near misses

* **MISTAKE-610 / THM-4651**
  (`01-canon/theorems/THM-4651-giuga-sum-at-p-squared-is-the-wilson-quotient-and-the-korselt-giuga-parity-divide.md`):
  from a primary-pseudoperfect source `N`, the prime-neighbor operation
  `N(N-1)` produces a Giuga number, while `N(N+1)` remains primary
  pseudoperfect. They are different signs.
* `05-knowledge/results/collatz_affine_guarded_lifts_20261004.md` repairs a
  composite-modulus CRT error: multiplier groups at different prime factors
  share one pair of exponents. At modulus 95 the multiplier group has order
  36, not `4*18=72`.
* **THM-4520**, `01-canon/theorems/THM-4520-collatz-level-operator-is-a-gauss-twisted-circulant-with-spectrum-on-the-half-circle.md`,
  expressly does not infer a spectral gap for the changing-level cocycle from
  the radius near 1/2 of each separate nonnormal operator. Its largest singular
  value is 1. This warning is load-bearing for the incoming analogy.

### Underused sidecars

* the exact derived subgroup of the finite affine group;
* the number of preimages of each carry residue, not only its support;
* primitive period and rotation stabilizer;
* the local factor `(prime, exponent, missing predecessor period)`;
* the global conductor, as distinct from separate congruences;
* the number of labels versus the number of target residues.

The concept board is **prime labels / common conductor / local defect index /
affine translation / integer occupancy / primitive necklaces**. The anchor is
the new paper's finite transfer, the niche is arbitrary-core completion, and
the wildcard is the sparse-reservoir obstruction obtained from an LRC lemma.

### What the incoming session changes

The new `HYP-9338-enumerative-geometry-apollonius-ford-descartes-and-carmichael-transfer-dictionary.md`
already records the finite product-one connection and the conditional
`(3/2)`-pseudoprime corollary. Do not publish these as new discoveries again.
Its assertion that the old compiler is “exactly Proposition 5.2” should be
read as the same mechanism: the old `p=dk+1,d|L` compiler has an auxiliary
coprime `k`; the new paper uses `p-1|L` directly. Their precise hypotheses
are not identical.

Its suggested Collatz-carry transfer is the new research target below.
The prior affine-group theorem and THM-2048 supply decisive information
that the new dictionary did not retain.

## 2. A complete one-prime compiler for an arbitrary core

The inherited one-prime theorem assumes that the source is already Carmichael.
That assumption is unnecessary for a finite complete classifier.

**Proposition 2.1.** Let `N>1` be squarefree and put
`lambda_N=lcm_{p|N}(p-1)`. A prime `q` not dividing `N` makes `Nq` Carmichael
if and only if

`q-1 | N-1` and `Nq=1 mod lambda_N`.                         (2.1)

Equivalently, if `gcd(N,lambda_N)=1`, enumerate the divisors `d|N-1`, put
`q=d+1`, and retain precisely the new primes satisfying
`q=N^{-1} mod lambda_N`. If that gcd is not 1, there are no children.

**Proof.** At an old prime, Korselt is exactly
`lambda_N | Nq-1`. At the new prime,
`Nq-1=N(q-1)+(N-1)`, so its condition is `q-1|N-1`.
Distinctness supplies squarefreeness and `Nq` is composite. Conversely these
are all the prime-divisor conditions. The old congruence is solvable only
when the stated gcd is 1. Also `q<=N`, so the divisor enumeration is complete.

Examples, including sources that initially fail:

* `N=33`, `lambda_N=10`: the only child is `q=17`, producing 561.
* `N=91`, `lambda_N=12`: the only children are `q=19,31`, producing 1729
  and 2821. The earlier hostile is a repairable source.
* `N=21`, `lambda_N=6`: the gcd obstruction excludes every child and,
  in fact, every Carmichael multiple.

**Hereditary obstruction.** If `N` divides a Carmichael number, then `N`
is squarefree and `gcd(N,lambda_N)=1`. Indeed `lambda_N` divides the larger
number minus one. This is only a necessary condition here.

The source prime support can therefore be viewed in the directed divisibility
graph `p -> r iff p | r-1`. Carmichael prime supports contain no edge of this
graph. This is an independent-set condition, not an antichain assertion about
the transitive closure. It is not sufficient: `{7,13}` contains no edge but
91 is not Carmichael. The smooth-predecessor bank in the paper makes the
condition automatic: every selected prime is larger than every prime dividing
any selected predecessor. The group identity then pays the remaining global
period obligation.

**Next hypothesis to advance a prescribed source.** A large bank must contain
enough subsets with product `N^{-1}` in a conductor containing both old and
new predecessor periods, with `N` a unit and that target in the actual
reachable subgroup/coset. An identity-product theorem alone does not supply
this target-membership assertion.

## 3. From Korselt defects to a finite-unit syndrome

This is an elementary quantitative extension of the marked Korselt/Giuga
interface in THM-4651, not a replacement for either condition.

For `N>=2`, factor `N=product p^{e_p}` and define

`d_p=(p-1)/gcd(N-1,p-1)`,
`I(N)=product_{p|N} p^{e_p-1} d_p`.

**Proposition 3.1.** The homomorphism

`Phi_N:(Z/N)^* -> (Z/N)^*, a -> a^{N-1}`

has kernel of size

`L(N)=product_{p|N} gcd(N-1,p-1)`,

and image of size `I(N)`. Thus the fraction of unit bases passing the Fermat
test is exactly `1/I(N)`. For composite `N`, `I(N)=1` if and only if `N` is
Carmichael.

**Proof.** Modulo `p`, cyclicity of `F_p^*` gives exactly
`gcd(N-1,p-1)` roots of `X^{N-1}-1`. At every such root the derivative is
a unit modulo `p`, because `p|N` and hence `p` does not divide `N-1`.
The elementary one-digit lifting equation therefore has exactly one solution
at every higher power of `p`. This also handles `p=2`. CRT multiplies the
local root counts. Divide the unit-group order by the kernel count, and use
the first isomorphism theorem. Index one forces every `e_p=1,d_p=1`, which
is Korselt's criterion, including the composite guard.

This identifies two different failure channels:

* `p^{e_p-1}` is the repeated-prime cost: every old root has a unique lift
  while the ambient unit group grows;
* `d_p` is the missing predecessor-period cost on the reduced field.

The scalar index alone loses their cause: 9 and 341 both have index 3, but 9
has a repeated-prime obstruction and 341 is squarefree with a period defect
at 31. Keep the prime-marked vector `(p,e_p,d_p)`.

Exact controls: `I(30)=8`, `I(91)=2`, `I(341)=3`,
`I(561)=I(1105)=1`. Although 341 passes base 2, only 100 of its 300 unit
bases pass. A rational base `u/v` coprime to `N` merely tests one element
`u*v^{-1}` of this kernel. This explains exactly why the paper's Carmichael
construction implies its rational-base pseudoprime corollary, while a large
set of rational-base pseudoprimes does not reverse that implication.

**Connection contract.** Source: the marked local conditions in THM-4651.
Target: the kernel/image of a group homomorphism, or a multiplicative syndrome
test. Map: `a -> a^{N-1}`, with CRT at every marked prime power. Preserved:
the all-unit Fermat predicate and its exact density. Lost by a scalar score:
which field period or repeated-prime direction failed. Required sidecar:
`(p,e_p,d_p)`. Cheapest test: 9 versus 341, then the 30/561 Giuga pair.
The exact next hypothesis for a reduced base test is that its chosen bases
generate the full unit group; passing an unproved generating set is not enough.

## 4. The affine-group bridge: why every scalar group character misses the carry

Fix `q>p>=1`, a binary word `w=(s_0,...,s_{q-1})` with `p` ones, and the
Terras branches

`f_0(x)=x/2`, `f_1(x)=(3x+1)/2`.

Write

`F_w=f_{s_{q-1}} ... f_{s_0}=(3^p x+B_w)/2^q`,
`D=2^q-3^p`, `m=|D|`.

The case `m=1` is the inherited unit-denominator family; suppose `m>1`.
Then `gcd(m,6)=1`.

**Inherited group theorem.** The exact result in
`05-knowledge/results/collatz_affine_guarded_lifts_20261004.md`, §2, is

`G_m=<f_0,f_1>=(Z/m,+) semidirect <2,3>`,
`[G_m,G_m]=(Z/m,+)`.

For the present forward-generator convention,
`[f_0,f_1](x)=x-1/4`, which generates all translations modulo `m`.
The quotient by translations is abelian, proving the derived-subgroup claim.

**Corollary 4.1 (scalar-character blindness).** Every one-dimensional complex
character of `G_m` takes the value 1 on every word of shape `(q,p)`, whether
or not `m|B_w`.

**Proof.** All such words have multiplier `3^p 2^{-q}=1 mod m`, so they
are translations `x -> x+2^{-q}B_w`. All one-dimensional characters kill
the commutator subgroup. The identity condition is exactly `m|B_w`, but
the entire translation subgroup has the same scalar-character value.

Hostile: for `(q,p)=(3,1)`, `m=5`, the three carries are `1,2,4`.
None is zero modulo 5, yet every one-dimensional group character evaluates
to 1 on all three words. This is complete invisibility, not merely a weak
estimate or insufficient numerical precision.

**Why this changes HYP-9338.** The paper works in an abelian multiplicative
group, whose scalar characters separate elements. The natural carry group
is metabelian; its abelianization deletes exactly the target-bearing
translation. One cannot import the paper's scalar-character energy to this
group without changing the representation. In particular THM-4520's radius
near 1/2 on different, 3-adic level spaces does not fill this hole.

### A faithful exact replacement

The carry is an additive residue. Define

`C_{j,k}(a;m)=sum_{|w|=j, #1(w)=k} exp(2 pi i a B_w/m)`.

All intermediate terms use the same fixed modulus `m`. Appending a final
letter gives `B_{w0}=B_w`, `B_{w1}=3B_w+2^j`, hence

`C_{j+1,k}(a;m)=C_{j,k}(a;m)`
`                  + exp(2 pi i a 2^j/m) C_{j,k-1}(3a;m)`.       (4.1)

Start with `C_{0,0}(a;m)=1` and zero outside `0<=k<=j`. The number of
zero-carry words is exactly

`Z(q,p)=m^{-1} sum_{a mod m} C_{q,p}(a;m)`.                       (4.2)

This is a finite fixed-composition transfer system. It retains ordered carry,
the number of ones, and frequency multiplication by 3. It is not a product
of independent scalar factors.

**Connection contract.** Source: the derived-group theorem, THM-4520's
non-normal operator boundary, and THM-4717's carry formula. Target: (4.1),
an additive-frequency transfer on the fixed clock modulus. Map: a word to
its affine translation, then additive Fourier transform. Preserved: exact
integer-cycle divisibility. Destroyed by abelianization: all translations.
Needed sidecars: ordered word, fixed composition, clock modulus, and actual
native parity/positivity. Cheapest test: `(q,p)=(3,1)`.
The exact next analytic obligation is a bound on the nonzero-frequency sum
in (4.2) at these growing clock moduli, with the required cardinality slice;
no currently inherited theorem supplies it.

## 5. LRC's quantized occupancy becomes a sparse-carry Fourier obstruction

The proof mechanism of **THM-2048**,
`01-canon/theorems/THM-2048-fiber-quantization-tax-for-lrc-discrepancy.md`,
is that integer-valued fiber occupancy cannot have arbitrary small variance.
Here it has an exact finite-group consumer.

**Proposition 5.1 (quantized Fourier mass).** Map `W>0` labelled objects into
a finite abelian group of order `m`, with fiber sizes `c_g`. Write
`W=a m+r`, `0<=r<m`, and set `C(chi)=sum_g c_g chi(g)`. Then

`sum_{chi!=1} |C(chi)|^2 >= r(m-r)`.                           (5.1)

Equality holds if and only if every `c_g` belongs to `{a,a+1}`.

**Proof.** The integer inequality `(c-a)(c-a-1)>=0`, summed over all
fibers, gives `sum_g c_g^2 >= m a^2+(2a+1)r`. Parseval and removal of the
principal coefficient give

`sum_{chi!=1}|C(chi)|^2=m sum_g c_g^2-W^2>=r(m-r)`.

Equality in the integer inequality is exactly the stated two-level occupancy.
This is the same two-neighbor variance mechanism as THM-2048; its LRC
containment hypotheses are not being imported into the carry problem.

When `W<m`, there are at most `W` occupied fibers. For the distribution of a
uniformly chosen label and Haar-uniform `U`, this gives

`TV(P,U)>=1-W/m`,
`max_{chi!=1}|C(chi)| >= sqrt(W(m-W)/(m-1))`.                    (5.2)

The first inequality follows just by testing the complement of the support;
the second follows from (5.1).

Apply this to the `W=binom(q,p)` words and `m=|2^q-3^p|`. In any regime
where `W/m -> 0`, the full carry distribution is maximally far from uniform
in total variation, and some unnormalized nonprincipal Fourier coefficient
has size at least `(1-o(1))*sqrt(W)`. Thus the “few expected hits” regime
cannot simultaneously support the strong full equidistribution suggested by
an unqualified transplant of the Carmichael proof.

This does not prohibit a random-model prediction for one specified residue.
It says a deterministic global equidistribution theorem is the wrong target.
The Carmichael proof's selected-subset reservoir has vastly more labelled
objects than group elements. The carry heuristic's sparse reservoir has the
opposite size ordering.

For scale illustration only, exact integer ratios give `W/m` about
0.02483 at `(50,32)`, 0.02778 at `(100,63)`, 0.0003385 at `(200,126)`,
and `3.4741e-9` at `(500,315)`. These examples do not establish an asymptotic
for every near-critical shape; no lower bound for every clock gap is assumed.

**Connection contract.** Source: THM-2048's integer occupancy variance lemma.
Target: exact Fourier mass and total-variation obstruction for carry fibers.
Map: `w -> B_w mod m`; preserved quantity: integral fiber multiplicity and
its variance. Lost by support-only summaries: collisions and second moments.
Sidecar: the full occupancy vector, or at least `W,m,sum c_g^2`. Cheapest
test: distinct labels in fewer than `m` residue classes. Next hypothesis:
prove a target-specific zero-fiber statement or define a genuine ensemble
over shapes; do not seek uniform small error in all residue fibers here.

### An integrality barrier to the conjecture's wording

An integer count `A_n` cannot satisfy `A_n=lambda_n(1+o(1))` with
`lambda_n>0` and `lambda_n->0`. Eventually a relative approximation would
put `A_n` strictly between 0 and 1. The same applies to the sum over a
specified shell of shapes if the total proposed main term tends to zero.

Accordingly HYP-9338's phrase “on average over shapes” must specify an
ensemble and normalization. A random expectation can tend to zero; a
deterministic relative asymptotic for a vanishing integer total cannot.
An eventual upper bound below 1 would be an actual exclusion theorem, much
stronger than an informal small expected count. None is proved here.

## 6. Correct necklace and native-cycle counting

For a left rotation beginning with letter `s`, the exact carry identity is

`2 B_rot = 3^s B_w+sD`.                                     (6.1)

Hence zero carry is rotation invariant. The identity also proves native
parity of an integral fixed point: if `B_w/D` is integral, successive
rotated fixed points are integral and
`2x_next=3^s x+s`; reducing modulo 2 gives `x=s mod 2`. Thus the formal
word is genuinely followed. Positivity is separate: with at least one 1,
all carries are positive, so the sign is the sign of `D`.

If `u` has length `j` and `k` ones, repeated `d` times, then

`B_{u^d}=B_u * (2^{jd}-3^{kd})/(2^j-3^k)`.

Consequently the repeated word has zero carry at its own clock if and only
if the primitive block does. Let `Z(q,p)` be the number of all labelled
zero-carry words. Burnside's lemma therefore gives the exact number of
zero-carry necklaces:

`N(q,p)=q^{-1} sum_{d|gcd(q,p)} phi(d) Z(q/d,p/d)`.             (6.2)

The number with primitive period exactly `q` is

`N_prim(q,p)=q^{-1} sum_{d|gcd(q,p)} mu(d) Z(q/d,p/d)`.         (6.3)

These are standard necklace identities with the required divisibility
predicate proved stable under repetition; they are not novel necklace
formulas. If `gcd(q,p)=1`, (6.2) reduces to `Z/q`. Otherwise replacing every
orbit size by `q` can be wrong. The `(4,2)` example supplies a minimal
nontrivial repeated-cycle control.

THM-4591's existing exact census already uses Lyndon words and minimal
period. Our small replay corroborates that interface; it does not extend
its period 301,993 universe. Through length 16 the only primitive mixed
shapes are `(2,1),(3,2),(11,7)`; repetitions appear at other shapes.

## 7. Fejer bridge and its precise boundary

**THM-2054**,
`01-canon/theorems/THM-2054-relative-fejer-whole-product-decorrelation.md`,
uses positive Fejer smoothing and a bounded-resonance lift to identify
constant terms on a line and a higher-dimensional torus. The paper uses
the same positive trigonometric polynomial in a different sampling measure:
the ambient prime counting measure.

The explicit map is `p -> xi(p)` on the unit circle. Its first `R-1`
moments are ambient character sums. If `ord(xi)>=R`, none is principal,
and the analytic input bounds the total Fejer mass. Nonnegativity then
bounds how many selected points can sit in a short arc, without estimating
the selected prime set's own character sums. This is the point where a
preserved ambient measure is valuable.

The alias guard is structural in both settings: THM-2054 needs all bounded
scalar relations to lift, whereas the paper needs no power `xi^j`,
`1<=j<R`, to become principal. Shared Fejer syntax supplies neither guard.
For a low-order character, a principal power makes ambient mass large;
that is why the paper has a distinct low-order/pruning case.

The exact next hypothesis for importing this into another repository lane
is an ambient moment bound at the actual retained frequencies and scale,
plus a target-bearing selection map. A fixed-modulus asymptotic or an
unrelated spectral radius does not satisfy that hypothesis.

The separate session worker on squarefree cycles owns the fourth bridge
to THM-4663's Cayley closed-walk traces. Its label-exclusion correction must
be retained: ordinary trace powers allow repeated primes, while Korselt
construction requires squarefree products of distinct labelled primes.

## 8. Independent fixed-divisor-of-factor-count extension

The exact-cardinality transfer worker has a stronger formulation. For a
short independent check, the paper's original Bernoulli proof itself
preserves the exponent-one conclusion after requiring `m|omega(N)` for
any fixed positive integer `m`, still conditional on (S),(A).

Replace `G` by `G x C_m` and label each prime by `(p,1)`. Prune characters
there. Extend a surviving character to `Psi=xi*zeta_m^j` on the ambient
product. If `ord(xi)>=R`, its first `R-1` powers still satisfy the same
ambient bounds: the extra phase has modulus 1. If `ord(xi)<R`, then
`ord(Psi)<=m ord(xi)<mR`, so the low-order energy lower bound merely loses
a factor `m^2`. With `D'=max(1,ceil(log_2(|G|m)))`, the uniform bound
`E>=N/(100 D' m^2 R^2)` suffices. For fixed `m`, `D'` has the same scale
and all exponential-error comparisons in the paper are unchanged.

Identity in the lifted group imposes both product one and cardinality zero
modulo `m`; unique prime factorization gives the claimed factor-count
condition. General nonzero cardinality classes are not automatic: if all
labels in `C_2` are its nonidentity element, every identity product has
even size. The exact-cardinality worker retains the coset obstruction and
the sharper `eta>2 delta+2 sigma` bookkeeping.

## 9. Reproduction, scope, and next moves

Run:

`python3 check_inherited_bridges.py`

`python3 -O check_inherited_bridges.py`

Both produce the identical JSON and pass **424,073 explicit checks**:

* every binary mixed word of length at most 16: 131,038 words;
* every modulus from 5 to 79 coprime to 6, independently generating the
  affine group and checking the exact commutator;
* all integers `2<=N<=1200`, testing every unit base against the local
  Fermat-kernel count;
* all squarefree cores `2<=N<500`, comparing every eligible prime `q<=N`
  against an independent factor-and-Korselt test;
* 144 fixed-composition residue distributions, direct enumeration versus
  the dynamic recurrence;
* arbitrary occupancy vectors in `{0,1,2,3}^m` for `2<=m<=6`, including
  equality conditions in the quantized Fourier bound;
* direct rotation orbits versus Burnside and Moebius counts, with native
  parity replays at every integral fixed point.

The script uses explicit exceptions rather than `assert`, so optimized mode
does not remove any check. No float certifies a theorem. Neither replay
turns the finite universes into an asymptotic result.

Next mathematical work should target one of three missing coordinates:

1. prescribed-seed target membership for a large Korselt bank;
2. an additive-frequency, fixed-composition estimate for one carry fiber;
3. an explicitly randomized or smoothed model in which a vanishing expected
   number of primitive cycles is a well-defined claim.

Another scalar-character test of the affine group, or a stronger claim of
full carry equidistribution in the sparse regime, cannot repair these gaps.
