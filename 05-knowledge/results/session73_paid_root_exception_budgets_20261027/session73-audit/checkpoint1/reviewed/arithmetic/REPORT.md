# Session 73: source-faithful arithmetic grammars and the two gap pairs

2026-10-27. **PROVED scoped deductions / CITED classical input / FINITE-EXACT controls.**
No novelty claim. No proof of prime-untouchable uniqueness, raw-H spectrum completeness,
native ROOT positivity, or a uniform tournament exception budget.

## Immediate findings

1. A finite bank of fixed prime-extension cores, **even supplemented by every
   squarefull source**, reaches relative-density-zero prime targets. Thus this
   particular finite-template route cannot establish that every prime except 2
   and 5 is touchable. Changing observation moduli does not add the missing cores.
2. An aliquot target fiber is maximally unstable under a common legal prime
   context: if `s(a)=s(b)` and prime `p` does not divide `ab`, then
   `s(ap)-s(bp)=a-b`. Every unequal source pair is separated. Target-mass
   pushforward therefore is not a closed state for future composition.
3. The exact positive carrier is `(n,s(n))`. Dropping divisor mass `delta` and
   extending by a coprime squarefree context `P` multiplies the error by
   `sigma(P)`, but the omitted-term count only by `2^omega(P)`. Even two missing
   terms can carry arbitrarily large error. A count budget is not a weighted budget.
4. Raw-H attained values form a multiplicative monoid; aliquot-touchable values
   do not: `4=s(9)` and `13=s(27)`, whereas `F(52)` is empty. This is an actual
   algebraic obstruction, rather than another numerical association of `{7,21}`
   with `{2,5}`.

## 0. Inheritance, exact pin, and what is being connected

Input repository commit: `92e7f12f3b95840cff92c5e6cf40b2ae8952cae4`.
`SOURCE_PINS.json` records exact repository paths and SHA-256 snapshots.

Closest mechanisms: Session66's guarded prime-extension and squarefull-complement
inverse search; Session69's complete-fiber bound and tournament firewall;
`weighted_renewal_moving_modulus_20261026`, R2/R4 and T1/T2, on preserved source
mass, precision versus dynamics, named ports, and paid exception cycles.
The near miss is Session66's noncongruence witness, which used a noncoprime
context; B below strengthens it to **every** shared fresh-prime context.
The least-used sidecar is source `n` inside an equal-target fiber.

Concept board: weighted positive grammar; context observability; source versus
target quotient; finite core versus growing core; exception count versus mass;
raw H versus automorphism-normalized tc. Anchor: a real gap-pair obstruction.
Niche: prime-target coverage of finite cores. Wildcard: exact positive omission debt.

The proposed bridge is a **proof-design transfer**, not an identification of
invariants: retain the source needed by composition, preserve its original mass,
and price omitted terms in their actual weights. No intrinsic tournament relation
on integers is manufactured.

## A. A finite core bank misses almost every prime target

Write `s(n)=sigma(n)-n`, `n>=1`. Let `A` be any **fixed finite** subset of
`{2,3,...}` and define

```
K_A = {s(ap): a in A, p prime, p does not divide a}
      union {s(n): n is squarefull}.
```

Here squarefull means every nonzero prime exponent is at least two.

**A1 (PROVED using the classical fixed-modulus PNT in progressions).**

`#{prime m<=X: m in K_A} / pi(X) -> 0` as `X->infinity`.

The quantifier fixes `A` **before** taking the limit. Nothing here permits
`A=A(X)` to grow, bounds a general recursive grammar, or says that targets
outside `K_A` are untouchable.

### A1a. Squarefull sources contribute only O(X^(2/3)) prime targets

**Elementary explicit bound:** their number is at most
`1 + floor(sqrt(X)) + 3 X^(2/3)` for `X>=2`.

If a squarefull source has total prime exponent `Omega(n)=2`, it is `p^2`,
with `s(n)=p+1`; this is prime only at `p=2`, giving target 3.
If `Omega(n)=3`, it is `p^3`; since `s(p^3)>=p^2`, there are at most
`floor(sqrt(X))` candidates with `s(n)<=X`.
Otherwise `Omega(n)>=4`. For its least prime divisor `p`, `n>=p^4`, while
the proper divisor `n/p` gives

`s(n)>=n/p>=n^(3/4)`.

Hence `n<=X^(4/3)`. Every squarefull integer has a unique expression
`u^2 v^3` with `v` squarefree; no `gcd(u,v)=1` condition is allowed.
The number at most `Y` is bounded by

`sqrt(Y) sum_(v>=1) v^(-3/2) <= 3 sqrt(Y)`.

Taking `Y=X^(4/3)` proves the bound. The prime-cube branch cannot be silently
folded into the last bound: `s(125)=31` but `125^3>31^4`.
This extends Session66's squarefull warning to the **prime-target** scale.

### A1b. Each fixed core contributes relative-density-zero prime targets

Put `S=s(a)>0`, `T=sigma(a)>0`. The coprimality guard gives
`s(ap)=Sp+T`. Fix finitely many primes `ell>2` not dividing `ST`.
For `m=Sp+T` prime, the primitive residue `m=T mod ell` forces `p=ell`,
which gives at most one exceptional target. The Chinese remainder theorem
therefore leaves a proportion

`product_ell (ell-2)/(ell-1)`

of the primitive residue classes modulo the product of these auxiliary primes.
For each **fixed** such modulus the prime number theorem in progressions gives
that same limiting proportion among prime targets. Thus the limsup for this
branch is at most the displayed product. Increasing the fixed auxiliary set
makes the product tend to zero, because reciprocal primes diverge and only
finitely many primes divide `ST`. A finite union over `a` and A1a finish A1.

**CITED classical input:** fixed primitive residue classes satisfy
`pi(X;q,r) ~ pi(X)/phi(q)`. It follows from equation (6), with partial summation,
in [Tao, *254A Notes 2*, 9 December 2014](https://terrytao.wordpress.com/2014/12/09/254a-notes-2-complex-analytic-multiplicative-number-theory/).
This is an iterated fixed-modulus argument, not a uniform-in-growing-q estimate.
No prime-pair conjecture is used. The independent audit confirmed A1a and this
fixed-q proof; its message is preserved in `AUDIT_MESSAGE.md`.

### Consequence and smallest finite controls

Let `c(m)` be the smallest admissible core `a>=2` in a prime-exponent-one
preimage of `m`, or infinity if there is no such preimage. For every fixed `C`,
the primes with `c(m)<=C`, together with all primes having any squarefull
preimage, form a relative-density-zero set. This is a required growth coordinate,
**not** evidence for many untouchable primes. Under distinct Goldbach, such
prime-exponent-one sources exist for every prime target `m>=11`, but their
smallest cores still escape every fixed bound in relative density.

The complete computation through target 300 exhibits touched primes missed by
the indicated bank **and all squarefull sources**: first target 19 for cores
`2..3`, 61 for cores `2..5`, and 149 for cores `2..10`. Their full fibers are
printed in the receipt. This is a finite control of the correct consequence
object, not the asymptotic proof.

## B. Source-fiber separation and the exact positive carrier

For coprime positive `a,b`, multiplication of divisor sums gives

`s(ab)=a s(b)+b s(a)+s(a)s(b)`.

Thus the nonnegative pair operation

`(a,u) odot (b,v)=(ab,av+bu+uv)`

represents coprime integer multiplication exactly on `(n,s(n))`. It becomes
coordinatewise multiplication after the invertible change `(n,s)->(n,n+s)`.
The arithmetic guard is indispensable: at `a=b=2` the formal operation predicts
5, whereas the true value is `s(4)=3`.

For a fresh prime `p`, the update is

```
(n',s')^t = [[p, 0], [1, p+1]] (n,s)^t,
n = s'-(p+1)s.
```

**B1 (PROVED).** If `s(a)=s(b)` and `p` is prime with `p` not dividing `ab`,
then `s(ap)-s(bp)=a-b`. Consequently all distinct sources in a common target
fiber are separated by every common admissible prime extension.
Example: `s(6)=s(25)=6`; at `p=7`, `s(42)=54` and `s(175)=73`.

**B2 (PROVED no-go).** There is no target-only transition function `F_p`
with `s(np)=F_p(s(n))` for all sources not divisible by `p` whenever that
domain contains two distinct equal-target sources. In particular no map to
tournaments can simultaneously (i) identify its Hamilton count with an
injective recoding of `s(n)` and (ii) send all coprime prime extensions to
ordered joins with a context tournament depending only on `p`.
Equal H values stay equal after a common ordered join, contradicting B1.
This obstructs that named operation-preserving bridge; it does not rule out
every possible relationship between the two subjects.

**B3 (PROVED source-mass no-go).** Pushforward of a source measure under `s`
is not a sufficient state for prime-extension response. The two priors
`delta_6` and `delta_25` have the same target law `delta_6`, but after the same
legal multiplication by 7 their target laws are `delta_54` and `delta_73`.
Preserving total mass on each old fiber does not repair the lost source.
Retain the named source and its supplied mass; do not replace these by a
uniform law on the fiber. This is the exact arithmetic analogue of the
pushforward/averaging distinction in weighted-renewal R2, not its ROOT conclusion.

The source `n` alone could of course encode all arithmetic data by recomputation.
No universal information-theoretic two-coordinate lower bound is claimed.
Also the commutative pair state forgets factor order and elapsed operation
history. Preserve any supplied word/clock separately; no Collatz or tournament
native clock is identified with the number of factor extensions.

## C. Weighted positive omissions have a precise debt

Retain the source/top divisor `n`, but omit a set `E` of proper divisors of
total weight `delta=sum_(d in E)d` from the divisor sum. For a squarefree
`P` coprime to `n`, unique divisor factorization transports precisely these
omissions to `{d e: d in E, e|P}`. Therefore

```
new omitted count = |E| 2^omega(P),
new omitted mass  = delta sigma(P).
```

Equivalently, evaluating the positive pair operation at `(n,s(n)-delta)`
instead of `(n,s(n))` changes the final output by exactly `sigma(P) delta`.
The new top-derived terms are retained, so there is no extra hidden omission.

For `E={1}` and one fresh prime `p`, the missing support has just two terms
`{1,p}` but mass `p+1`. As primes grow this defeats every uniform estimate
`mass <= C * number_of_exceptions` with constant independent of the source
weights. A valid repair is to retain the omitted weights or a proved upper
bound for them. This does not by itself refute a tournament edge-count budget:
it refutes **exporting** such a budget as a weighted arithmetic error bound
without the weight sidecar. It likewise supplies no native signed ROOT forcing.

## D. The gap firewall is algebraic, not just numerical

Ordered tournament join satisfies `H(T->U)=H(T)H(U)`, by unique concatenation
of the two Hamilton words. Hence the raw-H image is multiplicatively closed.
The positive aliquot image is not: 4 and 13 are attained, but 52 is not.
The full source bound `n<=(m-1)^2` proves that the divisor sieve through 2601
is already a complete test of `F(52)`; the delivered run uses 89401 and
certifies every fiber from 2 through 300. The only product of two touchable
values at most 52 that is untouchable is `4*13=52` (unordered factors >=2).

Recovered facts remain separate:

- `THM-1370-h-spectrum-omits-7-21-all-n.md` excludes raw H values 7 and 21
  at every order. Its proof is inherited, not independently re-proved here.
- 7 and 21 are aliquot-touchable: `s(8)=7`, `s(18)=21`.
- Of all 64 labelled tournaments of order four, 24 have H=5.
- An explicit order-eight tournament has H=63, verified by both subset DP and
  all 40320 permutations. Its mask is `137938847` in the receipt's convention.
  Thus the gap at 7 does not propagate up the whole `7*3^k` tower.
- `tc=H/|Aut|` remains a different target. A raw-H exclusion cannot be applied
  to tc without its automorphism factor.

Neither `{7,21}` being the entire odd-H gap set nor `{2,5}` being the entire
prime-untouchable set is proved by this packet.

## E. Literature status, with exact logical direction

Let G assert that every even integer at least 8 is the sum of two distinct
primes. Then `s(pq)=1+p+q` gives

`G => 5 is the only odd untouchable => {2,5} are the only prime untouchables`.

No reverse implication is asserted. A missing Goldbach witness is not an
empty aliquot fiber. The equal-prime substitution is invalid because
`s(p^2)=p+1`, not `2p+1`.

[Pomerance–Yang, *On Untouchable Numbers and Related Problems*, introduction](https://math.dartmouth.edu/~carlp/uupaper3.pdf)
is dated **14 June 2012**, irrespective of recent search crawl labels; it gives
the conditional route and distinguishes it from known small untouchables.
[Booker, *Finite connected components of the aliquot graph*, v2, 15 June 2017](https://arxiv.org/pdf/1610.07471)
states G and explains in footnote 2 why the published Goldbach verification
through `4*10^18` can be supplemented to enforce distinctness. These are cited
results, not computations rerun here. They exclude a second odd untouchable
through `4*10^18+1` as imported in Sessions66/69.

Primary sources and a targeted status search were checked on 2026-10-27.
No unconditional all-prime classification was located in that search. We do
not infer that a weaker prime-only assertion has the same status merely
because the stronger odd-target assertion remains conjectural in those sources.
The safe result delivered here is the conditional implication and the new
fixed-bank obstruction, not a claimed complete classification.

## F. Exact tests and independent audit

Run from this directory:

```
python3 check.py > normal.json
python3 -O check.py > optimized.json
cmp normal.json optimized.json
```

The matching receipts contain 13715 coprime carrier checks, 12104 fresh-prime
observability checks, all 5625 unordered collision pairs in complete target
fibers 2..300, 2574 exact local-affine residue tests, exact CRT survivor counts,
20 positive omission transports, the squarefull shape split through 89401,
all order-four tournaments, and the independently counted H=63 witness.
No floating-point arithmetic or removable `assert` is used.

Each of the four deliberately false controls must fail in both modes:

```
python3 check.py --negative-control target-only
python3 check.py --negative-control drop-guard
python3 check.py --negative-control omit-cube
python3 check.py --negative-control count-is-weight
```

`negative-controls.json` records all eight nonzero exits. The receipt labels
its finite universe and does not purport to prove an asymptotic by enumeration.
`SHA256SUMS` covers the final report, program, outputs, source pins, and audit.

The strongest next arithmetic target is a **growing-core coverage theorem**
with explicit primality/coprimality guards and residual-source coverage. A
successful finite-template proof would need an operation generating genuinely
unbounded cores, not only more residue precision on a fixed bank.
