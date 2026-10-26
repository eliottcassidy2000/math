# Session 66 independent aliquot audit

2026-10-26. **PASS**, with attribution and quantifier guards below. Auditor:
`/root/audit_aliquot_fiber_bounds`. Producer files were read or copied, never
edited. The assigned checkout at `2807ee83497f47e06311463c9e4fc5dc9dcf752e`
remained read-only; this audit does not authorize or perform publication.

## Verified scope and proof audit

Write `s(n)=sigma(n)-n` on positive integers, `s(1)=0`, and `F(m)={n:s(n)=m}`.

1. **PROVED, elementary.** If `m>1`, every member of `F(m)` is composite.
   For its least prime divisor `p`, either `n=p²`, giving `m=p+1`, or the
   three distinct proper divisors `1,p,n/p` give `m>=1+p+n/p>=1+2sqrt(n)`.
   Thus `n<=(m-1)²`, equality iff `n` is a prime square, and otherwise
   `n<=floor((m-1)²/4)`. This includes squares of composite roots: they are
   not exceptional prime squares. For `m=1` the fiber is all primes and
   the finite bound must not be applied. `F(0)={1}`; negative fibers empty.
2. **PROVED, elementary.** For even `n>=6`, the distinct proper divisors
   `1,2,n/2` give `n<=2m-6`. `n=4,m=3` is a required exception; `n=2,m=1`
   is also outside the stated range. For odd `m>=5`, no exception remains.
3. **PROVED, standard parity.** `sigma(n)` is odd iff `n` is a square or
   twice a square. Therefore an even source of an odd target has that form;
   an odd source of an even target is a square. With `n=a²`, `a>1` odd,
   prime `a` gives `m=a+1`. For composite `a`, its least prime `p<=sqrt(a)`
   gives the proper divisor `a²/p>=a^(3/2)`, plus other positive terms.
   Hence `s(a²)>a^(3/2)`, strictly. `a³<=X²` is a safe integer search
   cutoff, although equality is unnecessary. These arguments prove the
   complete even-target prefix census, not merely a bounded-source sample.
4. **PROVED, exact all-parity inversion.** If some prime has exponent one
   in `n`, write `n=ap`, with `a>=2`, `gcd(a,p)=1`. Multiplicativity of
   `sigma` gives `m=s(a)p+sigma(a)`. Since `a` and `1` are distinct proper
   divisors, `a<=m-1`. Enumerating these cores and retaining positive prime
   integral `(m-sigma(a))/s(a)` coprime to `a` is complete for this class.
   All other composite sources are squarefull. The unique representation
   is `n=u²v³` with `v` squarefree; **no coprimality condition on u,v**.
   This follows prime by prime from even exponents or odd exponents >=3.
   Under `B=(m-1)²`, the number of pairs is at most
   `sqrt(B) sum(v>=1,v^(-3/2)) <= 3sqrt(B)=3(m-1)`.
   Combined with at most `m-2` cores, this is an `O(m)` **candidate-count**
   bound, not an `O(m)` bit-runtime theorem. Factorization, primality,
   arithmetic cost, and preprocessing must be charged in a runtime claim.
5. **PROVED, classical identity rederived.** Set `S=s(a)>0,T=sigma(a)`
   for `a>=2`. For distinct primes `p<q` coprime to `a`,
   `s(apq)=S*p*q+T*(p+q+1)`. Completing the product gives
   `(Sp+T)(Sq+T)=S*m+a*T` iff `s(apq)=m`. Divisor enumeration must retain
   primality, distinctness, positivity, integrality, and coprimality guards.
   At `a=1`, `S=0` and this factorization degenerates: it does not solve
   Goldbach. This identity occurs in Booker's prior work (below).

## Density guard: two distinct parities must not be conflated

Let `E(X)` count **even targets** at most X reached from even sources, and
let `O(X)` count even targets reached from odd sources. The preceding proof
gives `O(X)<=pi(X-1)+floor(X^(2/3))=o(X)`. It does **not** prove that the
image of all odd sources has density zero; almost all odd targets are in
that image by the cited almost-all Goldbach theorem.

Let `U(X)` count all positive untouchables at most X. Using that cited
theorem to give `U_odd(X)=o(X)`, the correct reduction is

    U(X)=floor(X/2)-E(X)+o(X).

Consequently `U(X)/X` has a limit iff `E(X)/X` does, with limits summing to
1/2. One can instead count the entire image of even sources: its odd-target
part is `O(sqrt(X))` by the square/twice-square parity test and linear source
bound. This proves the same equivalence, not existence or a value of either
limit. A finite prefix fraction is not an asymptotic density.

For pointwise odd uniqueness the needed sufficient conjecture is: every
even integer >=8 is a sum of **two distinct** primes. If `m-1=p+q`, then
`s(pq)=m` only for `p!=q`; `s(p²)=p+1`, not `2p+1`. Targets 1,3,7 have
witnesses 2,4,8. The implication does not assert a converse or an equivalence
with ordinary Goldbach. `F(5)` is unconditionally empty by the finite bound.

## Primary-source and priority audit

- **Pomerance–Yang**, [On untouchable numbers and related problems](https://math.dartmouth.edu/~carlp/uupaper3.pdf),
  section 5, Proposition 5.2 and the ensuing paragraph (PDF page 8), already
  contains the parity split, prime-square handling, strict composite-root
  bound `a<X^(2/3)`, and near-linear prefix enumeration. These are not new
  Session 66 results. The PDF is dated June 2012; search-engine crawl dates
  must not be treated as publication dates. Introductory bounds in this
  older paper do not supersede later improvements.
- **Booker**, [Finite connected components of the aliquot graph](https://arxiv.org/pdf/1610.07471),
  section 2 (PDF page 4), gives the same two-prime-core hyperbola in variables
  normalized by `g=gcd(s(a),sigma(a))`. Theorem 2.2 proves the faster
  `O_epsilon(m^(1/2+epsilon))` inverse algorithm for even targets; the
  algorithmic theorem does not require its graph paper's Goldbach hypothesis.
  Thus Session 66's transparent all-parity implementation is not the first
  inverse algorithm nor a fastest-even-target claim.
- **Pollack–Pomerance**, [Some problems of Erdős on the sum-of-divisors function](https://math.dartmouth.edu/~carlp/reversal13.pdf),
  introduction, Conjecture 1.4 and section 3, distinguishes existence of the
  model constant from its conjectural equality to the untouchable density.
  The model sum is over **even** a. The quoted approximation near 0.17 is
  not an established density. The introduction records the almost-all
  distinct-prime Goldbach input and the lower-density bound greater than
  0.06. Its reference to Luca–Pomerance improves the older unknown-upper-
  bound status: positive lower density of even touchables gives upper
  density of untouchables strictly below 1/2.
- **Luca–Pomerance**, [The range of the sum-of-proper-divisors function](https://geodesic.mathdoc.fr/articles/10.4064/aa168-2-6/),
  Acta Arithmetica 168 (2015), 187–199, publisher-source abstract confirms
  positive proportion of even values. We did not independently audit that
  paper's proof. `<=1/2` remains a valid weaker bound, not strongest known.

Primary sources were opened on 2026-10-26. The Bristol-hosted Booker mirror
returned HTTP 403; the arXiv original was accessible, resolving that access
issue. No unproved recent search claim is imported. An absence of a newer
result in this bounded search is not a proof of historical exhaustiveness.

## Independent finite evidence

`audit.py` uses an independently authored direct proper-divisor sieve as
its complete-source oracle, a trial-factorization oracle, and a separately
implemented core/squarefull construction. It imports only a frozen copy of
the producer module for the final cross-comparison, not for the reference.

- Full fibers: every target 2..1000, all possible sources 1..998001.
- Bounds/parity: all 919641 composite sources in that range.
- Direct factorization cross-check: n=1..10000.
- Even-prefix endpoints: 0..40,99,100,121,500,999,1000.
- Guarded two-prime families: every a=2..50, m=2..1000 (48951 cases).
- Explicit hostile controls include n4, n121, n32, n27, repeated Goldbach
  primes, missing coprimality, and nonmultiplicative fiber equivalence.
- Both normal and optimized Python runs pass, with identical JSON output;
  checks use explicit exceptions rather than removable assertions.

Reproduce from the containing workspace:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 session66-audit-aliquot/audit.py > session66-audit-aliquot/normal.json
PYTHONDONTWRITEBYTECODE=1 python3 -O session66-audit-aliquot/audit.py > session66-audit-aliquot/optimized.json
cmp session66-audit-aliquot/normal.json session66-audit-aliquot/optimized.json
```

`producer-input-sha256.txt` pins the initially inspected code;
`final-producer-input-sha256.txt` pins the repaired density script and final
reports. `normal.json` includes
the complete fiber digest and all test universes. The final producer report
was subsequently reviewed: parity, historical priority, open-problem and
finite-census wording pass. The new density calculation and root divisor
bridge have separate `DENSITY_AUDIT.md` and `DIVISOR_AUDIT.md` reports.
Passing these pinned artifacts does not approve unseen later revisions.

## Blocker and authority record

Completed: local independent proof/code/source audit described above, final
producer wording, root divisor bridge, and new Chen–Zhao instantiation.
Denied: none. Pending mathematical/verification work: none.
Nonblocking startup limitation: `agents/processor.py --status --peek --limit8`
could not inspect the repository inbox because `.machine-id` is absent.
No registration/configuration write was attempted in the read-only checkout.
Collaboration messages were available and used for producer coordination.
Uncertain: novelty beyond the specific prior-art matches; not claimed.
No commit/push was attempted; parent assigned audit only, read-only checkout,
and no publication. No external mutation or new authorization is required
for the completed work.
