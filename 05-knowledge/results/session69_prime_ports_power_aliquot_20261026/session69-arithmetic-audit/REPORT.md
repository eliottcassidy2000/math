# Session69 independent arithmetic-interface audit

**ACCEPT in the scopes stated by the two source reports.** No mathematical
correction is required by this audit. The exact first-six statements, the
power-row normal forms, the complete small aliquot fibers, and both guarded
ABC comparisons survive independent proof review and fresh exact computation.
This is not a proof of the unrestricted tournament atlas, full metagraph
composition, global ABC, or completeness of the prime untouchables.

## Useful findings first

- Six is first in the **odd-root** square/later-doubling predicate. Dropping
  that qualifier makes four first. The predicate has density one half, so its
  first appearance at six does not make it rare in the sense of perfect numbers.
- Six is the first **finite** non-singleton aliquot target fiber, because
  `F(1)` is already infinite. Its complete fiber is `{6,25}`. This is a
  genuinely different property from the square/doubling onset.
- All-powers and repeated-squaring are different atlases. Their address
  multiplicities first disagree at 8; their earliest origins first disagree
  at 27. Neither is an atlas of primes: 15 is already a composite new seed.
- The proposed normalization checks out: every perfect number maps to
  `1+1=2`. To reconstruct the original source/target one must retain the
  removed gcd. The binary-term map from square-later integers first has
  `C>rad(ABC)` at `n=72`, giving `1+8=9` and radical 6.
- The reports correctly keep `{2,5}` as proved untouchable examples but do
  **not** prove their completeness among primes. They correctly keep raw
  Hamiltonian-path counts `H` separate from `tc=H/|Aut|`.

Fresh independent run: **3,173,300 active checks**, normal and optimized
outputs byte-identical. Every original receipt and manifest checked here
matches. Both the original aliquot negative control and our independent
negative control fail in normal and optimized modes, at their intended gates.

## 1. Inputs, inheritance, and review boundary

Reviewed inputs:

- `session69-power/REPORT.md`, `check.py`, both result JSON files, and
  `MANIFEST.json`, plus its pinned operation-fiber and Paley references;
- `session69-aliquot/REPORT.md`, `check.py`, original receipts, negative
  controls, `SOURCE_PINS.json`, and `source/aliquot_fibers.py`;
- read-only repository commit
  `61d7cca5d8df80169f939268269239c2f67ea724`, including the operating contract,
  START-HERE, protocol, guardrails, literature router, and the specific canon
  files named below.

Closest inherited mechanisms are the deleted operation-diagonal forests,
the certified aliquot source bound, and distinct-prime product witnesses.
Canonical hostiles are arbitrary root 8 at `4^2=2*8`, equal Goldbach primes,
and a prime dividing the extension core. The corrected near miss is radical
support overlap at core 4, extension prime 3. The indispensable sidecars are
the original occurrence address, the removed gcd, and the tournament
automorphism denominator.

No source repository was edited. This audit writes only its own directory.
Input hashes, rather than mutable path names alone, identify the reviewed
versions. The source checkout was clean at the pinned commit when checked.
The audit does not independently rerun the inherited all-order tournament
proof's entire finite census; it checks its mathematical scope and separately
recomputes the order-four and Paley-seven controls.

## 2. Proof review: powers and doubling

For `F_(k,r)(n)=2^k*n^r`, direct substitution gives
`F_(k,r) o F_(ell,s)=F_(k+r*ell,rs)`. This proves the proposed normal form for
interleaved operations, including the square-only restriction `r=2^j`.
No commutation of doubling and squaring is presumed.

For a nontrivial odd part, unique factorization gives the unique
`n=2^a*b^e`, with odd `b>1` not a proper perfect power. Every row containing
`n` has the same primitive odd base and therefore numerical base at least
`b`; row `b` itself contains `n`. This proves the earliest-origin statement.
In the square-only carrier, removing the largest common power of two from
the odd exponents instead gives the unique odd nonsquare seed. An odd cube
such as 27 must therefore be distinguished from an odd square.

For two all-power rows with the same primitive odd base, matching odd
exponents forces `r=(f/g)t`, `s=(e/g)t`, where `g=gcd(e,f)` and `t>=1`.
The least admissible dyadic exponent is
`t*max(af/g,ce/g)`. Every greater exponent occurs by doubling. This proves
the complete row-intersection formula, not merely its smallest witness.
Matching square heights yields the reported square-only formula by the same
argument. Pure powers of two are a separate component, including row one.

Address counts retain the base, power exponent, and doubling exponent.
For nontrivial odd part, an allowed exponent `r` gives precisely
`floor(A/r)+1` possible base dyadic valuations; pure-two bases require a
positive valuation and give `floor(A/r)`. Restricting `r` to powers of two
gives the square-only counts. Excluding base one is necessary to avoid
infinite indexed multiplicity for powers of two. No integer-value support
is silently substituted for address multiplicity.

For `n=2^a*u`, `u` odd, the unique odd root below `n^2` is `u^2`. Positive
doubling and later root are therefore exactly
`a>=1` and `u>2^a`. At fixed `a`, the first allowed odd `u` is `2^a+1`;
globally the first source is six. With arbitrary roots, every even `n>=4`
works via `n^2=2*(n^2/2)`, so four is the hostile minimum. The statement
about a later **odd root** is not a statement about a later **new power seed**:
nine is already in row three.

Density scope is correct. Among even sources up to `X`, failure requires
odd `u<=min(2^a,X/2^a)`. Splitting the sum at `a≈(log_2 X)/2` bounds the
number of failures by `O(sqrt(X)+log(X))`; all odd sources fail. Hence the
qualifying set has natural density `1/2`. This does not concern the density
of perfect numbers, twin-prime collisions, ABC exceptions, or prime rows.

## 3. Proof review: complete aliquot fibers and the first collision

For `m>1`, a source with `s(n)=m` is composite. Write `n=pq` with `p` its
least prime factor and `q>=p`. Proper divisor `q` and divisor one imply
`q<=m-1`, hence `n<=q^2<=(m-1)^2`. Equality requires `p=q=m-1` prime.
Except at prime squares, the distinct divisors `1,p,q` give
`p+q<=m-1` and consequently `4n<=(m-1)^2`. Thus the bounded sieve can
certify target fibers globally; it is not merely a bounded source scan.

This proves `F(2)=F(5)=empty` using the reported small cases. In particular,
the stronger non-prime-square bound and prime-square alternative make the
proof for 5 exhaustive. The complete initial fibers establish that six is
the first finite non-singleton target; `F(1)` is the infinite prime fiber.
The corrected transcription is `F(3)={4}`, `F(4)={9}`.

If odd `p,p+2` are prime, the distinct proper divisors of `2p` give
`s(2p)=p+3`, while `s((p+2)^2)=p+3`. The witnesses are distinct. The first
one equals its target iff `2p=p+3`, hence iff `p=3`. This proves uniqueness
of the self-preimage **within this family**, not uniqueness of perfect
numbers or a general description of every collision fiber.

For prime `p`, `s(p^3)=1+p+p^2=Phi_3(p)`. For composite base 4 the omitted
proper divisors matter: `s(64)=63`, not `Phi_3(4)=21`. The polynomial
identity `Phi_3(x^2)=Phi_3(x)Phi_6(x)` remains true, but it cannot transport
an image-complement statement between unrelated invariants.

## 4. Proof review: both ABC maps and every guard

### Source normalization

For `n>1`, set `m=s(n)>0` and `g=gcd(n,m)`. The triple
`(n/g,m/g,(n+m)/g)` is positive and pairwise coprime. Together with `g`,
it reconstructs `(n,m)` exactly. Without `g`, every perfect number maps
to `(1,1,2)`; conversely that primitive triple in this map implies `n=m`.
Equal aliquot targets do not force equal primitive triples: sources 6 and
25 over target 6 give `(1,1,2)` and `(25,6,31)`.

### Square-later normalization

The positive difference `u^2-n=u(u-2^a)` yields
`(u-2^a)+2^a=u`. Since `u` is odd, this triple is primitive. Conversely,
for a positive primitive triple with distinguished term `B=2^a>=2` and
odd `C`, source `n=BC` has exactly the required odd-root property.
The reverse construction retains the distinguished binary term; forgetting
which summand is distinguished would change the carrier. The least source
with `C>rad(ABC)` is 72, namely `1+8=9`, radical 6. This is a bounded
minimality check, not an ABC proof or an ABC counterexample.

### Fixed-core branch normalization

For fixed `a>=2` and prime `p` with `p` not dividing `a`, multiplicativity
of sigma gives `s(ap)=S*p+T`, where `S=s(a)>0`, `T=sigma(a)`. Writing
`d=gcd(S,T)`, `U=S/d`, `V=T/d` and `h=gcd(p,V)`, coprimality of `U,V`
gives exactly `gcd(Sp,T)=d*h`. Thus the reported primitive branch triple
is correct. If `h=1`, its radical contains `p` and
`C=Up+V<=(U+V)p`; if `h=p`, `C=U+V/p<=U+V`. Both imply
`C<=(U+V)rad(ABC)`.

The generic radical identity must divide by `gcd(p,U)`. At `a=4,p=3`,
the triple `9+7=16` has radical 42, not 126. The producer already repaired
this case; its final formula is correct. The cancellation case `h=p`
is not excluded from the theorem. The core-dependent constant is genuinely
unbounded: choosing prime core `a=q` gives `U+V=q+2`. A fixed-core or
finite-bank result therefore does not supply global ABC's uniform constant.

## 5. Primary-literature and tournament scope audit

[Pomerance–Yang, author manuscript dated 14 June 2012](https://math.dartmouth.edu/~carlp/uupaper3.pdf),
page 1, explicitly supplies the distinct-prime Goldbach implication for odd
aliquot targets and treats uniqueness of odd untouchable 5 conjecturally.
The recent crawler date is not the manuscript date. This supports the
conditional route, not an unconditional classification of prime untouchables.

[Booker, arXiv:1610.07471v2, 15 June 2017](https://arxiv.org/pdf/1610.07471),
Hypothesis G and footnote 2, distinguishes ordinary Goldbach from the needed
distinct-prime version. The published computation gives a smaller summand
at most 9781, reducing the equal-prime issue to primes `5<=p<=9781`.
Our independent audit verifies all **1,204** cases in this supplement.
The cited large run through `4*10^18` was not repeated; exclusion of further
odd untouchables through `4*10^18+1` remains a cited consequence plus these
elementary small targets, not a new run by this audit.

The two known prime examples in the retrieved literature are 2 and 5.
For this packet the complete all-prime assertion is **UNPROVED / CONDITIONAL**.
The stronger odd-uniqueness assertion being conjectural does not logically
prove that the weaker prime-only assertion is globally open. This audit does
not infer a new global openness claim merely from failure to locate a theorem.

[Scholze–Stix, Section 1.1, July 2018](https://www.math.uni-bonn.de/people/scholze/WhyABCisStillaConjecture.pdf)
states a uniform epsilon-dependent ABC bound and records a proof dispute.
[Lichtman's May 2025 expository note](https://arxiv.org/abs/2505.13991)
states an almost-all result, not the all-triples conclusion. These support
using ABC as a conditional comparison here, without turning a local
radical estimate or a literature search into a resolution of that dispute.

The current pinned theorem
`THM-1370-h-spectrum-omits-7-21-all-n.md` concerns **raw** directed
Hamiltonian vertex-word counts. Its proof uses strong-component
multiplicativity, insertion, the strong-order floor, and a finite census.
It explicitly records `H=63` as attained, so no `7*3^k` gap tower follows.
`THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md`
keeps completeness of the raw-H spectrum separate. The status in
`THM-4741-presentation-counts-single-view-classes-are-k1-c3-stacks-rotation-inequality-and-forbidden-fibre-sizes-7-and-21.md`
does not promote finite non-rigid `tc` exclusions to unrestricted ones.

The automorphism action on Hamiltonian vertex words is free, so
`tc=H/|Aut|`. Raw `H!=7,21` does not exclude `H=7|Aut|` or `21|Aut|`.
Independent subset-DP gives Paley-seven `H=189`, and the independently
replayed complete automorphism enumeration gives `|Aut|=21`, hence `tc=9`.
Neither 7 (order) nor 21 (symmetry order) is its raw path count.
Independent subset-DP over all 64 labelled order-four tournaments finds
exactly 24 with `H=5`. Alongside `s(8)=7` and `s(18)=21`, this is a decisive
hostile to identifying the two image complements.

## 6. Reproduction and exact finite coverage

`python3 check_independent.py` and `python3 -O check_independent.py` are
self-contained, stdlib-only, and import no producer module. The independent
paths use exact integer-root tests for power origins, explicit address
enumeration, and a smallest-prime-factor/multiplicative sigma sieve for
aliquot values, unlike the producer's ordinary proper-divisor sieve.
`--negative-control` changes the expected canonical first source from six
to four; both interpreter modes must fail at that named gate.

| Check | Independent universe/result |
|---|---|
| Power origins and address counts | every integer through 20,000 |
| Row intersections | 32,580 pairs; bases through 180; value cap 2,000,000 |
| Square-later criterion | every integer through 1,000,000 |
| Primitive binary-term bijection | every represented source through 100,000 |
| Complete aliquot fibers | targets 2..1000; certified source cap 998,001 |
| Primitive aliquot source map | sources 2..20,000, including perfect iff control |
| Fixed-core branches | cores 2..300; primes through 499; 27,827 cases |
| Branch edge cases | 242 radical-support overlaps; 536 cancellations |
| Twin-prime family | 24 pairs with smaller prime below 500 |
| Booker finite supplement | all 1,204 prime cases 5..9781 |
| Native tournament paths | subset-DP Paley-seven; all 64 labelled order-four tournaments |

The complete fresh fiber digest is separately compared with the inherited
inverse algorithm for **every** target 2..1000; this comparison is recorded
in `verification.json`, not inferred from the displayed small table.

Square-later counts `(X,qualifying,even failures)` are
`(100,37,13)`, `(1000,454,46)`, `(10000,4859,141)`,
`(100000,49550,450)`, `(1000000,498513,1487)`. They agree with the proof's
density boundary but do not replace that proof.

## Disposition

**Completed / high certainty:** stated elementary proofs, explicit finite
controls, both always-active negative controls, original normal/-O receipts,
source hashes, and literature-quantifier review. No correction requested.

**Not established:** global prime-untouchable completeness, uniform ABC,
unrestricted tournament-prime classification, full metagraph composition,
or a target-preserving correspondence between aliquot gaps and H/tc gaps.
The required next mathematical action is a map preserving the actual
target predicates and composition data, not another equality of small labels.

**Operational blocker:** none for this audit. No shared repository edit or
push was attempted; this standalone audit is returned to the parent for
integration under its existing session authorization.
