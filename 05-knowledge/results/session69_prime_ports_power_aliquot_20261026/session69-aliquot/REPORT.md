# Session 69: aliquot fibers, prime gaps, and a guarded radical bridge

**Status:** elementary proofs and FINITE-EXACT verification; cited literature distinguished below. No proof of ABC, prime-untouchable uniqueness, a complete tournament spectrum, or a full unrestricted prime atlas is claimed. No literature-priority claim.

## Immediate conclusions

1. The user's association is the inverse fiber `F(m)={n>=1: sigma(n)-n=m}`. The corrected entry is `F(4)={9}`: `s(9)=4`, not 3. Complete fibers begin

| target m | complete F(m) |
|---:|:---|
|0|{1}|
|1|all primes (infinite)|
|2|empty|
|3|{4}|
|4|{9}|
|5|empty|
|6|{6,25}|
|7|{8}|
|8|{10,49}|
|9|{15}|
|10|{14}|
|11|{21}|
|12|{121}|
|13|{27,35}|
|14|{22,169}|
|21|{18,51,91}|

2. **6 is the first finite non-singleton target fiber**, not the first non-singleton fiber without qualification: `F(1)` is already infinite. It is also the least positive perfect number, the first `n` with `s(n)=n`. Here the domain is positive source integers; no value `s(0)` is assigned.
3. The first collisions have a uniform explanation: if `p,p+2` are twin primes with `p` odd, then
   `s(2p)=p+3=s((p+2)^2)`.
   Thus `p=3,5,11,17,...` creates the pairs in targets `6,8,14,20,...`. Within this family `2p=p+3` holds **iff p=3**, explaining the special self-preimage at 6. The statement supplies two distinct witnesses, not a complete fiber in general, and does not assume infinitely many twin primes.
4. **2 and 5 are proved untouchable.** The all-prime assertion that they are the only prime untouchables is not established here. Distinct-prime Goldbach implies it; the present independent verification establishes it through `10^6`. The cited much larger verification is separate from our own run.
5. **7 and 21 are both touchable**: `s(8)=7`, `s(18)=21`. Conversely, 5 is attained as a tournament Hamiltonian-path count. Thus aliquot emptiness and raw-H omission cannot be identified merely by matching small numbers.

## 1. Inheritance, guards, and source integrity

Read-only input: `/workspace/scratch/07d5532f2fce/session68-repo`, exact commit `61d7cca5d8df80169f939268269239c2f67ea724`. `AGENTS.md`, startup packet, START-HERE, relevant frontier, ACTIVE-GUARDRAILS, RESEARCH-PROTOCOL, CORE-PAPERS, relevant corrections and source canon were inspected. The targeted message utility reported `.machine-id not found`; no registration was invented and no source checkout was changed. The sparse working checkout omits `CORE-PAPERS-ABC-IUT.md` and HYP-9352, but both exist at the pinned commit and were recovered read-only through `git show HEAD:path`; their current status updates were inspected. Working-checkout absence is not repository absence.

Closest mechanism: Session66's complete inverse algorithm and guarded identity, in `05-knowledge/results/session66_root27_aliquot_graph_interfaces_20261026/session66-aliquot/REPORT.md` and `scripts/aliquot_fibers.py`. The latter is copied unchanged under `source/`, SHA-256 `90e8584d4df789dcb1b35c04914d94d470d646c4bc5a50e5c7f64b36ad3ff124`.

Canonical hostile: allowing a prime already dividing the core gives the false `s(4)=5`; in fact `s(4)=3`. Corrected near miss: ordinary Goldbach need not supply distinct primes. Least-used sidecar: prime-exponent-one versus squarefull source shape, together with the gcd removed before an ABC comparison.

Concept board: inverse fibers versus forward orbits; prime-extension guards; divisor exponent chains; radical versus valuation; raw H versus normalized tc; source-cutoff versus certified target-cutoff. Anchor: exact fibers and prime-only quantifiers. Niche: explain the first collision by twin-prime branches. Wildcard: an honest ABC connection and its fixed-core limitation.

## 2. Why the fibers really are complete

For `m>1` and `s(n)=m`, `n` is composite. Write `n=pq`, where `p` is its least prime divisor and `q>=p`. The distinct proper divisors `1,q` give `q<=m-1`, hence `n<=q^2<=(m-1)^2`. Equality requires `p=q=m-1` prime. Unless `n` is a prime square, the divisors `1,p,q` are distinct, so `p+q<=m-1` and `4n<=(m-1)^2`.

This proves finiteness of every fiber except `F(1)`, and permits independent complete checking of targets `2..1000` by a divisor sieve through `998001`. The copied Session66 algorithm uses a smaller candidate universe: either `n=ap`, `p` prime and `p` not dividing `a`, whence

`m=s(a)p+sigma(a)`,

or `n=u^2 v^3` with `v` squarefree. The latter representation is unique but **does not require gcd(u,v)=1**: `32=2^2*2^3` is a required hostile. A fixed source scan without the bound cannot certify empty fibers.

For 2 the global bound gives `n<=1`, impossible. For 5 a non-prime-square source would satisfy `n<=4`, yet the only composite there is the prime square 4; a prime-square source would need prime `p=m-1=4`, impossible. This proves the two empty fibers independently of conjectures.

## 3. Prime untouchables: exact quantifier and source critique

Let `G` assert: every even integer at least 8 is a sum of two **distinct** primes. For any odd `m>=9`, `m-1=p+q` then gives `m=s(pq)`. The remaining odd targets `1,3,7` have sources `2,4,8`; target 5 is empty. Hence

`G => 5 is the only odd untouchable => {2,5} are the only prime untouchables.`

The second implication uses the fact that 2 is the only even prime. **Neither implication has been reversed.** In particular, the open status of the stronger odd-target statement by itself does not logically establish the status of the weaker prime-only statement. Primary sources retrieved here establish the conditional route; no unconditional theorem settling the all-prime assertion was located. Safe research status for this packet is **UNPROVED / CONDITIONAL**, with 2 and 5 the only known examples, rather than promoting the user's global assertion to a theorem.

[Pomerance–Yang, introduction, author manuscript dated 14 June 2012](https://math.dartmouth.edu/~carlp/uupaper3.pdf) states the distinctness hypothesis and conjectural odd uniqueness. The search engine's recent crawl/publication wrapper is not the paper's date. [Booker, arXiv:1610.07471v2, 15 June 2017, Hypothesis G and footnote 2](https://arxiv.org/pdf/1610.07471) explains the verification through `4*10^18`: the published Goldbach computation found a smaller prime at most 9781, with the possible equal-prime cases reduced to a finite supplement. Consequently the cited result excludes another odd—and therefore prime—untouchable through `4*10^18+1`. This large run was not repeated here. The present code independently checks all 78,498 primes at most `10^6`: 78,496 receive witnesses, while 2 and 5 are excluded by the proofs above.

Ordinary Goldbach alone does not make the displayed proof work: at equal primes `s(p^2)=p+1`, not `2p+1`. For example `3+3+1=7`, but `s(9)=4`; 7 is instead touched by 8. A failed distinct-Goldbach witness search is not an untouchability certificate, because other source shapes remain.

## 4. The 1–2–4 decomposition: a real bridge and a decisive boundary

For a **prime** `p`, the proper divisors of `p^3` are exactly `1,p,p^2`, so

`s(p^3)=1+p+p^2=Phi_3(p)`.

Thus `7=1+2+4=s(8)` and `13=1+3+9=s(27)`. But

`21=1+4+16=Phi_3(4)` **does not** imply `s(64)=21`:

4 is composite and the omitted divisors `2,8,32` give `s(64)=63`. This is the missing divisor-closure coordinate in the tempting polynomial analogy.

The exact identity `Phi_3(x^2)=Phi_3(x)Phi_6(x)` specializes at `x=2` to `21=7*3`. It is an algebraic identity, not a proof that an invariant with a gap at 7 must have a gap at 21, much less at all `7*3^k`.

## 5. ABC: two precise maps, and what they lose

### 5.1 Every aliquot source supplies an additive triple

For `n>1`, put `m=s(n)` and `g=gcd(n,m)`. Then

`A=n/g`, `B=m/g`, `C=sigma(n)/g`

are positive pairwise-coprime integers with `A+B=C`. This is the proper ABC carrier. The tuple `(A,B,C,g)` retains the source and target; discarding `g` loses scale. Every perfect number maps to the **same** primitive triple `(1,1,2)`, including 6 and 28. Thus the exceptional behavior of 6 in the original aliquot graph is not an ABC exception. Conversely the equal target `s(6)=s(25)=6` produces primitive triples `(1,1,2)` and `(25,6,31)`, so a fiber does not have one canonical primitive ABC value. For the second triple `rad(25*6*31)=930`, much larger than `C=31`; this is not an unusually small radical. The chain coincidence `6^2=4*9` is a multiplicative identity, not an ABC triple. In contrast `8+1=9` is a genuine primitive additive triple with radical 6. Its having `C>rad(ABC)` does not contradict ABC, whose exponent and constant are essential.

The standard conjecture asks that for each `epsilon>0` there be a **uniform** `K_epsilon` with `C<=K_epsilon rad(ABC)^(1+epsilon)` for every primitive positive triple. We use this only as a conjectural comparison, not a proved input. [Scholze–Stix, Section 1.1](https://www.math.uni-bonn.de/people/scholze/WhyABCisStillaConjecture.pdf) supplies this formulation and documents a dispute about a proposed proof; the repository's THM-3833 (`THM-3833-abc-conditional-cube-radical-and-hyperbolic-power-finiteness.md`) likewise keeps ABC conditional. The more recent [Lichtman expository note, arXiv:2505.13991, 20 May 2025](https://arxiv.org/abs/2505.13991) treats ABC as a conjecture and explains an almost-all result; an almost-all theorem is not the required all-triples theorem. No global literature-consensus claim is inferred merely from a search result.

### 5.2 Fixed-core prime-extension families satisfy an elementary radical bound

Fix `a>=2`. Set `S=s(a)>0`, `T=sigma(a)`, `d=gcd(S,T)`, `U=S/d`, `V=T/d`. For prime `p` not dividing `a`, the inverse branch has target `m=Sp+T`. Let `h=gcd(p,V)`, which is either 1 or p. Its **branch** additive triple, distinct from the source triple above, is

`A=Up/h`, `B=V/h`, `C=(Up+V)/h`.

Because `gcd(U,V)=1`, it is primitive and `gcd(Sp,T)=dh`.

**Proposition (elementary, unconditional).** Every such prime branch satisfies

`C <= (U+V) rad(ABC)`.

If `h=1`, `p|A` and hence `rad(ABC)>=p`, while `C=Up+V<=(U+V)p`. If `h=p`, then `C=U+V/p<=U+V` and `rad(ABC)>=1`. This proves the claim, including all exceptional divisor primes. In the noncancelled case the exact radical expression is

`rad(ABC) = p*rad(U V C)/gcd(p,U)`.

The denominator is essential: `a=4,p=3` gives `U=3,V=7`, primitive `9+7=16`; multiplying `rad(UVC)=42` by p would incorrectly give 126 instead of 42. The initial scratch checker caught this overlap before delivery; the repaired expression is checked explicitly.

This is a genuine elementary ABC-shaped result on each fixed core, **not progress proving uniform ABC**: the constant `U+V` depends on the core and is unbounded as the core varies. A finite bank of cores has an elementary uniform bound obtained by taking the largest of these constants, but cannot represent the unrestricted problem. This echoes Session66's separate finite-core density obstruction; the unbounded core is a real missing coordinate in both settings.

## 6. Tournament firewall

`H(T)` counts directed Hamiltonian vertex words. `tc(T)=H(T)/|Aut(T)|` counts automorphism orbits: the action is free because an automorphism fixing an ordered word fixes every vertex.

- `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md` proves raw H never equals 7 or 21. Their being **all** positive odd omissions is still open (`THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md`).
- `THM-4741-presentation-counts-single-view-classes-are-k1-c3-stacks-rotation-inequality-and-forbidden-fibre-sizes-7-and-21.md` and `HYP-9352-presentation-count-laws-strong-spectrum-forbidden-fibre-sizes-and-the-three-move-law-at-n-11.md` retain the unrestricted non-rigid tc exclusions as open. `H!=7,21` does not exclude `H=7|Aut|,21|Aut|`.
- The raw-H gap does not lift to a whole `7*3^k` tower: the inherited theorem records H=63 as attained.
- The new code brute-counts all 64 labelled four-vertex tournaments and finds 24 with H=5, directly defeating an identification of prime untouchability with raw-H gaps.

The useful shared *method* is the guarded image problem: distinguish a certified empty fiber from failure to find a witness, preserve composition data, and do not collapse quotient coordinates. No predicate-preserving map from aliquot fibers to Hamiltonian-path fibers has been constructed. Similar small-number patterns remain an analogy without one.

## 7. Reproducible finite verification

Run `python3 check.py` and `python3 -O check.py`. Both use active exceptions, not removable assertions, and produce identical JSON. The complete-target comparison is between the copied factorization-based Session66 algorithm and an independent ordinary divisor sieve.

The final receipt records: targets `2..1000`, proved common source bound 998001; all prime targets through `10^6`; 8,785 guarded ABC branches (`a=2..200`, primes at most 199), including 341 cancellation cases and 146 support-overlap cases; twin-prime collisions with `p<=997`; exact 1–2–4/cyclotomic boundaries; and all labelled order-four tournaments.

`--negative-control` intentionally asserts `F(6)={6}`. Normal and optimized modes both must fail with the target6-completeness error. Final source and receipt hashes are in `SHA256SUMS`. No results in this packet are a full all-order computation, an independent rerun of the large published Goldbach verification, or a solution to the open claims.
