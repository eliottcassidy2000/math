# Independent binary-carrier audit

Status: PASS for the claims and scopes below. This audit does not prove Chan–Pak or all-target F2 coverage.

Source snapshot: `session219-stage` at acbad6c8636cf8cbb675b4c947fa19e33ac0df23; HYP-9478 setting. Worker draft reviewed: `session220-join-bank/PROOF.md` and `session220-three-core/verify.py`. No worker files modified.

## 1. Carrier identity and envelope

For f(0)=0, f(2n)=f(n), f(2n+1)=f(n)+n+(n mod 2), and the HYP-9478 g, induction gives

    f(v)+g(v)=v+1−(v mod 2).

In particular 1<=g(v)<=v+1. An independent bit-pair formula is

    f(v)=sum_{i<j, bit_i=bit_j=1} [2^(j-i-1)+1_{j=i+1}].

For K,h>0,S>=0, every maximizer of Kp−hf(p)−S(p mod 2) on 0<=p<2^m belongs to {0} union {2^m−2^j:0<=j<m}. Here is a direct proof of the rounding step. Write an internal-zero word as p=(2U)2^k+V with U>=1 and 0<V<2^k. Concatenation yields

    f(p)=f(U)+f(V)+2U rev_k(V)>=f(U)+2U.

Rounding to p'=(2U+1)2^k gives f(p')=f(U)+U+(U mod 2)<=f(p), while p'>p and p' is even. Thus the objective strictly increases. Repetition removes all internal zeroes before a later one. If a nonzero candidate has fewer than m bits, multiplying by 2 increases its linear term, preserves f, and does not increase its parity penalty. Hence the remaining word consists of top-aligned ones followed by zeroes.

Sharp hypotheses/control: S>=0 matters. For m=3,K=1,h=10,S=−100 the unique maximizer is p=1, outside the candidate list. The m=0 domain consists only of 0. Extremal values never imply all intermediate odd values occur: the worker's explicit padding support starts 517,561,... and misses 519 inside its hull.

## 2. Full graph nonregularity, including a finite unlabelled bank

Take any fixed finite nonempty integer bank with a,c>0, b>=0, and F_y(v)=a g(v)−b e(v)+c(v+1). For v_m=(10)^m=2(4^m−1)/3, direct recurrence gives

    9F_y(v_m)=(4a+6c)4^m+6am+5a+3c−9b.

Suppose the exact relation containing all (v,F_y(v)) were synchronous regular in standard positional binary. Restrict its input to alternating words. On a sufficiently long accepted path, two equal states at even input positions beyond the finite leading padding give an even-length pump. Pumped input words remain alternating with exponent m+t*j, t>0. The pumped output has numerical value C*4^(t*j)+D, with rational constants C,D. For any fixed bank row its required value is A_y*4^(t*j)+(2a_y*t/3)*j+E_y. Equality holds for only finitely many j: a nonzero exponential plus nonzero affine term has at most two real zeroes, or is a nonconstant affine term. A finite union of rows cannot supply all pumping repetitions. This proves nonregularity also when the label y is omitted.

Encoding scope: MSB-first strings with usual leading-zero equalization, or explicit finite padding symbols with the regular canonical-format restriction. The loop can be chosen wholly in actual input digits. LSBF with the corresponding end padding follows by reversing all tracks. Arbitrary common leading zeroes can first be restricted to minimal padding. Finite labels do not help. This says nothing about regularity of the projection onto outputs alone.

## 3. Stronger selector theorem reviewed

The join-bank worker proves every regular partial selector of valid pairs has bounded input runs, hence only polylogarithmically many distinct outputs <=X. The proof in the reviewed draft is valid:

- For any mixed digit period w, the g increment on the parity-stabilized affine orbit has a strictly positive linear drift D.
- At the negative virtual fixed prefix, every cyclic zero run of length l contributes more than 2^(−l), proving positivity with no asymptotic assumption.
- A fixed suffix preserves the coefficient of g and hence preserves D.
- A mixed loop on an accepting path contradicts the pumped-output exponential-plus-constant form.
- Every trimmed SCC must have only one nonempty input digit type; the condensation DAG then bounds the total run count. Epsilon edges, if allowed, add no runs.
- Since F_y(v)>=c(v+1)+a−b, accepted positive outputs <=X have v=O(X); bounded-run inputs of that size are polynomially many in log X.

This also extends to an unlabelled finite bank by the finite-root argument above; one cannot merely infer nonregularity of a union from nonregularity of its members.

A direct extension to standard asynchronous rational relations is available, without invoking a general synchronization theorem. On a mixed input loop of length k, v_m grows like 2^(km). The corresponding output loop has some length l and either its output remains constant or grows like 2^(lm). The bound F_y(v)=Theta_Y(v), following from 1<=g(v)<=v+1, excludes a constant output and forces l=k. The same positive-drift argument applies. This requires the standard two-tape finite-state relation on literal positional strings; it does not address richer numerical encodings or extra unbounded coordinates.

Positive control: the all-ones input language has g(v)=1,e(v)=0, so F_y(v)=a+c(v+1) is an affine relation on that language and can be recognized by finite-state arithmetic. Thus the theorem does not ban all infinite regular selectors. If a=0, ordinary affine/parity graphs can be synchronous regular: the positive drift coefficient hypothesis is essential.

## Reproduction

Run `python3 session220-audit/audit.py` from the common workspace. Independent code verifies 16,384 recurrence/bit-pair/carrier identities, 1,375 generic envelope parameter domains (all maximizers checked), 297 alternating formulas, and 640 pumped-word numeric identities. `RESULT.json` records exact counts and the negative-S control. Finite verification supports the code; the arguments above provide all-length proofs.

No claim of literature novelty is made. Composition novelty and any additional worker results not discussed here require separate review.

## 4. Full draft review: three-core maxima and cut minima

Reviewed `session220-three-core/PROOF.md` Sections 3–4. PASS: at fixed suffix the prefix coefficient K is positive by the stated physical inequalities; the positive-prefix constraint is preserved by internal rounding. At fixed prefix the suffix rounding increases R, decreases its reversal, weakly decreases f, and preserves the first bit. Thus any global maximizer must satisfy both one-run restrictions simultaneously. No unjustified exchange of maxima or optimization order is needed.

For minima, every within-word quadratic coefficient arises from −h f, and every cross-word coefficient from −q P rev(R). All are nonpositive. Fixing arbitrary bits only changes constants and linear terms. The displayed doubling identity, undirected capacities, source/sink signs, and additive constant are correct. `audit_cylinders.py` independently enumerates 2,268 arbitrary fixed-bit domains across twelve core states and compares their literal minima with the worker cut evaluator; all pass. This supplements the worker's unconstrained rectangle tests.

A minor boundary defect was found and sent for correction: the low-bit cylinder formula (5) must impose t>=1; for t=0 its S' should be S rather than zero. This does not affect Sections3–4 or the ordinary positive-length cylinder application. Freeze the corrected draft before publication.

## 5. Join-bank final expanded draft

SHA256 reviewed: `7669a26ac38aaa9d406d8a20e27fa9fcecf4148dba9153d37fe817cd36a6297b`.
PASS, including new Sections5–6. Uniform binary normalization is excluded for the precisely stated identity and physical odd a coefficients. Coefficient comparison is justified by linear independence of g from the affine/parity coordinates; this does not exclude piecewise or target-dependent descent. The rational/asynchronous and unlabelled-bank proof explicitly stabilizes parity after the first repetition, allows zero/empty output words, uses growth before comparing equal-base exponential forms, and treats finite row switching correctly. Leading-zero input padding does not invalidate the argument: a mixed word still has positive integer value and hence exponential growth on repetition. The SCC argument bounds actual runs even with epsilon transitions.

### Corrected three-core freeze

Final reviewed proof SHA256: `cb1e5a996378a4606851fecd25a8688cb3ebcded295de916208b27fb5cef1b27`. The t=0 issue is repaired explicitly. Final status: PASS.

Sections(3a)–(3b) also PASS. Writing K2^j+h2^(m−j)=2sqrt(Kh)2^(m/2)cosh((ln2)(j−j*)) proves nearest-integer minimization and the tie at a half-integer. The adjacent-candidate threshold K2^(2j+1)=h2^m agrees. Under 1<=j*<=m−1 the best interior candidate beats j=0 even without its penalty and beats zero. The exact coefficient oscillates with period two in m because j* increases by one when m increases by two. Its stated bounds 2sqrt(Kh) and (3/sqrt2)sqrt(Kh) are correct. These are sharp hull calculations, not an occupancy result.

## 6. Final checkpoint revisions and portability

Latest join-bank proof SHA256 `e028d8b063953f7f0bcf8b48b941b0eff3dd49f8da0a48988dfae3904c05056f`: PASS, including Section7. For an ordinary two-track rational witness relation, precomposing the v decoder inverse and postcomposing the N decoder produces a rational relation on literal binary strings, to which Section6 applies. Fixed block radix conversion is a special case. The streaming-register updates are exactly the defining g recurrence; the use of unbounded integer registers is explicitly outside the finite-state obstruction and gives evaluation only, not inverse coverage.

Latest extrema verifier SHA256 `69130a4542b639c1c733dbfcf999fe1e00fcc8390faa8032831da7398660ccda`: PASS after the capacity-bound repair. The former arbitrary 10^100 augment ceiling has been replaced by sum(cap[source]), which bounds the first edge of every residual source-to-sink path. Taking the minimum along that path now yields its true bottleneck, so BFS saturates an edge and the usual Edmonds–Karp augment-count bound applies independently of capacity magnitudes. An added independent regression checks twelve scaled instances with coefficients multiplied by 10^150, all against exhaustive minima.

The cylinder checker now first resolves `../extrema/verify.py`, the portable published layout, and only falls back to the original raw-workspace sibling when that file is absent. No absolute workspace path is required. Both source verifiers were rerun successfully after the revisions. Earlier hashes in this audit record the lineage of actually reviewed drafts; this section identifies the current freeze.

The extrema worker subsequently added a 400-bit-prefix regression with two free bits, checked against all four assignments. This changes only the regression note/code and is also PASS. Current extrema freeze supersedes the preceding hashes: proof `db8070a5d4a0929058ddf43f0889a94009e360c01ba80930cd4d6cfbe240c94a`; verifier `d5950aa44ed4313a8057c0be0f8ecc289beef71a2265cf0d06659e5825be1f6a`.

Portable-layout smoke test: copied the cylinder checker and the extrema verifier/core states into an isolated temporary directory containing only `audit/` and `extrema/`, with no raw-workspace fallback available. All 2,268 cylinder checks and twelve large-capacity checks passed. Receipt: `PORTABLE_CYLINDERS.json`. Temporary test copies were removed afterward.
