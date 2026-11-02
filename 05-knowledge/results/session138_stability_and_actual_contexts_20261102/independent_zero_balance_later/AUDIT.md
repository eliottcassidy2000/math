# Independent audit: exact cyclic sign chambers of the three Q4 modules

**Verdict: ACCEPT, all stated universal quantifiers.** Final reviewed source:
`zero_balance_later/NO_ZERO_THEOREM.md`, SHA256
`6adadf59d2628d4905e835ac82b83b4723e89f76d7829043f4930c7db227fd6b`.
This audit is separate from the frozen checkpoint-one profile-dichotomy audit.

## Accepted scope

Let M be the actual twelve-vertex substitution of Q4 with masks `(8,1,3,6)`
by four directed triangles, and let M_s be its stated one-vertex completion
at cut s. For `A_a=TT_a[C3]`, all positive integers a,b satisfy:

- s=867 or995: the response in `C3[M,A_a,A_b]` is strictly positive if
  `min(a,b)=1`, and strictly negative if a,b>=2. There is no zero anywhere
  in this unbounded integer-parameter family.
- s=996: for **arbitrary nonempty tournaments** A,B, the response is at least
  `135618 H(A)H(B)>0`. The constant is sharp at A=B=K1.
- For867/995, one side C3 and the other **any** nonempty tournament gives a
  strictly positive response.
- For867/995, the response is strictly negative whenever both actual side
  polynomials have the sufficient form `P_A=p^2 R_A`, `P_B=p^2 R_B`, where
  `p=m^3+2m` and R_A,R_B are nonzero polynomials with nonnegative monomial
  coefficients. In particular, ordinal sides containing at least two C3
  components satisfy this sufficient condition.

These are actual one-vertex replacements retaining the outer incidences,
not formal capacity networks. The theorem closes these three seeds in the
specified contexts; **arbitrary strong-quotient scalar zero rigidity remains
OPEN**. The factorization chamber is sufficient, not necessary, and the
response of867/995 for completely arbitrary side tournaments is not classified.

## 1. The tridiagonal kernel is the exact context, including endpoints

An outer C3 run word is forced once its starting phase and length are fixed.
For r occurrences of M, the side multiplicities are `(r,r)` in three phases,
and the six neighboring pairs listed in the manuscript in one phase each.
These are linear Hamilton-path words, so equal multiplicities really have
three different starts; no cyclic-rotation divisor is appropriate.

The local path components are ordered by their successive appearances in the
word. Thus multiplying the ordered cover counts introduces no missing r!
or additional matching factor. Empty local profiles are zero because all
three modules are nonempty. The singletons and r=1 boundary are included.

Collecting the coefficient of `f_A(i)f_B(j)` therefore gives exactly

```
K_ii = d_(i-1)+3d_i+d_(i+1),
K_(i,i+1) = K_(i+1,i) = d_i+d_(i+1),
```

with all other entries zero and d outside its actual support zero. I rebuilt
all507 entries of the three thirteen-by-thirteen kernels by directly counting
forced run words, rather than importing the producer's kernel function. The
formula also agrees with direct Hamilton-subset counts on all1,331 ordered
triples of labelled tournament modules of orders1..3.

The fixed module profiles are actual: all three twelve-vertex differences
were recomputed using the earlier independent anchored-partition kernel and
match the frozen input. The original Hamilton count is241137. The underlying
same-deletion-card construction and strong original quotient were already
independently checked in checkpoint one; no abstract profile is substituted
for an actual tournament here.

## 2. The likelihood-ratio lemma is valid—and has the right hypotheses

The crucial claim is **not** that every product of positive ordered-cover
polynomials preserves arbitrary likelihood-ratio orders. Its precise input
is a positive, log-concave binomial coefficient sequence v, multiplied by a
polynomial R with nonnegative **monomial** coefficients.

From `m binom(m,r)=r binom(m,r)+(r+1)binom(m,r+1)`, multiplication by m acts as

```
(Lv)_r = r(v_r+v_(r-1)).
```

The manuscript's expansion of the adjacent log-concavity difference for
`v_r+v_(r-1)` is correct. The two ordinary adjacent minors are nonnegative;
the remaining cross-product inequality follows from monotonicity of adjacent
ratios. With zero extension there are no internal holes, and the boundary
cases remain nonnegative. Multiplying termwise by the positive log-concave
sequence r preserves log-concavity. L therefore preserves log-concavity and
extends the positive support from1..N to1..N+1.

On the original support,

```
(Lv)_r/v_r = r*(1+v_(r-1)/v_r)
```

is nondecreasing. Iteration is legitimate: every intermediate L^j v is
log-concave, and restricting each successive ratio to the original common
support gives a product of positive nondecreasing functions. Thus
`(L^e v)_r/v_r` is nondecreasing. Taking the nonnegative sum corresponding to
R proves the stated ratio claim. A nonzero R also gives strictly positive
coefficients at every original support coordinate. The added coordinates
are nonnegative.

Importantly, the proof does **not** require the final sum `R(L)v` to be
log-concave. It only needs each iterate L^e v to be log-concave in order to
prove the ratio comparison to the fixed base. This distinction resolves the
potential multiplication/likelihood-ratio gap.

The log-concavity hypothesis cannot simply be dropped. For the positive but
non-log-concave sequence `v=(2,1,2)`, multiplication by m gives
`Lv=(2,6,9,8)`. On the original support the ratios are `(1,6,9/2)`, which
are not nondecreasing. This exact hostile is retained in the independent
receipt; positivity alone is not being used as a substitute for log-concavity.

For the actual base

```
v=f_(TT2[C3])=(9,126,684,1656,1800,720),
```

all adjacent log-concavity inequalities hold. Since
`P_(A_a)=p^2*p^(a-2)` for a>=2, the multiplier has nonnegative monomial
coefficients and the lemma gives the required comparison to v for every a.
No log-concavity of the later side profile is assumed. The same reasoning
works for the more general multipliers R_A,R_B in the final corollary.

## 3. One crossing is enough, but the tail condition must be retained

The867 kernel has nonnegative rows1..3, one mixed row4, and nonpositive rows5
onward. The995 kernel has nonnegative rows1..4, mixed row5, and nonpositive
rows6 onward. Multiplying either kernel by a nonnegative actual profile
therefore gives a vector whose signs are nonnegative followed by nonpositive;
the single mixed row can lie on either side of the crossing but cannot create
a second crossing. Its first coordinate is positive on the relevant supports.

In particular the last positive coordinate c lies at most at4 or5, hence
inside the base's support1..6, and every coordinate beyond6 is nonpositive.
This last condition is necessary: extra positive tail contributions would
invalidate the comparison. It is genuinely supplied by the recorded row
signs rather than assumed from a graph-size cutoff.

For a negative base pairing, choose `lambda=w_c/v_c>0`. At positive h_r the
ratio w_r/v_r is at most lambda; at negative h_r it is at least lambda.
Consequently the common-support pairing is at most lambda times the negative
base pairing, and the extra tail can only decrease it. Strict negativity
survives. A vector with no positive coordinates is handled directly.

The exact base pairings are independently confirmed:

```
v^T K_867 v = -373209154588134,
v^T K_995 v =  -93517858646046.
```

Apply the comparison once to `h=Kv`. Fix the first enlarged side, use K's
symmetry to identify the new negative base pairing, then apply it to
`h=K w_A` for the other side. The same row-sign pattern applies to any
nonnegative w_A, so this second step does not smuggle in log-concavity of
w_A. This proves the entire negative quadrant and the broader factorization
chamber, not merely the tested parameter box.

For the ordinal-component corollary, colored path-cover polynomials multiply
under ordinal sum. Every remaining actual tournament factor has nonnegative
monomial coefficients by the colored OCF. Two C3 components therefore supply
p^2 and leave a permitted multiplier. There is no necessity claim for this
particular factorization.

## 4. Thin sides, arbitrary-side positivity, and failure boundaries

With one C3 side the profile is `(3,6,6)` and vanishes after coordinate3.
Only d_1 through d_4 can contribute. All four are strictly positive for867
and995, and the d_1 coefficient contains a positive multiple of the other
side's H. This proves positivity for **any** nonempty other side, not just
A_b, and covers a=1 or b=1 without invoking the negative-quadrant argument.

For996 the whole tridiagonal matrix is entrywise nonnegative and
`K_11=135618`. Keeping that one entry proves the universal lower bound
`135618 H(A)H(B)`. At A=B=K1 no other profile coordinate survives, so equality
holds. The independent checker also directly counts these actual order14
contexts and obtains135618, rather than merely evaluating the kernel formula.

Nonemptiness is material. Under the natural degenerate extension with one
outer side empty, the context becomes an ordinal join and the996 response
is its local H difference times the remaining H, namely `-3726 H(A)<0`.
That degenerate case is outside the theorem. The proof does not claim a
positive arbitrary-side response after removing this hypothesis.

The integer parameters do not interpolate continuously. Positive thin-side
responses and strictly negative responses throughout a,b>=2 leave no integer
zero between them. Nothing here asserts that every mixed-profile seed has a
one-crossing kernel, that every side lies in the negative chamber, or that
arbitrary outside-attachment changes retain this fixed context functional.

## 5. Independent exact receipts

`audit_zero_balance.py` imports only the already-audited independent
`../independent/audit_profile_dichotomy.py`, whose hash is pinned. It imports
no producer algorithm from `zero_balance_later/`. The actual module inputs
are read as asserted data and independently evaluated. Literal run-word
counts rebuild K. Side profiles are computed by the positive monomial/Stirling
expansion of `(m^3+2m)^a`, separately compared with multiplication by L^3+2L;
this is not the producer's alternating finite-difference evaluation.

Each normal/optimized run completes:

- 3 actual twelve-vertex profile comparisons;
- 507 literal run-word kernel entries;
- 1,331 direct actual three-module context counts;
- 64 positive Stirling side profiles and their recurrence checks;
- 852 positive log-concave integer sequences (length1..5, entries1..6),
  5,112 log-concavity/likelihood-ratio lift checks, and2,556 nonnegative
  polynomial-multiplier checks;
- the explicit non-log-concave failed-assumption hostile;
- 12,288 finite sign controls over1<=a,b<=64 for the three seeds;
- a direct sharp996 singleton-context count.

The finite parameter and abstract-sequence controls are **diagnostics**, not
the source of the all-order quantifiers. Those come from the audited proofs
in Sections1–4. Normal and `-O` outputs are byte-identical, SHA256
`2510edaf1ee4c59cf1d96619ab71c5bea57376a539e464423c0f2e2f8fa5d770`.
All gates raise explicitly; optimization does not remove them.

Reproduce from this audit directory, retaining the sibling producer packet
and frozen independent-kernel directory:

```sh
python3 -B audit_zero_balance.py > replay.normal.json
python3 -B -O audit_zero_balance.py > replay.optimized.json
cmp replay.normal.json normal.json
cmp replay.optimized.json optimized.json
sha256sum -c SHA256SUMS
```

Reviewed source versions are in `SOURCE_PINS.json`. The producer's final
curated manifest also verifies. Six producer payload files plus its two
manifest files form an eight-file later packet; this audit's six files are
separate. No frozen checkpoint-one file, source worktree or Git state was
changed. No verification blocker remains.
