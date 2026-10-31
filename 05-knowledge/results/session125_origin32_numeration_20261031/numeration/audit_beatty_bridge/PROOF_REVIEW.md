# Independent audit of the exact Beatty supplement

## Verdict

**PASS.** No correction is required to `BEATTY_BRIDGE.md`, `beatty_bridge.py`,
or their stated mathematical scope. The supplement correctly identifies
canonical Chung–Graham append-zero with

`P(n)=ceil(beta*n)`, where `beta=(3+sqrt(5))/2`.

It derives exact irrational-slope division, the Free/Locked phase threshold,
balanced numeral branching and the `1 -> 3 -> 8` orbit. These are statements
about the canonical numeral system, not square-sum graph Hamiltonicity or
ninja-pair existence.

## Universal proof review

### Strict fractional bound

For `A_i=F_(2i+2)`, the exact identity
`A_(i+1)-beta*A_i=beta^(-i-1)` gives
`P-beta*n=rho=sum_i d_i beta^(-i-1)`.
Reversal preserves legality, since each forbidden block `2 1* 2` is a
palindrome. The two-state fractional induction in the supplement is valid
at the empty word and all four permitted transition types. In particular,
a Free initial state gives `rho<1`, not merely `rho<=1`.

For a finite nonzero word, positivity of the digits gives `rho>0`.
Thus `P` is the unique integer strictly above `beta*n` and less than
`beta*n+1`. At zero, `P=rho=0` directly. The limiting infinite word
`2,1,1,...` has fractional value one, but is not a counterexample to the
finite statement. This is exactly where strictness matters.

### State threshold and first gap

A Locked word ends with `2 1^r`; a Free word ends with `0 1^r`, or is all
ones (including the empty-word case). After reversal, compare with the
infinite fractional all-ones value
`tau=1/(beta-1)=(sqrt(5)-1)/2`.
The common initial `r` ones only multiply the residual comparison by a
positive factor. A next digit 2 gives residual at least `2/beta>tau`;
a next digit 0 gives residual strictly below `1/beta<tau`.
Finite all-ones words lie strictly below the infinite sum.
This proves the two directions and excludes equality, without assuming the
desired numerical gap or relying on equal tree counts.

Since
`P(n+1)-P(n)=ceil(beta-rho)`, the gap is 3 for Free and 2 for Locked.
This includes `n=0`: `rho=0<tau`, `P(1)-P(0)=3`.

### Division, termination and zero

For integer `N>=0`, the choice `p=floor(N/beta)` is equivalent to
`ceil(beta*p)<=N<ceil(beta*(p+1))`. The right equivalence uses integer `N`:
`N<x` if and only if `N<ceil(x)`, including integral `x`.
Therefore `d=N-P(p)` is exactly one of the permitted digits after the
canonical word of `p`. These half-open integer intervals are disjoint and
cover all nonnegative integers, proving uniqueness. For `N>0`, `p<N`, so
repeated division terminates; `(p,d)=(0,0)` is the special zero case.
Leading zero prefixes and the empty representation preserve the arithmetic.

### Integer square-root implementation

For `n>0`, let `q=floor(n*sqrt(5))` and write `n*sqrt(5)=q+theta`, with
`0<theta<1`. Then

- `ceil((3n+q+theta)/2) = floor((3n+q)/2)+1`;
- `floor((3n-q-theta)/2) = floor((3n-q-1)/2)`.

These yield exactly the two code formulas. They are valid regardless of
parity of the integer numerators. The floor formula would give `-1` if
blindly evaluated at zero; the explicit zero return is therefore necessary
and present. The same zero exception is correctly present in the Wythoff
formula `P(n)=floor(beta*n)+1=n+floor(phi*n)+1`, which applies only to
positive `n`.

### Mechanical counts and scope

The identity `P(n)=2n+ceil(tau*n)` gives the Free indicator by an exact
ceiling difference. On any nonnegative integer interval its sum telescopes.
The ceiling error lies in `[0,1)` (zero at `n=0`, strictly between zero and
one otherwise), so the interval-count error has absolute value strictly
less than one. Hence counts on equal-length intervals differ by at most
one. The irrational limiting frequency rules out eventual periodicity.
The supplement limits this conclusion to numeral branching; no graph
occurrence sequence is inferred to be mechanical or Sturmian.

The lattice-strip formula is also exact: for a legal next digit `d`,
`(n,P)->(P+d,3P-n+3d)` gives `rho'=(rho+d)/beta`.
Digit 0 puts `rho'` below `1/beta<tau`; digit 1 gives
`rho'-tau=(rho-tau)/beta`, preserving the threshold side. Digit 2 stays
below one exactly when `rho<tau`, and then puts the new state above `tau`.
Thus the window recovers the legality automaton rather than introducing an
independent hidden degree of freedom.

### Shift orbit and two noncommutations

A zero appended to a canonical word remains canonical. Therefore applying
append-zero twice shifts each weight twice, giving
`P(P(n))=3P(n)-n`. The orbit of 1 is exactly the weight sequence beginning
`1,3,8`, and is the same sequence used in the earlier `D(1)` proof.

Both cautions are correct and necessary:

1. `P(P(1))=8`, whereas `ceil(beta^2)=7`; nested ceilings cannot be dropped.
2. Raw words `22` and `100` both have value 8, but their raw shifts have
   values 22 and 21. The Beatty formula requires normalizing first. The
   one-place relation `3A_0=A_1` does not persist unchanged after shifting.

### Norm, Diophantine triples and exact orbit classification

The identity `1/beta=3-beta` gives the useful exact relation
`parent(n)=3n-P(n)` for every `n>=0`, including zero.
For positive `n`, write `s=P(n)`, `t=P(s)=3s-n`, and define
`C=n^2-3ns+s^2`. Substitution of `s=beta*n+rho` gives
`C=sqrt(5)*n*rho+rho^2>0`. It is an integer, invariant under
`(n,s)->(s,3s-n)`. Direct expansion gives

`ns+C=(s-n)^2`, `nt+C=s^2`, and `st+C=(t-s)^2`.

Thus every positive shift orbit has consecutive `D(C)` triples. The
associated discriminant `D=2s-3n` obeys
`D^2-5n^2=4C` and `D congruent n (mod 2)`.
The quadratic integer `alpha=(D+n*sqrt(5))/2=P-n/beta` has conjugate
`rho=P-beta*n`, norm `C`, and transforms to `beta*alpha` under the
canonical shift. The parity condition places it in the quadratic integer
ring, so the unit interpretation is exact; no completeness claim about
arbitrary Pell solutions is being inferred.
For a legal appended digit `d`, weighted shifting gives
`P(s+d)=3s-n+3d`; expanding the quadratic form yields

`C(s+d)=C(n)+d*(2s-3n)+d^2`.

For `n>0`, the integer `2s-3n>sqrt(5)*n` is at least 3. A later nonzero
digit therefore increases the norm by at least four. The first nonzero
digits 1 and 2 have norms 1 and 4 respectively. Consequently `C=1` holds
exactly on the `1` shift orbit, and `C=4` exactly on the `2` shift orbit.
This classification is for positive canonical integers; zero is excluded.

For `n>0`, removing all trailing zero digits determines one primitive
positive word and a unique depth. Appending that many zero digits reconstructs
`n`. Equivalently, invert `P` only while the division remainder is zero.
The implementation returns `(0,0)` at zero as a convention: because `P(0)=0`,
zero does not have an intrinsically unique shift depth.

The norm is **not** a complete orbit invariant: 5 and 7 have norm 11, while
their distinct primitive words are `12` and `21`. Also, the general quotient
parent does not preserve the canonical norm. For instance `parent(5)=1`,
but `C(5)=11` and `C(1)=1`. The backward quadratic-form transformation acts
on the pair `(parent(n),n)`; this is the canonical pair
`(parent(n),P(parent(n)))` only when the remainder is zero.

## Independent executable checks

`independent_beatty_audit.py` imports the producer implementation only as
the unit under test. Its oracles use exact rational arithmetic in
`Q(sqrt(5))`, sign comparison by squaring, independently generated full
Fibonacci weights, and an integer-interval inversion. The oracle never calls
`isqrt`, never uses floating point, and contains no `assert` checks.

It checks:

- 250,001 consecutive integer shifts, parents, divisions, encodings and
  strict phase inequalities;
- every legal word through length 10, evaluating fractional expansions in
  the exact quadratic field and checking both Free/Locked upper bounds;
- random large inputs at 512, 4096 and 16384 bits;
- 10,000 random balanced intervals and zero/empty-word boundaries;
- 100,000 positive norms, `D(C)` triples, primitive-core decompositions,
  norm-one/four classifications and every legal next-digit norm increment.

The producer suite was separately replayed in normal and optimized Python;
both final replays passed **1,983,501 explicit checks**.
Both independent modes pass **4,676,438 explicit checks** against the frozen
code SHA256 `125bf43db526570cdd75d65ad21eed9ff7b8d637ea4a5e90ac5618bcff9cbb57`.
The two independent reports agree apart from the optimization-mode flag.
Reports in this audit directory record the actual check counts and outcomes.
All finite tests supplement the universal proof review above.

Run from the containing `numeration` directory:

```sh
python audit_beatty_bridge/independent_beatty_audit.py \
  audit_beatty_bridge/independent.normal.json
python -O audit_beatty_bridge/independent_beatty_audit.py \
  audit_beatty_bridge/independent.optimized.json
```

No producer source, source worktree or remote repository was edited.
`producer_input.sha256` records the supplement snapshot audited.

## Citation check and novelty boundary

Independently inspected
[Rob Burns, *Shifting Zeckendorf and Chung-Graham representations*,
Section 3](https://arxiv.org/html/2507.06162v1) on 2026-10-31.
Its Fibonacci identity (5), treatment of shifts and Lemma 3 support the
supplement's characterization of the source. The convention retains zero
odd Fibonacci slots, so its two-place shift corresponds to one compressed
even-Fibonacci place here. Its introduction discusses Wythoff and generalized
Beatty descriptions. No claim of historical novelty follows from these
self-contained derivations or from the finite audit.

## Portable package checks

The audit script locates its producer module relative to its own file, not
the current working directory. `producer_input.sha256` uses paths relative
to this audit directory; verify it here with `sha256sum -c producer_input.sha256`.
`portability_verification.json` records normal and optimized replays after
relocating this package and invoking it from an unrelated working directory.
The original scratch audit is preserved; these copies differ only in package
paths, reproduction instructions and regenerated packaging metadata.
