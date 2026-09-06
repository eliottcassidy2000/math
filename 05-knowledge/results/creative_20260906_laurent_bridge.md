# Coupled coefficient windows and the minimal derivative cone

**Status: PROVED relative to THM-4440's stated signed coefficient input;
INDEPENDENTLY AUDITED + FINITE-EXACT controls.**
No theorem ID is reserved. The new result is a constraint-preserving observer
cone with exactly `min(k,n-k)` generators. It gives an explicit sufficient
certificate for signed quadratic responses at a vanished coefficient.
General trinomial doubled-row separation and the sharp Laurent first-return
bound remain **OPEN**; this note claims no new closure of either.

The [independent referee](creative_20260906_laurent_audit.md) accepts the
written proof and verifies fresh carriers and cone coordinates by a separate
exact method. Root also reviewed the complete proof, including the inherited
THM-4440 scope and strict endpoint conditions, and accepted the result.

## Inheritance and the live objects

The closest proved mechanism is
[THM-4440, signed duplication SOS and real-rooted Laurent return](../../01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md):
for a real-rooted real polynomial `H`,
`ord_0(H)<k<deg(H)` and `H_k=0` imply `[u^(2k)]H^2<0`.
The existing [virtual Hadamard construction](nc2_hadamard_transport_overnight_hexagon_sep05.md)
uses exactly this implication on a common carrier
`H(u)=u^q F(u)G(rho/u)`.

The canonical hostile is an individually positive opposite-coefficient
product, retained explicitly below. The corrected near miss is the
[three-group midpoint note](overnight7_20260906_laurent_midpoint_transport.md),
Section 4: separate contiguous replacements of the two path factors remain
real-rooted but lose the required common vanishing coefficient.
[MISTAKE-547](../../01-canon/MISTAKES.md) supplies the observer warning;
its p-adic metric results are not mathematical dependencies of this theorem.
The least-used sidecar here is the **paired displacement from the selected
coefficient**, retained through a coupled derivative operation.

The board has five objects: the actual doubled row; its virtual common
carrier; coefficient-preserving windows; the centered-square kernel;
and the nonnegative certificate cone. The decisive connection is between
the last three. It does not identify actual endpoint doubling with squaring
the auxiliary factors. Targeted repository searches for mixed derivatives,
polar derivatives, derivative duplication, real-rooted coefficient windows,
and Euler cones found no existing statement of the cone reduction below.
This is a repository novelty check, not a literature priority claim.

## 1. Every coupled window preserves the zero and gives a strict sign

Let `H(u)=sum_(j=0)^n h_j u^j` be a real polynomial, all of whose roots
are real, with exact degree `n>=2` and `h_0 h_n != 0`. Fix an integer
`0<k<n` and suppose `h_k=0`. Write `D=d/du` and
`rev_d P=u^d P(1/u)`, with the displayed ambient degree retained even
if an endpoint coefficient vanishes. For integers

```text
0<=r<k,             0<=s<n-k,
H_(r,s)=D^r rev_(n-s) D^s rev_n H,
W_(r,s)(j)=(j)_r (2k-j)_r (n-j)_s (n-2k+j)_s,
```

where `(a)_b` is a falling factorial, a coefficient computation gives

```text
H_(r,s)(u)=sum_j (j)_r (n-j)_s h_j u^(j-r),
[u^(k-r)] H_(r,s)=(k)_r (n-k)_s h_k=0,                 (1)
Q_(r,s)=sum_j W_(r,s)(j) h_j h_(2k-j)
       =[u^(2k-2r)] H_(r,s)^2 < 0.                    (2)
```

Coefficients with indices outside `0..n` are zero. Falling factorials with
nonnegative integer first argument smaller than their order are zero.
Thus `(2)` includes the actual endpoint deletions, rather than assigning
new coefficients outside the window.

Here is the complete strictness argument. Reversal, differentiation, and
reversal again preserve real-rootedness; a zero root after differentiation
is removed by the next reversal and introduces no nonreal root. The second
reversal uses ambient degree `n-s` because `D^s rev_n H` has that exact
degree, as `h_0!=0`. The final derivative is also real-rooted.

Moreover, `h_(k-1)` and `h_(k+1)` are both nonzero. To see this without
assuming a strict Newton inequality, a derivative of a real-rooted
polynomial cannot have a multiple root at zero unless the original
polynomial already has a root at zero of higher multiplicity. If `P(0)!=0`,
then at a critical point the logarithmic derivative identity gives

```text
P''(0)/P(0) - (P'(0)/P(0))^2 = -sum_i 1/r_i^2 < 0,
```

where the nonzero real roots are `r_i`; hence `P'(0)=P''(0)=0` is
impossible. If `P(0)=0`, the multiplicity statement follows by factoring
its zero root. Induction back through derivatives proves that two
consecutive zero coefficients of `H` are impossible when `h_0!=0`.
Apply this at the pairs `(k-1,k)` and `(k,k+1)`. Applying the same
logarithmic derivative identity to `H^(k-1)`, whose value at zero is now
nonzero and whose derivative at zero vanishes, also gives
`h_(k-1) h_(k+1)<0`.

Both adjacent coefficients survive the window: `r<=k-1` and
`s<=n-k-1`. Consequently
`ord_0(H_(r,s))<k-r<deg(H_(r,s))`. Formula `(2)` now follows directly
from the precisely scoped coefficient implication in THM-4440.
No simplicity of the roots of `H` is required.

At `r=k` or `s=n-k` the selected coefficient reaches an endpoint.
Its squared coefficient is zero, so the strict parameter ranges are sharp.
Factoring a zero power of `u` first extends the statement to a carrier with
zero roots, using its actual lower and upper nonzero endpoints.

## 2. The two-sided cone has only `min(k,n-k)` extreme rays

Reflect `H` if needed and write `M=k<=N=n-k`. On the actual paired
support put `d=|j-k|`, so `d=1,...,M`; `d=0` contributes zero because
`h_k=0`. The paired kernel depends only on `X=d^2`:

```text
W_(r,s)(X)=product_(a=0)^(r-1) ((M-a)^2-X)
             *product_(b=0)^(s-1) ((N-b)^2-X).
```

Define the `M` basic kernels

```text
B_t(X)=product_(a=0)^(t-1) ((M-a)^2-X),   0<=t<M.       (3)
```

**Cone reduction.** On `X in {1^2,...,M^2}`, every `W_(r,s)` with
`r<M,s<N` is a nonzero nonnegative linear combination of the `B_t`.
Conversely every `B_t` is itself the allowed window `(r,s)=(t,0)`.
The `B_t` are linearly independent on this grid. Thus the cone generated
by all two-sided windows is simplicial with exactly these `M` extreme
rays. Reflection gives the count `min(k,n-k)` without the ordering choice.

For a proof, multiplication in this basis obeys the exact identity

```text
(c-X) B_t = B_(t+1) + (c-(M-t)^2) B_t.                 (4)
```

Start with `B_r`, then multiply the `s` factors with
`c=(N-b)^2`, in order `b=0,...,s-1`. Discard `B_M`, which vanishes on
the entire paired grid. Before step `b`, every nonzero basis coefficient
has index

```text
t >= max(r, M-N+b).                                   (5)
```

This is true initially since `N>=M`. Under `(5)`,
`N-b>=M-t>0`, so the second coefficient on the right of `(4)` is
nonnegative. If equality holds it vanishes, and the least remaining index
increases by one, exactly as required by `(5)` at the next step. If the
inequality is strict the index can remain unchanged. This proves the
invariant and coefficient nonnegativity inductively. The resulting vector
is nonzero because `W_(r,s)(1)>0` in the strict ranges.

Finally `B_t((M-i)^2)=0` when `t>i`, and its diagonal value for `t=i`
is strictly positive. The evaluation matrix is triangular and invertible.
Its independent generators are therefore precisely the extreme rays of
their cone.

This is a reduction of observer complexity, not scalar atom separation.
The complete grid of `k(n-k)` derivative windows supplies only
`min(k,n-k)` distinct generating inequalities at a vanished coefficient.

## 3. Exact certificate interface and a hostile boundary

For any real weight function `w(d)`, `d=1,...,M`, define the signed
response

```text
R_w(H)=2 sum_(d=1)^M w(d) h_(k-d) h_(k+d).
```

Recover its unique coordinates in `(3)` by rational back-substitution:

```text
c_t = [w(M-t)-sum_(i<t) c_i B_i((M-t)^2)]
      / B_t((M-t)^2),              t=0,...,M-1.        (6)
```

If every `c_t>=0` and at least one is positive, then `R_w(H)<0` by
`(2)`. For rational input weights, `(6)` is an exact rational certificate
of membership in the entire coupled derivative cone. Membership is
necessary and sufficient **for this cone**, and only sufficient for
universal negativity. No converse about all possible signed responses is
asserted.

Nonnegative weights alone do not suffice. The real-rooted carrier

```text
H=(1-3u)^2(1+u)(1+3u/5)
 =1-(22/5)u+0u^2+(54/5)u^3+(27/5)u^4
```

has `k=M=2`. Its outer product `h_0 h_4=27/5` is positive. The
nonnegative kernel `w(1)=0,w(2)=1` gives response `54/5>0` and has
cone coordinates `(1,-1/3)`. Thus it lies outside the cone for an explicit
structural reason. The four allowed two-sided window responses are

```text
Q_(0,0)=-2106/25,
Q_(0,1)=Q_(1,0)=-7128/25,
Q_(1,1)=-21384/25.
```

Preserving the selected zero alone also fails: for `H=1-u^2,k=1`,
`(uD-1)H=-1-u^2` still has zero linear coefficient but squared quadratic
coefficient `+2`. The first failed implication is dropping real-rootedness
of the transformed common carrier. The repaired operation is the coupled
window in `(1)`. This also explains why this theorem does not authorize
arbitrary independent replacements of the auxiliary path factors.

## 4. Laurent / Hadamard connection and remaining obligation

Given the existing common carrier `H(u)=u^q F(u)G(rho/u)` at a real
nonzero root of `(F star G)(rho)`, remove its actual zero-root factor,
then apply `(1)--(6)` at the shifted selected coefficient. This map
preserves the one common zero and retains the original opposite-coefficient
products as explicit weighted summands. It forgets the separate ownership of the real roots by `F` and `G`; the
needed sidecar is their full original factorization and the literal
coefficient map to the actual doubled row.

For the existing path control on support `(-3,1,9)`, one has
`F=(1+u)^4`, `G=4+u`, `rho=-1`, and

```text
H=u F(u)G(-1/u)=(1+u)^4(4u-1)
 =-1+0u+10u^2+20u^3+15u^4+4u^5.
```

The allowed windows `r=0,s=0,1,2,3` have squared selected coefficients
`-20,-300,-2400,-7200`. Here `M=1`, so the four observations all lie
on the same single generating ray. Additional derivatives have added no
independent sign certificate, a concrete instance of the cone reduction.

The new precise next test is to express a proposed grouped midpoint
response as `R_w(H)`, with an explicit equality preserving all carries,
and apply `(6)`. No such expression for the general actual-minus-virtual
row is established here. This stopping boundary is stricter than treating
separate real-rooted replacement factors as sufficient: it requires the
same carrier and a verifiable paired kernel.

## 5. Reproducibility and scope of the computation

The [standalone checker](../../04-computation/creative_20260906_laurent_bridge.py)
imports no repository producer. Its real-root universe has exact degrees
`3..8`, every multiset of `n-1` factors from `(-3,-1,1,2)`, and every
interior `k` for which the last factor can be tuned to a nonzero rational
value making `h_k=0`. It retains 1,752 indexed carriers. For every allowed
`r,s`, literal reversed differentiation and ordinary convolution are
compared against the independent weighted-pair formula, followed by its
actual strict sign, centered-square identity, and endpoint controls.

The independent cone universe is `1<=M<=8`, `M<=N<=20`, and every
`0<=r<M,0<=s<N`. The constructive recurrence is checked against all
paired grid values and a separately coded triangular inverse. The
mixed-sign, real-rootedness, arbitrary-Euler, and endpoint hostiles above
are retained. These finite controls challenge the proof, not its
unbounded justification.

```bash
python3 -B 04-computation/creative_20260906_laurent_bridge.py
python3 -B -O 04-computation/creative_20260906_laurent_bridge.py
```

Matching output: [creative_20260906_laurent_bridge.out](creative_20260906_laurent_bridge.out).
The producer passes **492,428 explicit gates**, including 7,014 cone rows.
Raw LF SHA256 values are:

```text
source d258c2b9b450e2c6bc00bd4e303b6ecc6cea01d1983ab8c6d815928948398547
output 3ce7d347909de389021f7bc9d49b00b501ccb447d4a1cf4b29abd95fa7871946
semantic 6158b4bc952a64a39ec5ccb3490744357cc57dae9d1b1c0a81f7edc7b8ae964e
```
