# Independent hostile audit: balanced nonsplit factor/Stieltjes classification

## Verdict

**ACCEPT AFTER TWO LOCAL REPAIRS.**

The factor/Stieltjes normal form, moment indexing, exact square-contact
defect, complete `e<=1` classification, and chord-family monodromy formula
all pass independent symbolic and finite exact checks.  Two statements need
repair before canonization:

1. the converse must explicitly retain the displayed balanced degree
   constraint; and
2. the Faber-flux roadmap must treat THM-2245's degree-fourteen carrier as a
   retrospective control, not a current target.  THM-2247 already closes
   degree fourteen, and the inherited degree-twenty-two nonsplit branch is
   also closed.  The live nonsplit polynomial exact-prefix train starts at
   degree at least twenty six.

Neither repair changes the proved balanced classification.

## Audited inputs and replay

LF-normalized SHA-256:

```text
029f872ec550e9b340c5fd572dcce374bf18b4372685ed37f7199f8ec4c07732
  .scratch/nonsplit_balanced_passport_next/BALANCED-PASSPORT-NEXT-FRONTIER.md

3a11ba7b3f8b49f83e324b03ef64eff0d3dae474f4334b5977f2875d17190a46
  .scratch/nonsplit_balanced_passport_next/balanced_passport_linear_ode_census.py

c6041980d5f27bc32556184a7c2fa4d09c04e33cfec5050f4c827fb0ccc00428
  .scratch/nonsplit_balanced_passport_next/balanced_passport_linear_ode_census.out
```

The source script was rerun under ordinary Python and `python3 -O`.
Both fresh transcripts agree byte-for-byte with each other and with the
stored transcript.  The script has zero AST `assert` nodes.  It reports
`1485` exact gates and `FAILED CHECKS: NONE`.

Independent audit artifacts:

```text
e9292c9731e17dd70ed674ae4936fbf6b29373a15c7a6bb9b9417c2eb417fb88
  .scratch/nonsplit_balanced_passport_hostile_audit/audit.py

c791f6846d4abfe52f4885724c02dd6bc6d0c1b1502c28cddd50a947680d5f0c
  .scratch/nonsplit_balanced_passport_hostile_audit/audit.out

c791f6846d4abfe52f4885724c02dd6bc6d0c1b1502c28cddd50a947680d5f0c
  .scratch/nonsplit_balanced_passport_hostile_audit/audit.opt.out
```

The independent script has zero AST `assert` nodes, compiles, and gives
identical ordinary and optimized transcripts.

## 1. Factor identity and direction audit

Under the balanced hypotheses, with monic pairwise-disjoint squarefree
`S,E,T`, pole divisor `D=prod_j (x-beta_j)^(p_j)`, and

```text
V = v S D T^2,
G = g E/(D T),
F = mu S E^2/D,
lambda = vg,
mu = vg^2,
M = S E T,
```

direct substitution gives

```text
2 S T E' + E T S' - E S sum_j p_j T_j = C,
C = 2 kappa/lambda.
```

The reverse calculation is also valid once the balanced degree constraints
are retained.  The zero/pole divisor argument is bidirectional: no hidden
zero of `G` or extra finite support is introduced.

### Minimal converse witness if balance is omitted

The sentence beginning “Conversely, pairwise-disjoint squarefree
`S,E,T`...” is false if read as making the constant identity alone
sufficient:

```text
S=E=1, T=D=x, p_1=1, C=-1,
kappa=1, v=1, g=-2.
```

Then

```text
V=x^3, G=-2/x^2, F=4/x,
2VG'+V'G=2,
```

and the factor identity is the nonzero constant `-1`.  But `F(infinity)=0`
and the finite zero and pole degrees are `0` and `1`, so this is not a
balanced response.

**Required repair.**  Replace the converse sentence by:

> Conversely, within the displayed balanced degree constraints
> `N=sum_j p_j=s+2e` (and hence in the present nonconstant chamber
> `h>=1`, `r>=1`), pairwise-disjoint squarefree `S,E,T`, positive `p_j`,
> and a nonzero constant identity construct a balanced square-potential
> solution.

Equivalently, explicitly say the converse inherits all hypotheses at the
start of the theorem.  The witness does not refute the theorem in that
intended scope.

## 2. Moments and exact defect

From

```text
F'/F = C/M,
deg(M)=r+1,
```

and the expansion of the signed residue Cauchy transform at infinity, the
indexing is exactly

```text
m_0=...=m_(r-1)=0,
m_r=C.
```

The `m_0` convention includes the degree contribution (`p_0=degree`) in
the Newton-sum implementation.  Thirty-one independently generated
symbolic packets reproduce the claimed vanishing range and first active
moment.

Also,

```text
V = v M^2 (mu/F),
F/mu = 1-(C/r)x^(-r)+O(x^(-r-1)),
```

so

```text
deg(V-vM^2)=r+2,
LC(V-vM^2)=vC/r=2 kappa v/(r lambda).
```

There is no off-by-one error.  Division by `r` is legitimate in the
balanced nonconstant chamber: `h>=1`, `N>=1`, and the balanced constraints
force `r>=1`.  All thirty-one independent symbolic packets reproduce both
the degree and leading coefficient.

## 3. Completeness of the `e<=1` chamber

The inequality `h<=e+1` makes the three listed cases exhaustive:

```text
e=0:       h=1, p=(N), cyclic family S=y^N+a;
e=1,h=1:   unique B_N=x^N-Nx+N-1 family;
e=1,h=2:   reciprocal-monomial chord family indexed by
            1<=d<=floor(N/2), modulo d<->N-d.
```

The `e=h=1` family has `S_N` monodromy; `N=2` is split and every `N>=3`
member is genuinely nonsplit.  The chord family has exactly
`floor(N/2)` affine/target-equivalence classes.

For the chord family, with `g=gcd(N,d)` and `m=N/g`, conjugates of the
chord transposition generate `S_m^g`, while the ambient `N`-cycle rotates
the `g` blocks.  Thus

```text
|Mon| = g (m!)^g,
Mon=S_N iff gcd(N,d)=1.
```

Twenty-five independent exact controls through `N=10` reproduce the pole
partition, transitivity, order, and full-symmetric criterion.

The first passport collision is also real.  At `N=5`, fixing the third
cycle `(0 1 2 3)` with fixed point `4`, the zero involutions

```text
(1 2)(3 4),   (1 4)(2 3)
```

both yield pole type `(4,1)` and monodromy order `20`, but are not conjugate
under the centralizer of the fixed third inertia.  Hence the bounded census
correctly identifies `e=2` as the first point at which passport data alone
lose a Nielsen/dessin coordinate.

## 4. THM-2760 and THM-2245 scope

The THM-2760 comparison is correctly labeled **methodological analogy**.
No Jacobi substitution, root transport, or import of its simple-root theorem
is proved here.

The THM-2245 interface admits a sharper exact formulation.  In the actual
nonsplit response setup,

```text
q=A/U,
R_Q=qK,
G=F_resp'/(2kappa)=R_Q/U,
U^2=V.
```

Therefore

```text
A K = V G = lambda M.                               (A)
```

THM-2245's Keller one-form

```text
A(2T K'+K T')=2kappa T
```

then gives

```text
(T K^2)'/(T K^2)=2kappa/(AK)=C/M,                  (B)
V T K^2=(AK)^2=lambda^2 M^2.                       (C)
```

Thus, in a live Faber computation, `P=AK` is the direct carrier.  Before
introducing an unknown `M` in the linear Padé gate, compute `P` and test:

1. polynomiality in the inherited polynomial branch;
2. `deg(P)=r+1`;
3. squarefreeness and support factorization `P=lambda S E T`;
4. logarithmic residues in the alphabet `(+1,+2,-p_j)`;
5. balance `sum p_j=s+2e`; and
6. the forced degree-`r+2` square-contact defect.

In the forward THM-2245 branch, polynomiality of `A,V` is inherited and is
not an additional conjecture.  When reconstructing a source backwards from
an abstract spectral solution, polynomial `A,V` remains a separate sidecar.

### Required roadmap repair

Section 8 should not call the THM-2245 degree-fourteen carrier “current.”
THM-2247 already proves that branch empty; later canon closes the inherited
nonsplit degree-twenty-two branch.  Use degree fourteen only as a
retrospective exact control for `(A)`--`(C)`.

The live program is:

```text
derive K=R_Q/q in the degree-26 nonsplit chart;
form P=AK;
test P against polynomiality, factor support, residue alphabet, balance,
and the forced square-contact defect;
only surviving rows proceed to spectral normalization.
```

This remains an **OPEN program**, not a closure result.

## Canonization recommendation

After the two repairs above, canonize with status split explicitly:

- **PROVED:** balanced factor/Stieltjes equivalence; Padé/log-derivative
  identity; moment annihilation and first active moment; exact square-contact
  theorem; complete all-degree `e<=1` classification and chord monodromy.
- **FINITE-EXACT:** the dessin census through `N=8` and its first duplicate
  passport; the additional symbolic and chord controls in the transcript.
- **ANALOGY ONLY:** relation to THM-2760.
- **OPEN:** degree-26 `K`, `AK`, and residue/defect elimination program.

Do not claim chart entry, polynomial reconstruction from an arbitrary
spectral point, nonsplit closure in degree at least twenty six, `JC(2)`, or
`DC(2)`.  With those boundaries stated, the result is strong enough for
canonization and gives a genuinely cheaper live carrier than the proposed
unknown-`M` search.
