# A prefix-changing weight-22 transport closes the row-fifteen boundary

**Status: PROVED exact finite-row transport; exact certificate passed;
[independent audit PASS](planar_jc48_sep06_weight14_audit.md).** The weight-24 payer of THM-4438 can be replaced
by a compensated weight-22 packet on its entire rational boundary `G_m`.
The construction changes the earlier row solution. It does not contradict
the proved least-weight statement for perturbations that preserve that
earlier solution. No full `B_2` Keller pair, termination, chart entry, or
planar Jacobian conjecture conclusion is asserted.

## 1. Inheritance and the precise obligation

The inherited chart is **THM-4308**,
[source-normal bracket and Hasse truncation](../../01-canon/theorems/THM-4308-source-normal-bracket-hasse-truncation-through-row-eight.md).
The complete weight-14 face is restored by **THM-4390**,
[row-nine face absorption](../../01-canon/theorems/THM-4390-source-normal-weight-fourteen-row-nine-face-absorption.md),
and **THM-4395**,
[row-ten affine absorption](../../01-canon/theorems/THM-4395-source-normal-weight-fourteen-row-ten-global-affine-absorption.md).
The selected late source continues through rows 11, 12 and 13 by
**THM-4399**,
[row-eleven response](../../01-canon/theorems/THM-4399-source-normal-row-eleven-late-weight-eighteen-response-absorption.md),
**THM-4403**,
[row-twelve affine continuation](../../01-canon/theorems/THM-4403-source-normal-two-channel-weight-eighteen-row-twelve-affine-continuation.md),
and **THM-4410**,
[row-thirteen continuation](../../01-canon/theorems/THM-4410-source-normal-least-weight-twenty-row-thirteen-affine-continuation.md).

The closest successful repair is **THM-4426**,
[weight-eighteen memory repair](../../01-canon/theorems/THM-4426-source-normal-row-fourteen-weight-eighteen-memory-repair.md).
Its two restored coefficients `z=[p^9]H` and `h=[p^6y^2]H` change earlier
rows and remove the obstruction left by the frozen-prefix response of
**THM-4415**,
[same-row rank obstruction](../../01-canon/theorems/THM-4415-source-normal-row-fourteen-same-row-response-rank-obstruction.md).
On `Phi=eta=0, alpha11=1, c51=1087/135`, simultaneous row-fourteen bracket
and projected depth give the rational curve `Q(z,h)=0`, isomorphic to `G_m`.

**THM-4438**,
[row-fifteen relative response](../../01-canon/theorems/THM-4438-jc2-row-fifteen-relative-response-on-boundary-gm.md),
already clears row-fifteen projected depth on that curve. Its no-response
debt is the **bracket quartic `P`**, while its depth debt is the already
imposed conic **`Q=0`**. These are different polynomials. The remaining
question treated here is whether a lower-weight source change can replace
its weight-24 bracket payer while retaining projected depth.

The canonical hostile is THM-4438's irreducible quartic after pulling `P`
back to `G_m`: every rational parameter fails without a new response.
The corrected near miss is to interpret “weight 24 is least in the
valuation-15 packet” as a lower bound when earlier rows may change.
The least-used sidecar is the kernel of the row-13 response: THM-4410 gives
the ratio `7/10` between `p^5y^4` and `py^6`. Here that kernel is transported
two rows further instead of discarded.

The five-concept board is: first source valuation; response-neutral
combinations; earlier-row memory; exact depth representatives; and the
separate quartic/conic obligations. The connection sends a finite row-15
solution to another by an explicit additive translation. It preserves the
Keller and `P-G` equations to their stated orders, degree caps, projected
depth, and rows zero through eleven. It changes source rows thirteen and
fourteen and the unknown rows beginning at twelve. It does not preserve
later equations. Targeted statement searches in the source-normal canon
and JC2 reports found no existing compensated transport with the constants
below. No external priority claim is made.

## 2. Universal finite-row translation on the declared low jet

Work over any characteristic-zero field. Put

```
p=t(1+x^2t), y=xtp, u=x^2t,
P(A,C)=C^2-A^3+(3/4)A+1/4.
```

Fix the first four rows of `(A,C)` to THM-4308's values with
`Phi=0, Delta=896/15`, so `K=-32/5`. In particular these are the first
four rows of every THM-4438 boundary solution. All later background
coefficients may be arbitrary: the translation identity below does not
depend on them.

Define the lower packet and the nonzero rational response constant

```
K22 = p^5 y^4 - (7/10)p y^6 - (508/135)p^2 y^6,
beta = 235202/27945,
T = K22 + beta p^3 y^6.                              (1)
```

There are explicit polynomials `dA,dC`, supported in rows 12 through 15,
such that, for every scalar `w`,

```
A' = A+w dA, C' = C+w dC, H' = H+w T                (2)
```

preserves

```
J_(x,t)(A,C)=1 mod t^15,
P(A,C)=-u/2+H(p,y) mod t^16,
pi_15(A) in pi_15(P_2), pi_15(C) in pi_15(P_3).       (3)
```

It also preserves the inherited degree caps and the entire unknown prefix
through row eleven. Its inverse is the same formula with `-w`, and the
translations compose by addition of their parameters. This is an exact
finite additive action, not only a tangent-space calculation.

Here are all its nonzero rows, with `dA=sum a_n t^n` and `dC=sum c_n t^n`:

```
a12 = 13x^6/30,
c12 = -13x^5(x^2+2)/40,

a13 = x^4(441x^4+490x^2-60)/135,
c13 = -x^5(441x^4+1606x^2+608)/180,

a14 = x^4(978075x^6+2573190x^4-304888x^2+147200)/93150,
c14 = -x^5(978075x^6+5746500x^4+4183508x^2-270848)/124200,

a15 = 2x^2(3912300x^10-8802675x^8+13143708x^6
             -20841156x^4+44222160x^2-92124320)/419175,
c15 = -2x(195615x^12+391230x^10-9212432)/27945.       (4)
```

For the depth assertions, the source and transcript give explicit
representatives `F2 in P_2`, `F3 in P_3` with

```
pi_15(F2)=dA, pi_15(F3)=dC.                           (5)
```

Their fourteen and thirteen summands respectively have the form
`c x^a u^b p^d y^e`, with `a+b<=2` or `a+b<=3`, and all have valuation at
least twelve. This is an explicit span witness using THM-4308's depth
module, not a test using degree caps alone. The complete rational
coefficients are retained in the short certificate's `WA,WC` tables.

**Proof.** Direct substitution of (4) gives, coefficient by coefficient,

```
J(A,dC)+J(dA,C)=0 mod t^15,
2C dC - 3A^2 dA +(3/4)dA = T mod t^16.               (6)
```

Only background rows zero through three enter (6): every omitted
background term has valuation at least four and every correction starts
at twelve. Differentiation in `t` lowers valuation by at most one.
Thus a missing background term cannot affect a Keller coefficient below
15 or a `P-G` coefficient below 16.

The nonlinear Keller variation is `w^2 J(dA,dC)`, of valuation at least
23. The nonlinear `P` variation is

```
w^2((dC)^2-3A(dA)^2)-w^3(dA)^3,
```

of valuation at least 24. Consequently (6) gives the full identities (3)
for arbitrary `w`, not merely to first order in `w`. The explicit
representatives (5) and linearity of projected depth give its preservation.
The rows in (4) satisfy `deg a_n<=n+1`, `deg c_n<=n+2`. This proves every
assertion about (2).

## 3. The weight-22 section on the entire inherited boundary

Write `H_pre(s)` for the partial source in THM-4426 on its boundary
parameter `s in k*`, before adding the THM-4438 response. Its largest
retained weight is 22. THM-4438 supplies row-15 solutions with source

```
H_old(s)=H_pre(s)+j(s)p^3y^6,
j(s)=-145 P(z(s),h(s))/(24 C),
C=10852621164972710686787843667734315747451565056000000000000000.
                                                               (7)
```

Choose the translation parameter

```
w(s)=-27945 j(s)/235202
    =1350675 P(z(s),h(s))/(1881616 C).                 (8)
```

The weight-24 coefficient in `H_old+wT` is exactly zero. The new source is

```
H_new(s)=H_pre(s)+w(s)K22.                             (9)
```

Its weight is at most **22**, and its precise changed coefficients are

```
[p^5y^4]H_new = w(s),
kappa20_new = kappa20_old -(7/10)w(s),
rho22_new   = rho22_old   -(508/135)w(s).              (10)
```

Every other coefficient of `H_pre` remains unchanged. In particular the
coordinates `z,h,c51` remain on the inherited conic `Q=0`; the earlier
bracket solution changes together with (10). Equations (2) carry all the
inherited row-15 solutions to solutions for (9), with exact projected
depth. The inherited `A^10` terminal fibres are carried isomorphically.

The constant `beta` never vanishes, and no new source coordinate is
inverted. These statements hold over every characteristic-zero field and
every `s in k*`. The only parameter poles are those already present in
the inherited `G_m` chart. The zeros of `w(s)` cause no rank singularity:
the translation is then simply the identity. Over the rationals, THM-4438's
irreducibility hostile shows `w(s)!=0` for every rational `s!=0`. Over an
extension containing a root of its quartic, `w=0` is allowed.

This constructs a new partial source of weight at most 22. It does not
classify all sources with that weight bound or prove that weight 22 is
minimal once arbitrary earlier coefficients may vary.

## 4. Failure controls and the mechanism

The initial cheap probe was the already bracket-neutral combination
`p^5y^4-(7/10)py^6` at row thirteen. The ratio `7/10` was inherited from
THM-4410, not rediscovered as a new single-row fact. Propagating its actual
unknown-row change gives the additional coefficient `-508/135` at row
fourteen. After this compensation, projected depth remains unchanged,
while row fifteen sees the nonzero response `beta`.

The exact failure sequence for the displayed correction (4) is:

* Omitting both compensators leaves row-thirteen discrepancy
  `-(7/10)x^6`.
* Including only the first leaves row-fourteen discrepancy
  `-(508/135)x^6`.
* Including both leaves row-fifteen discrepancy `beta x^6`; it is precisely
  the response used to replace the older payer.

A separate triangular coefficient calculation recovers the unique three
coefficients `(-7/10,-508/135,beta)` for this displayed transport. That
uniqueness is not a classification of every possible lower-weight
deformation. The proof uses every literal row coefficient and explicit
depth representatives; it never infers solvability from one scalar
observable without those checks.

The first failed implication in the earlier shortcut would be from a
fixed-prefix rank obstruction to an obstruction for moving prefixes. The
strongest survivor is THM-4438's original least-weight theorem in its own
valuation-15 packet. The new operation uses an earlier response kernel
as a later response direction. Its applicability here rests on the exact
valuation gap, which removes nonlinear dependence on the late source.

## 5. Reproduction and stopping boundary

The [standalone source](../../04-computation/planar_jc48_sep06_weight14.py)
and [output](planar_jc48_sep06_weight14.out) use exact SymPy arithmetic,
with no research-producer imports. All **169 checks** remain active under
optimization. The certificate retains all eight correction rows, both
complete depth representatives, every coefficient in (6), the finite
parameter identities, failure controls, and the exact replacement sign.

```
python3 -B 04-computation/planar_jc48_sep06_weight14.py
python3 -B -O 04-computation/planar_jc48_sep06_weight14.py
```

The proof's computation is a fixed polynomial identity certificate; it is
not a parameter census or an extension of a fixed-prefix response table.
The boundary existence input remains THM-4438. The certificate proves the
new transport directly and does not re-prove that input by sampling it.

Later row equations, infinite compatible lifts, polynomial termination,
full depth membership of a Keller solution, and chart entry remain open.
The next useful probes are the earlier response kernels at valuations
eleven and twelve, and the first later row at which the present translation
sees unknown background coefficients. Neither is treated by this result.
There is no inference to `JC(2)` or `DC(2)`.

Normal and optimized outputs are byte-identical. Raw SHA256 pins are

```
source a4a140ab538620e7885d8b77758c02e16a5c3a8de7e53e3daac51f01bda321f7
output 74d12cf42c23d9402fdb892faab12fd4ce7309d43561b3dc63ee92e9c6b2e132
semantic 2ee3438f27b78d55dc9275c443d019c13554bdf34be31224de414d736d4f9aba
```

Source and output are frozen. The independent analytic audit and a separate
93-gate Fraction reconstruction pass; see the linked audit.
