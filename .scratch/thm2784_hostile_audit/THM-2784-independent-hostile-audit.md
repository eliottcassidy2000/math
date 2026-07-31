# THM-2784 independent hostile audit

## Verdict

**PASS WITH TWO MINIMAL PASSPORT WORDING REPAIRS.**

The function-field equivalence, finite DVR classification, both infinity
chambers, global divisor identities, Riemann--Hurwitz count, squarefree
classification, hostile family, and transfer to the all-degree Faber
response layer are correct.  The scope boundary is also correct: this is a
response-layer sieve inside an inherited polynomial exact-prefix chart, not
chart entry or nonsplit `JC(2)` closure.

Two phrases in the balanced-passport interpretation overstate the degenerate
degree-one case:

1. `2^e 1^s` can be the identity when `e=0`, so call the zero permutation
   "an involution or identity", or simply state its cycle type.
2. When `r=1`, the distinguished third cycle is a fixed point, not a
   nontrivial cycle.  Say "one distinguished cycle of length `r`, the only
   possibly nontrivial cycle" or "at most one nontrivial cycle".

No equation, classification, or consequence needs repair.

## 1. Both directions of the square-potential equivalence

Let

```text
K=C(x),       L=K(U),       U^2=V,       kappa!=0.
```

If `R` is anti-invariant and `R'=kappa/U`, then `F=R^2` is fixed by the
deck involution and

```text
F'=2RR'=2kappa R/U,
V(F')^2=4kappa^2F.                                    (1)
```

The derivative condition makes `R`, hence `F`, nonconstant.

Conversely, for nonconstant `F in K` satisfying `(1)`, put

```text
R=UF'/(2kappa).                                       (2)
```

It is anti-invariant and `(1)` gives `R^2=F`.  Differentiating either this
identity or `(1)` gives `R'=kappa/U`.  A division-free verification is

```text
V'(F')^2+2VF'F''=4kappa^2F',
V'F'+2VF''=4kappa^2,

U R'=(V'F'+2VF'')/(4kappa)=kappa.                     (3)
```

Here `F'!=0` uses characteristic zero.  Formula `(2)` selects the sign whose
derivative is `+kappa/U`; its negative has the opposite derivative.

The nonconstant hypothesis is essential: `F=0` satisfies the displayed
quadratic equation but does not produce the required response.

### Sharper base linearization

Put

```text
G=F'/(2kappa).
```

Then the square-potential equation is equivalently

```text
F=VG^2,
2VG'+V'G=2kappa.                                      (4)
```

The first equation gives the exact squareclass statement

```text
[F]=[V] in K*/K*2,          K(sqrt(F))=K(sqrt(V))=L.  (5)
```

Thus squaring the response loses the marked sign but not the quadratic
field.  Equation `(4)` is a useful linear rational ODE for future
classification of polynomial potentials.

## 2. Complete finite DVR audit

At a finite place write

```text
m=ord(V)>=0,        n=ord(F),        t=ord(F').
```

The valuation equation is

```text
m+2t=n.                                                (6)
```

If `n!=0`, characteristic zero gives `t=n-1`, hence

```text
n=2-m.                                                 (7)
```

If `n=0`, then `F` is a regular nonzero DVR unit and `t>=0`; `(6)` forces
`m=t=0`.  The exhaustive classification is therefore

```text
m=0:      n=2 at a zero; n=0,t=0 at a unit; no pole;
m=1:      n=1;
m=2:      impossible;
m>=3:     n=2-m, a pole of order m-2.                  (8)
```

In particular, every finite critical point of the rational map `F` lies
over zero or infinity.  The independent script exhausts all integer cases
`0<=m<=10`, `-10<=n<=10`, including every unit derivative order used in
the proof.

Equation `(5)` supplies a parity cross-check:

```text
ord(F)=ord(V) mod 2
```

at every place, exactly as `(7)` predicts.

## 3. Infinity signs and global ledgers

Let

```text
D=deg(V),     n_inf=ord_inf(F)=P-Z.
```

The sign convention is correct: with `t=1/x`,

```text
d/dx=-t^2 d/dt.
```

Hence

```text
ord_inf(F')=n_inf+1                   if n_inf!=0;     (9)
ord_inf(F')=r+1
  if n_inf=0 and r=ord_inf(F-F_inf).                 (10)
```

Since `ord_inf(V)=-D`, valuation of `(1)` gives:

### Unbalanced chamber

```text
n_inf=D-2,
s+h+e=1.                                               (11)
```

Because `V` is a nonconstant nonsquare polynomial over an algebraically
closed field, the sole finite marked point is either one simple root or one
high root.  Thus

```text
V=c(x-a),  or  V=c(x-a)^m with odd m>=3.              (12)
```

The oddness in the second case is forced by genuine nonsplitting.  The
response potential is a scalar multiple of `(x-a)^(2-m)`.  The
Riemann--Hurwitz contributions are the matching order-`m-2` zero/pole at
infinity/the finite point.

### Balanced chamber

```text
P=Z=N,
D=2r+2,
r=s+e+h-1.                                             (13)
```

The three partitions are exactly

```text
2^e 1^s,             (m_1-2,...,m_h-2),
(r,1^(N-r)).                                           (14)
```

Their sums are all `N`, and

```text
e+(N-h)+(r-1)=2N-2.                                   (15)
```

This is the genus-zero Riemann--Hurwitz identity.  All other finite
preimages of `F_inf` are unramified because the unit case in `(8)` forces
`F'` to be a unit.

For an actual rational map, `(14)` also gives useful inequalities omitted
from the main statement:

```text
N-r=e-h+1>=0,             hence h<=e+1;               (16)
m_i-2<=N for every high root.                         (17)
```

The exact permutation wording should be:

```text
the zero permutation has cycle type 2^e 1^s
  (an involution or identity);
the third permutation has one distinguished r-cycle,
  its only possibly nontrivial cycle.                 (18)
```

The independent integer ledger checks 1,568 infinity packets, including 96
balanced and eight unbalanced solutions of the displayed valuation
identities.  An independent rational atlas checks that `r<=N`, all partition
sums, and `(15)` for every realized map.

## 4. Squarefree classification

If `V` is squarefree, `(8)` gives no finite pole of `F`, so `F` is a
nonconstant polynomial.  If `n=deg(F)`, comparison of leading degrees in
`(1)` gives

```text
D+2(n-1)=n,          n=2-D.                           (19)
```

Thus `D=n=1`.  For `V=cx+b`, the common simple zero and leading coefficient
comparison give uniquely

```text
F=(4kappa^2/c^2)V,
R=(2kappa/c)U.                                         (20)
```

The residue argument in degree two and the holomorphic-differential
argument in degrees at least three have the correct infinity orders.  They
are geometric explanations of the cheaper degree proof, not extra
assumptions.

## 5. Hostile family

For

```text
F_n=4(1-x^(-n)),
V_n=x^(n+2)(x^n-1)/n^2,
```

one has

```text
F_n'=4n x^(-n-1),
V_n(F_n')^2=4F_n.                                     (21)
```

Every `V_n` is nonsquare because the roots of `x^n-1` are simple.  Its
radical has degree `n+1`, while the multiplicity `n+2` at zero supplies the
pole of `F_n`.  This exactly refutes replacing full squarefreeness with a
condition on the radical degree.  The independent controls verify
`1<=n<=9`, reaching radical degree ten.

The family is balanced with passport

```text
over 0:        1^n;
over infinity:(n);
over F_inf:    (n),
```

which also independently confirms `(13)--(15)`.

## 6. Exact transfer to the Faber response layer

The transfer requires the following inherited objects, all of which are
present in the theorem's stated chart:

1. THM-2202 makes the canonical quadratic approximate root polynomial:

   ```text
   P=H^2+L,  H=Vz^2+Bz+C_0,  L=Az+E.
   ```

2. Genuine nonsplitting gives `L=K(U)`, `U^2=V`, and

   ```text
   w=Uz+B/(2U),        dw/dz=U.
   ```

   Hence

   ```text
   Jac_(x,z)(P,Q)=U Jac_(x,w)(P,Q),
   Jac_(x,w)(P,Q)=kappa/U.                             (22)
   ```

3. The polynomial mate `Q` and every target shear `H(P)` are deck-fixed.
   In the unique full Faber gauge, triangularity by distinct leading
   `w`-degrees kills every odd seed.

4. For a surviving even seed, Laurent parity is

   ```text
   sigma(c_j)=(-1)^j c_j.
   ```

   Therefore `Phi` and `R=4c_3+p c_1` are anti-invariant, while `Psi` is
   invariant.

5. THM-2129's exact identity is

   ```text
   Jac(P,Q)=(w^2+p/4)Phi_Q'+w Psi_Q'+R_Q'.             (23)
   ```

   Comparing the `w^2,w,1` coefficients with `(22)` yields

   ```text
   Phi_Q'=Psi_Q'=0,          R_Q'=kappa/U.             (24)
   ```

6. The constant field of `C(x,U)` is `C`.  Since constants are deck-fixed
   and `Phi_Q` is anti-invariant,

   ```text
   Phi_Q=0,             Psi_Q in C.                    (25)
   ```

The ingredients in steps 2--6 are all-degree; the low-degree spectral
atlases in THM-2214 are not being extrapolated.  Conversely, this transfer
does not produce the polynomial exact-prefix chart for an arbitrary Keller
pair and does not use the square-potential equation to close the surviving
linear case.

## 7. Exact controls

The candidate companion at the time of audit passes normal, optimized, and
stored-output replay, with no Python `assert` node:

```text
script  7fcbe78730b20ffe06dcc16bf640600ca5810631789e9b6d74eda113f7b053d1
output  07922e330700c495d0f2f541727f5a4383b17dd37c03151ea6658c111f5b5fee
gates   316
```

The independent companion uses the linearization `(4)`, a different
four-point rational atlas, and an independent integer-ledger enumeration.
It checks:

```text
36 polynomial-potential atlas rows, 20 nonsquare;
9 hostile-family rows;
96 realized local-factor rows;
48 realized rational maps;
275 exhaustive finite local integer cases;
1,568 infinity ledgers;
12 squarefree degree rows;
30 reduced Faber rows and their deck characters;
3,185 exact gates;
0 Python assert nodes.
```

Artifacts:

```text
.scratch/thm2784_hostile_audit/thm2784_independent_hostile_audit.py
.scratch/thm2784_hostile_audit/thm2784_independent_hostile_audit.out
```

LF-normalized hashes:

```text
script 211cead903f60526ae5bf8fef6c8cdaf657478c6bcad30a09da308523d75246e
output 3dd561b26a565c5cf5803830bc638a42a0855c8fd88d4db5c0091a6f3aded14a
```

Normal and optimized executions byte-match the stored independent output.
After the two wording repairs in `(18)`, THM-2784 is suitable for
`PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED` promotion.  Its
response-layer and chart-entry scope boundaries should remain unchanged.
