---
id: THM-2428
title: "Degree-twenty-two B-axis trigonal ramification closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  In the open first-flux chart of the genuine nonsplit polynomial
  exact-square-prefix degree-twenty-two branch, the last one-sparse
  coefficient axis B is empty. The natural weighted quotient gives an
  absolutely irreducible cubic cover of the v-line. Its discriminant is
  one squared quadratic times a squarefree degree-nine polynomial; the
  squared factor is exactly the intersection with the excluded first-flux
  wall. The nine remaining finite places are simple ramification. By
  Riemann--Hurwitz the normalization has genus at least three, so no
  nonconstant rational trajectory exists. This closes the B-axis, not
  mixed strata, degree twenty two, JC(2), or DC(2).
source: klein-2026-07-26-degree-twenty-two-b-axis-trigonal-triage
depends_on:
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
related:
  - THM-2406-degree-eighteen-H4-weighted-pole-deep-wall-collapse
  - THM-2423-degree-twenty-two-W-axis-genus-two-and-origin-cusp-closure
  - THM-2425-degree-twenty-two-CDE-axis-hyperelliptic-closure
script: 04-computation/jc2_degree22_b_axis_trigonal_ramification_thm2428.py
output: 05-knowledge/results/jc2_degree22_b_axis_trigonal_ramification_thm2428.out
script_sha256: 33e1f4c70ddab133b7d5502d6e75c1c2309f658af5c578c774d55e3eb8aff6d1
output_sha256: 0035079a9fa9e0ad8c7b586768a40e0790478e5ec04b3b75e9a0bca07e491ef9
hash_basis: working-tree bytes (LF)
---

# THM-2428 -- the degree-twenty-two B-axis is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2411 closes the first-flux pole divisor. In its complementary chart,
THM-2423 closes the `W`-axis and THM-2425 closes the hyperelliptic `C,D,E`
axes. This theorem treats the remaining one-sparse axis. Its exact
conclusion is

```text
genuine degree-22 trajectory,
mathcal A!=0,
B!=0,                  C=D=E=W=0
    => contradiction.                                           (1)
```

The mechanism is a trigonal branch count. Nine visible finite ramification
places force an even tenth contribution by Riemann--Hurwitz; locating that
last contribution at infinity is unnecessary.

## 1. The weighted B-axis quotient

Use the normalized degree-twenty-two coordinates of THM-2411,

```text
y=11s,                  u=dT,                  Z=T^2,

wt(y,u,Z,B,C,D,E,W)=(1,2,3,2,3,4,5,6).                 (2)
```

First suppose `y!=0` in `C(x)` and put

```text
v=u/y^2,                 zeta=Z/y^3,              p=B/y^2.
                                                                  (3)
```

After dividing the normalized fluxes `N_1,N_2` by `y^5,y^6`, respectively,
they become

```text
f_1
 =1331(616p-1089v+63)zeta
  +4[(-745360v+6160)p+922383v^2-25410v+63]
 =0,                                                        (4)

f_2
 =15944049zeta^2-162339408zeta v+2236080zeta
  -1190488992v^3+147581280v^2-1219680v+672
  +(65591680zeta+1443016960v^2-71554560v+98560)p
 =0.                                                        (5)
```

The open first-flux chart is exactly

```text
616p-1089v+63!=0.                                        (6)
```

Thus (4) reconstructs `zeta` uniquely. Eliminating it from (4)--(5)
introduces no component on (6) and gives

```text
Res_zeta(f_1,f_2)=198414832 R_B(v,p),                    (7)

R_B=A_3(v)p^3+A_2(v)p^2+A_1(v)p+A_0(v),                 (8)
```

where

```text
A_3=333921280(14641v^2+1694v-19),

A_2=-325248(65547757v^3+4172685v^2-305525v+2643),

A_1=1584(18649222647v^4+563356398v^3
          -136161300v^2+3108490v-17451),

A_0=-81(155624547606v^5+3215383215v^4-1700698560v^3
         +58124770v^2-855470v+2583).                    (9)
```

The four coefficients in (9) have gcd one, so no vertical value of `v`
makes (8) an identity in `p`.

## 2. Arithmetic and geometric irreducibility

The integer polynomial `R_B` is primitive. Modulo `13`, its four
coefficients as polynomials in `v` still have gcd one. At `v=1`,

```text
R_B(1,p)=4p^3+9p^2+4p+1             in F_13[p].          (10)
```

Direct evaluation at the thirteen residues shows that (10) has no root,
so it is irreducible. If `R_B mod 13` factored, coefficient primitivity
would force both factors to have positive `p`-degree. Since the leading
coefficient at `v=1` is `4`, both degrees would survive specialization,
contradicting (10). Hence `R_B mod 13`, and therefore `R_B` over
`Q(v)`, is irreducible.

It is also geometrically irreducible. Indeed, rational irreducibility makes
the Galois group of `Qbar/Q` transitive on the absolute factors. Since the
total `p`-degree is three, any nontrivial absolute factorization would have
three conjugate linear factors. The cubic discriminant would then be a
square in `Qbar(v)`. Section 3 gives a nonsquare discriminant, a
contradiction.

Thus the normalization `mathcal C_B` of (8) is one geometrically integral
curve, and projection to `v` is a degree-three morphism

```text
pi: mathcal C_B -> P^1_v.                                (11)
```

## 3. Nine finite branch needles force positive genus

The exact cubic discriminant is

```text
Disc_p(R_B)
 =-729915048896495616 Q_2(v)^2 K_9(v),                  (12)
```

where

```text
Q_2=131769v^2-20570v+189,                               (13)

K_9
 =6210933086199321858048v^9
  -636900851142634892319v^8
  -150815920271976966600v^7
  +6846896917609271116v^6
  +1075176441098315688v^5
  +42156959680552310v^4
  +515709064108744v^3
  -56346612564v^2-48979998760v-659045583.              (14)
```

Exact Euclidean gcds give

```text
gcd(K_9,K_9')=gcd(K_9,Q_2)=gcd(K_9,A_3)=1.              (15)
```

Hence the nine roots of `K_9` are distinct; none is a degree-drop point
of the cubic or a root of the squared factor. At a root `alpha`, divide by
the unit `A_3(alpha)` and work over `C[[v-alpha]]`. The local polynomial
discriminant has valuation exactly one. The order-to-normalization index
enters its discriminant with even valuation, so the order is already
maximal. Tameness in characteristic zero then forces local type `(2,1)`:
one simple ramification point. These nine finite places contribute at least
nine to the ramification divisor.

Let `g_B` be the genus of `mathcal C_B` and `R` the total ramification,
including infinity. Riemann--Hurwitz gives

```text
2g_B-2=-6+R,             equivalently R=2g_B+4.          (16)
```

Thus `R` is even. Since the nine finite simple contributions already give
`R>=9`, in fact

```text
R>=10,                         g_B>=3.                    (17)
```

The parity-forced tenth contribution is the missing global completion
sidecar. It may occur at infinity or at another place not used in the lower
bound; its location and ramification type are immaterial for exclusion.

The squared factor in (12) is also understood rather than discarded.
Substitution of the first-flux wall into (8) gives

```text
R_B(v,(1089v-63)/616)=81Q_2(v)^2/7.                     (18)
```

Thus its two roots are precisely the double intersections with the
already-excluded chart boundary (6). They are not nine lost branch points.

## 4. Excluding the rational trajectory

The pair `(v,p)` in (3) gives a rational map from `P^1` to
`mathcal C_B`. Since (17) gives positive genus, that map is constant.
Thus `v,p in C`. Because `B in C*`, equation (3) makes `y`, then `u`,
constant. First-flux reconstruction makes `zeta`, hence `Z,T,q`, constant.
The genuine deck fixes the constant field but sends `q` to `-q`, contrary
to `q^2=T!=0`.

It remains to treat `y=0`. Here

```text
mathcal K=0,

N_1=1331(616B-1089u)Z.                                  (19)
```

The open chart makes the parenthesized coefficient nonzero, so `Z=0`,
again contradicting `T!=0`. This proves (1).

## 5. Scope and next obstruction

Together, THM-2423, THM-2425, and this theorem empty all five one-sparse
weighted axes `B,C,D,E,W` in the complementary `mathcal A!=0`
degree-twenty-two chart. Candidate THM-2429 subsequently closes the first
support-two plane `C,W`; the other mixed parameter strata remain open. The
next exact target is therefore another support-two plane: retain first-flux
reconstruction, compute the discriminant divisor over the weighted parameter
line, and keep the third flux and full polynomial mate sidecar whenever the
two-flux curve has a rational component.

Branches outside the inherited reduction, split/even short edges, and
integral order raising remain open. Nothing here proves degree twenty two,
`JC(2)`, or `DC(2)`.

No tournament, knot, additive graph, carry process, or Kakeya set is
identified with the trigonal curve. The transferable operation is exact:
nine visible local branch obligations plus the global parity sidecar force
a tenth. Equation (18) separately identifies the double wall contacts, so
branch count and boundary contact are not blended.

## 6. Exact verification

Run

```bash
python3 04-computation/jc2_degree22_b_axis_trigonal_ramification_thm2428.py
python3 -O 04-computation/jc2_degree22_b_axis_trigonal_ramification_thm2428.py
```

The companion checks (4)--(9), the modulo-thirteen primitivity and
specialization certificate (10), the discriminant and all gcds
(12)--(15), the exact wall identity (18), and the `y=0` boundary (19).
All truth-bearing checks use explicit exceptions and remain active under
optimized Python.

Normal, optimized, and stored transcripts byte-match after LF
normalization. The declared hashes are over the working-tree bytes.

## 7. Independent hostile audit

An independent audit reran the companion normally and under `-O`,
byte-compared both transcripts with the stored output, and verified both
declared hashes. Direct first-flux substitution independently gave

```text
f_2(zeta=-4K/[1331 mathcal A])=112R_B/mathcal A^2,
```

which reproduces the resultant scalar in (7). Modulo `13`, it independently
found coefficient gcd one and the root-free specialization
`4p^3-4p^2+4p+1`; this is the same residue polynomial as (10), since
`9=-4` in `F_13`. It checked the Gauss-primitivity and nonzero-leading-term
steps, then rechecked that any absolute factorization would split the cubic
linearly and make its discriminant a square.

The auditor independently recomputed the discriminant as
`-Res_p(R_B,partial_p R_B)/A_3`, recovered (12), and verified all three gcds
in (15). It specifically attacked normalization at every `K_9` root: the
valuation-one polynomial discriminant, even index-square correction, and
tame local classification preserve all nine simple ramification
contributions. Riemann--Hurwitz then gives (17) exactly as stated.

Finally, the fixed-constant-field and `y=0` deck closures were checked
against THM-2411. No mathematical, typing, boundary, scope, or
reproducibility defect remains. **QED.**
