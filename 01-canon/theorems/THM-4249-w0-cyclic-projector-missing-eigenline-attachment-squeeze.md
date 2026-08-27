---
id: THM-4249
title: "W=0 cyclic-projector and torsion-envelope attachment squeeze"
status: >
  PROVED RELATIVE TO THM-4230/4241 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. On the W=0 gate U*Z*(U+Z) nonzero, every degree-34/42 map that
  collapses the twelve attachments has zero u-eigenline coordinate. The
  residual full-Hom shells have 176 and 132 source-target symmetry classes,
  and THM-4230's finite marked-ratio sets S_34,S_42 lie in exact CM-torsion
  envelopes of cardinality 55 and 34 after an exact good-reduction exclusion
  of the common ratio 1/3. Their emptiness, W=0, M=12, seam entry, JC(2),
  and DC(2) remain OPEN.
source: codex-planar-jacobian-breakthrough-20260826
depends_on:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4241-w0-hidden-cm-saturation-and-visible-hidden-index-four-gluing
related:
  - THM-4218-exact-weight-ten-hidden-elliptic-tail-degree-three-planar-jacobian-exclusion
  - THM-4247-w0-involution-degree-twelve-attachment-exclusion
mistake_firewall:
  - MISTAKE-521
  - MISTAKE-522
scripts:
  - 04-computation/jc23_w0_cyclic_projector_squeeze_thm4249.py
  - 04-computation/jc23_w0_cyclic_projector_squeeze_thm4249_independent_audit.py
  - 04-computation/jc23_w0_ratio_one_third_hidden_projection_exclusion_thm4249.py
  - 04-computation/jc23_w0_ratio_one_third_hidden_projection_exclusion_independent_audit_thm4249.py
outputs:
  - 05-knowledge/results/jc23_w0_cyclic_projector_squeeze_thm4249.out
  - 05-knowledge/results/jc23_w0_cyclic_projector_squeeze_thm4249_independent_audit.out
  - 05-knowledge/results/jc23_w0_ratio_one_third_hidden_projection_exclusion_thm4249.out
  - 05-knowledge/results/jc23_w0_ratio_one_third_hidden_projection_exclusion_independent_audit_thm4249.out
script_sha256:
  - 64cefef1ab610cdeab05eeaaeff25ae03bb2c69095f734e86d69b92fdccfea10
  - 1c1ae0d47f5218af5978cb840c0f6f9c564a6df338a7b650700cbca774e5e3c4
  - 408acb2025e40fbebcfd90ee7ecf5725a8c1cc961ee969d439f49ae713fb4d07
  - 0903c89f249c3428e17a1d15733fd03ec8f27a329d01f5b2bdfe6da09cdbb11e
output_sha256:
  - 10b26a5e47fcf75594b3956dbbb96d8458b9169c6cef47d95548624891f58a64
  - ec91ec63a3a2e58670d8d5b40d84027110b66c1c70b0f30d55953ac6021a4704
  - 23188f7b34dac946af984971353eafab7698865ae342f36bda7a40e389f7e42a
  - e321e56720110d780e25bae4701e9712464bda6115ed8df13a65dbe0a0576722
hash_basis: raw LF bytes
audit: >
  PASS/ACCEPT. The primary exact path verifies the C12 action, integral
  projectors, low hidden shells, full theta shells, spectral sieve, orbit
  census, boundary resultant, and envelope counts. A standard-library-only
  clean-room path reconstructs all 53,280 raw shell vectors, every projector
  value, the complete residual and orbit quotients, principal-ideal nesting,
  exact torsion-kernel unions and intersections, the raw 55/35 free
  mu_6-orbit counts, and the raw 1,644 map-ratio incidence frontier. Two
  further independent paths reconstruct all 3,168 degree-42 residual vectors,
  prove the uniform norm-three hidden-projector obstruction, and exclude the
  132 incidences over the common ratio 1/3, leaving the final 55/34 envelopes
  and 1,512-incidence frontier. Normal, optimized, and fixed-hash-seed streams
  byte-match their frozen outputs.
---

# THM-4249 -- `W=0` cyclic projectors and torsion envelopes

**PROVED RELATIVE TO THM-4230/4241 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED. THE FINITE SETS `S_34,S_42` ARE SQUEEZED, NOT PROVED EMPTY.**

## 1. Statement and inheritance

Retain the `W=0` gate and normalization of THM-4230/4241:

```text
C_0: x^6+y^4=1,                 E_0: Y^2=X^3+1,
O=Z[omega],                     omega^2+omega+1=0,
Q_j=tau^j Q_0,                  0<=j<12,
U*Z*(U+Z)!=0.                                           (1)
```

For `e in {34,42}`, let `S_e` be THM-4230's finite set of marked ratios
`R=U/Z` for which some degree-`e` map `C_0->E_0` sends all twelve `Q_j` to
one target point.

> **Theorem.** In the normalized full `O`-basis `[u,f,g,h]` of THM-4241,
> every degree-`34` or degree-`42` attachment-collapse candidate has zero
> `u` coefficient. After all three cyclic-projector constraints and the
> inherited pure-hidden exclusion, the residual contains
>
> ```text
> degree 34: 4224 vectors = 176 (mu_6 x <tau>)-orbits,
> degree 42: 3168 vectors = 132 (mu_6 x <tau>)-orbits,  (2)
> ```
>
> every orbit having size `24`.

> **Torsion-envelope addendum.** There are explicitly defined finite sets
> `R_34,R_42` of marked ratios, obtained from the CM torsion kernels below,
> such that
>
> ```text
> S_34 subset R_34,      |R_34|=55,
> S_42 subset R_42,      |R_42|=34.                    (3)
> ```

The ratio `1/3` is the unique admissible ratio coming from the common
nontrivial `3`-torsion orbit in the raw degree-`42` torsion envelope. The
good-reduction argument in Section 7 proves `1/3 notin S_42` and removes it
from the final set `R_42` in `(3)`.

The inheritance quartet is:

- closest proved mechanism: THM-4230's finite marked-ratio reduction;
- independent same-geometry control: THM-4247's involution/denominator
  exclusion of the degree-twelve hidden projection;
- canonical hostile: THM-4241's degree-`34/42` vectors in the index-four
  visible-hidden glue;
- corrected near miss: MISTAKE-521, because rational eigenspaces do not
  determine the integral Hom lattice; and
- least-used sidecar: the integral `C_12` group-ring action, whose
  non-idempotent projectors retain exactly the denominator information lost
  by rational spectral projection.

## 2. Integral cyclic projectors

Use THM-4241's normalized basis

```text
2h=v+omega^2 f+g,                                      (4)
```

where the visible maps are

```text
u=(-x^2,y^2),              v=(-x^-2,i y^2/x^3).       (5)
```

Precomposition by `tau` acts exactly on the induced Hom lattice as

```text
Tu=-omega u,       Tf=g,       Tg=-omega f,
Th=omega^2 h-omega f,          Tv=omega^2 v.           (6)
```

Write

```text
m=a_u u+b f+c g+d h,
ell=b f+c g+d(omega^2 f+g)/2.                          (7)
```

Direct multiplication in `O[T]` gives three integral spectral operators:

```text
P_u=-(T-omega^2)(T^2+omega),
P_v=-omega^2(T+omega)(T^2+omega),
P_L=-omega^2(T^2-omega^2)(T^2-omega),                 (8)

P_u(m)=a_u u,             P_v(m)=d v,
P_L(m)=2ell,              2m=2P_u(m)+P_v(m)+P_L(m).   (9)
```

These are integral operators, not rational idempotents. This is why they
remain valid on the index-four glue.

Put `D_j=[Q_(j+1)-Q_j]` in `J(C_0)`, with indices modulo twelve. Attachment
collapse for `m` is exactly `m(D_j)=O` for all `j`. Since

```text
(T^k m)(D_j)=m(tau_*^k D_j)=m(D_(j+k)),                (10)
```

every `P(T)m` also kills every `D_j`. Thus all three maps in `(9)` collapse
whenever `m` does. Equation `(10)` is a homomorphism statement and is
unchanged by translating a curve-map representative.

The Hermitian form of THM-4241 splits `(7)` as

```text
q(m)=4N(a_u)+N(d)+3K,             q(2ell)=12K,         (11)
```

where `K` is a nonnegative integer on the glued lattice.

## 3. Visible fibre bounds

On the attachment orbit,

```text
u(Q_j)=[-omega]^j u(Q_0),       v(Q_j)=[omega^2]^j v(Q_0). (12)
```

The gate makes `Q_j` pairwise distinct. It also makes the first six
`u(Q_j)` distinct and nonzero. If `a_u u` collapses, its common value is
fixed by `[-omega]`; because `1+omega` is a unit, that value is `O`.
Therefore `ker[a_u]` contains those six points and `O`, so

```text
a_u!=0  =>  N(a_u)>=7.                                (13)
```

If `d v` collapses and `d!=0`, one fibre of a degree-`4N(d)` map contains
the twelve distinct `Q_j`. Hence

```text
d!=0  =>  N(d)>=3.                                    (14)
```

These implications are necessary, not sufficient.

## 4. Exact low-hidden-shell obstruction

Let `H=2ell=A f+B g`. In the explicit normalization inherited from
THM-4230,

```text
T^6H=-H,                     T^8H=omega H.             (15)
```

If `H` collapses, its common value is killed by both `2` and `1-omega`.
These elements are comaximal in `O`, so the common value is `O`.

The exact low-shell census is

| `q(H)` | vectors | unit-`T` orbits | cyclic determinant |
|---:|---:|---:|:---|
| `6` | `24` | `2` | a unit |
| `12` | `24` | `2` | `2` times a unit |
| `24` | `24` | `2` | `4` times a unit |

Moreover, the degree-`24` shell is exactly twice the degree-`6` shell.
For `q(H)=12`, adjugating the matrix with rows `(H,TH)` shows

```text
[2]f(Q_0)=[2]g(Q_0)=O.                                (16)
```

For `q(H)=24`, first write `H=2r` with `q(r)=6`; the unit cyclic
determinant of `(r,Tr)` gives `(16)` again. This divisibility refinement is
essential: using only the determinant of `(H,TH)` would yield an unnecessarily
weak four-torsion statement.

To exclude `(16)`, let `rho` denote the quartic parameter called `a` in
THM-4230, Section 6.1:

```text
rho^4-2rho^3-2rho+1=0.                                 (17)
```

Up to nonzero constants, the `Y`-numerators of `f=f_rho` and `g=Tf` at an
attachment are

```text
y(t^2+rho^3),                  y(1+rho^3 t^2).         (18)
```

The gate makes both target values affine and gives `y!=0`. Thus `(16)`
forces both expressions in `(18)` to vanish, whence

```text
t^2=-rho^3,                 1+rho^3t^2=0,
rho^6=1.                                                (19)
```

But the exact resultant is

```text
Res_rho(rho^4-2rho^3-2rho+1,rho^6-1)=-108!=0.         (20)
```

Therefore a collapsing projector can have neither `K=1` nor `K=2`.
THM-4247 independently reaches the `K=1` exclusion by bounding and factoring
the reciprocal hidden denominator; the resultant route here is what puts that
row and the doubled `K=2` row under one projector argument.

Combining `(11)--(14)` with the exact theta shells leaves, for `a_u!=0`,
only

```text
e=34: (N(a_u),N(d),K)=(7,3,1),(7,0,2),
e=42: (N(a_u),N(d),K)=(9,3,1),(9,0,2).                (21)
```

Equation `(20)` removes all four profiles. Hence every collapse candidate
has

```text
a_u=0.                                                 (22)
```

## 5. Complete residual census

Starting from THM-4241's full theta shells, the exact sieve is

| degree | raw shell | `a_u=0` | projector bounds | pure-hidden removal | final |
|---:|---:|---:|---:|---:|---:|
| `34` | `36288` | `5184` | `4224` | `0` | `4224` |
| `42` | `16992` | `3600` | `3360` | `192` | `3168` |

The last `192` are precisely THM-4230's pure-hidden degree-`42` vectors,
already excluded from attachment collapse. The final `(N(d),K)` profiles
are

```text
e=34:
(4,10):864, (7,9):1248, (13,7):768,
(16,6):576, (19,5):576, (25,3):192;

e=42:
(3,13):672, (9,11):576, (12,10):864,
(21,7):768, (27,5):288.                               (23)
```

Quotienting by target units and the source `T` action gives `(2)`. The
counts are over the full glued Hom lattice, not the rational eigenspace or
the visible-hidden direct sum.

## 6. Exact CM-torsion envelopes

Put

```text
pi=omega^2-1,                 N(pi)=3,
P=v(Q_0).                                               (24)
```

Because `dv` collapses and `v(Q_1)=[omega^2]v(Q_0)`, every residual
candidate satisfies

```text
[d*pi]P=O.                                             (25)
```

Every residual principal ideal `(d)` divides a maximal one in the following
lists. After multiplication by `pi`, the maximal annihilator norms are

```text
e=34: 21,21,39,39,48,57,57,75,
e=42: 36,63,63,81.                                    (26)
```

Let `T_e` be the union of the corresponding kernels `E_0[d*pi]`. Exact
principal-ideal enumeration gives

```text
|T_34|=336,       every pairwise intersection=E_0[pi],
|T_42|=216,       every pairwise intersection=E_0[3]. (27)
```

The target unit group `mu_6` acts on either union. Its orbit-size profiles
are

```text
T_34: {size 1:1, size 2:1, size 3:1, size 6:55},
T_42: {size 1:1, size 2:1, size 3:1, size 6:35}.       (28)
```

For `P=(X,Y)` from `(24)`, direct substitution of `(1)` and `(5)` gives

```text
R=U/Z=-1/(X^3+1).                                     (29)
```

The function `X^3` is the exact quotient by `mu_6`. The three short orbits
in `(28)` are inadmissible:

- the origin is nonaffine and corresponds to the vanished endpoint;
- `E_0[pi]-{O}` has `X=0`, hence `R=-1`; and
- `E_0[2]-{O}` has `X^3=-1`, hence `Z=0`.

Every other orbit has size six and produces one distinct admissible ratio.
Let `R_34^tor,R_42^tor` be their images under `(29)`. Thus

```text
|R_34^tor|=55,                    |R_42^tor|=35.       (30)
```

For degree `42`, all four maximal kernels meet in `E_0[3]`. Its six points
outside `E_0[pi]` form one admissible unit orbit. Since

```text
psi_3(X)=3X(X^3+4),
```

that orbit has `X^3=-4` and `(29)` gives `R=1/3`.

Retaining the actual ideal `(d)` for each of the `176+132` map orbits,
rather than replacing it by a maximal ideal, leaves exactly

```text
degree 34: 864 map-ratio incidences,
degree 42: 780 map-ratio incidences,        total 1644. (31)
```

For any fixed `A^6,B^4`, the 24 root choices split into two `C_12` orbits.
The involution `(x,y)->(-x,-y)` exchanges them and permutes the full Hom
shell. Thus one canonical node orbit suffices for every incidence in `(31)`.

## 7. Exact exclusion of the common ratio `1/3`

At `R=1/3`, the point `P` lies in `E_0[3]-E_0[pi]`. Every degree-`42`
residual norm in `(23)` is divisible by three, hence `pi|d`; because
`pi^2` is associate to `3`, one has `E_0[3] subset E_0[d*pi]`. Thus this
ratio occurs as one necessary incidence for **all** `3168` residual vectors,
or all `132` size-`24` map orbits. Their profiles and hidden projectors are

```text
(N(d),K)=(3,13),(9,11),(12,10),(21,7),(27,5),
H=2ell=A f+B g,                    q(H)=12K.            (32)
```

If the original map collapses, so do `H` and `TH`; Section 4 makes their
common value `O`. Since

```text
TH=-omega B f+A g,
delta=A^2+omega B^2,                                  (33)
```

adjugating the two rows gives

```text
[delta]f(Q_0)=[delta]g(Q_0)=O.                        (34)
```

The norm-three element `pi=omega^2-1` cannot divide `delta`. Indeed, modulo
`pi` one has `omega=1` and `O/(pi)=F_3`, so `pi|delta` would give

```text
A^2+B^2=0 mod 3.
```

Because `-1` is not a square in `F_3`, this forces `pi|A` and `pi|B`. Then
`H=pi H'` for an integral hidden vector and

```text
q(H')=12K/N(pi)=4K.                                   (35)
```

Every listed `K` is nonzero modulo three, so `4K` is not divisible by six.
This contradicts THM-4230's fact that every degree in `O f+O g` is divisible
by six. Since `pi` is the unique Eisenstein prime above three,

```text
3 does not divide N(delta).                            (36)
```

Now specialize the coefficient field to the exact good prime

```text
q=397,
(zeta_12,rho,s,x,y)=(157,161,27,15,28).               (37)
```

Here `rho` is the quartic root from `(17)`, `s` is the sixth-root scale in
THM-4230's explicit map, and all cyclotomic, quartic-pair, scale, source,
and separability relations hold. The chosen attachment has

```text
x^6=1/4,             y^4=3/4,             t=379,
R=x^6/y^4=1/3.                                       (38)
```

Every displayed denominator is nonzero. Direct evaluation of THM-4230's full
explicit hidden-map formula `(34e)` gives

```text
F=f(Q_0)=(340,181) in E_0(F_397),
[6]F=(0,396),        [9]F=(35,0),        [18]F=O.    (39)
```

Thus `F` has exact additive order `18`. If `(34)` held in characteristic
zero, good specialization would give `[delta]F=O`. Applying
`[conjugate(delta)]` would then give `[N(delta)]F=O`, so `18|N(delta)`,
contradicting `(36)`.

The exact certificate also reconstructs all `132` orbit representatives as a
hostile control. In this reduction `g(Q_0)=[15]F`, the full `O`-annihilator
  of `F` is `(6+12omega)` of norm `108`, and none of the `132` coefficients
`A+15B` lies in that ideal. The 24 radical choices are the two `C_12` node
orbits exchanged by `(x,y)->(-x,-y)`; precomposition permutes the full lane,
so `(37)--(39)` test a complete canonical orbit rather than one arbitrary
root choice.

Therefore

```text
1/3 notin S_42,
R_34=R_34^tor,        R_42=R_42^tor-{1/3},            (40)
```

which proves `(3)`. Removing the `132` one-ratio map orbits from `(31)` leaves
the exact decisive frontier

```text
degree 34: 864 incidences,
degree 42: 648 incidences,                  total 1512. (41)
```

The remaining mixed attachment evaluations are undone.

## 8. Preserved data, sharp boundary, and nonclaims

The operation

```text
full Hom lattice -> integral T-projectors -> residual ideals
                 -> CM torsion kernels -> marked ratios              (42)
```

preserves the actual integral glue, attachment equality, source cyclic
symmetry, target units, and the marked ratio. It discards the hidden
coordinates after they have supplied necessary inequalities; consequently
membership in `R_e` is only necessary and can admit false positives.

The gate factors are sharp for this proof. `U=0` or `Z=0` destroys an
attachment endpoint, while `U+Z=0` merges the two six-point fibres and gives
the excluded ratio `-1`. No continuity claim crosses these walls.

The theorem does **not** prove `S_34` or `S_42` empty, exhibit a member of
either set, close `W=0` or exact `M=12`, classify the hidden-Hom locus away
from `W=0`, prove entry into the reduced seam, or prove `JC(2)`/`DC(2)`.

## 9. Reproduction

From the repository root:

```bash
python3 -B 04-computation/jc23_w0_cyclic_projector_squeeze_thm4249.py
python3 -B -O 04-computation/jc23_w0_cyclic_projector_squeeze_thm4249.py
PYTHONHASHSEED=4249 python3 -B \
  04-computation/jc23_w0_cyclic_projector_squeeze_thm4249.py

python3 -B \
  04-computation/jc23_w0_cyclic_projector_squeeze_thm4249_independent_audit.py
python3 -B -O \
  04-computation/jc23_w0_cyclic_projector_squeeze_thm4249_independent_audit.py
PYTHONHASHSEED=4249 python3 -B \
  04-computation/jc23_w0_cyclic_projector_squeeze_thm4249_independent_audit.py

python3 -B \
  04-computation/jc23_w0_ratio_one_third_hidden_projection_exclusion_thm4249.py
python3 -B -O \
  04-computation/jc23_w0_ratio_one_third_hidden_projection_exclusion_thm4249.py
PYTHONHASHSEED=4249 python3 -B \
  04-computation/jc23_w0_ratio_one_third_hidden_projection_exclusion_thm4249.py

python3 -B \
  04-computation/jc23_w0_ratio_one_third_hidden_projection_exclusion_independent_audit_thm4249.py
python3 -B -O \
  04-computation/jc23_w0_ratio_one_third_hidden_projection_exclusion_independent_audit_thm4249.py
PYTHONHASHSEED=4249 python3 -B \
  04-computation/jc23_w0_ratio_one_third_hidden_projection_exclusion_independent_audit_thm4249.py
```

The primary path uses exact SymPy/Eisenstein arithmetic. The clean-room path
uses only the Python standard library, imports neither the primary nor
SymPy, and reconstructs the full shells, projectors, sieve, ideal divisibility,
torsion kernels, intersections, unit orbits, and raw incidence frontier
independently. The two ratio-one-third paths use distinct shell/group
implementations and certify the final 34-ratio, 1,512-incidence frontier.
All frozen outputs report necessary envelopes and surviving open tests rather
than treating a projector certificate as an equivalence.

**QED.**
