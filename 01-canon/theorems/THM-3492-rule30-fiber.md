---
id: THM-3492
title: "Rule 30 slack-phase corner fiber product and pointed-carrier ambiguity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  The complete inward slack polynomial and the calibrated terminal phase
  profile meet canonically in one center bit.  Their exact common refinement
  is a corner fiber product, but every bivariate lift has a large mixed
  kernel, even after imposing the physical cyclic-arc image.  Each pointed
  Pascal-face coordinate gives an explicit section carrying the monic inward
  face, and depth six gives a physical two-section hostile.  No Rule 30 prize
  consequence is claimed.
source: root-rule30-next-targets-20260816
audit: >
  An independent hostile audit rederived both exact sequences and their
  dimensions, rank-checked additional non-power phase sets, reconstructed the
  depth-five and depth-six controls from a separate local-rule/Green path, and
  verified the pointed-carrier routing.  Ordinary and optimized runs equal the
  stored transcript byte-for-byte: ACCEPT.
depends_on:
  - THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum
  - THM-3481-rule30-cyclic-arc-norm-rank-and-marked-innovation-spectrum
  - THM-3488-rule30-inward-slack-monicity-and-parity-cartier-ramification
  - THM-3489-rule30-packed-restart-and-pointed-pascal-face
related:
  - THM-3476-rule30-depth-four-transverse-jet-barrier-and-slack-pascal-atlas
  - THM-2538-anchored-transverse-gain-and-common-ancestry-arrival-boundary
script: 04-computation/rule30_slack_phase_fiber_thm3492.py
output: 05-knowledge/results/rule30_slack_phase_fiber_thm3492.out
script_sha256: f27374a9660c55b402f4acd31ddd4bc88fddd44200e7571a298b6606c4dbb57f
output_sha256: 3101f5286e4d6949d3fa2bdc33632d8a3e9e60a7605ed9a4392e05388333e81f
hash_basis: raw bytes
---

# THM-3492 -- Rule 30 slack-phase corner fiber product and pointed-carrier ambiguity

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-3488 retains the complete transport-slack polynomial and proves that its
top ballistic face is unique.  THM-3489 retains the complete terminal phase
profile and proves that its physical point is a sparse, transverse Pascal
face.  Those are two different coordinates.  This theorem gives their exact
minimal common object and then proves why the two one-variable records do not
select a unique joint slack-by-phase coupling.

## 1. Inheritance and physical corner data

All vector spaces and polynomials are over `F_2`.  The inheritance pass is:

1. closest mechanism: THM-3488's ballistic monicity together with
   THM-3489's pointed Pascal inversion;
2. canonical hostile: THM-3471's nonzero marked slack current which dies at
   `q=1`;
3. corrected near miss: retain both axes and their common marked value, but
   do not infer an unproved joint event law; and
4. least-used sidecar: the physical arc-image ideal and the choice of a
   pointed Hasse carrier inside it.

Fix a Rule 30 target/depth `t>=5`.  Retain THM-3471 and THM-3488's notation

```text
B_t=H_1(t-2),
C_t(q)=[z^t]R^in(z,q),
D_t(q)=B_t+C_t(q),
N=t-4.                                                (1)
```

THM-3471's three-strip circuit and THM-3488's definition of `C_t` give

```text
D_t(1)=B_t+C_t(1)=c_t.                               (2)
```

Since `N>=1`, adding the scalar `B_t` cannot alter the top coefficient, and
THM-3488 gives

```text
deg D_t=N,       [q^N]D_t=1.                         (3)
```

Let `p=P_t` be the edge period.  This theorem treats a **proper terminal
arc**:

```text
s=t mod p,       0<s<p,       h_*=p-s.               (4)
```

Put `V_p={f:Z/pZ -> F_2}` and let

```text
F_t=A_t Q_t in V_p                                  (5)
```

be THM-3489's complete terminal phase profile.  The light-cone calibration
is exactly

```text
F_t(h_*)=c_t.                                        (6)
```

Equations (2) and (6), rather than an analogy between the two Pascal
transforms, are the common corner.

## 2. The universal corner fiber product

For `N>=1`, write

```text
P_N=F_2[q]_(degree <=N),
e_q(D)=D(1),       e_h(F)=F(h_*),
K_(p,N)=V_p x_(F_2) P_N
       ={(F,D):e_h(F)=e_q(D)}.                       (7)
```

A bivariate table is an element of `P_N tensor V_p`; write it as
`Z(q,h)`.  Its two axis restrictions are

```text
Res(Z)=(Z(1,h),Z(q,h_*)).                            (8)
```

### Theorem 2.1 (exact corner sequence)

For every finite phase set of size `p>=1`, every marked point `h_*`, and
every `N>=1`, restriction gives the exact sequence

```text
boxed:
0 -> (q+1)P_(N-1) tensor ker(e_h)
  -> P_N tensor V_p
  -Res-> K_(p,N) -> 0.                               (9)
```

In particular,

```text
dim K_(p,N)=p+N,
dim ker(Res)=N(p-1).                                 (10)
```

This is the same universal tensor-marginal sequence as THM-2538's
transportation kernel.  For vector spaces `U,V` with nonzero functionals
`e_U,e_V`, contraction of a joint tensor onto its two marginals has

```text
0 -> ker(e_U) tensor ker(e_V)
  -> U tensor V
  -> V x_(F_2) U -> 0,                               (10a)
```

where the last map is
`T -> ((e_U tensor id_V)T,(id_U tensor e_V)T)`.
Here `U=P_N`, `V=V_p`, `ker(e_q)=(q+1)P_(N-1)`, and
`ker(e_h)` is the zero-marked phase space, so `(10a)` is exactly `(9)`.
The mixed-Haar checkerboard in THM-2538 is the smallest categorical instance;
the present theorem is its polynomial-by-phase instance.

### Proof

The two restrictions of any `Z` agree at `(q,h)=(1,h_*)`, so the image lies
in the fiber product.  Conversely, let `(F,D)` be compatible and put

```text
c=F(h_*)=D(1),
Z_0(q,h)=F(h)+D(q)+c.                                (11)
```

Then `Z_0(1,h)=F(h)` and `Z_0(q,h_*)=D(q)`, proving surjectivity.  If both
restrictions of `Z` vanish, coefficientwise polynomial division gives a
unique

```text
Z=(q+1)H,       H in P_(N-1) tensor V_p.             (12)
```

Marked restriction now gives `(q+1)H(q,h_*)=0` in the polynomial ring, so
`H(q,h_*)=0`.  This is precisely the left term of (9).  Its dimension is
`N(p-1)`, and (10) follows.  `square`

Thus the fiber pair is canonical, while a bivariate lift is not.  Choosing
any profile `g` with `g(h_*)=1` gives the explicit section

```text
sigma_g(F,D)=F+(D+c)g.                               (13)
```

For two such carriers,

```text
sigma_g(F,D)+sigma_(g')(F,D)=(D+c)(g+g') in ker Res. (14)
```

Equation (14) is the basic mixed-correlation ambiguity.

## 3. The physical arc-image constraint

The lower bound `ceil(t/2)<=P_t` gives `t<=2p`.  Under (4), either `t=s` or
`t=p+s`, and hence

```text
nu_2(t)=nu_2(s),
ell=2^nu_2(s)-1.                                     (15)
```

Identify `V_p` with

```text
F_2[X]/(X^p-1)=F_2[Y]/(Y^p),       Y=X+1.            (16)
```

THM-3481 gives the exact physical image

```text
I_t=im(A_t)=Y^ell V_p,
dim I_t=p-ell.                                       (17)
```

THM-3489 proves that marked evaluation is nonzero on `I_t`.  More generally,
for any subspace `I<=V_p` on which `e_h` is nonzero, put

```text
K_(I,N)=I x_(F_2) P_N.                               (18)
```

### Theorem 3.1 (arc-constrained exact sequence)

Restriction gives

```text
boxed:
0 -> (q+1)P_(N-1) tensor (I intersect ker(e_h))
  -> P_N tensor I
  -Res-> K_(I,N) -> 0.                               (19)
```

Consequently, for the physical arc image,

```text
dim ker(Res|_(P_N tensor I_t))
 =N(p-ell-1).                                        (20)
```

### Proof

Choose `g in I` with `g(h_*)=1`.  Formula (13) now lies coefficientwise in
`I` and proves surjectivity onto (18).  The kernel proof in Theorem 2.1
applies without change.  Since a nonzero functional on `I` has a
codimension-one kernel, its dimension is `dim(I)-1`; (20) follows from
(17).  `square`

The constraint reduces the ambiguity from `N(p-1)` to
`N(p-ell-1)`.  It does not normally remove it.

## 4. Pointed-coordinate routing of the monic face

THM-3489's profile face at the physical point is

```text
J_(p,s)={p-s+a : a subseteq_bits (s-1)}.             (21)
```

For every `j_0 in J_(p,s)`, define

```text
g_(j_0)(h)=[X^h](X+1)^j_0=[X^h]Y^j_0.               (22)
```

Every face index satisfies

```text
j_0>=p-s>=2^nu_2(s)=ell+1,                           (23)
```

so `g_(j_0) in I_t`.  Its phase Hasse moments and marked value are

```text
M_j(g_(j_0))=1 if j=j_0, and 0 otherwise,
g_(j_0)(h_*)=binom(j_0,h_*)=1.                       (24)
```

The second equality is Lucas' theorem: every `j_0` in (21) contains the bits
of `h_*=p-s`.

### Theorem 4.1 (all pointed carriers are exact sections)

For every `j_0 in J_(p,s)`, the bivariate polynomial

```text
Z_(j_0)(q,h)
 =F_t(h)+(D_t(q)+c_t)g_(j_0)(h)                     (25)
```

has all of the following properties:

```text
Z_(j_0)(1,h)=F_t(h),
Z_(j_0)(q,h_*)=D_t(q),
[q^e]Z_(j_0) in I_t for every e,

M_j^h(Z_(j_0))
 =M_j(F_t)+(D_t+c_t)*[j=j_0].                        (26)
```

In particular, the unique inward top face is preserved through marked
evaluation and is carried by exactly one chosen phase-Hasse coordinate:

```text
[q^N]M_j^h(Z_(j_0))=[j=j_0].                         (27)
```

### Proof

Use (2), (6), and (24) in (25).  The first two lines of (26) follow by
substitution.  Both `F_t` and `g_(j_0)` lie in `I_t`, proving the coefficient
constraint.  Taking phase Hasse moments and applying the singleton table in
(24) proves the last line.  Equation (3) proves (27).  `square`

If `s>1`, then the face has at least two indices.  Any distinct
`j_0,j_1 in J_(p,s)` give two arc-constrained bivariate lifts with identical
axes but different mixed Hasse routing.  Their nonzero difference is

```text
(D_t+c_t)(g_(j_0)+g_(j_1)),                          (28)
```

an element of the kernel in (19).  If `s=1`, the pointed face is a singleton,
but (20) is still positive for the physical `N>=1`; the full mixed lift is
still not unique.  Thus a sparse pointed face does not canonically choose a
joint coupling.

## 5. Cheapest nontrivially arc-constrained physical hostile: depth six

Depth five is the absolute first two-carrier hostile, but there `s=5` is odd
and hence `ell=0`, so its arc image is the whole phase space.  Depth six is
the first physical hostile for which the proper ideal `I_t=Y^ell V_p`
actually constrains every coefficient profile.

At `t=6`, exact Rule 30 evolution gives

```text
p=P_6=8,       s=6,       h_*=2,       ell=1,
N=2,           c_6=B_6=0,
C_6(q)=D_6(q)=1+q^2.                                 (29)
```

The terminal profile and its Hasse table are

```text
F_6=(1,0,0,0,1,0,1,1),
(M_0,...,M_7)(F_6)=(0,1,0,1,1,1,0,1).               (30)
```

The pointed face is

```text
J_(8,6)={2,3,6,7},
M_2(F_6)+M_3(F_6)+M_6(F_6)+M_7(F_6)=0=c_6.          (31)
```

Choose `g_2=Y^2` and `g_3=Y^3`.  Both lie in `I_6=YV_8`, both equal one at
phase two, and their Hasse supports are the singletons `{2}` and `{3}`.
Therefore

```text
Z_2=F_6+(1+q^2)g_2,
Z_3=F_6+(1+q^2)g_3                                  (32)
```

have the same physical phase axis `F_6`, the same marked slack axis
`1+q^2`, and all coefficient profiles in `YV_8`.  Nevertheless their top
slack faces occur in different phase-Hasse coordinates.  The constrained
mixed kernel already has dimension

```text
N(p-ell-1)=2(8-1-1)=12.                              (33)
```

This is a physical marginal pair and a physical image constraint.  The
displayed `Z_2,Z_3` are algebraic lifts of those data; neither is asserted to
be the raw Rule 30 event-by-event joint distribution.

## 6. Source, target, map, and loss ledger

| item | exact content |
|---|---|
| source | `P_N tensor V_p`, or `P_N tensor I_t` after the arc constraint |
| target | the corner pair `(F_t,D_t)` in `V_p x_(F_2) P_N` |
| map | `Z -> (Z(1,h),Z(q,h_*))` |
| preserved | full terminal phase profile, center-completed slack polynomial, common center, and monic top slack coefficient |
| destroyed | every element of the mixed kernel, hence slack-by-phase correlation and the phase carrier of the top face |
| sufficient algebraic sidecar | a section, equivalently a chosen carrier plus a mixed-kernel coordinate |
| needed physical sidecar | a raw space-time common refinement or an independently proved event-level coupling/chain homotopy |

The two Pascal structures are not identified.  The transport Hasse order is
an exponent in `q+1`; the phase Hasse order is an exponent in `X+1`.  The
fiber product records their one proved common scalar and no more.

## 7. Exact verifier and validity boundary

Reproduce the companion in ordinary and optimized mode:

```bash
python3 04-computation/rule30_slack_phase_fiber_thm3492.py
python3 -O 04-computation/rule30_slack_phase_fiber_thm3492.py
```

It uses explicit `check` gates.  The abstract universe checks (9), explicit
kernel bases, restriction ranks, and separated lifts for

```text
p in {2,4,8,16,32,64},
every marked phase h_*,
1<=N<=4,                                             (34)
```

for `504` exact parameter triples.  The physical-ideal universe checks direct
arc matrices, image ideals, (19), all pointed carriers, and two-carrier
hostiles for

```text
p in {2,4,8,16,32,64},
1<=s<p,
1<=N<=4,                                             (35)
```

for `480` exact parameter triples.  Bit-packed linear maps check whole finite
vector spaces, not sampled vectors.  An independent local Rule 30 evolution,
Green coefficient recursion, and inward-source reconstruction check every
entry in the depth-six hostile (29)--(33).  Universal scope comes from the
proofs above; the finite universes audit their implementations.

The specialization `q=1` still destroys the monic inward face: already
`D_6(1)=0`.  The arc rank and pointed face do not determine any mixed
slack-by-phase statistic because (19) measures the remaining ambiguity
exactly.  No infinite sequence of physical sections is constructed, and no
distributional or computational bound follows from choosing the sparse
section (25).

Accordingly this theorem proves no Rule 30 center nonperiodicity, balance,
density, unpredictability, or random-access complexity lower bound.  The next
exact target must supply a physical coupling sidecar -- for example the raw
space-time event refinement -- and prove that it selects or controls a class
in the mixed kernel, rather than treating an algebraic section as physical.
