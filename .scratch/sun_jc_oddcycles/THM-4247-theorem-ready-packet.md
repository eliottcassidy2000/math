# THM-4247 theorem-ready scratch packet

## Proposed title

`W=0 involution-projection theta refinement and degree-twelve attachment-wall exclusion`

## Proposed status

`PROVED RELATIVE TO THM-4230/4241 + VERIFIED-EXACT` after the two scratch
certificates are combined into one maintained certificate and the uniform
degree-12 denominator lemma below is written into the proof.  The theorem is
not a planar-JC closure and does not empty `S_34` or `S_42`.

## Exact statement

Use the notation of THM-4230/4241:

```text
C0: x^6+y^4=1,              E0:Y^2=X^3+1,
M=Hom(J(C0),E0),            iota:(x,y)->(x,-y),
V=M^(iota=+1),              L=M^(iota=-1).
```

For `m in M`, put

```text
v=m+m o iota,               ell=m-m o iota.           (1)
```

Then:

1. `v in V`, `ell in L`, and

   ```text
   q(v)+q(ell)=4q(m).                                (2)
   ```

2. The exact projection-degree histograms for all raw Hom vectors of the two
   response degrees are

   ```text
   q(m)=34:
   q(ell):count =
   12:1536, 24:2304, 36:5952, 60:5760, 72:8064,
   84:5376, 108:4992, 120:1728, 132:576.

   q(m)=42:
   q(ell):count =
   12:672, 24:288, 36:2304, 60:288, 72:3456,
   84:3072, 108:3744, 120:1728, 132:576,
   156:672, 168:192.
   ```

   In particular every hidden projection degree is a positive multiple of
   `12`.

3. The hidden degree-12 shell consists of exactly `24` maps in two orbits of
   size `12` under postcomposition by `mu_6` and precomposition by `T`.  No
   degree-12 hidden map sends all twelve admissible gate-interior attachments
   to `O`.

4. If a degree `34` or `42` full map collapses all twelve attachments, then
   its hidden projection does send all twelve attachments to `O`, and its
   visible projection sends them to one common value.  Consequently:

   - the `q(ell)=12` rows are impossible, removing `1,536` degree-34 vectors
     and `672` degree-42 vectors;
   - the degree-34 row `(q(v),q(ell))=(4,132)` is impossible by the fibre
     degree bound, removing `576` vectors;
   - the degree-42 row `(q(v),q(ell))=(0,168)` is the `192`-vector pure-hidden
     shell already excluded by THM-4230.

   The exact remaining raw-vector histograms are therefore

   ```text
   degree 34:
   {24:2304, 36:5952, 60:5760, 72:8064,
    84:5376, 108:4992, 120:1728}, total 34,176;

   degree 42:
   {24:288, 36:2304, 60:288, 72:3456, 84:3072,
    108:3744, 120:1728, 132:576, 156:672}, total 16,128.
   ```

These are raw lattice-vector counts, not quotient-orbit counts.  They do not
enumerate marked ratios and do not prove either finite set `S_34,S_42` empty.

## Proof

### A. Integral involution projections and the bivariate theta census

THM-4241 proves `L=M^-`, gives the saturated visible lattice, and gives an
`O=Z[omega]` basis `[u,f,g,h]` of `M` with

```text
2h=v0+P,                    P=omega^2 f+g,
q(u)=q(v0)=4,
H_L=[[6,-4-2omega],[-4-2omega^2,6]].                 (3)
```

Write a full vector as

```text
m=a u+b f+c g+d h,              a,b,c,d in O.
```

Since `iota` fixes `u,v0` and negates `f,g`,

```text
ell=(2b+omega^2d)f+(2c+d)g,                          (4)
v=2a u+d v0.                                         (5)
```

The eigenspaces are Rosati-orthogonal, proving `(2)`, and `(5)` gives

```text
q(v)=16N(a)+4N(d).                                    (6)
```

Conversely, `(4)` recovers `b,c` uniquely whenever the hidden coefficients
lie in the residue class `(omega^2d,d) mod 2L`.  Enumerating the bounded
positive-definite shells in `(3)`, grouped by this exact residue class and
then by `(6)`, is therefore an exact bijective enumeration of `M`, not an
upper bound.  The scratch certificate reproduces THM-4241's univariate theta
counts `36,288` and `16,992` before reporting the displayed refinements.  Its
coordinate boxes do not touch their asserted boundaries.

Modulo four, `(2)` and `q(v) in 4Z` give `4|q(ell)`.  THM-4241 gives
`6|q(ell)`, hence `12|q(ell)`.  A zero hidden projection would put `m in V`,
whose degree is divisible by four, impossible for `34` or `42`.  This also
proves positivity without relying on the census.

### B. Collapse implies the two projected collapse conditions

Choose `Q_*` among the `iota`-fixed points `y=0`, and use curve maps

```text
phi_m(Q)=m([Q-Q_*]).
```

The twelve attachments obey `iota Q_j=Q_(j+6)`.  If all
`phi_m(Q_j)=P`, then

```text
phi_ell(Q_j)=phi_m(Q_j)-phi_m(iota Q_j)=O,
phi_v(Q_j)=phi_m(Q_j)+phi_m(iota Q_j)=2P.             (7)
```

This is a necessary condition only.  The converse can lose a two-torsion
difference, exactly as predicted by the `F_4` glue sidecar.

### C. Uniform denominator bound for hidden degree 12

This paragraph must be written explicitly because THM-4230 states the
degree-bound argument in its degree-42 application even though its proof is
uniform.

For every nonzero hidden map `r in L`, the coefficient pair `(X_r/x,Y_r/y)`
lies in `t^-1K[t^2]`, and the group law preserves oddness.  Both points above
`t=0` map to `O` under `f,g`, hence under `r`.  Therefore the reduced
denominator `d_r(t)` of `X_r/x` is odd and has odd degree.  A finite root of
multiplicity `e` contributes at least `6e` to the `X`-pole divisor, except
that `t=0` contributes `6e-2`; the `t=+/-i` fibres cost more.  Thus

```text
6 deg(d_r)-2 <= 2 deg(r).                              (8)
```

At `deg(r)=12`, `(8)` gives `deg(d_r)<=4`, and oddness gives

```text
deg(d_r)<=3.                                          (9)
```

### D. The complete degree-12 hidden shell

Exact enumeration of `(3)` gives `24` degree-12 vectors.  The actions

```text
(a,b) -> epsilon(a,b), epsilon in mu_6,
(a,b) -> (-omega b,a)  (precomposition T)
```

split them into two disjoint free orbits of size `12`.  Convenient orbit
representatives are

```text
r0=(1-omega)f+g,                 r1=f+omega g.         (10)
```

Post-units fix the target origin and do not alter poles.  The `T` action
sends `t` to `-1/t`, preserving the wall pair `{+1,-1}`.  It is therefore
enough to audit `(10)`.

For `r1`, an equivalent representative `omega^2f+g` gives, by direct
coefficient-form elliptic addition and the relations

```text
z^4-z^2+1=0,
p^2-(1+2z-z^3)p+1=0,
```

the reduced denominator

```text
d_r1(t)=C1 t(t^2-1),                                  (11)
```

where the absolute resultant norm of `C1` is `65,536`.  The raw numerator is
a nonzero scale-square multiple; after deleting that harmless nonzero scale,
its absolute resultant norms at `t=0,+1,-1` are

```text
16, 2,176,782,336, 2,176,782,336,                     (12)
```

so no factor in `(11)` cancels.

For `r0`, direct addition first forms `f+(-omega f)=(1-omega)f` and then adds
`g`.  After the same algebraic reductions the raw denominator factors as

```text
C0 t(t^2-1) R8(t^2),                                  (13)
```

with nonzero leading-coefficient resultant `18,339,659,776`.  The raw
numerator values at `0,+1,-1` have nonzero absolute resultant norms

```text
4,477,456,
40,290,721,869,103,654,477,234,176,
40,290,721,869,103,654,477,234,176.                  (14)
```

Thus the reduced denominator contains `t(t^2-1)`.  Bound `(9)` forces all of
`R8` to cancel and proves

```text
d_r0(t)=C0' t(t^2-1).                                 (15)
```

As hostile controls, four exact good reductions

```text
(q,z,p,s)=(313,29,135,21),(349,24,246,28),
          (373,69,297,33),(397,157,161,27)
```

give reduced denominator degree `3` and monic reciprocal gcd `t^2-1` for
both orbit representatives.

For a hidden map, all twelve attachments map to `O` only if `r` and `Tr`
have a simultaneous pole on the base, equivalently if

```text
gcd(d_r(t),t^deg(d_r)d_r(-1/t)) != 1.                 (16)
```

Equations `(11),(15)` make the gcd exactly `t^2-1`.  But THM-4230's node
coordinate gives

```text
Z/U=((t^2-1)/(2t))^2.                                 (17)
```

Hence its roots `t=+/-1` lie on `Z=0`, excluded by the gate.  This proves
the degree-12 exclusion.

### E. Exact row deletions

By `(7)`, the entire `q(ell)=12` rows are impossible.  Their census sizes are
exactly `1,536` and `672`.

For degree `34`, `(2)` pairs `q(ell)=132` with `q(v)=4`.  A nonconstant
degree-four curve map has fibre divisor of total degree four and cannot
contain the twelve distinct gate-interior attachments, so all `576` vectors
in that row are impossible.

For degree `42`, `q(ell)=168` pairs with `q(v)=0`, hence `m` itself is pure
hidden of degree `42`.  The row has `192` vectors, matching THM-4230's exact
pure-hidden shell; THM-4230 excludes their attachment collapse.  Subtracting
these disjoint rows gives the displayed remaining histograms.

## Proof-gap audit

### Closed in this packet

- **Integrality of projections:** supplied by THM-4241 and formulas `(4),(5)`.
- **Basepoint translations:** closed by choosing an `iota`-fixed `Q_*`.
- **Degree-12 quantifier:** the uniform pole/oddness proof is restated in
  Section C; do not cite only THM-4230's degree-42 wording.
- **Both shell orbits:** audited separately; post-unit and `T` preservation is
  explicit.
- **Algebraic noncancellation:** numerator values at every proposed wall root
  have nonzero absolute resultants.
- **Gate legality:** the only reciprocal roots give `Z=0` by `(17)`.
- **Count subtraction:** deleted rows are disjoint and their remainder sums
  are asserted by the certificate.

### Remaining presentation/audit work before promotion

1. The two symbolic orbit calculations currently live in two scratch scripts.
   They should be combined into one maintained script so the proof does not
   depend on manually joining outputs.
2. The orbit-0 scratch probe now explicitly asserts the coefficient-form
   homogeneity `A -> s^2 A`, `B -> s^3 B`; retain that assertion when the
   scripts are consolidated.
3. The four finite-field embeddings are hostile controls, not the logical
   characteristic-zero proof.  The maintained theorem should present
   `(11)--(15)` and the resultant nonvanishing as load-bearing.
4. A genuinely independent implementation of the orbit-0 symbolic group law
   would strengthen audit status but is not mathematically required once the
   displayed identities/resultants are recorded and replayed.

No unresolved mathematical gap was found in the scoped theorem after
restating the uniform denominator lemma.  The packet does **not** support any
claim beyond the exact row deletions.

## Proposed maintained artifacts

Use one theorem certificate plus one compact output:

```text
04-computation/jc23_w0_involution_projection_degree12_attachment_thm4247.py
05-knowledge/results/jc23_w0_involution_projection_degree12_attachment_thm4247.out
```

Recommended output fields:

```text
hidden_gram=[[6,-4-2omega],[-4-2omega^2,6]]
full_theta_34=36288
full_theta_42=16992
projection_hist_34={...}
projection_hist_42={...}
projection_identity=qv+qell=4q:PASS
hidden_projection_divisible_by_12=PASS
degree12_vectors=24
degree12_unit_T_orbits=2,12,12
orbit0_denominator=t*(t^2-1)
orbit1_denominator=t*(t^2-1)
orbit0_resultants_nonzero=PASS
orbit1_resultants_nonzero=PASS
good_reductions=313,349,373,397
good_reduction_denominator_degrees=3,3
good_reduction_reciprocal_gcd=t^2-1,t^2-1
gate_wall_of_common_roots=Z=0
eliminated_34={qell12:1536,qv4:576}
remaining_hist_34={24:2304,36:5952,60:5760,72:8064,84:5376,108:4992,120:1728}
remaining_34=34176
eliminated_42={qell12:672,pure_hidden:192}
remaining_hist_42={24:288,36:2304,60:288,72:3456,84:3072,108:3744,120:1728,132:576,156:672}
remaining_42=16128
S34_S42=FINITE_UNENUMERATED_OPEN
result=PASS
```

Reproduction sources in this scratch packet:

```text
.scratch/sun_jc_oddcycles/jc_w0_projection_histogram.py
.scratch/sun_jc_oddcycles/jc_w0_hidden_degree12_attachment_audit.py
.scratch/sun_jc_oddcycles/symbolic_degree12_orbit0_probe.py
```
