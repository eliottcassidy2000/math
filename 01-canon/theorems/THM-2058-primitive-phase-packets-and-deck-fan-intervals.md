---
id: THM-2058
title: "Primitive phase packets and one-dimensional deck/fan intervals"
status: >
  PROVED from THM-2053 and settled lower-dimensional LRC; all new arguments
  are elementary. A transverse denominator-N safe packet is the disjoint
  union of its reduced phase-order packets. Counts obey divisor summation,
  Mobius inversion, an exact rational Ehrhart/Beatty shift law, and a primitive-
  phase discrepancy bound. Across a longitudinal fiber the entire labelled
  packet transports by the unit M. After fixing a bad N and a THM-2055 hull
  owner, positivity, owner, and determinant constraints cut M to one explicit
  interval of coprime integers, minus finitely many collision points. Non-hull
  adjacent-pair representatives remain essential deck sidecars. This is an
  exact finite carrier theorem, not LRC(14).
source: codex-2026-07-21-LRC-primitive-packets
depends_on:
  - THM-2053
  - THM-2055
  - THM-2047
  - THM-1065-doubling-family-mod-six-characterization
related:
  - THM-1002
  - THM-1605-tnc-proved-monodromy-transitivity
  - THM-685
  - THM-2041
  - HYP-3036
  - HYP-3456
  - HYP-8846
  - HYP-8871
script: 04-computation/lrc14_primitive_phase_packet_referee_codex_20260721.py
output: 05-knowledge/results/lrc14_primitive_phase_packet_referee_codex_20260721.out
script_sha256: ad40c2702303288a6315d96dcf783a22fcbb202f1308d4cf4f84e86126f459c3
output_sha256: abacef8c372ad141806d9ce7f70b416b1fe3f5ee062f876dcb6dda648efbf813
---

# THM-2058 -- primitive phase packets and deck/fan intervals

## 1. Exact reduced-phase decomposition

Fix a finite nonempty template of nonzero integers

```text
m=(m_1,...,m_s)
```

and put `delta=1/14`. For an integer `x`, write `|x|_N` for its least absolute
residue modulo `N`. Define the complete safe packet and its size by

```text
S_N(m)={x mod N:14 min_k |x m_k|_N>=N},
L_N(m)=|S_N(m)|.                                      (1)
```

Let `U_d` be the units modulo `d`, with the convention `U_1={0}`, and define

```text
S_d^prim(m)=S_d(m) intersect U_d,
p_d(m)=|S_d^prim(m)|.                                 (2)
```

Then there is a disjoint, label-preserving decomposition in `Z/NZ`:

```text
S_N(m)=disjoint_union_(d|N) (N/d) S_d^prim(m).         (3)
```

In particular,

```text
L_N(m)=sum_(d|N) p_d(m),                              (4)
p_N(m)=sum_(d|N) mu_Mob(N/d)L_d(m),                   (5)
S_N(m)!=empty iff p_d(m)>0 for some d|N.              (6)
```

Here `mu_Mob` is the ordinary Mobius function. Equation (6) strengthens the
divisor monotonicity in THM-2053: the exact obstruction is absence of primitive
safe phase packets on every divisor, not merely failure at `N` itself.

### Proof

Given `x mod N`, put

```text
g=gcd(x,N),       d=N/g,       x=(N/d)a.
```

Then `a mod d` is primitive and `d` is uniquely determined. For every `k`,

```text
|x m_k|_N/N=|a m_k|_d/d.                              (7)
```

Thus `x` is safe exactly when `a` belongs to `S_d^prim(m)`. This proves the
set-level disjoint union (3), including the normalized residue vector and its
active runner labels. Counting proves (4), ordinary divisor-lattice Mobius
inversion proves (5), and (6) follows. QED.

## 2. Canonical transport across the longitudinal fiber

Keep the THM-2053 transverse coordinates

```text
v_k=a_kN+m_kM,       gcd(N,M)=1.                      (8)
```

At phase `ell/N`,

```text
ell v_k == ell m_k M (mod N).                         (9)
```

Consequently the actual denominator-`N` packet of the row is

```text
S_N(v)=M^(-1)S_N(m).                                 (10)
```

This is more than count invariance. The map `x |-> M^(-1)x` preserves reduced
denominator and transports the complete labelled residue/slack vector, active
runner indices, endpoint signs, and weak-versus-strict status. Therefore (3)
is a lawful whole-fiber carrier.

A single untransported phase is not lawful. For

```text
m=(1,-1,2,...,12),       N=27,
```

the primitive safe packet is `{2,25}`. For `M=2`, the actual packet is
`{1,26}`. Packet nonemptiness and all labels survive, but the phase `2` does
not. This is the exact LRC analogue of the repo's rule that a tied Gaussian
face must be preserved as a whole rather than atom by atom.

The analogy stops there. The packets in (2) are exact **phase-order** strata,
as in HYP-3036, not THM-2041's exact-frequency character projectors. Unit
Frobenius merely permutes a `0/1` safe indicator and supplies no safe seed.

### Stabilizer and orbit norm

The whole-packet rule has a second exact consequence imported from the
orbit-product mechanism in `THM-1605-tnc-proved-monodromy-transitivity.md`.
Put

```text
G=U_N,       A=S_N^prim(m),       H={u in G:uA=A},
B_M=M^(-1)A  (M in G).                                  (10a)
```

Then

```text
B_M=B_(M') iff M'H=MH,       #{B_M:M in G}=phi(N)/|H|.  (10b)
```

Moreover every primitive phase belongs to exactly `|A|/|H|` distinct packet
images. Consequently, for commuting indeterminates `(z_x)_(x in G)`,

```text
product_(B in G.A) product_(x in B) z_x
  = (product_(x in G) z_x)^(|A|/|H|).                  (10c)
```

Indeed (10b) is orbit--stabilizer. The unit action on primitive phases is
transitive, so double-counting packet--phase incidences makes their number
constant and gives the exponent in (10c). This is a useful exact norm identity,
but unlike the Gaussian/TNC application it has no nonconstant `ct` term to
oppose a constant Vieta product; it creates no safe seed.

There is also an unavoidable quotient loss. Since `A=-A`, one has `-1 in H`
for `N>2`, so the unlabelled packet forgets at least the orientation
`M <-> -M`. At `N=27` above, `H={1,26}`: the two orientations have identical
safe phase locations, while their signed residue vectors at a fixed phase are
negatives. Signed endpoint or owner arguments must therefore retain an
orientation/residue sidecar. The same orbit norm applies separately on every
divisor packet in (3).

## 3. Exact Ehrhart/Beatty law and primitive discrepancy

Let

```text
G(m)={theta in R/Z:min_k ||m_k theta||>=1/14}.         (11)
```

Because every `m_k` is nonzero, its boundary lies on the rational walls

```text
theta=(14r+1)/(14|m_k|) or (14r-1)/(14|m_k|).         (12)
```

Write all its components as disjoint closed circular intervals
`I_1,...,I_K`, allowing singleton intervals. Choose a cut outside `G(m)` when
the set is nonempty and proper, and let `Q` be a common denominator of all
endpoints (take `Q=1` when it is empty). Put

```text
rho=measure(G(m)).
```

Componentwise floor counting gives the exact shift law

```text
L_(N+Q)(m)=L_N(m)+Q rho.                              (13)
```

Here `Q rho` is an integer. Equivalently, for `1<=s<=Q` and `h>=0`,

```text
L_(Qh+s)(m)=h Q rho+L_s(m).                           (14)
```

Thus `L_N` is an Ehrhart quasipolynomial of degree at most one. When `rho>0`--
in particular for the repeated-speed THM-2053 templates, whose maximum is at
least `1/13`--it has degree one and on each residue class its zero set is an
initial segment. When `rho=0`, singleton contributions are periodic. This is
the exact count-level version of the HYP-3456 floor-clock transfer used
qualitatively in THM-2053.

There is also a uniform discrepancy bound:

```text
|L_N(m)-N rho|<=K.                                    (15)
```

Indeed, the number of `N`-grid points in one closed interval differs from its
length times `N` by at most one. Combining (5) and (15) yields

```text
p_N(m)>=rho phi(N)-K 2^(omega(N)).                    (16)
```

Therefore the explicit inequality on the right certifies a primitive safe
phase of exact order `N`. For mere safe existence, THM-2053's widest-component
cutoff can be sharper; (16) is useful when the primitive order or owner-labelled
packet itself is required.

## 4. A bounded primitive packet seed

Suppose the distinct absolute values among the `m_k` number at most twelve and
put `R=max_k|m_k|`. Settled lower-dimensional LRC supplies a phase of height at
least `1/13`. The arity-free pair-sum maximizer theorem THM-2047, Section 2,
then supplies a reduced maximizer `a/q` with

```text
q<=2R,       q divides x+y for two absolute template speeds x,y.
```

It follows that

```text
p_q(m)>0,       D_N(m)>=1/13 for every N divisible by q.  (17)
```

Thus every template has a bounded primitive good divisor; every bad modulus
must avoid at least one explicit divisor `q<=2R`. This is a seed theorem, not a
claim that primitive safe support is multiplicative.

## 5. The deck/fan intersection is one-dimensional

Fix one adjacent exposed pair from THM-2053. In a saturated coefficient-lattice
basis write

```text
B=[w_0 eta] in GL_2(Z),
c_k=a_k w_0+m_k eta.
```

If `d` is the primitive parameter covector and

```text
(N,M)^T=B^T d,
```

then

```text
d=B^(-T)(N,M)^T,
v_k=a_kN+m_kM,
gcd(N,M)=1.                                           (18)
```

Now fix `N>0`, a signed hull owner `p` from THM-2055, and let `w_p` be its
quarter-turned
linear functional, so the determinant norm in the owner cone is `d dot w_p`.
After substitution, failure of the exact THM-2053 gate is

```text
Q_B(N,M)<91 L_(p,B)(N,M),                             (19)
Q_B(N,M)=||B^(-T)(N,M)^T||^2,
L_(p,B)(N,M)=w_p dot B^(-T)(N,M)^T.
```

For fixed `N`, the left side minus the right side is a convex quadratic in
`M`. Its strict negative set is one bounded open interval or is empty. The
THM-2055 owner cone and all positivity constraints

```text
a_kN+m_kM>0                                           (20)
```

are linear half-lines in `M`. With a fixed half-open convention on owner ties,
their intersection is therefore one explicit interval `I_(N,p)`, possibly
empty, open, closed, or half-open.

The remaining primitive candidates before speed collisions are exactly

```text
{M in Z intersect I_(N,p):gcd(M,N)=1}.                (21)
```

Their number is the exact divisor/floor sum

```text
C_(N,p)=sum_(e|N) mu_Mob(e) #{r in Z:e r in I_(N,p)}. (22)
```

The quadratic endpoints need no floating-point evaluation: integer membership
is decided by the sign of the integer quadratic in (19), or equivalently by
integer-square comparisons at its algebraic roots. Distinctness removes only
explicit collision points

```text
(a_i-a_j)N+(m_i-m_j)M=0,                              (23)
```

at most one per pair unless two columns are identical, in which case that
chamber cannot contain a distinct row. Hence no two-dimensional disk scan is
needed: enumerate bad `N`, hull owners, and coprime integers in one interval,
then delete at most `binom(13,2)` walls.

## 6. Representative fibers are a load-bearing sidecar

Adjacent normalized **values** can have labelled fibers. If `P,Q` are adjacent
values with fibers `I,J`, then every `(i,j) in I x J` is a valid exposed pair.
Write

```text
c_i+c_j=g_(ij) w^0_(ij),
N_(ij)=phi(w^0_(ij))=(v_i+v_j)/g_(ij).                (24)
```

The identity-component deck for this pair samples a reduced phase `a/q`
exactly when

```text
q divides N_(ij).                                     (25)
```

Indeed, `a/q=ell/N_(ij) mod 1` with `gcd(a,q)=1` is equivalent to (25).
Therefore one arbitrary representative per normalized value is unsafe. A
divisor-covering subfamily can suffice, but its lawful sidecar is

```text
(i,j,g_(ij),N_(ij),m^(ij)).                           (26)
```

If `q|g_(ij)N_(ij)` but `q` does not divide `N_(ij)`, the pair-sum phase lies
in a nonidentity component of `ker(c_i+c_j)`; one must retain the component or
longitudinal residues, or the resolved phase itself.

This is precisely where THM-2055 hull deletion must not be overextended. It is
valid for the determinant support norm, but not for choosing deck
representatives.

### Exact one-tail witness

Take

```text
c_k=(k,0) (k!=12),       c_12=(12,1),
phi_w(x,y)=x+(w-12)y,
S_w={1,...,11,13,w}.
```

All unchanged columns have normalized value `(1,0)`. Pairing `c_r` with
`c_12` and taking `eta=(1,0)` gives

```text
w^0_r=(r+12,1),
m^(r)=(1,2,...,11,-r,13),
N_r=w+r,       M=1.                                   (27)
```

For the strict row `S_38`, retaining only the hull representative `r=13` gives

```text
D_51(1,...,11,-13,13)=3/51=1/17<1/14.                (28)
```

It misses the strict exit. The non-hull representative `r=10` gives

```text
D_48(1,...,11,-10,13)=4/48=1/12.                     (29)
```

At `ell=4`, this is the actual phase `1/12`, and the core
`{1,...,11}` gives the matching pigeonhole upper bound. Thus

```text
M(S_38)=1/12.
```

For the Goddyn--Wong boundary row `S_24`, whose exact height `1/14` is the
`n=13` case of `THM-1065-doubling-family-mod-six-characterization.md`, the
same hull choice gives `D_37=2/37` and `r=10`
gives `D_34=2/34`; both remain below `1/14`. These values give a concrete,
exact failure of a hull-only deck quotient.

More positively, every reduced safe phase `a/q` with `q<=12` on this plane is
captured by choosing

```text
r == -w (mod q),       1<=r<=q-1.                     (30)
```

Indeed `q|N_r`. For a uniform direct `q=12` capture, representatives covering
all eleven nonzero residues modulo `12` are necessary and sufficient. A safe
phase cannot have `w==0 (mod q)`, because the tail would sit at zero, so the
representative in (30) is always available. This is a finite labelled sweep
rather than retention of every raw column.

## 7. Sharp limits on unrelated transfers

The archive analogies become useful only after stating what they preserve.

1. **CRT is not available.** For `m=(1,-1,2,4,7)`, every primitive phase is
   safe at `3` and at `5`, so `p_3=2` and `p_5=4`, but `p_15=0`. Indeed
   `+/-{1,2,4,7}=U_15`; for every unit `a`, one listed `m` has `am=+/-1`, below
   the integer threshold `ceil(15/14)=2`.
2. **Positive spanning-tree faces do not create a seed.** For
   `m=(1,-1,2,...,12)` at `N=15`, THM-2053 gives `D_15=1/15<1/14`, while all
   danger sets contain phase `0`. Thus every pair has nonempty full-deck
   danger-set overlap; every resulting raw overlap edge, and hence every
   positive spanning-tree monomial, is nonzero although the safe deck is empty.
3. **Tiling-cube mask walks do not preserve the orbit image.** The map from a
   cyclic phase to its runner-danger mask is not translation by a fixed cube
   mask, and its image can omit the empty mask. THM-346's congruence mechanism
   therefore does not apply.
4. **Fixed characteristic-seven jets remain local.** THM-2043's complete
   period-14 parity-Hasse packet can be held fixed while a strict exit moves to
   denominator `113`. For `7^k u_k==1 (mod 113)`, put

   ```text
   w_k=24+84*7^k u_k.
   ```

   Then `w_k==24 (mod 14)`, its first `k` lift digits agree with the boundary
   row, and `w_k==108 (mod 113)`. At `t=47/113`, every row `S_(w_k)` has minimum
   `9/113>1/14`. The base row `S_108` is in the same strict THM-2055 one-tail
   owner sector as `S_24` and still fails the determinant gate. The non-hull
   representative `r=5` restores the exit because `113|w_k+5` and
   `D_113(1,...,11,-5,13)=9/113`.

The surviving synthesis is exact: the signed hull is the determinant sidecar,
the labelled primitive phase packet is the residue sidecar, and (22) counts the
hidden rows left by their quotient. A zero count closes a cell and a singleton
is rigid; a larger count is only an atlas-priority statistic. Pair-sum,
relative-Fejer, or owner-labelled Euler certificates are still required to
discharge the surviving rows. QED.

## 8. Exact referee

The stored referee checks 720 set-level packet decompositions and 720 Mobius
inversions, the exact Beatty shift law including the zero-measure `S_24`
singleton regression, unit transport with labelled residue vectors, the
primitive orbit--stabilizer norm, the CRT and spanning-tree no-gos, all six
one-tail deck values, six deck/fan interval counts, and six depths of the
same-sector Hasse family. It
passes both

```text
python 04-computation/lrc14_primitive_phase_packet_referee_codex_20260721.py
python -O 04-computation/lrc14_primitive_phase_packet_referee_codex_20260721.py
```

and byte-matches the stored output. The SHA-256 hashes are frozen above.
