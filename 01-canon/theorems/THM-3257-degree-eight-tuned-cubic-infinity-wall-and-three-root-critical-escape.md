---
id: THM-3257
title: "Degree-eight tuned cubic infinity wall and three-root critical escape"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In both THM-3237 cubic accessory families, enlarge
  the clutch to B_(s,t)=1+s*x^8+t*x^9.  The degree-eight coefficient cannot
  create the generic infinity face, but on the primitive degree-nine wall it
  linearly tunes the next strict-transform coefficient.  One exact value
  s=s_* cancels that coefficient in both families; a nonzero cubic carry
  remains, the saturated critical resultant drops from degree 55 to 52, and
  exactly three reciprocal roots escape with a cubic-root law and one C3
  local monodromy cycle.  Exact good reductions leave 52 squarefree,
  boundary-disjoint critical points.  This is a critical-resultant
  bifurcation, not a Keller map or a conclusion about JC(2).
source: root/creative-synthesis/2026-08-03
audit: >
  The exact companion pins the promoted THM-3237 companion; independently
  expands the universal two-parameter top jet; derives the response-identity
  jet relations, unique tuning coefficient, cubic carry, and local law;
  reconstructs both characteristic-zero accessory families; proves the
  required normalized jet factors nonzero by exact cubic norms; and obtains
  exact degree-52 residuals.  Separate good reductions at (113,85) and
  (101,64) certify squarefreeness and boundary-disjointness.  Normal and
  optimized executions match the frozen transcript; the source has no
  assertion node or floating literal.  An independent direct-resultant and
  algebraic-field implementation rederived every response jet, norm, exact
  quotient, degree and cubic law; separate finite-field reconstruction
  matched both good reductions and their gcd gates.
depends_on:
  - THM-3237-degree-nine-jacobian-infinity-wall-and-square-root-escape
related:
  - THM-3059-quartic-twojet-even-jelonek-c3-escape-counterexample
  - THM-3068-c3-escape-inverse-pole-ledger-and-reciprocal-cofactor-shift
  - THM-3233-fraction-free-exceptional-quotient-tail-renormalization
script: 04-computation/jc_degree8_tuned_cubic_infinity_wall_thm3257.py
output: 05-knowledge/results/jc_degree8_tuned_cubic_infinity_wall_thm3257.out
script_sha256: f7c9b6d92204ab0af0271311bfbcfa11fb51014a1e3d875d271ff406898da06d
output_sha256: 1ade6b505bf92f7cb1a17395a3fcff22ee9ac9692ea1ac449b2bee812831ab1b
hash_basis: LF-normalized bytes
---

# THM-3257 -- a dormant degree-eight term tunes a cubic infinity escape

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3237 proves that the monomial clutch

~~~text
B_t=1+t*x^9
~~~

first changes the infinity face at degree nine.  At its primitive wall, two
top coefficients vanish and a quadratic reciprocal initial form makes two
critical roots escape.  The next question is whether terms that were inert
for the generic face remain inert after restriction to that wall.

They do not.  The degree-eight term becomes a coordinate on the exceptional
divisor.  It can cancel the quadratic carry and expose a cubic one.

## 1. Two-parameter clutch and normalized jets

Retain either characteristic-zero cubic accessory field, response pair
(V,A), monic degree-44 passport boundary, and saturated critical resultant
of THM-3237.  Write

~~~text
V=v*x^16+v15*x^15+v14*x^14+v13*x^13+...,
A=a*x^8+a7*x^7+a6*x^6+a5*x^5+... .                    (1)
~~~

The inherited normalization and response identity give

~~~text
v=a^2,                         a=2/Gamma,
a7=a*v15/(2v),
a6=a*(4v*v14-v15^2)/(8v^2),
a5=a*(8v^2*v13-4v*v14*v15+v15^3)/(16v^3).              (2)
~~~

Put

~~~text
r1=v15/v,             r2=v14/v,             r3=v13/v,
Kappa=4r2-3r1^2,
Rho=r1^4+16r1*r3-16r2^2-8r1.                           (3)
~~~

Now use the two-parameter clutch

~~~text
B_(s,t)=1+s*x^8+t*x^9                                  (4)
~~~

in THM-3237's polynomial

~~~text
P_(s,t)(x,z)=(V(x)z^2+B_(s,t)(x)z)^2+A(x)z+x.           (5)
~~~

Let K_(s,t) be the universal saturated y-resultant before division by the
passport boundary and let H_(s,t) be the quotient.  The boundary is monic,
so whenever the higher coefficients vanish the first surviving coefficient
of K is also the leading coefficient of H.

## 2. The codimension-two cubic wall

Exact universal expansion of the top four layers gives

~~~text
[x^99]K=-16*t^4*v^3*(a*t-v).                            (6)
~~~

Thus the primitive degree-nine wall remains

~~~text
t=t_*=v/a=a=2/Gamma.                                   (7)
~~~

After substituting (2) and (7), the next layer vanishes for every s:

~~~text
[x^98]K_(s,t_*)=0.                                     (8)
~~~

The degree-eight coefficient first appears one layer later:

~~~text
[x^97]K_(s,t_*)=-2*a^11*(a*Kappa+4s*r1).               (9)
~~~

In both accessory fields r1 is nonzero.  Hence there is one and only one
tuning value

~~~text
s_*= -a*Kappa/(4r1).                                   (10)
~~~

At (s_*,t_*), the next layer is

~~~text
[x^96]K_(s_*,t_*)=-a^12*Rho/r1.                        (11)
~~~

The exact norms below prove that neither r1 nor Rho vanishes:

| passport | Norm(r1) | Norm(Kappa) | Norm(Rho) |
|---|---:|---:|---:|
| (4,1,1,1) | -1024/175 | 23328/30625 | 58451308942336/937890625 |
| (3,2,1,1) | -69632/525 | -2147586048/30625 | 122694886932611072/8441015625 |

For reference,

~~~text
Norm(s_*)=-250047/32768000             in (4,1,1,1),
Norm(s_*)=-173456171875/35651584       in (3,2,1,1).     (12)
~~~

Equations (6)--(11) therefore imply the exact degree ledger

~~~text
deg H_(s,t)=55              when t!=0,t_*;
deg H_(s,t_*)=53            when s!=s_*;
deg H_(s_*,t_*)=52.                                     (13)
~~~

The first line is a generic statement in this two-parameter family; special
lower-degree cancellations away from the displayed top face are not being
classified.

## 3. Cubic-root escape and local inertia

Fix s=s_*, set

~~~text
delta=t-t_*,
epsilon=Gamma*t-2=Gamma*delta,
w=1/x,
Hhat(t,w)=w^55 H_(s_*,t)(1/w).                          (14)
~~~

Differentiating (6) at the wall and using (2) gives

~~~text
partial_t [x^55]H_(s_*,t)|_(t=t_*)=-16a^11.             (15)
~~~

Together with (11), the completed local expansion is

~~~text
Hhat=
 -16a^11*delta
 -a^12*(Rho/r1)*w^3
 +O(delta^2,delta*w,delta*w^2,w^4).                     (16)
~~~

Its Newton initial form has one compact edge, and its three Puiseux branches
satisfy

~~~text
w^3 ~ -16*r1*delta/(a*Rho)
    = -8*r1*epsilon/Rho.                                (17)
~~~

Therefore exactly three units of critical-resultant intersection
multiplicity move to x=infinity.  After any complex embedding, a meridian
delta -> exp(2*pi*i)delta
cycles the three cube-root branches, so the local permutation on the
escaping critical roots is one 3-cycle.

This is the strict-transform mechanism in its cleanest form:

~~~text
degree-eight term inert on generic face
  -> live coordinate after two wall cancellations
  -> tuned third cancellation
  -> cubic reciprocal initial form.                    (18)
~~~

## 4. Squarefree residual controls

The exact characteristic-zero quotient has degree 52 in both accessory
fields.  Two good reductions give:

| passport | (p,u) | (t_*,s_*) | deg H | deg gcd(B,g) | deg gcd(H,g) | deg gcd(H,H') | [x^52]H |
|---|---:|---:|---:|---:|---:|---:|---:|
| (4,1,1,1) | (113,85) | (85,65) | 52 | 0 | 0 | 0 | 64 |
| (3,2,1,1) | (101,64) | (89,73) | 52 | 0 | 0 | 0 | 91 |

Here g=S*T is the inherited owner boundary.  The accessory reductions are
simple, all denominators are units, and the displayed degrees are preserved.
Thus the tuned characteristic-zero residuals are squarefree and disjoint
from the owner boundary.  As in THM-3237, the constant leading
y-coefficient of the first gradient equation converts them into 52
reduced transverse critical points.

## 5. Connection contract and scope

The typed connection to THM-3233's exceptional PRS tower is procedural:

~~~text
source:      top reciprocal jet of the critical resultant;
operation:   restrict to a wall, divide exceptional multiplicity, retest;
target:      the next nonzero normal coefficient and its Newton edge;
preserved:   vanishing order, reciprocal root multiplicity, local monodromy;
destroyed:   polynomial-map fibre labels, Keller cofactors, global chart;
sidecar:     branchwise cofactor units and a second polynomial coordinate.  (19)
~~~

The 3-cycle in Section 3 is even.  It resembles the ternary infinity
inertia in THM-3059, but it acts on roots of a **critical resultant**, not on
the sheets of a generically finite polynomial map.  It therefore supplies
neither the odd infinity component sought in the current Keller lane nor a
counterexample to it.

Most decisively, 52 Morse critical points survive at the tuned wall.
Consequently (5) has no constant-Jacobian mate in the current chart.
Nothing here constructs a Keller map, proves or refutes JC(2), or classifies
arbitrary nonlinear B,C_0,E_0 deformations.

The durable gain is a sharper search rule: a term below the first live
generic degree can become the controlling coordinate after restriction to
an exceptional wall.  Face-minimality must therefore be recomputed on each
strict transform, not inherited from the ambient family.

## 6. Exact reproduction

Run

~~~text
python3 04-computation/jc_degree8_tuned_cubic_infinity_wall_thm3257.py
python3 -O 04-computation/jc_degree8_tuned_cubic_infinity_wall_thm3257.py
~~~

and compare LF-normalized bytes with the declared output.  The companion
pins THM-3237's exact engine, performs a fresh symbolic top-jet derivation,
reconstructs both accessory cases, and checks the independent finite-field
controls without floating point, randomness, or optimization-sensitive
assertions.

QED.
