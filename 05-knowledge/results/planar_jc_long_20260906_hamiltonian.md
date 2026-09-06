# Depth-lowering Hamiltonians and actual completed source automorphisms

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENT ANALYTIC AUDIT.**
For every polynomial `R(p,y)` in characteristic zero, the Hamiltonian
`S=c+p^2(p^3-y^2)R` lowers each positive source depth by one and increases
the cusp divisibility of depth-zero functions. Consequently its exponential
at every scalar parameter defines an automorphism of the actual `D`-adic
completion of `B_2`, with inverse obtained by negating the parameter.
It preserves the completed depth filtration, the completed affine source
family, and the regular source two-form. These are completed automorphisms;
no polynomial termination or Keller pair is asserted.

## 1. Inheritance and the topology sidecar

The supplier is the independently audited universal carrier theorem in
[planar_jc48_sep06_hamiltonian.md](planar_jc48_sep06_hamiltonian.md),
read first from incoming commit `d0208173b2`. It proves precisely

```
{-u/2+H,S} in K[p,y] for every H in K[p,y]
iff S in K+p^2(p^3-y^2)K[p,y],
```

and proves that every nonconstant fixed-source carrier is not locally
nilpotent on the polynomial source. Those statements are inherited.
The new point is depth lowering and integration at an actual scalar time
in a specified complete ring, not merely a formal series in an additional
time variable.

The actual depth modules and logarithmic chart are supplied by
[THM-3989 / cusp-log-laurent-conductor-and-nondividing-depth-reduction](../../01-canon/theorems/THM-3989-cusp-log-laurent-conductor-and-nondividing-depth-reduction.md)
and [THM-4308 / source-normal-bracket-hasse-truncation-through-row-eight](../../01-canon/theorems/THM-4308-source-normal-bracket-hasse-truncation-through-row-eight.md).
The global regular source form is supplied by
[THM-3973 / exact-volume-simple-cubic-determinantal-affine-plane-completion](../../01-canon/theorems/THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion.md).
The index `2` in `B_2` is the completion height; `Spec B_2` is a surface.
It is not a change in the ambient dimension of the Jacobian conjecture.

The canonical hostile is the earlier five-term `S0`: its finite transport
does not retain the affine source form at all orders. The smaller
fixed-source control `S=D`, valid at `H=0`, already fails the universal
carrier. The corrected near miss is to infer a polynomial flow from a
Hamiltonian derivation, or to identify two different completions. The
least-used sidecar is the uniform modulus estimate on the complete ring.

The live board has five entries: actual depth generators; cusp divisibility;
quotient-by-quotient exponential; logarithmic source volume; and polynomial
specialization. Targeted searches for depth-lowering Hamiltonians,
`D`-adic exponentials, and their actual source modules found the inherited
time-formal carrier but no version of the completed consumer below.
No literature priority claim is made.

## 2. Every generator image, with signs fixed

Work over a characteristic-zero field `K`, with

```
u=x^2 t, p=t(1+x^2 t), y=x t p,
B_2=K[x,u,p,y] subset K[x,t], A=K[p,y], D=p^3-y^2,
P_d=sum_(a+b<=d) x^a u^b A,
{f,g}=f_x g_t-f_t g_x,
delta(f)={f,S},               S=c+p^2 D R(p,y).
```

The symbol `D` here is a function. It is distinct from the boundary divisor
also denoted `D` in THM-3973. Define

```
A_R=(2D+3p^3)R+pD R_p,
B_R=p(-2yR+D R_y),
C_R=4yR+2pyR_p+3p^3R_y,
E_R=10R+2pR_p+3yR_y.
```

The complete polynomial generator formulas are

```
delta p=-D B_R,                  delta y=D A_R,
delta D=-D^2 C_R,
delta u=p^3 y E_R,
delta x=p(5p^3+4y^2)R
        +p^2(p^3+y^2)R_p+py(2p^3+y^2)R_y.             (1)
```

In particular the four source generator images lie in `A`, the first two
in `D A`. Thus `delta` restricts to a derivation of the actual ring `B_2`.

There is a stronger uniform factorization **in `B_2`**. The relations
`p^3=D(1+u)` and `y^2=Du` turn the last two equations into

```
delta x=D[p(5+9u)R+p^2(1+2u)R_p+py(2+3u)R_y],
delta u=D(1+u)y E_R.                                  (1a)
```

The bracketed coefficient has depth at most one. This factorization
retains the actual source ring; divisibility in `A` and in `B_2` must not
be confused.

**Derivation.** The inherited brackets are
`{p,y}=-D/p`, `{u,p}=2yp^2/D`, `{u,y}=3py^2/D`.
Since `S_p=p A_R` and `S_y=p B_R`, the first two equations follow.
Differentiate `D=p^3-y^2` to obtain the third. The last two follow either
by literal source differentiation, or from `u=y^2/D`, `x=yp/D` and the
quotient rule. All denominators cancel exactly as displayed. This also
fixes the sign of `delta D`.

## 3. Depth lowering and a sharp iterated estimate

The relations

```
Dx=yp,                   Du=y^2                         (2)
```

show `D P_d subset P_(d-1)` for `d>=1`. Moreover (1) gives
`delta A subset D A` and `delta(D^k A) subset D^(k+1)A` for `k>=0`.
For a generator `x^a u^b f(p,y)`, differentiating an `x` or `u` factor
lowers `a+b`; differentiating `f` supplies a factor `D`, which (2)
absorbs into one positive-depth factor. If `a+b=0`, the derivative is
already in `A`. Therefore

```
delta P_d subset P_(d-1)                    (d>=1),
delta P_0 subset D P_0,
delta(D^k P_0) subset D^(k+1)P_0            (k>=0),
delta^(d+k) P_d subset D^k P_0              (d,k>=0).   (3)
```

Independently, (1a) and the first two equations of (1) give the uniform
containments

```
delta P_d subset D P_d,
delta(D^N P_d) subset D^(N+1)P_d,
delta(D^N B_2) subset D^(N+1)B_2               (d,N>=0). (3a)
```

Indeed differentiating a positive-depth generator removes one factor and
replaces it by `D` times a factor of depth at most one. Differentiating
its coefficient in `A` also supplies `D`. The Leibniz rule and
`delta D in D^2 A` give the remaining two claims. In particular
`delta^N(B_2) subset D^N B_2` uniformly in the input.

These are statements about the full actual modules, not a restriction to
one displayed monomial basis or a finite source jet.

There is no uniform extra power of `D` hidden in (3). For `R=1`,
`delta p=2pyD` and `delta D=-4yD^2` have exactly the displayed orders in
`A`. The logarithmic leading-symbol calculation in Section 5 also shows
`delta(x^d)` has depth exactly `d-1` for every `d>=1`.

## 4. Scalar-time exponentials on the full D-adic completion

Let

```
widehat B=lim_N B_2/D^N B_2.
```

The ring `B_2` is a finitely generated domain, hence Noetherian. The natural
map into this completion is injective: in the source chart
`D=t^3(1+x^2t)^2`, so a nonzero polynomial cannot belong to every
`D^N B_2`. The uniform estimate (3a) gives

```
delta(D^N B_2) subset D^(N+1) B_2.                    (4)
```

Consequently the induced derivation on the entire quotient obeys

```
delta_N^N=0 on B_2/D^N B_2.                           (5)
```

For every scalar `a in K`, define there

```
E_(a,N)(f)=sum_(j=0)^(N-1) a^j delta_N^j(f)/j!.
```

This is a uniform finite sum. Leibniz's formula makes it an algebra
homomorphism, and the finite binomial cancellation gives
`E_(a,N)E_(b,N)=E_(a+b,N)`. In particular its inverse is `E_(-a,N)`.
The quotient maps intertwine all these operators. Taking inverse limits
therefore gives **continuous inverse automorphisms**

```
E_a:widehat B -> widehat B,
E_a E_b=E_(a+b),                  E_a^(-1)=E_(-a).      (6)
```

This proves the convergence assertion for arbitrary elements of the full
completion, not just polynomials. The derivation extends continuously
with `delta(D^N widehat B) subset D^(N+1)widehat B`, so the tail beginning
at index `N` lies in `D^N widehat B` for every input simultaneously.
Modulo `D^N`, one may choose any polynomial representative and apply the
same finite sum (5). Changing representatives does not change the answer
by (4). Every `E_a` preserves each kernel of reduction modulo `D^N`;
equivalently it preserves `D^N widehat B`. Moreover `E_a` is the identity
modulo `D`, and `E_a-id` raises `D`-adic order by at least one.

There is a useful unit control. Iterating the third equation of (1) gives

```
delta^j D in D^(j+1)A,
E_a(D)=D U_a,
U_a=1+sum_(j>=1) a^j(delta^j D/D)/j! in 1+D widehat B.
```

The geometric series makes `U_a` a unit in the `D`-adic completion, in
agreement with exact preservation of the `D`-adic ideals.

Define the completed depth filtration intrinsically by closures

```
overline P_d=closure of P_d in widehat B.
```

Every finite partial sum of `E_a(f)` for `f in P_d` lies in `P_d`, and
its difference from `f` lies in `P_(d-1)` when `d>=1`. Continuity and the
inverse in (6) imply

```
E_a(overline P_d)=overline P_d,
(E_a-id)(overline P_d) subset overline P_(d-1) (d>=1). (7)
```

Thus the action fixes the associated graded of this completed filtration.
It also preserves the completed affine source family:

```
E_a(-u/2+overline P_0)=-u/2+overline P_0,              (8)
```

because `delta u in A` and every subsequent derivative remains in `A`.
The updated source coefficient in (8) is completed; it is not asserted
to be a polynomial `H(p,y)`.

## 5. A separate logarithmic convergence and volume statement

There is first a direct statement in the literal source chart:

```
S-c=t^5(1+x^2t)^4 R(t(1+x^2t),x t^2(1+x^2t)),
S_t in t^4 K[x,t],                 S_x in t^6 K[x,t].
```

The coefficient of `t^5` is the scalar `R(0,0)`, which explains the
second bound even when that coefficient is nonzero. Consequently
`delta(t^j K[x][[t]]) subset t^(j+4)K[x][[t]]` for `j>=0`. The same
scalar exponential therefore defines a formal source automorphism of
`K[x][[t]]`, with `E_a(x)=x+O(t^4)` and `E_a(t)=t+O(t^6)`.
The inclusion `D^N B_2 subset t^(3N)K[x,t]` induces a continuous map from
`widehat B` to this source completion, and the two flows intertwine under
it. No injectivity or equality of these complete rings is asserted.

Use the exact THM-3989 chart

```
s=y/p=xt, tau=D/p^2=t,
p=s^2+tau, y=s(s^2+tau), x=s/tau, u=s^2/tau,
omega=dx wedge dt=ds wedge dtau/tau.
```

Here the Hamiltonian is

```
S-c=tau F(s,tau),
F=(s^2+tau)^4 R(s^2+tau,s(s^2+tau)).
```

The bracket is `tau(f_s g_tau-f_tau g_s)`. Consequently

```
delta s=tau(F+tau F_tau),
delta tau=-tau^2 F_s,
delta(tau^j K[s][[tau]]) subset tau^(j+1)K[s][[tau]]
                                                    for every integer j. (9)
```

Thus `exp(a delta)` is also a well-defined continuous automorphism of
`K[s]((tau))`, preserving `K[s][[tau]]`, with

```
E_a(s)=s+O(tau),
E_a(tau)=tau(1+O(tau)).                               (10)
```

The parenthesized factor is a unit. The two exponentials agree term by
term on every finite-depth input from `B_2`. For such an input the
Laurent series keeps its exact negative leading coefficient and pole
depth: each positive iterate has strictly greater `tau` valuation.
This is stronger than mere bounded pole depth for these inputs.

For `R=1`, the leading action on `s^a tau^j` is

```
delta(s^a tau^j)=(a-8j)s^(a+7)tau^(j+1)
                  +higher tau powers.               (11)
```

In particular `delta(x^d)` has leading coefficient `9d s^(d+7)` at
power `tau^(-d+1)`. This proves the sharp depth assertion after (3).

The source form is preserved exactly. In the logarithmic chart its Lie
derivative vanishes because

```
partial_s(delta s/tau)+partial_tau(delta tau/tau)
 =partial_s(F+tau F_tau)-partial_tau(tau F_s)=0.
```

Equivalently in the source chart `i_delta omega=dS`, so
`Lie_delta omega=d(dS)=0`. The vector field is regular on `Spec B_2` by
(1), and `omega` is the globally regular two-form of THM-3973. The identity
on the dense source chart is therefore an identity of global regular
forms. It persists under every automorphism (6) in the completed module
of regular differentials. To justify passage to completion explicitly,
`d(D^N f)=N D^(N-1)f dD+D^N df` loses at most one order. Thus the
convergence estimates above also give convergence after differentiation;
the formal-time pullback identity may be evaluated at any scalar `a`.

In logarithmic coordinates the same statement is

```
d E_a(s) wedge d E_a(tau)/E_a(tau)=ds wedge dtau/tau.
```

One must retain the denominator: this is not the assertion that the
ordinary `(s,tau)` coordinate Jacobian equals one.

## 6. Topology and polynomial-scope hostiles

The two complete rings must not be identified. For every `N>=1`,

```
D^N x^(2N) in D^N B_2,
D^N x^(2N)=s^(2N)(s^2+tau)^(2N)tau^(-N).
```

This sequence tends to zero `D`-adically but has logarithmic valuation
`-N`. Hence the embedding of `B_2` in the Laurent chart is not continuous
for the full `D`-adic topology, and it does not automatically extend to
a map from the entire `widehat B` to `K[s]((tau))`. Claims (7) and (9)
are stated separately for exactly this reason. No arbitrary infinite
element of `widehat B` is assigned an unproved Laurent pole depth.

The fixed-source positive control `S=D,H=0` has
`{-u/2,D}=-3py`, but `delta p=2yD/p` is not in `A`. Thus it cannot replace
the universal carrier in the depth-zero assertion. The earlier `S0`
also remains outside that carrier. These are actual generator controls.

For nonzero `R`, the inherited factorial-closure argument shows that
`delta` is not locally nilpotent on the polynomial source. This forbids
a polynomial additive-group action having this infinitesimal generator.
It is fully compatible with (6), whose finite quotients do have such
actions. Non-local-nilpotence alone does **not** prove that every individual
scalar specialization of this completed exponential is nonpolynomial;
that stronger claim is not made. Nor does (6) supply a polynomial `H`,
termination of source coordinates, or a planar Keller pair.

## 7. Connection and next obligation

The source of the new connection is the universal affine-source carrier;
the target is the complete actual depth module and its `D`-adic quotients.
The map is Hamiltonian differentiation, followed by quotientwise
exponentiation. It retains depth, common source form in completion, and
the differential. It forgets polynomial termination. The needed sidecar
is the chosen topology and, for a Laurent interpretation, finite starting
depth. The cheapest positive check is (1), and the topology hostile in
Section 6 is decisive against a stronger completion identification.

This supplies an all-order completed repair operation carrying the actual
depth and differential conditions together. The remaining useful test is
whether a required finite source response lies in this carrier's image;
the earlier finite transport does not acquire that property automatically.
Even a positive response test would leave the polynomial specialization
and termination obligations open.

## 8. Exact controls

The standalone [source](../../04-computation/planar_jc_long_20260906_hamiltonian.py)
imports no repository producer. The [frozen output](planar_jc_long_20260906_hamiltonian.out)
passes **528 always-active exact gates**, with byte-identical normal and
optimized output. Universal statements rest on the displayed proofs;
the finite universe is:

- Eight carriers `R=1,p,y,p^2,py,y^2,D,1+p+y+py`, with every generator image
  checked by literal `(x,t)` differentiation, both uniform `D`
  factorizations, the source valuation bounds, and divergence zero.
- The complete 120-element bank given by `R in {1,p+y}`, `a+b<=3`,
  `c+e<=2`, and `x^a u^b p^c y^e`, with actual lower-depth representatives
  and their Laurent valuations checked separately.
- Two carriers, four core iterate seeds, and five derivative levels,
  checking divisibility in `K[p,y]` itself.
- Scalar flows for `R=1,p` through logarithmic order six, with composition
  inverse, parameter-two group law, Hamiltonian invariant, full
  logarithmic volume equation, and unit coordinate checked exactly.
- Sharp depth symbols for `x^d`, `1<=d<=6`, and the topology hostile for
  `1<=N<=7`, together with the actual fixed-source control `S=D`.

```bash
python3 -B 04-computation/planar_jc_long_20260906_hamiltonian.py > /tmp/hamiltonian_normal.out
python3 -B -O 04-computation/planar_jc_long_20260906_hamiltonian.py > /tmp/hamiltonian_optimized.out
cmp /tmp/hamiltonian_normal.out /tmp/hamiltonian_optimized.out
cmp /tmp/hamiltonian_normal.out 05-knowledge/results/planar_jc_long_20260906_hamiltonian.out
```

Raw LF-byte SHA-256 manifest:

```
source 009d982dfdc4ae88dc5fd2888057e2616df5684ff352209e6c20d5e8c31b3894
output ff7d5d116428a6e652581091cc148e796e8884c192bb6aff6dca3c85f4b01f35
```

The [independent analytic audit](planar_jc_long_20260906_hamiltonian_analytic_audit.md)
accepts the all-`R` generator formulas, both filtration mechanisms,
uniform complete-ring convergence, inverse/group law, source form,
and the topology distinction. No theorem identifier has been allocated.
