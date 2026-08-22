# Degree six retunes the quartic wall to an even quintic escape

**Status:** VERIFIED-EXACT WILDCARD ARTIFACT. This is not a canonical theorem,
does not reserve an ID, and does not change the planar-JC ledger.

## Inheritance pass

- Closest proved mechanism: THM-3263,
  `01-canon/theorems/THM-3263-degree-seven-retuned-quartic-infinity-wall-and-odd-critical-resultant-inertia.md`,
  retunes the degree-seven coefficient after the primitive wall is imposed and
  produces a quartic reciprocal Newton edge.
- Corrected near miss: adding `q*x^6` does not directly tune THM-3263's
  `x^95` carry. After the old `s` retuning it first reopens `x^96`; `r` must
  be retuned as a function of `q` before `q` reaches `x^95`.
- Canonical hostile: THM-3068,
  `01-canon/theorems/THM-3068-c3-escape-inverse-pole-ledger-and-reciprocal-cofactor-shift.md`,
  shows why an exact reciprocal critical-root cycle is not yet polynomial-map
  inertia.
- Least-used relevant sidecar: THM-3064,
  `01-canon/theorems/THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary.md`,
  supplies the pointed branchwise cofactor datum that this scalar construction
  still lacks.

The live concept board was: reciprocal resultant jet; nested strict-transform
walls; residual critical locus; local inertia sign; branchwise Keller cofactor.

## Universal exact jet

Write

```text
V=a^2*x^16*(1+r1/x+r2/x^2+r3/x^3+r4/x^4+r5/x^5+...),
B=1+q*x^6+r*x^7+s*x^8+t*x^9
```

and retain THM-3263's invariants

```text
Kappa=4*r2-3*r1^2,
Lambda=4*r2-r1^2,
Rho=r1^4+16*r1*r3-16*r2^2-8*r1,

Sigma=
 128-128*r3+256*r2*r4-256*r3^2-64*r1^2*r4
 +128*r1*r2*r3-64*r2^3-32*r1^3*r3
 +48*r1^2*r2^2-12*r1^4*r2+r1^6.
```

The next normalized invariant is

```text
Theta=r1^3-4*r1*r2+8*r3+8.
```

On the primitive wall `t=a`, the coefficients of degrees 99 and 98 vanish,
and

```text
c97=-2*a^11*(a*Kappa+4*r1*s-8*r).                       (1)
```

The first strict transform is therefore

```text
s=(8*r-a*Kappa)/(4*r1).                                 (2)
```

Along (2), the dormant degree-six parameter appears one step earlier than the
quartic carry:

```text
c96=-a^11*(a*Rho+8*r*Lambda-32*q*r1)/r1.                (3)
```

When `Lambda!=0`, the second retuning is

```text
r(q)=(32*q*r1-a*Rho)/(8*Lambda),
s(q)=(8*r(q)-a*Kappa)/(4*r1).                            (4)
```

Now the next carry is linear in `q`:

```text
c95=-3*a^11*(a*Sigma+64*q*Theta)/(8*Lambda).             (5)
```

Thus, when `Theta!=0`, there is a unique third tuning

```text
q_*=-a*Sigma/(64*Theta),
r_*=r(q_*),
s_*=s(q_*).                                               (6)
```

At (6), the first surviving coefficient is

```text
c94=-a^12*Psi/(16*Theta^2).                              (7)
```

Here `Psi` is the exact 58-term integral polynomial in `r1,...,r5` extracted
by the companion; it has total degree 11 and canonical-expression digest
`6b1814756894c97cfe6a31a3bd4ae51bc1dc591879211ec5f67d65f2a960e860`.
Setting `q=0` in (4)--(5) exactly recovers THM-3263's `r_*` and quartic carry.

## Both accessory fields

All denominators and `Psi` are nonzero in both inherited characteristic-zero
accessory fields. The exact rational norms are:

| passport | `Norm(Theta)` | `Norm(Sigma)` | `Norm(Psi)` |
|---|---:|---:|---:|
| `(4,1,1,1)` | `-37200002176/5359375` | `-157881139404422797328384/28722900390625` | `-572298646716642082095739982919499776/301716116790771484375` |
| `(3,2,1,1)` | `249238609408/48234375` | `64559323245631953174528/28722900390625` | `-4402683590765317259172542063793800401125376/1697153156948089599609375` |

The tuned-parameter norms and exact residual digests are:

| passport | `Norm(q_*)` | `Norm(r_*)` | `Norm(s_*)` | residual digest |
|---|---:|---:|---:|---|
| `(4,1,1,1)` | `2107940627729338801/186000010880000000` | `-3018470836755969/1860000108800000` | `1568249062443/1162500068000` | `7926bec47cfbfad38a1de91dfeebf84b824da3c780f478240bf66bc6e9981d26` |
| `(3,2,1,1)` | `-8978749567656192625/9187932097216512` | `8285100359018853125/484519856689152` | `-4924346007190859375/161506618896384` | `631e99ab7332025f3a338fd15919afa90cd9e542e2ed3c314dc0a4b4867015ec` |

After exact division by the inherited degree-44 passport boundary, each
residual has degree 50. Independent good reductions give:

| passport | `(p,u)` | `(t_*,s_*,r_*,q_*)` | `deg H` | `gcd(B,g)` | `gcd(H,g)` | `gcd(H,H')` | `[x^50]H` |
|---|---:|---:|---:|---:|---:|---:|---:|
| `(4,1,1,1)` | `(113,85)` | `(85,1,43,54)` | `50` | `0` | `0` | `0` | `97` |
| `(3,2,1,1)` | `(101,64)` | `(89,75,80,100)` | `50` | `0` | `0` | `0` | `56` |

Here `g=S*T` is the reduced owner boundary. Therefore both
characteristic-zero residuals are squarefree and boundary-disjoint; 50
reduced critical points remain.

## Even quintic local inertia

Fix `(q,r,s)=(q_*,r_*,s_*)`, put

```text
delta=t-a,              epsilon=Gamma*t-2=2*delta/a,
w=1/x,                  Hhat=w^55*H(1/w).
```

The transverse leading derivative is still `-16*a^11`, while (7) is the
`w^5` coefficient. Hence the initial form is

```text
Hhat=-16*a^11*delta-a^12*Psi*w^5/(16*Theta^2)
     +O(delta^2,delta*w,delta*w^2,delta*w^3,delta*w^4,w^6),

w^5~-256*Theta^2*delta/(a*Psi)
   =-128*Theta^2*epsilon/Psi.                            (8)
```

Equation (8) has five Puiseux branches and one local 5-cycle. Its sign is
even. This extends the inherited cubic/odd-quartic pattern by one exact nested
wall, but only for escaping roots of the critical resultant.

## Connection contract and stopping boundary

```text
source:      THM-3263's retuned quartic critical-resultant wall;
map:         add q*x^6, retune s, retune r(q), then take the unique q_*;
target:      a nonzero quintic strict-transform coefficient and 5-cycle;
preserved:   reciprocal vanishing order, escaping multiplicity, local inertia;
destroyed:   inverse-sheet labels, affine/Keller realization, physical Jacobian;
sidecar:     a genuine second polynomial coordinate and pointed branchwise
             primitive-element cofactors with inverse-different units;
hostiles:    the reopened c96 term, THM-3068's disconnected Laurent model,
             and the 50 surviving critical points;
cheapest test: attempt any cofactor comparison only after a lawful polynomial
               inverse cover supplies branch labels.
```

The artifact does not imply a polynomial inverse cover and does not advance
JC(2). The even 5-cycle cannot be retyped as Jelonek inertia without the
missing cover and cofactor data.

## Reproduction and hashes

```text
python3 04-computation/jc_degree6_retuned_quintic_infinity_wall_wildcard_20260803.py
python3 -O 04-computation/jc_degree6_retuned_quintic_infinity_wall_wildcard_20260803.py
```

Both executions byte-match the frozen output.

```text
script_sha256 = 8b0c633366869aa97c207e9ba3535f9b393d5b8441b6cbda5419638401b7a1a0
output_sha256 = c7b08d089213a6b5737880a537af8c72c4febd9a061172ddf0b1989a285fd578
```

The companion pins 16 exact script/output artifacts, has no assertion node,
floating literal, randomness, or discovery cache, and derives the reciprocal
jet through `x^94` symbolically before evaluating either accessory field.
