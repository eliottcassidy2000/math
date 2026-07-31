# THM-2780 marked-`D3` torsor descent report

## Verdict

**READY: the canonical candidate at `cf898c7e34`, with its Smith-null
strengthening at `6500f6c74d`, passes independent hostile audit.  Promote to
`PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED` after installing
or explicitly routing the two sharpenings below.**

The current canonical claims were checked independently:

- the four weightings and their `S4/A4` stabilizers;
- the regular `V4` action and field-of-definition statement;
- local fixed-colouring equivalence with the zero Kummer row;
- the unique fixed root-line pair and mod-two null for each nonzero row;
- mod-two invertibility of every determinant-one/three frame;
- the orbit-symmetrized `2,2,2,2` / `1,1,1,3` stopping boundary;
- the local/global section distinction.

No false implication was found.  In particular, the canonical warning that
zero rows need not give a global marked section is correct.  The additions
in this report sharpen which tournament information is lost before orbit
symmetrization and prove that even the complete residue packet does not
identify the unramified torsor.

The proposed finite object has one genuine positive interpretation and two
sharp losses.

1. The four absolute `1,2,3` determinant colourings are
   `W(D3)`-equivariantly the four even sign states.  They are a regular
   `V4`-set and therefore the correct abstract fibre for the quartic Kummer
   torsor.
2. The directed tournament is chiral only with a retained chamber.  Switch
   `a12` and relabel `(a23 a13)(b13 b23)` to obtain the weighted converse.
   Hence the chamber-free weighted switching class is self-converse and
   loses the quartic sign character.
3. A realized Kummer plane descends quasi-etale exactly when every
   normalized divisor row is zero, but those rows do not determine the
   unramified torsor.  Two distinct simultaneous-`S3` standard unit planes
   on one four-torus have identical zero rows, identical common sum, and the
   same marked finite frame.

The exact minimal missing coordinate is the embedded `Q`-equivariant Kummer
class in `H^1_et(R_reg,mu_2)`, equivalently the actual squareclass plane and
its unit/class-group realization.  A marked `h` is either merely an abstract
fibre label or, if evaluated on the `d_i`, already a point upstairs and not
base data.

## Positive descent theorem

For the six unoriented root lines

```text
L={R(e_i-e_j),R(e_i+e_j)}
```

and the four even states

```text
Omega={delta in {+1,-1}^3:product(delta)=1},
```

put

```text
kappa_delta({l,m})=|det(alpha,beta,delta)|.
```

The exact classification is

```text
weight 3 <=> both lines lie in delta^perp;
weight 2 <=> the lines are orthogonal;
weight 1 otherwise.
```

The weight-three triangle recovers `delta`, so `delta -> kappa_delta` is
injective.  Determinant covariance makes it `W(D3)`-equivariant, and the
diagonal `V4` acts simply transitively.

For a product-oriented rank-two Kummer packet

```text
s_i^2=tau_i,  s_1s_2s_3=c,
```

the map

```text
kappa_delta -> e_1/4+(delta.s)/2
```

is `V4`-equivariant and reconstructs the quartic roots.  By Galois descent it
is an isomorphism from the twist of the colouring set to the quartic
four-root scheme.  Rank two prevents root collisions.

## Exact affine gate

At a height-one divisor `D`,

```text
nu_D=(v_D(tau_1),v_D(tau_2),v_D(tau_3)) mod 2
```

is an even word.  Tame inertia acts on colourings by the fixed-point-free
translation

```text
delta -> (-1)^nu_D delta.
```

Therefore the normalization is quasi-etale exactly when every row is
`000`.  The unordered family of four colourings survives a nonzero
translation, so static weights do not reveal affine placement of a `V4`
element.  THM-2769's `110` row is the sharp ramified control.

## Exact unramified twin hostile

On

```text
R=C[x_1^+-,x_2^+-,x_3^+-,y_1^+-,y_2^+-,y_3^+-]/
  (x_1x_2x_3-1,y_1y_2y_3-1),
```

with simultaneous `S3` index permutation, take

```text
W_x=<[x_1],[x_2],[x_3]>,
W_y=<[y_1],[y_2],[y_3]>.
```

Both are irreducible standard planes, all their classes are units, and all
divisor rows vanish.  But

```text
R*/R*2=<[x_1],[x_2],[y_1],[y_2]>=F2^4
```

and the two planes meet only in zero.  Thus the corresponding connected
`V4` covers are nonisomorphic over the identity of `R`, even modulo
`Aut(V4)`.  Setting `e_1=0` gives two explicit inverse quartics with the same
finite sidecars and different Kummer classes.

This hostile proves that even the **complete** divisor-row collection is a
ramification test, not a reconstruction coordinate.

## Exact controls

Scratch artifacts:

```text
.scratch/d3_torsor_descent/THM-2780-proof-candidate.md
.scratch/d3_torsor_descent/thm2780_weight_colouring_descent.py
.scratch/d3_torsor_descent/thm2780_weight_colouring_descent.out
```

The companion has `26` explicit `require` call sites and zero Python
`assert` nodes.  It checks:

- all `4` state colourings and their `1^9 2^3 3^3` spectra;
- all `24*4*15=1440` Weyl covariance gates;
- the regular/fixed-point-free `V4` translations;
- all `6!=720` relabelings for retained-chamber chirality;
- the explicit weighted switching self-converse witness;
- exact symbolic half-Hadamard quartic reconstruction;
- two distinct `S3`-stable standard planes in the four-torus unit lattice.

Normal and optimized executions byte-match the stored output.

LF-normalized hashes:

```text
proof   4c767c16bb61be226ded790d34f132a2a4b23c5d9820007bcdfa3ebb8657ba6b
script  456e87aa2411798ac6678d03cb8e69b4e741a944cb1f6d1e80c29d9d8d736385
output  c04b424af35ac8dccc48e518f851d85c4f370f30c71e8df8dc0459ed52d16279
```

Independent replay of the current canonical companion also passes normal,
optimized, and stored-output comparison.  Its exact audit coordinates are

```text
canonical theorem
  14e8ade6cb4f608731da34d62bf14f79a62be222490e195448ac1a267423f262
canonical script
  22cbd2e0097e9a48f5695deae2bbd881142cedafcd5048d8001d11bb1be599f1
canonical output
  2ab9c0ceca0993263f62c3637b58d13f0101af7c1c3d810b3ff8a1565642fe76
canonical AST
  32 require call sites; 0 assert nodes
```

## Promotion boundary

When promoting the existing `THM-2780` candidate:

- record this independent audit and promote the status;
- retain the current canonical computation/results artifacts and hashes;
- add the chamber-switching self-converse witness, or route it explicitly as
  the reason directed sign is absent from the chamber-free target;
- add the simultaneous-`S3` unit-plane twin, or route it explicitly as the
  sharp proof that rows plus finite sidecars do not reconstruct `W`;
- state "`all rows zero iff quasi-etale`" only for a **given** Kummer plane;
- do not call a marked state a global section of a connected torsor;
- do not claim chamber-free chirality or a sign-valued affine invariant;
- do not infer an `A4/S4` Keller exclusion, `JC(2)`, `DC(2)`, Graceful Tree,
  or `LRC(14)`.
