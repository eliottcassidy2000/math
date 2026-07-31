# THM-2779 merged-candidate independent hostile audit

**Verdict: PASS.  Three wording sharpenings were applied during audit; none
remain open.**

Audited candidate:

```text
01-canon/theorems/
THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate.md

Final post-repair LF SHA-256:
21ac7ec9b19b8ed0abe4763e0b7e13ebc1e5eb776168c8e0088540f29daabccb
```

The abstract all-prime addendum is correct.  The p=13 endpoint/odometer and
coefficient-square consequences are also correct at their stated
coefficient/algebraic scope.  No repair to a formula, count, quantifier, or
logical implication is required.

## 1. Independent controls

The fresh companion `audit.py` imports no THM-2779 code.  It:

- exhaustively transports the merged group law
  `(x,y,z)(x',y',z')=(x+x',y+y',z+z'-yx')` to the tertiary convention by
  `(A,B,C)=(-y,x,z)` for `p=2,3,5,7,11,13`;
- checks center, derived group, determinant commutator, faithful `p^2`
  action, normalized-frame fibres, and exact-degree stabilizer conjugacy
  classes for those six primes;
- exhausts every unit `V` in `F_p[epsilon]/epsilon^p` for `p=2,3,5`
  (`2`, `18`, and `2,500` units), proving the `p`-decoder fibre and the
  unique zero-section for every chosen `u`-coefficient;
- independently reconstructs `S_beta`, `K_beta`, all thirteen local chart
  gauges, the socle-qualified inverse, and its `+3N_13` full repair;
- checks the endpoint action, the exact group coordinates of its three
  generators, all `169` odometer conjugacies, the `156/13` successor split,
  and the class-`1`/class-`7` kernel renormalization;
- verifies the displayed dual-field Hadamard witnesses, their nonzero
  coordinates and Pluecker identity, and the exact
  `2,184*169^2=62,377,224` edge-product inference;
- checks the conditional flat-zero and six-grade scalar-line dilation gates.

It also records the useful hostile boundary omitted from the theorem because
it is not needed there:

```text
faithful transitive degree-p^2 H_p-sets:
  p+1 classes for odd p; 2 for p=2;

all transitive degree-p^2 H_p-sets:
  p+2 classes for odd p; 3 for p=2,
```

the extra class being the nonfaithful coset action with central stabilizer.
Thus the word **faithful** in the theorem's `p+1` refinement is load-bearing.

## 2. Mathematical audit

### All-prime Heisenberg theorem

The cocycle, inverse, commutator, center, and derived-subgroup formulas are
correct for every prime.  The same proof includes `p=2`, where the order
census is `1^1,2^5,4^2`, hence `H_2 isomorphic D_8`.

For an action on fewer than `p^2` points, every orbit has size `1` or `p`.
Every `p`-point image is a `p`-subgroup of `S_p`, whose order is at most
`p` because `v_p(p!)=1`; hence it is cyclic and kills the derived center.
The affine action

```text
(x,y,z)(r,w)=(r+x,w+z-yr)
```

is faithful and transitive on `p^2` points.  At exact degree a faithful
action must be transitive; its order-`p` stabilizer is core-free exactly
when noncentral.  The odd-prime projected-line classification and the
quadratic-refinement exception at `p=2` are correct.

### Decoder socle torsor

For `S=epsilon V` with `V(0)!=0`, `V` is a unit and

```text
Ann(epsilon)=F_p epsilon^(p-1)=F_p N_p.
```

Therefore all solutions of `SK=epsilon` are
`V^(-1)+F_pN_p`, multiplication by `S` has rank `p-1`, and every chosen
`u`-coefficient gives a unique zero-normalized section because all
coefficients of `N_p` equal one.  The `p=13` values, local chart direction,
uniform `-N_7` invoice, and `K_beta+3N_13=V_beta^(-1)` specialization all
check exactly.

### Endpoint-origin action and group-law conventions

For `R=ws+vt`, the operations

```text
X_end(v,w)=(v,w+v),
Y_end(v,w)=(v+1,w),
Z_end(v,w)=(v,w+1)
```

satisfy `[X_end,Y_end]=Z_end` and give a faithful `H_13` action.
Relative to the merged convention and its affine action, their exact
coordinates are

```text
X_end=rho(0,-1,0),
Y_end=rho(1,0,0),
Z_end=rho(0,0,1).
```

Relative to the tertiary companion convention the global isomorphism is
exactly

```text
(A,B,C)=(-y,x,z).
```

There is no sign conflict.  There is only a generator-name reuse between
Sections 4 and 5.

### Universal coefficient-square inference

The secondary companion checks all `28,392` nonzero source edges and all
`28,392` nonzero target edges in each certified field.  For a normalized
frame, both steps are nonzero.  Hence each of

```text
(P0+P1)(Q0+Q1), (P0+P1)(Q0-Q1),
(P0-P1)(Q0+Q1), (P0-P1)(Q0-Q1)
```

is nonzero for every origin pair, and the Pluecker identity is the
rank-one factorization identity.  This proves all `62,377,224` gates by
edge-product reduction; it is not a literal enumeration of that many
squares.  Nonzero finite-field image implies nonzero characteristic-zero
cyclotomic coefficient because the imported Lucas/order checks certify ring
homomorphisms.  The theorem correctly stops at common endpoint-array
coefficients and does not infer one Boolean physical ancestry.

### Odometer and THM-2657 normalization

Under `iota(v,w)=v+13w`,

```text
X_end(n)=14n,   Y_end=low-digit +1 without carry,   Z_end(n)=n+13.
```

Physical `+1` equals `Y_end` on `156` addresses and `Z_end Y_end` on the
`13` carry addresses.  This is exact.

The standard two-digit section has extension class `1`; THM-2657's fixed
root normalization has class `7`.  Rescaling the kernel generator by `7`
sends class `1` to class `7`; the inverse rescaling by `2` sends `7` to
`1`.  They are isomorphic after this unit choice, not equal with both fixed
kernel generators.

### Dilation gate

Under the theorem's explicit hypotheses

```text
D_target=id,   D_root=0,
```

a linear natural map obeys `0*pi=pi`, hence vanishes.  Retaining all six
root grades changes the target: dilation transports one digit copy to the
next, and a chain map is determined by one common scalar, giving thirteen
maps and twelve nonzero maps.  This is a conditional algebraic survivor,
not a constructed physical intertwiner.

## 3. Wording sharpenings applied

1. Near endpoint formula `(31)`, the generator coordinates now state
   `X_end=rho(0,-1,0)`, `Y_end=rho(1,0,0)`,
   `Z_end=rho(0,0,1)`.  Sections 4 and 5 currently reuse `X,Y` after a
   sign/swap; the relations are correct, but the identification is implicit.
2. “Multiplication by `7` identifies them” was replaced by the directed
   statement:
   “rescaling the kernel generator by `7` sends the standard class `1` to
   THM-2657's class `7` (inverse rescaling by `2`)”.  This prevents equality
   from being read with fixed kernel normalization.
3. “Over every faithful carrier address” was replaced by “over every carrier
   address”.
   Faithfulness belongs to the action, not to an individual address.

The theorem also preserves **faithful** in “minimal faithful transitive
sets”; deleting it would make the `p+1` count false by one.  There are no
remaining wording repairs.

## 4. Replay and scope

All canonical companion normal/optimized/stored triples byte-match:

```text
primary
  script  4c6a58c80ddd4be0fd9bdd297b310df054bbc08996eb223d519d3cce6b8ed13a
  output  f7c96259777a3ab4a5e46cac8666181ae77a3be2e440cee8785997507706791a

secondary (full dual-field endpoint reconstruction)
  script  004e06c617f9305e2f0bc30871926e3faa7843f47dcf63af1fd8a892e63101e4
  output  a89c00a3830ee9ff282cc5e4557d41293af5d6f0e7feabd5d3c7e7808591e754

tertiary
  script  ef6e9f9bcb4f11152d291342a11ae215245d1d19b96c49940a01ba9ea850cbd9
  output  1feb463864015035ab8d7fcfcddf9cfe8b0ec0a3ed36481f2f66d6a9149182e6

tracked second independent audit
  script  5019f87b24500a5a13825d3be01908ea983a08b360a384fd614107f476201f46
  output  1b9ad37b35e92a14dd90d0db8c1d0cf225761e2c37ca8e2fe2120bd0f64c47d4
```

The scratch reimplementation also matches under normal and optimized Python,
has no `assert` node, and has LF hashes:

```text
scratch script  a38d7d0113cf1d06b49b90a13e673f09441bad8f6598017a3b1b81662f22face
scratch output  4e6d89f08ab80cb398cc1012ba45eee85d42a6bd8ec3fef1a996b5b03c461618
```

No common physical ancestry, root-deck intertwiner, current equivariance,
semantic transition, row exclusion, or LRC(14) conclusion was inferred.
