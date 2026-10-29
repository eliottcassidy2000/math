# Independent second-checkpoint audit

Status: **VERIFIED** exact certificates and proof review. Source remains pinned
to `5d6eb91dd76eee86b91e04d856daa32c32b3d61d`; producer additions are the
explicitly named Session97 local files. No producer modules are imported.

## Fruit three-isogeny and local reconstruction

`check_fruit_isogeny_independent.py` uses a separate left-to-right binary
multiplication routine to verify:

- `3H=phi(G)=(-225,210)` for `G=(-4,28)`, `H=(-196,196)`;
- ten instances of `phi(nG)=3nH`;
- `35H` lies strictly in the positive real cubic chamber;
- its rational orientation discriminant is a square;
- its q has even numerator and odd denominator, reducing the reconstruction
  cubic to `t³+t²+1` modulo 2. This has no root, so the monic polynomial has
  no Q2 root, and in particular no rational lift.

The all-order claim `nH in phi(E(Q)) iff 3|n` follows from the quotient-group
argument: `3H` is in the image, while H is not by the q=0 mod2 obstruction.
The class of H consequently has exact order three. No full cokernel-size,
rank, or Sha assertion is required.

**Independent removal of a citation dependency.** Doubling G gives
`x(2G)=676/49`. On the integral short Weierstrass model obtained by
`X=9x+327`, `Y=27y`, the equation is

    Y²=X³−302643X+63998478,

and `X(2G)=22107/49` is not integral. Nagell–Lutz therefore proves that G
has infinite order. The finite-kernel isogeny and `3H=phi(G)` imply H has
infinite order. The point 6H lies in the identity component of E'(R), a
circle; its infinite cyclic subgroup is dense. The coset `35H+6ZH` thus
hits the open positive chamber infinitely often. Every index is 2 modulo 3,
so none lifts rationally. This is a valid infinite-family conclusion; it
is not a statement about effective hitting time or Sha.

## Genus-three nonmatching kernels and Fano geometry

`check_genus_three_independent.py` directly enumerates isotropic 3-spaces in
the six-dimensional quotient of even subsets of eight points by complement.
It obtains 135 Lagrangians: 105 contain a pair class; 30 do not and lift to
doubly-even extended Hamming codes. Puncturing at the named eighth point
produces seven triples with every pair on exactly one line, verified for
all thirty codes. The permutation fixing that point and cycling the other
seven fixes exactly two nonmatching Lagrangians. Affine maps of F7 preserving
either Fano structure number 21, with multipliers exactly {1,2,4}.

The arithmetic witness `x^7+x+1` over F5 was independently proved irreducible
by ruling out all 155 monic factors of degrees one, two and three. Every
proper factorization of a degree-seven polynomial has such a factor. This
uses a different certificate from the producer's Frobenius-power test.
The hyperelliptic curve therefore has Frobenius type (1,7) on its branch
points, giving an actual arithmetic source for the two stable kernels.
No assertion of supersingularity, target hyperellipticity, or tournament
Hamilton-path exclusion follows.

## Remaining scope

All independent scripts use explicit exceptions, and normal/-O outputs agree.
The ranks certificate is separately covered by
`INDEPENDENCE_CERTIFICATE_AUDIT.md`. The first-return controller and Hasse
point-count scripts were frozen for the preceding publication checkpoint.
No external theorem about all isogeny target types or global rank distribution
has been certified by these finite computations.
