# Averaging pairings repairs an energy split, not its named guards

## A universal four-child identity

For four vectors `x0,...,x3` in a real or complex Hilbert space put
`P=||Σxi||²`, `C=Σ||xi||²`. The three pairings of four children have middle
energies `Mχ`, indexed by nonzero characters of `F2²`. Writing
`Wχ=Σ_i (−1)^(χ·i)xi`, the parallelogram identity gives

    2Mχ=P+||Wχ||²,
    Σ_(χ≠0)||Wχ||²=4C−P,
    average(Mχ)=(P+2C)/3.

This proves the advertised exact repair without selecting a supposedly
always-best pairing. Apply it in the direct sum over named terminal ports:
charges at different ports cannot first be cancelled against one another.
Averaging intermediate energies preserves every source and the terminal
sum, but does not preserve every chronological question about a split.

## The integral reconstruction obstruction is exactly parity

Write the six pair sums in edge order `01,02,03,12,13,23` as `a,b,c,d,e,f`.
They reconstruct the children as

    u=(a+b−d)/2,       (x0,x1,x2,x3)=(u,a−u,b−u,c−u).

Over the integers their image consists exactly of

    a+f=b+e=c+d,       a+b−d even.

Over Eisenstein integers the parity condition holds on both integral
coordinates. These conditions are sufficient by the displayed formula.
The incidence matrix has Smith factors `1,1,1,2`: it has minors of sizes
one, two and three equal to one; every full minor is even since the columns
sum to zero modulo two; and a full minor has determinant two. Thus modulo
two a common shift is invisible, whereas in every other characteristic
the map is injective. Its Gram matrix `2I+J` is singular modulo three,
but the original map is still injective there. A Gram degeneracy must not
be confused with a kernel of the observation map.

## Named bits are an additional requirement

On the full cube `F2^d`, suppose observed linear bits have row span `W`
and future named guards have span `G`. Every required guard is recoverable
exactly when `G⊆W`, and the minimum number of additional independent
linear observations is `dim(W+G)−dim(W)`. The proof is kernel duality:
a bit constant on observation fibres must annihilate their common kernel.
On a restricted domain this criterion need not be necessary. For example,
on `{00,11}`, observing the second bit determines the first despite their
linearly independent coefficient rows.

For actual ROOT seeds `33,41,49,57`, the branch bits run through all four
patterns. Their successors are `25,31,9,43`, with clocks `2,2,4,2`.
Keeping only the second bit collides `33` and `41`; keeping only the xor
collides `33` and `57`. Their successors differ modulo four. Hence an
alternative pairing can be a perfectly valid energy tree while failing
to retain the bits needed by a later native chronological operation.

Positive-part billing also depends on order. For scalar children
`(−2,−1,0,1)`, all three binary trees have signed total energy increment
`−2`; their sums of positive node increments are respectively `0,4,2`
for the pairings `02|13`, `01|23`, `03|12`. Thus equality of signed totals
is not equality of paid invoices. Even marked energy observations do not
reconstruct arbitrary charges: moving a unit charge between leaf positions
can leave all relevant squared norms unchanged.

The checker verifies these identities and hostile examples exactly. The
full-cube guard theorem is algebraic, not a claim that every native state
space is a full independent bit cube.
