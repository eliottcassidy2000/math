# Independent audit of the sparse-transport raw-roof formula

**Status: PASS, analytic audit.** Auditor: moments lane, 2026-09-05.
Scope: Proposition 4 of
[the LRC sparse-transport report](synthesis_20260905_lrc_sparse_transport.md),
read against the literal sheet definitions in THM-4409 and the raw-address
and owner-gate proofs in THM-4386. This is an independent proof reading,
not an additional claimed finite census. Root already replayed the script.

The factor three and clipping convention are correct under the stated
primitive ternary-unit assumptions. The following explicit formulas make
the two delicate steps reviewable.

For a raw component choose `y` in `(0,1)`, with nearest integers `n_i` and
`|w_i*y-n_i|<3/14`. Write `o_i=-w_i^(-1)*n_i mod 3`. The raw gate makes
these three residues distinct. For each `j=0,1,2`, set

```text
x=(y+j)/3,  pi_j(i)=o_i-j mod 3.
```

There is a unique integer `k_i` with

```text
3*k_i=n_i+w_i*(j+pi_j(i)).
```

Thus `|w_i*(x+pi_j(i)/3)-k_i|<1/14`. Conversely every physical component
with a sheet permutation arises in this way, using the unique third of
`[0,1)` containing it. At a cut `x=j/3`, distinct sheet teeth cannot have
a common point: only sheet `-j mod 3` contains that cut. Hence a pair
component cannot cross these cuts, and no component is split by this
correspondence. Each raw edge contributes three physical edges with full
lengths divided by three. Their edge minima sum to the raw `K_ij`.

The same three preimages realize the cyclic permutations of the selected
owner ordering. The opposite owner ordering is carried by other raw
carriers, including the reflected `-C`; the sum already includes them.
There is no additional factor two. Primitivity is necessary for the cited
one-to-one raw-address theorem and has not been dropped.

Only a third-sheet tooth that crosses the physical boundary zero can be
clipped in a potentially relevant edge minimum. It is centered at zero
and has full length `2h`, where `h=lambda/w_k`. Let the contacting pair
component be `I`, and choose its nonzero-sheet constituent of speed `a`.
On the positive side of zero,

```text
left(I)>=(1/3-lambda)/a,
|I|<=2lambda/a,
h>left(I)                  (positive-length contact),
(1/3-lambda)/(2lambda)=11/6>1.
```

Therefore `h>|I|`, so
`min(|I|,h)=|I|=min(|I|,2h)`. The negative side follows by reflection.
This proves that clipping preserves every relevant minimum, even when the
third tooth itself crosses zero. The published report's shorter explanation
is valid; these inequalities are an optional clarification.

Finally, a contact edge has a unique positive-length triple intersection,
and interval convexity gives the raw overlap roof. The gate `L(C)>0` is
essential: it ensures all pairwise overlaps and prevents noncontact raw
vectors from entering the one-pair sum. Thus the network-to-roof identity
does not replace the true carrier set by a weaker pairwise support.

I also checked the interval-star proof, the shortest-parent intersection
argument, and the bounded-weight proof: their inequalities and strict
endpoint convention are consistent. No mathematical repair was found.
