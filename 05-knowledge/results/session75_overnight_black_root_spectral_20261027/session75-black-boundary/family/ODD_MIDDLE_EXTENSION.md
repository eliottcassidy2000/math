# Odd middle-vertex extension: exact counterexamples at every order n>=10

**Session75, 2026-10-27. PROOF CANDIDATE for independent audit.**
This supplements `ROOTED_CONVERSE_COUNTEREXAMPLE.md`; the definitions,
rooted-orbit mechanism, and even-family proof there are unchanged.

## 1. Construction and claim

For `a>=1`, use the same rooted `U=C3(T_a,C3,{r})`, its converse copy `V`,
and `m=a+4`. Introduce one additional vertex `w`. Start with

`V -> w -> U`, with every `V -> U` arc,

then reverse exactly two arcs:

`r_U -> w`, `w -> r_V`.

In particular there is **no direct bridge U->V**: all direct interblock arcs
remain `V->U`. Call the result `J_a^odd`, of order `2a+9=2m+1`.

**Claim.** It is strong, has minimum semidegree2, and has exactly two minimum
outgoing shores, `U` and `U union {w}`, each with cut value1. Its black distance
to nonstrong is exactly2. Its two lower-bound cones coincide and consist of
exactly nine blue Hamilton words.

Together with the even-family theorem, this supplies counterexamples to
`d_black(nonstrong)=lambda` at **every order n>=10**. Every example has
`lambda=1`, `delta=2`, and exact black nonstrong distance2. This does not prove
or refute the separate all-order core-distance proposal.

## 2. Strongness, scores, and the complete minimum-cut family

Strongness follows from strong `U,V` and the directed cycle through the three
parts `r_U -> w -> r_V -> U`. The root-root pair of the blocks has orientation
`r_V -> r_U`, as do all other direct `V -> U` pairs.

The full scores on `U` remain exactly the even-family scores: those of `U`
internally plus1 at its root. They have minimum2. Vertex `w` has score `m`.
The involution swapping corresponding `U,V` vertices and fixing `w` is an
anti-automorphism. Thus the `V` full scores are the complements of the `U`
scores around `2m`; minimum semidegree is2.

The shores `U` and `U+w` each have exactly one outgoing arc, respectively
`r_U -> w` and `w -> r_V`. To exclude every other cut:

- A shore with nonempty proper intersections with both strong blocks has at
  least one internal outgoing arc in each, hence cut value at least2.
- A proper nonempty `X subset U` with `w` absent has outgoing count
  `c_U(X)+1_(r in X)`, at least2 by the rooted cut margin and strongness proved
  in the even-family packet.
- A shore `X+w` with nonempty proper `X subset U` has an internal `U` outgoing
  arc and the arc `w -> r_V`, hence at least2.
- A nonempty subset of `V`, with or without `w`, has at least `m>=5` outgoing
  arcs into `U` before any other arcs are counted.
- The cut of `{w}` and that of its complement both have value `m`.
- Taking the complement after the block-exchange anti-automorphism preserves
  directed cut value. It reduces shores containing all of one block and
  part of the other to the one-block cases just checked.

These cases exhaust the nonempty proper shores. The only cuts of value1 are
therefore the two stated shores; no arbitrary-shore filter is hidden.

## 3. Why both geodesic cones are the same nine blue words

A Hamilton word starting in `U` and ending outside `U` can leave `U` only via
`r_U -> w`. Net transition counting forces exactly one such departure and no
return, so it must first exhaust `U`, then visit `w`, then exhaust `V`.
In particular it cannot end at `w`, because the vertices of `V` would be
unvisited.

For the other minimum shore `U+w`, its only outgoing arc is `w -> r_V`.
A decreasing Hamilton word must exhaust `U+w` before crossing this arc;
its first part ends at `w`, so it exhausts `U` ending at `r_U`, then visits `w`.
Consequently the same form is forced:

`alpha_U, w, beta_V`,

with `alpha` ending at `r_U` and `beta` starting at `r_V`. Each has three
choices, giving exactly nine words. Their positional reversal fixes `w` and
anti-isomorphically exchanges the two blocks by the same root-fixing word-orbit
argument as in the even case. Its cross-pair equations, including the two
exceptional root-w arcs, hold exactly. Thus every one is blue.

The cut-current lower bound now proves that no black one-step route can reach
nonstrongness: any terminal zero cut would have to be one of these initial
minimum shores, and its complete decreasing cone is blue.

## 4. Two literal black steps

Keep labels `V=0,...,m-1`, `U=m,...,2m-1`, and put `w=2m`,
`r=a+3=m-1`. The first word is the even word with `w` inserted after its
whole `U` segment:

`P_odd=(a+1,a,a-1,...,0, m,m+1,...,2m-1, w,r,a+2)`.

It cuts the directed Hamilton cycle

`alpha_U, w, reverse(alpha)_V`

at the internal `V` arc `(a+2)->(a+1)`. Its endpoints are both `V` vertices
of full score `2m-2`; their score sum differs from `2m=n-1`, proving blackness.
Both initial minimum shores exclude both endpoints, so both cuts stay1.
All other cuts were at least2 and can fall by at most1. The intermediate
therefore remains strong.

After this first reversal the induced `U` and `V` tournaments are exactly
those in the even construction. The only direct `U->V` arc is now `m->0`.
Moreover all arcs between `w` and the two blocks have become uniform:

`V -> w -> U`.

The second word is the even second word with `w` appended:

`Q_odd=(2m-1,2m-2,...,m, 0,1,...,a-1,a,a+2,r,a+1,w)`.

All arcs before the final `w` are directed by the verification in the even
packet; the final arc is directed because every `V` vertex now beats `w`.
This second word is automatically black from its endpoint scores. Its first
vertex `r_U` still has full outdegree `a+1` (it was interior in the first
Hamilton word); its last vertex `w` still has score `m` for the same reason.
Their sum is

`a+1+m=2a+5=n-4`,

not `n-1`, as a blue positional reversal would require. It starts in `U` and
ends at `w`, so the outgoing cut of `U` falls from1 to0. The terminal is
nonstrong. This proves exact distance2 for every member of the odd family.

Only the cut `U` is reduced at this second step: `U+w` contains both endpoints
and stays at value1. There is no claim that both initial bottlenecks must be
removed simultaneously.

## 5. Exact controls and reproduction

The first member, order11, has rows

`(2032,2025,2019,2021,1006,448,640,768,576,1056,496)`.

Its two minimum shores are mask992 (`U`) and mask2016 (`U+w`). There are1215
Hamilton words,1116 black. Every black one-step neighbor is strong. Exactly
nine words decrease either initial minimum shore, and all nine are blue.
The checker explicitly enumerates all2046 nonempty proper cuts.

Its two black words are

`P=(2,1,0,5,6,7,8,9,10,4,3)`;

`Q=(9,8,7,6,5,0,1,3,4,2,10)`.

Reproduce from this directory:

```sh
python3 -B check_odd.py > odd-normal.json
python3 -B -O check_odd.py > odd-optimized.json
cmp odd-normal.json odd-optimized.json
```

The checker imports only the co-delivered stdlib verifier `check_family.py`.
It checks `a=1,...,20` (orders11 through49): actual route directedness and
blackness, the exact anti-equations of all nine decreasing words, strongness
of the actual intermediate, terminal nonstrongness, complete minimum-cut
sizes/labels from the exact sorted-score cut formula, and literal `K` words
with colours at both ends. At order11 it additionally exhausts the whole
Hamilton-word and cut universes. Normal and optimized receipts are identical.
The all-order proof is structural; these examples are finite controls.
