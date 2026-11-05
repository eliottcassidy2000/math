# A primitive bounded-span interface with a sparse, exactly decodable image

Session 197. This packet gives a concrete branching construction, a backward reduction for its full interface state, and an all-order obstruction to scalar completeness. In particular, even primitive observable continuation content is not enough for arbitrary-target reachability.

## 1. Three actual cells and a physical seed

Use the distinct-terminal one-seam construction of `../ribbon-terminals/PROOF.md`. The source root is its last vertex, each cell exit is its last vertex, and all cross edges point forward except the cell entry pointing to the old root.

Adjacency rows encode outneighbors by zero-based bit positions. The three cells are:

- A: rows `[6,4,8,3]`, entry 2, exit 3. A Hamiltonian cycle is `(0,1,2,3)`.
- B: rows `[54,60,56,1,8,24]`, entry 2, exit 5. A Hamiltonian cycle is `(0,1,2,5,4,3)`.
- C: rows `[56,33,3,22,6,28]`, entry 0, exit 5. A Hamiltonian cycle is `(0,3,4,1,5,2)`.

Each listed order has backward span at most three, and each entry is in the first three positions. Consequently arbitrary sequences of these cells stay at span three by literal concatenation. In the oldest-first width-three mask convention, their append words are respectively `0013`, `001146`, and `446637`.

The seed is the directed triangle with rows `[4,1,2]`, rooted at vertex 2; its mask word is `044`. Its returning profile is `(h,e,R)=(3,1,3)`. Thus the entire construction is the explicit word language

`044 (0013 | 001146 | 446637)*`.

Every output has a Hamiltonian cycle, hence is strong. Indeed cut an existing cycle at the old root to obtain an old path starting there, and cut the cell's displayed cycle immediately after its entry to obtain a cell path ending there. Join the two paths by a forward cross edge and close by the sole reversed seam. Neither endpoint of the forward join is the exceptional old-root/cell-entry pair. This inductively constructs a cycle at every step, and allows the next source path to start at any chosen exit.

## 2. Exact matrices and preserved primitive charge

The physical transfer theorem gives the following column-state matrices on `(h,e,R)`:

```
A = [[5,-1,7],     B = [[17,-4,28],    C = [[41,-7,49],
     [1, 0,1],         [4,-1,6],          [5,-1,6],
     [7,-1,14]],       [28,-6,51]],       [75,-12,93]].
```

Their gadget path constants are recomputed from the literal graphs in the verifier. Each matrix preserves the row `ell=(-2,1,1)`. Since the seed has charge `ell*(h,e,R)=-2`, every reachable state satisfies

`R=2h-e-2`.

For any rooted tournament, `e<=h`, and `R>=h`: split each Hamiltonian path immediately before the root to inject its paths into the covers counted by `R`. In this family the charge therefore strengthens the physical bounds to

`h>=3`, `1<=e<=h-2`, `R>=h`.

The inequality `e>=1` follows from the constructed Hamiltonian cycle. The value `h` stays odd directly from the matrices: `e` and `R` have the same parity, each first row has odd first coefficient, and its other two coefficients sum to an even integer.

Any common divisor of `h,e,R` divides `-2h+e+R=-2`; because `h` is odd, it is one. Thus every actual returning state is primitive.

### Even the observable continuation content is one

Subscripts denote appended cells in chronological order, so `h_BA` means B then A. Direct matrix multiplication gives

```
h_A  = 5h-e+7R,
h_AA = 73h-12e+132R,
h_BA = 277h-61e+491R.
```

It follows that

`h_AA-12h_A-13h=48R`,

`h_BA-61h_A+28h=64R`.

A common divisor of `h,h_A,h_AA,h_BA` is odd and divides both `48R` and `64R`, hence divides `R`. The first equation for `h_A` then makes it divide `e`, so the charge forces that divisor to be one. Therefore

`gcd(h,h_A,h_AA,h_BA)=1`

at every reachable state. In particular the gcd of all future output counts, including the current output, is one. This is stronger than merely having a primitive state vector.

## 3. Exact backward reduction for the pair state

Eliminate `R` with the charge. The maps on `(h,e)` are

```
A: (19h-8e-14,    3h-e-2),
B: (73h-32e-56,  16h-7e-12),
C: (139h-56e-98, 17h-7e-12).
```

Their two-dimensional determinants are respectively `5,1,-21`; in particular B is a unimodular affine map. Nevertheless their positive images occupy disjoint narrow regions.

Write an output pair as `(u,v)`. For any input satisfying the displayed physical bounds, direct subtraction gives:

- after A, `1/7 < v/u < 1/5`;
- after B, `1/5 < v/u < 1/4`;
- after C, `1/9 < v/u < 1/7`.

For example after A, `7v-u=2h+e>0` and `u-5v=4h-3e-4>=h+2>0`. For B the two corresponding positive differences are `5v-u=7h-3e-4>=4h+2` and `u-4v=9h-4e-8>=5h`. For C they are `9v-u=14h-7e-10>=7h+4` and `u-7v=20h-7e-14>=13h`.

Thus the last cell is uniquely determined by the pair. Its inverse is:

```
A: h=(-u+8v+2)/5,       e=(-3u+19v-4)/5;
B: h=-7u+32v-8,         e=-16u+73v-20;
C: h=(u-8v+2)/3,        e=(17u-139v-2)/21.
```

A complete pair-membership test is therefore: accept the seed `(3,1)`; otherwise determine the unique interval, apply its inverse, reject nonintegral values or failure of `h>=3`, `h` odd, `1<=e<=h-2`, and continue. If no interval applies, reject. Each accepted inverse step strictly decreases `h`, as shown below, so this terminates. Reversing its successful steps constructs exactly the original cell word. This proves both soundness and completeness for the pair-state language, and uniqueness of its word.

This is a backward reduction on a two-coordinate interface. A scalar target `h` does not specify `e`; pair recognition is not a scalar all-target compiler.

## 4. Strong sparsity despite primitive content

Using `R>=h` and `e<=h` in the physical first rows gives

`h_A>=11h`, `h_B>=41h`, `h_C>=83h`.

Consequently a word of `d` cells has count at least `3*11^d`. Among counts at most `X>=3`, its word length is at most `D=floor(log_11(X/3))`, and the number of distinct scalar outputs is at most

`1+3+...+3^D <= (3/2)*(X/3)^(log_11 3)`.

Since `log_11 3<1`, this class has zero asymptotic density, even among odd integers. The same argument applies after adding any fixed finite set of physical seeds with the same charge and odd count. Collisions in scalar outputs only make this upper bound smaller.

There is also a simple residue obstruction for the stated seed: each pair map preserves `h mod 4` when `h` is odd, so every output here is `3 mod 4`. The first one-cell counts are 35, 131, and 263.

Thus neither a conserved primitive interface, primitive observable content, a branching physical grammar, nor a unimodular inverse branch establishes scalar completeness. These cells give all four features while their image remains sparse. No obstruction is claimed for all bounded-span tournaments, other cell sets, infinitely many seeds, or a different state-dependent construction.

## 5. Verification and dependency

Run `python session197/ribbon-primitive/verify.py`. The checker imports only the graph-level helpers in `../ribbon-terminals/verify.py` (SHA-256 `9ada1cc1cfaf33030a89122dcf4ed65103aa23f5a2ad01e5dd5c288cde198a3e`). It recomputes all three gadget matrices, checks each displayed Hamiltonian cycle and append word, recounts every literal zero-, one-, and two-cell graph, checks the decoder on every cell word through length eight, and compares it against all admissible pair states with `h<=501`.

The proof of the all-order distinct-terminal transfer is in the separately frozen `../ribbon-terminals/PROOF.md` (SHA-256 `45383a2eeb947587a6ce35c2367053c81c2db1813da990cc19aa4a4a1e93f91e`). The present scalar sparsity and decoder conclusions concern only the three explicit cells above.
