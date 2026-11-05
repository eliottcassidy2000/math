# A sign-changing primitive interface already at span two

Session 197. This addendum gives a physical change of charge, an exact backward decoder, and a sparse image. It strengthens the fixed-charge example without asserting target completeness.

## Two literal cells

Use the reviewed distinct-terminal transfer on `(h,e,R)`. Start from the triangle `[4,1,2]` rooted at its last vertex, with profile `(3,1,3)` and charge `q=-2h+e+R=-2`.

The two cells, in their displayed orders, are:

- A: adjacency rows `[8,5,1,6]`, entry 0, exit 3. This is the earlier A cell reordered as `(2,0,1,3)`.
- D: adjacency rows `[2,0,3]`, entry 1, exit 2. This gadget is transitive.

Both have backward span at most two, their entries lie in their first two positions, and their exits are last. Every cell attaches by the single reversed edge from its entry to the old last root. In the oldest-first width-two word convention the complete language is

`022 (2213 | 013)*`.

All actual outputs are strong. The seed has a Hamiltonian cycle. The A cell has the Hamiltonian path `(3,1,2,0)` ending at its entry; D has `(2,0,1)`. Cut the old Hamiltonian cycle at the source root, join its path to the chosen cell path using a forward cross edge, and close by the reversed seam. This constructs a Hamiltonian cycle after either cell even though D itself is not strong.

Their exact matrices are

```
A = [[5,-1,7],     D = [[1, 0,4],
     [1, 0,1],         [0, 0,1],
     [7,-1,14]],       [4,-1,6]].
```

With `ell=(-2,1,1)`, one has `ell A=ell` and `ell D=-ell`. Thus `q` always lies in `{-2,+2}`. Moreover

```
D^2 = [[17,-4,28],
       [4,-1,6],
       [28,-6,51]],
```

exactly the earlier B matrix. Its unimodular return is two physical sign changes, not evidence that its scalar image fills an interval.

## Primitive state and primitive observable content

Every rooted physical state has `e<=h` and `R>=h`. Here strongness also gives `e>=1`. The first rows preserve oddness because `e` and `R` have the same parity, the first coefficient is odd, and the remaining coefficient sum is even.

Since `-2h+e+R=q=+/-2`, the state gcd is one. More strongly, use the already literal continuations A, AA, and DDA. Because D squared equals B, their count rows are

`h_A=5h-e+7R`,

`h_AA=73h-12e+132R`,

`h_DDA=277h-61e+491R`.

Therefore

`h_AA-12h_A-13h=48R`,

`h_DDA-61h_A+28h=64R`.

A common divisor of these four counts is odd, divides R, then e, and finally q. It must be one. Thus primitive observable continuation content persists throughout this span-two language.

## Backward reduction with the charge recorded

Eliminating `R=2h-e+q` gives

```
A: (h,e,q) -> (19h-8e+7q, 3h-e+q, q),
D: (h,e,q) -> (9h-4e+4q, 2h-e+q, -q).
```

The physical input conditions are `h>=3` odd, `q=+/-2`, `1<=e<=h`, and `e<=h+q` (the latter is `R>=h`). For an output `(u,v,q_out)`:

- A gives `1/7 < v/u < 1/5`;
- D gives `1/5 <= v/u < 1/3`.

For A, `7v-u=2h+e>0`, while `u-5v=4h-3e+2q>0` under either charge and the displayed input bounds. For D, `5v-u=h-e+q=R-h>=0`, and `u-3v=3h-e+q>0`. The equality endpoint `v/u=1/5` belongs only to D; it occurs at the seed's first D step.

The unique inverse is consequently

```
A: h=(-u+8v-q_out)/5,
   e=(-3u+19v+2q_out)/5,
   q=q_out;

D: h=u-4v,
   e=2u-9v-q_out,
   q=-q_out.
```

The complete state-language decoder accepts `(3,1,-2)` and otherwise identifies the interval, applies the inverse, rejects nonintegrality or failure of the physical input bounds, and repeats. Its successful predecessor has smaller h because the forward maps obey `h_A>=11h` and `h_D>=5h`. Reversing a successful reduction constructs the graph word, so the decoder is both necessary and sufficient and the state representation is unique.

This is a decoder for the interface state, not for a scalar target without its other coordinates.

## Reachability limitation

All outputs still satisfy `h=3 mod 4`: D changes h by `4R`; A changes it modulo four by `h -> 3h+2`, which fixes each odd residue. The first D outputs are 15 and 131, with profiles `(15,3,29)` and `(131,29,231)` and charges +2 and -2.

A length-d cell word has `h>=3*5^d`. Hence the number of scalar outputs at most X is at most

`1+2+...+2^floor(log_5(X/3)) <= 2*(X/3)^(log_5 2)`.

It has zero density among odd targets. Adding a fixed finite collection of compatible odd physical seeds does not change that density conclusion. The argument concerns these two cells, not all span-two tournaments or all sign-changing interfaces. In particular it does not establish a new global forbidden count.

## Checks and dependency

Run `python session197/ribbon-sign-flip/verify.py`. It recomputes the matrices from both actual gadgets, compares every graph word through three cells against an independent subset-DP profile, checks all interface words through length ten against the exact decoder, and tests the growth, parity, charge, and observable-gcd identities.

The pinned graph helper is `../ribbon-terminals/verify.py`, SHA-256 `9ada1cc1cfaf33030a89122dcf4ed65103aa23f5a2ad01e5dd5c288cde198a3e`. The reviewed physical matrix theorem is `../ribbon-terminals/PROOF.md`, SHA-256 `45383a2eeb947587a6ce35c2367053c81c2db1813da990cc19aa4a4a1e93f91e`.
