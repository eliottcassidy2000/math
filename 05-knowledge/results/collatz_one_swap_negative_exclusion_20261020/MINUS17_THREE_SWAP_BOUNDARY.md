# The -17 cycle is exactly three adjacent swaps from balance

Date: 20 October 2026. Status: **PROVED by an explicit path and an integer-flow lower bound; FINITE-EXACT independent full graph check**.

## 1. Distance definition

At shape (q,p)=(11,7), let

    c=01011011011

be the zero-intercept lower mechanical word, and let

    u=11110111000

be the chronological parity word of the actual shortcut cycle starting at -17. Define the distance between their rotation orbits to be the minimum number of cyclic adjacent swaps of unequal bits required to pass from any rotation of c to any rotation of u.

The cycle graph's swap operation commutes with rotation. Therefore one may fix c and minimize over all 11 rotations of u. This does not discard a better path starting from a different balanced rotation.

**Theorem.** This orbit distance is exactly 3, not 2.

## 2. Explicit upper bound

The following path uses three allowed swaps, with zero-based endpoint pairs shown:

    01011011011
      -- swap (0,1) --> 10011011011
      -- swap (4,5) --> 10010111011
      -- swap (3,4) --> 10001111011.

The final word is the left rotation of u by seven positions. Its carry is 12649 and its clock is 2^11-3^7=-139, so its fixed point is -91, a member of the same actual negative cycle through -17.

This proves distance at most 3. All three edges swap unequal bits; no inert move or nonadjacent interchange is counted.

## 3. Integer-flow lower bound

For any fixed source c and same-weight target v on a q-cycle, put

    b_i=sum_(j=0)^i (c_j-v_j),    0<=i<q.

Take any path of adjacent swaps. Let f_i be its signed net number of one-tokens crossing the directed edge i->i+1 modulo q. Conservation at each vertex gives

    f_i-f_(i-1)=c_i-v_i.

Hence f_i=t+b_i for one integer t, and the total number of swaps is at least

    sum_i |f_i| = sum_i |t+b_i|
                 >= min_(t in Z) sum_i |t+b_i|.        (1)

The minimum is attained by t equal to the negative of an integer median of the b_i. This is a lower bound for every swap path, regardless of how many opposite-direction moves cancel. No general attainability theorem for the flow minimum is needed here.

For the left rotations v=rot_r(u), r=0,...,10, the median formula gives respectively

    6, 5, 5, 6, 4, 6, 5, 3, 5, 4, 4.

Thus every target rotation requires at least 3 swaps. The displayed path attains 3, proving the theorem.

The exact prefix vectors, medians and lower bounds are recorded in `minus17-distance.json`. A separate exhaustive breadth-first search over all binomial(11,7)=330 binary words gives the same distance for each target rotation. Its layer sizes are

    1,8,29,60,82,90,43,14,3

at distances 0 through 8.

## 4. Correct failure boundary

The -17 cycle is a concrete counterexample to extending the nontrivial all-sign one-swap exclusion to an unrestricted three-edit family. It is not a two-edit counterexample. Hamming distance 4 only supplied a lower bound of 2 and was insufficient to identify this boundary.

Nothing here resolves the primitive two-edit negative family or the nontrivial positive two-edit family. As always, several sequential swaps may undo each other, so edit count and minimum orbit distance must not be conflated.

Reproduce:

    python3 verify_minus17_swap_distance.py > minus17-distance.replay.json
    python3 -O verify_minus17_swap_distance.py > minus17-distance.optimized.replay.json
    cmp minus17-distance.json minus17-distance.replay.json
    cmp minus17-distance.json minus17-distance.optimized.replay.json
