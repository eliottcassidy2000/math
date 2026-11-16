# Exact 27→63 rewrites, coalescence, and the cost of a usable port

- [The tournament result](TOURNAMENT.md): no 27-path tournament
  embeds in a 63-path witness. The obstruction persists for every
  27-path witness padded by any number j of directed triangles,
  when the target is 63·3^j. A fixed small core rewrite bypasses
  this at every j, and a separate strong run from the original
  27-path seed reaches 137781 using just one old reversal and
  one new vertex at each step. The scope is explicit: that strong
  run is finite, while the padded branch is an exact infinite family.
- [Cost by source](CONTAINMENT.md): all five strong 27-path classes
  admit sharp one-reversal repairs with at most two new vertices.
  [The three-triangle source](FACTORED_EDIT.md) needs two reversals
  even with unlimited extra vertices, plus at least three vertices
  at that optimum; with no extra vertices it needs five reversals.
  The odd-cycle-pack shape changes from (1+z)^3 to 1+31z.
- [Square-sum ports, 8/9, 1472 and 1729](PORTS.md): the same
  Collatz coalescence family has infinitely many quadratic-cost
  two-edge connections and infinitely many no-overhead ones.
  The pair 1472,1729 cannot be connected by three edges with
  the same outer square, but distinct outer squares connect it
  at the smallest possible height. The shortest route between
  the original 8/9 ports emerges at height 41.
- [Golden integer carry](GOLDEN.md): consecutive Fibonacci
  parameters produce exact primitive square-sum rotations.
  One complete low-digit turn advances the high digit by one.
  The first path is the familiar 8→9 path on 1…15.

These constructions preserve different interfaces. The prescribed
counts in one infinite tower are now covered by an explicit rewrite,
but neither that proof nor the numerical analogies give a selector
for arbitrary odd targets. Chan–Pak remains open.

Reproduction: Python 3 check_moves.py, check_factored.py, check_ports.py;
C++ edit.cpp and factored_edit.cpp. The optional complete-class
replays check_27_63.py, rewires27.py and check_padded_27.py take the
order-8 CSV from the [census](../actual_hamilton_multiplicities_20261021/session24-multiplicities/).
