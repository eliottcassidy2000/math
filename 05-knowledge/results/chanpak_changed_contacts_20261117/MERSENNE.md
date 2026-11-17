# The sevenfold Mersenne ray

**Exactly:** `7(2^k−1)` is attainable by a tournament for every
integer `k ≥ 3`. For `k=1,2` the counts are the known holes 7
and 21. This is an exact infinite family, not an assertion of
infinitely many prime quotients or of sevenfold Fermat coverage.

Use the left scaffold with outgoing bit masks

    (4, 5, 72, 19, 39, 79, 155, 63).

Its terminal port is vertex **5**. Replace vertex **0** by a
nonempty transitive block of size m, keeping uniform contacts
to the other seven vertices. This is the only growing module.
Write a for the full path count and b for the weighted end-at-5
read: choose a path ending at 5 and an optional complementary
path on all omitted vertices. These are ordered covers, with
the marked first path nonempty.

Vertex 0 beats only vertex 2 in the scaffold. In a full
Hamiltonian path every clone run either uses vertex 2 immediately
after it or ends the full path, so there are at most two.
In a weighted end-at-5 cover every clone run either uses vertex
2 immediately after it or ends the optional second path, again
at most two. The first path ends at old vertex 5 and cannot end
in the clone module. Ordering the visits fixes the run slots.
Every one-slot allocation has one increasing ordering; every
two-slot allocation has `2^m−2` choices. Thus a and b are affine
in `2^m` at all heights. The literal m=1 and m=2 checks are

    m = 1: (a,b) = (129,190)
    m = 2: (a,b) = (257,382),

giving exactly `a=64·2^m+1` and `b=96·2^m−2`.
A terminal frame for port 5 exists at m=1 and survives by
replacing its occurrence of vertex 0 with the consecutive
transitive block.

Now use an actual right transitive tournament TT3, with its
source named u. Orient all arcs right-to-left except for the
single seam `5→u`. In the terminal-profile convention the
right reader is its converse at terminal u, with tuple
`(d,e,f)=(1,4,0)`. A direct proof prevents a source/sink
convention error:

- A Hamiltonian path not using the seam must traverse the
  whole right TT3, in its natural order, then the entire left.
  The right path ends at its sink, which beats every left
  vertex, so there are a choices.
- If the path uses `5→u`, any subset of the other two
  right vertices can be in an initial right run (four choices).
  Then read the first left run ending at 5, the seam, the
  remaining right vertices starting with u, and the optional
  second left run. These left choices are exactly b.
  Right-to-left links are valid because the second left run
  avoids the already used terminal port 5.

So the assembled count is

    a+4b = 448·2^m−7 = 7(2^(m+6)−1),  m≥1.

These witnesses are strong: concatenate a terminal-left frame
with the full right frame starting at u, and close it by the
arc from the right sink to the initial left vertex.

For k=3,4,5,6 the respective counts 49,105,217,441 are handled
by literal 49 and 217 witnesses and ordinal products 3·35 and
49·3². All are checked in the standalone replay. The strong
tail construction covers every remaining height.

Run `python replay_mersenne.py` in this directory, or pass its path.
The independent literal subset path count checks the two
profile base cases, two out-of-base controls, the assembled
seam, spanning cycles via reachability, and the four small
counts. No archive census is necessary for this positive
assertion.
