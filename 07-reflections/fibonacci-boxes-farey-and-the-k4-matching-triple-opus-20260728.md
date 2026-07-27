# Fibonacci boxes, reduced fractions, and the K4-matching form of Pythagoras

**opus-2026-07-28. Provenance note / wildcard synthesis, exact checks run
inline (transcribed below), continuing the size-4/V4 thread.** Owner's
prompt: relate the tree of primitive ("primal") Fibonacci-box triples to
the reduced fractions in (0,1), tournaments of size 4 and 6, and the
Fibonacci sequence.

## The verified skeleton (all exact, this session)

1. **Fibonacci window = primitive triple.** Every primitive Pythagorean
   triple is a 4-term generalized-Fibonacci window
   `e = (n-m, m, n, n+m)` (consecutive-sum law holds; coprime,
   opposite parity). The K4 on the window's 4 entries has 6 pairwise
   products (the "size 4 and 6") and exactly 3 perfect matchings, and
   **the three matchings ARE the triple** [PROVED symbolically]:

   ```text
   (03)(12):  e0*e3 = n^2 - m^2 = a,     e1*e2 = mn = b/2
   (02)(13):  e0*e2 + e1*e3 = n^2 + m^2 = c     (crossing, as a SUM)
   (01)(23):  e2*e3 - e0*e1 = n^2 + m^2 = c     (adjacent, as a DIFFERENCE)
   a^2 + b^2 - c^2 = 0  identically.
   ```

   Pythagoras is the statement that the two c-matchings agree with the
   (a, b)-matching quadratically -- a K4-matching identity. The 3
   matchings are the 3 nonzero elements of V4 (kernel of S4 acting on
   them = exactly {id + three double transpositions}, verified), so the
   triple `(a, b, c)` is literally a V4∖{0}-labeling and the previous
   note's `Aut(V4) = S3 = S4/V4` permutes the roles. The octahedral
   picture: the 6 products are the vertices of the line graph
   L(K4) = octahedron; the 3 matchings are its 3 diagonals.

2. **Reduced fractions.** `(m, n) -> m/n` puts the triple tree inside
   the reduced fractions of `(0,1)`; the ternary child moves
   `(m,n) -> (n, 2n-m), (n, 2n+m), (m, 2m+n)` generate, from `(1,2)`
   (the `(3,4,5)` root), every coprime opposite-parity fraction exactly
   once [verified: depth-6 injectivity 364/364; completeness for all
   331 fractions with denominator <= 40]. The Farey/Stern-Brocot
   mediant is Fibonacci addition of boxes; the triple tree is the
   ternary refinement of the binary fraction tree living on the
   opposite-parity stratum.

3. **The Fibonacci sequence is the golden path.** Consecutive true-
   Fibonacci windows give triples with hypotenuse `F_{2k+1}`:
   `(3,4,5), (5,12,13), (16,30,34), (39,80,89), (105,208,233), ...`
   [verified]. Note `13 = F_7` (the pair (7,13) as index/value) --
   FLAGGED ONLY: per MISTAKE-228, Fibonacci labels get no LRC role
   without a map.

## The typed connections (board, not claims)

- **Binary vs ternary = the two free factors.** The fraction tree is
  binary (Stern-Brocot; the repo's dyadic/parity towers, LEM-020),
  the triple tree ternary (the repo's C3-towers, the Keller 3-adic
  monoid, THM-2450). `PSL_2(Z) = Z_2 * Z_3`: the repo's two recurring
  tower grammars are the two free factors of the modular group, and
  the triple/fraction pair realizes both on ONE object. (Bold; flagged
  as a frame, not a theorem.)
- **Farey side.** The repo's LRC machinery already lives on reduced
  fractions (THM-2056's Kelvin-Farey certificate, the W_q Farey
  flanks). The box coordinates `(n-m, m, n, n+m)` are exactly a
  mediant-pair packet: a candidate reparametrization of Farey flank
  data with a built-in V4/S3 action -- worth checking against
  THM-2055/2056's signed-hull owner fan (does the fan's sector
  structure respect the three matchings?).
- **Resolvent continuity.** Same V4 as the resolvent-transfer note:
  S4 on 4 box entries, S4/V4 = S3 on the three matchings = on the
  triple's components. The Keller graph-quartic's resolvent cubic and
  the Pythagorean triple's three matchings are the same quotient seen
  in two problems.

## Cheapest next tests

(i) Do the three ternary tree moves act on the matchings/V4 labels in
a computable S3-pattern (i.e., is the Barning tree a V4-torsor tree)?
(ii) The Farey-flank question above (one exact session against
THM-2056's fan data). (iii) The 2-adic/3-adic free-factor frame: test
whether the repo's dyadic-tower objects and C3-tower objects ever
co-occur as the two factors of one modular-group object (the
triple/fraction pair is the model case).

No LRC/JC role is assigned to any of this (MISTAKE-228 discipline);
it is a carrier identification with three exact verified laws.
