---
id: THM-4551
title: "The boundary words of the Fibonacci-monodromy monotile are leaper flights. For the paper's boundary substitution sigma(a) = cac, sigma(c) = cacac on the theta graph (sigma = iota rho = Inn(cac) q^3), sigma^n(a) and sigma^n(c) are the rook-step cutting sequences of the straight segments joining the centres of two squares offset by the Fibonacci leapers (F_{3n-1}, F_{3n}) and (F_{3n}, F_{3n+1}); one level replaces each rook step by the knight's flight cac and the zebra's flight cacac. sigma = E phi phi phi~ E is the unique palindromic member of the seven Sturmian morphisms with matrix Q^3, its fixed point is the centred rounding word of slope 1/phi, and a centre-to-centre flight meets a lattice corner iff the leaper keeps its colour iff its two-speed lonely-runner loneliness is exactly 1/2; every boundary word is the flight of a colour-changing leaper"
status: >
  PROVED (palindrome induction on Christoffel classes; Sturmian-morphism facts CITED: Lothaire ch. 2,
  Mignosi-Seebold, Seebold's conjugate chain) + FINITE-EXACT through n = 8 (words of length 121393)
  + INDEPENDENTLY AUDITED (2026-10-06; the lane's first route to the fixed-point intercept was unjustified
  and is replaced, MISTAKE-564).
  Convention: a = horizontal rook step, c = vertical rook step; cut(p,q) = sequence of rook steps
  crossed by the segment from the centre of a square to the centre of the square offset by (p,q).
  (1) sigma^n(a) = cut(F_{3n-1}, F_{3n}), sigma^n(c) = cut(F_{3n}, F_{3n+1}) (columns of Q^{3n}).
  (2) sigma = E o phi o phi o phi~ o E = q' o q' o q (q' the mirror of q); among the 7 Sturmian morphisms
      with matrix Q^3 (Seebold's conjugate chain from q'^3 to q^3) sigma is the middle one and the only one
      with both images palindromes.
  (3) x = lim sigma^n(c) has x(k) = round((k+1)/phi) - round(k/phi): the cutting sequence of the ray of slope
      phi from a square centre (same factors as, but not equal to, the Fibonacci word); its palindromic
      prefixes are exactly the words of odd Fibonacci length, i.e. the sigma^n(a), sigma^n(c).
  (4) For coprime (p,q): cut(p,q) meets a lattice corner iff p, q both odd iff the (p,q)-leaper keeps its
      colour iff delta(p,q) = 1/2 (THM-4550 (D)); Q permutes the nonzero classes of (Z/2)^2 cyclically and
      -Q^3 = I mod 2, so every boundary word is the flight of a colour-changing leaper and never meets a
      corner; the colour-preserving Fibonacci leapers (n = 1 mod 3) are the even-length prefixes F_{3k} of x,
      where x resolves the corner tie alternately ca, ac.
  (5) Cat map: on the F_n x F_n torus board Q^n = F_{n-1} I, Q^{2n} = (-1)^n I (mod F_n) (Cassini);
      pi(F_n) = 2n (n even >= 4), 4n (n odd >= 5); on Z_8^2 orbit lengths 1,3,6,6,12,12,12,12 and
      |Fix(Q^j)| = prod gcd(d_i, 8) over the invariant factors of O_j; no ring is preserved.
  NOT CLAIMED: anything about the Euclidean fractal boundary curve beyond the paper; the statement is in
  the cycle lattice of the theta graph (a = AB, c = CB).
source: mac-mini-2026-10-06-chessboard, Fibonacci lane (owner's prompt pairing PingYou Ltd, "A chiral aperiodic polygon with Fibonacci monodromy", with the chessboard)
depends_on: []
related:
  - 05-knowledge/results/chessboard_weave_20261006.md (section 4)
  - 05-knowledge/results/collatz_lucas_monotile_discrepancy_20260930.md (the same paper, Collatz side)
  - 05-knowledge/results/collatz_circulant_20260930_circulants_lucas_cubic_monotile.md
  - 01-canon/theorems/THM-4550-chessboard-lonely-runner-dictionary-and-the-knight-weave-law.md
scripts:
  - 04-computation/experiments/chessboard_weave_20261006_fib_sturmian.py
  - 04-computation/experiments/chessboard_weave_20261006_fib_catmap.py
  - 04-computation/experiments/chessboard_weave_20261006_fib_lonely.py
---

# THM-4551 -- the monotile's boundary words are leaper flights

Statement and proof: [`chessboard_weave_20261006.md`](../../05-knowledge/results/chessboard_weave_20261006.md),
section 4 (Theorem 4.1, Corollary 4.2, Proposition 4.3) and the audit of section 10.

**Proof sketch.** `sigma = E∘phi∘phi∘phi~∘E` is a Sturmian morphism whose two images `cac`, `cacac` are
palindromes, so it commutes with reversal. For coprime `(p,q)` with `p + q` odd, `cut(p,q)` is the unique
palindrome of the Christoffel class with Parikh vector `(p,q)` (the centre-to-centre segment is symmetric about
its midpoint, which is not a lattice point). Sturmian morphisms carry Christoffel classes to Christoffel classes,
`sigma` acts on Parikh vectors by `Q^3`, and `Q^3 = I mod 2` keeps `p + q` odd; so `sigma(cut(p,q)) =
cut(Q^3(p,q))`, and induction from `cut(1,0) = a`, `cut(0,1) = c` gives the theorem. The words
`sigma^n(c)` are nested prefixes whose segments converge to the ray of slope `phi` from a square centre, which
meets no lattice point; hence the fixed point is that ray's cutting sequence, the centred rounding word.

**Reading.** One level of the supertile hierarchy replaces every rook step by a knight's or a zebra's flight;
the flights that would meet a lattice corner are exactly the colour-preserving leapers, which are exactly the
perfectly lonely two-speed riders (THM-4550 (D)), and the boundary never meets one.
