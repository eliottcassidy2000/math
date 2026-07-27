# The size-4 tournament: V4/XOR as the common carrier, and the resolvent transfer

**opus-2026-07-28. Provenance note / synthesis, not a truth source.**
Owner's reframe: "XOR and the previous concepts are all just a
tournament of size 4 essentially." Typed out, this is one object seen
five ways, and it yields one new concrete attack (the resolvent
transfer, probe in flight).

## The object

The four states of a trierment edge -- {none, ->, <-, both} -- form
the Klein four-group `V_4 = F_2^2` under XOR, coordinates
`(a, b) = ([i->j], [j->i])`. Everything this week's threads touched is
a face of it:

1. **The Gram/current split** (THM-2501/2346/2504; the RXTX analogy)
   is V_4's two coordinates: the orientation bit `a XOR b` and the
   symmetric bit. Trierments = V_4-valued edge labelings; tournaments
   = the odd coset {->, <-}; ties/missing = the even coset.
2. **THM-2504's no-go in miniature**: XOR-translations of V_4 swap
   every pair, so the four states admit no invariant tournament;
   Hamming weight orients 5 of 6 state-pairs and forces exactly one
   tie -- the orientation pair {->, <-} itself. The affine seven-cube
   no-go is this at F_2^7.
3. **Aut(V_4) = GL_2(F_2) = S_3** acting on the three nonzero states:
   the S_3 of the D5 dictionary's geometric side. The collision
   anatomy "1 fixed + 2 swapped" is an involution of this S_3 acting
   on {->, <-, tie}. The char-3 trierment algebra (previous
   reflection) is the action of this S_3 on V_4 \ {0}.
4. **The Smith no-go / rotation repair** (mac-mini S143/S144): char 2
   degenerates the split because V_4 is 2-torsion; the lawful forcing
   moved to the odd primes 7, 13 -- away from the carrier's own
   characteristic.
5. **G1's forbidden degree is |V_4| = 4**, and quartic Galois theory
   factors through `S_4 / V_4 = S_3` -- the RESOLVENT map -- with
   `disc(quartic) = disc(resolvent cubic)` identically. The V_4
   kernel is exactly where a degree-4 wild map would have to hide its
   monodromy, and the quotient lands in the SOLVED grade-3 world.

## The new attack: the resolvent transfer (probe in flight)

If a degree-4 point-cap Keller map existed (G1), its graph quartic
`N(t; P)` over the target would have a resolvent cubic `R(w; P)` with
the SAME discriminant. The solved grade-3 anatomy (THM-2473/2546:
depressed escape eliminant, disc = -4 (square)^2 L, Jelonek lead,
cuspidal law) therefore potentially constrains grade 4 THROUGH the
V_4 quotient: if the graded Keller system forces R into the grade-3
normal form (or forces any identity violating a proved grade-3 law),
boxes of G1 die. Probe: construct N and R symbolically on THM-2465's
smallest live box, impose the graded Keller ideal, and test for
forced structure (depressed-resolvent law; lead(R) vs lead(N) = ell;
the disc identity chain). Either outcome is progress: a forced
identity is a new necessary condition on G1; no constraint kills the
idea cleanly.

## Losses / honesty

- "Tournament of size 4" is here a carrier statement (V_4 with its
  orientation no-go), not a literal 4-vertex tournament theorem; the
  four 4-vertex tournament isomorphism classes are not used.
- The resolvent transfer needs the Keller system to actually reach R;
  the probe decides. Nothing above excludes G1 or changes any ledger.
