# Cauldron Schur Removal Game (S618)

This is a live stub for codex-2026-06-03-S618.

The user-defined cauldron game starts with natural numbers placed one at a
time into `k` cauldrons. A cauldron boils when it contains values satisfying an
additive equation, with the base rule `A+B=C`. The first reduction is classical:
if one boil ends the game, the `k=3`, two-term, distinct-value version is just
the finite Schur coloring problem, so the expected safe prefix is `S(3)=13`.

The removal version is not the same problem. Once a cauldron boils and is
removed, the player can treat that color as a sacrificed resource. This turns
the question from "does a global sum-free coloring exist on `[1,n]`?" into
"how long can active sum-free resources be staged before every cauldron has
been spent?"

Pending computation: `04-computation/cauldron_game_s618.py`.

## Assumption Challenge

Do not assume the vertices of the Tournament Analysis are the cauldrons. Other
candidate vertices are active-state signatures, residue classes, subset-sum
obligations, first-boil witnesses, removal order types, finite-sum arities, and
proof routes. The quotient should preserve the relevant survival predicate and
state honestly what witness information it forgets.

