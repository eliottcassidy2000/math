---
source: codex-2026-06-03-S599
status: exploratory extension of HYP-2142; bounded-window certificate extraction
tags: [LRC, Cprime, two-block, determinant, Helly, CRT, endpoint-owners, n14]
---

# Two-block determinant Helly witnesses

The repo already had the right algebraic object before this session:

```text
det[[u_a,r_a],[u_b,r_b]] = u_a r_b - u_b r_a
                         = w n u_a u_b length(C_i).
```

S595 made that the large-owner two-block.  S581 made it a bounded CRT automaton
in `w`.  S598 put it after aggregate cap filters.  The missing intermediate
question was: if the bounded automaton is empty, how small is the contradiction?

S599 turns that into a Helly audit.  Each component of `G(S')` contributes a
bounded language in

```text
1 <= w <= floor((n-1) max(S') / n).
```

The cover attempt is the intersection of those languages.  Instead of
materializing the full CRT modulus, the script extracts a minimal empty
subfamily up to size four.

## What popped up

The determinant/two-block idea is now connected to four repo layers:

- S574/S581: endpoint-owner congruence windows become component languages in `w`.
- S595/HYP-2142: those languages are the owner/slack 2x2 determinant rows.
- S596/HYP-2143: Collatz has the same bounded two-block Diophantine shape.
- S598/HYP-2140: cap filters remove aggregate overloads before the determinant layer.

There is also a useful contrast with HYP-2042.  Raw LRC clearance does not seem
to have a fixed small Helly number.  The determinant residual might.  That is
not a contradiction because the vertices have changed: not runners, but
post-owner proof obligations in a dominance-bounded one-multiplier language.

## Computed signal

In the `BC_only` regime, after Bprime and Lemma C:

```text
n=6:  singleton 286, pair 12, live 0
n=8:  singleton 196, pair 2,  live 0
n=10: singleton 89,  pair 0,  live 0
n=12: singleton 49,  pair 0,  live 0
n=14: singleton 15,  pair 0,  live 0
```

So the sampled n=14 determinant residual did not need a global CRT modulus at
all: every bounded automaton case that reached S599 collapsed to a singleton
component wall.

The first pair witness is small but conceptually important:

```text
n=6 row=(1,12,13,14,15)
C0 allowed={1}
C1 allowed={2,4}
```

After the full S598 + owner stack, the first pair witness is all-large-owner:

```text
n=6 row=(2,6,8,10,13)
C0 allowed={1}
C1 allowed={3,7,10}
```

That second example is the cleaner future lemma target: two large-owner
determinant rows force disjoint bounded multiplier windows.

## Rebase integration

Origin added monad-compute S598 during push: a widened exhaustive Cprime box plus
n=9 check, 6.24M configs and 0 tight.  Good news, and a different kind of good
news.  Monad makes the finite empirical wall taller; S599 tries to name the
small determinant certificate that would explain why the wall is there.

## New proof shape

The proof program should try this order:

1. Singleton determinant wall.
2. Pair determinant Helly wall.
3. Triple or bounded-size determinant Helly wall.
4. Factored CRT/ZDD only for rows surviving the small witnesses.

The taste of the thing is good: S599 replaces "the automaton is empty" with
"which two-block rows make it empty?"  That is much closer to a proof one could
write by hand.

## Tournament Analysis

Vertices are certificate sizes, not runners or arc centers:

```text
singleton_empty, pair_empty, triple_empty, quad_empty,
high_order_empty, bounded_live, preempted_gate
```

Observable:

```text
(certificate_rank, sampled_route_count, name)
```

Switch/gauge: smaller empty determinant subfamilies beat larger residual
burdens.  The fingerprint is transitive, with no directed 3-cycles, singleton
SCCs, and one Hamiltonian path.

## Assumption challenge

I considered runners, gaps, cap centers, components, endpoint owners, owner
pairs, residue classes, prime powers, and proof obligations as vertices.  The
chosen quotient keeps exactly the predicate needed for Cprime after dominance:
whether a bounded `w` can cover all components.  It discards the full CRT lcm
and phase order beyond that bounded window.

Challenged assumption: the next proof does not have to begin with the global
CRT automaton.  The small empty determinant subfamily may be the actual
human-scale certificate.
