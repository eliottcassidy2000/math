# Cycle-preserving deletion and a second all-order obstruction proof

This is an extension of the frozen `PROOF.md` packet. It uses exactly the same
two-apex orientation-word model and the exact 32-state cycle automaton, with
no appeal to the unverified general low-cycle gluing metadata.

## Theorem: every sufficiently large low-cycle word can lose a base vertex

1. If `C(w)=1` and `|w|>=3`, some position can be deleted while preserving `C=1`.
2. If `C(w)=2` and `|w|>=4`, some position can be deleted while preserving `C=2`.

The complete list of deletion-irreducible words is

- cycle count 1: `2`, `30`;
- cycle count 2: `120`, `212`, `230`, `302`, `310`, `320`, `321`.

Thus every two-apex tournament with two Hamiltonian cycles and order at least
six has a base vertex whose deletion preserves exactly two Hamiltonian cycles.
The named roots remain fixed. This theorem is restricted to the two-apex
class; no assertion is made about arbitrary low-cycle tournaments.

## Exact finite proof

Let `delta` be the 32-state cycle transition, `s_0` its start, and `o` its output.
While reading a word `w`, retain the triple `(s,D,l)`, where

- `s=delta(s_0,w)`;
- `D` is the set of states reached by deleting exactly one letter from `w`;
- `l=min(4,|w|)`.

The initial triple is `(s_0,empty,0)`. On reading letter `a`, the exact update is

```
s' = delta(s,a)
D' = {s} union {delta(d,a): d in D}
l' = min(4,l+1).
```

The first term of `D'` deletes the new letter; the second deletes an earlier
letter. This proves the stated invariant by induction.

The complete reachable graph has **1,066 states**. For each state whose output
`c=o(s)` is 1 or 2, inspect whether any `d in D` has output `c`. The only failures
have length at most 2 for `c=1` and at most 3 for `c=2`. Since lengths at least
4 are saturated to 4, this check includes every arbitrarily long word. Direct
enumeration of the words up to length 3 then gives exactly the displayed list.

`deletion_certificate.py` constructs the complete graph and checks every state;
`deletion_certificate.json` records the state count and exceptional words. This
is a computer-assisted finite proof of the all-order assertion, not an
extrapolation from words of bounded length.

## Independent second proof that 67 needs leaving the low-cycle subclass

The following finite minima are independently evaluated by adjacency/subset-DP
on every word with one or two cycles and base size at most seven:

| Base size | Minimum H with C=1 | Minimum H with C=2 |
| --- | --- | --- |
| 1 | 3 | none |
| 2 | 5 | none |
| 3 | 9 | 13 |
| 4 | 15 | 23 |
| 5 | 25 | 41 |
| 6 | 45 | 69 |
| 7 | 75 | 115 |

There are 1,817 graph checks. The path count 67 is absent at base size at most
six. Any word of larger base size with `C=1` or `C=2` can, by repeated
cycle-preserving deletion, be reduced to base size seven. The induced-extension
lemma in `PROOF.md` then gives `H>=75` (indeed `H>=115` for `C=2`). Therefore no
word of any length has both `H=67` and `C` in `{1,2}`.

This argument is logically separate from the bounded-image saturation search:
it combines a cycle-deletion language certificate, a finite independent graph
census, and the elementary injection proving H monotonicity.

`lowcycle_minima.py` checks the table, all path counts in the finite census,
and the positive comparison witness `w=12301`, for which `H=67,C=4`.

## Constructive use and remaining question

Reversing the deletion theorem builds every two-apex low-cycle tournament from
one of nine tiny rooted word cores by base-vertex insertions that preserve its
cycle count. The automaton decides exactly which insertion is permitted, and
the 54-coordinate path transfer computes its count. This is a recursive
physical construction for the whole defined class, not a proof that every
admissible target integer has a construction in it.
