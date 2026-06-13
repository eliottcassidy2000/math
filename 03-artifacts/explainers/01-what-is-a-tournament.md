# What is a Tournament?

## The Idea

Imagine a round-robin sports tournament. Every team plays every other team exactly once. No ties allowed — one team always wins. After all the games, every pair of teams has a definite result: A beat B, or B beat A.

That's a **tournament** in the mathematical sense — just the collection of "who beat whom" results for every pair of players.

We don't care about scores or margins. We only care about the arrow: who beat whom.

---

## Drawing It

Mathematicians draw tournaments as **directed graphs** — a set of dots (called vertices or nodes) connected by arrows.

- Each dot = a player
- Each arrow = a match result (the arrow points FROM the winner TO the loser)

For 3 players (A, B, C), there's only one possible structure (up to relabeling) where everyone beats someone:

```
A → B → C → A
```

That's a cycle. A beats B, B beats C, C beats A. No one is the clear champion — it's circular.

But there's also a version where one player dominates completely:

```
A → B
A → C
B → C
```

A beats everyone, B beats C. This is called a **transitive tournament** — there's a clear ranking.

---

## How Many Are There?

With 3 players: 2 different tournaments (the cycle and the ranking).
With 4 players: 4 different tournaments.
With 5 players: 12 different tournaments.
With 7 players: 456 different tournaments.

The number explodes rapidly. With 10 players there are already over 9 million structurally distinct tournaments. With 20 players, more than there are atoms in the observable universe.

This explosion is one reason the subject is rich: there's an enormous variety of possible structures.

---

## What Makes Two Tournaments "The Same"?

Two tournaments are considered the same (isomorphic) if one can be obtained from the other just by relabeling the players. The teams might be called different names, but the pattern of wins and losses is identical.

For example: "Alice beats Bob beats Carol beats Alice" is the same tournament as "X beats Y beats Z beats X" — just with different labels.

When we count "how many tournaments are there," we count the distinct patterns, ignoring the labels.

---

## Why Tournaments?

Tournaments turn up naturally anywhere there's a pairwise comparison with no ties:

- Round-robin sports competitions
- Election outcomes (candidate A is preferred over B by more than half the voters)
- Dominance hierarchies in animal groups
- Ranking systems where every pair is compared
- Any situation where you have n items and want to rank them by pairwise votes

The math of tournaments also appears unexpectedly in coding theory, statistical physics, and topology — connections we'll explore in later write-ups.

---

## Key Words

- **Tournament**: a complete directed graph — every pair of vertices has exactly one arrow between them
- **Vertex**: a player or team (a dot in the diagram)
- **Arc**: a match result — an arrow from winner to loser
- **Isomorphism**: two tournaments are "the same" if one is just a relabeled copy of the other
- **Transitive tournament**: there's a complete ranking — one player beats everyone, another beats everyone except the first, and so on all the way down
