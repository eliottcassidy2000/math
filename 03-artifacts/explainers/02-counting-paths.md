# The Core Question: Counting Paths

## What is a Hamiltonian Path?

Once you have a tournament — a set of players with arrows showing who beat whom — a natural question is:

**Can you arrange all the players in a line, left to right, so that every player beat the player to their right?**

A line like this — where you visit every player exactly once, always moving in the direction of the arrow — is called a **Hamiltonian path**. (Named after the 19th-century mathematician William Rowan Hamilton, who loved these kinds of traversal puzzles.)

In the transitive tournament (clear ranking), there's exactly one such path: the ranking order itself. Player 1 beat everyone, Player 2 beat everyone except Player 1, and so on.

In the cyclic tournament (A beats B beats C beats A), it's less obvious. Let's check:
- A → B → C? Does C beat anyone after C? No one's left, so this works! ✓
- B → C → A? Same — works! ✓
- C → A → B? Works! ✓

So the 3-player cycle has exactly **3** Hamiltonian paths.

---

## The Big Discovery: It's Always Odd

Here is a beautiful, surprising fact that is NOT obvious:

**No matter what tournament you have — no matter how many players, no matter how complex the win/loss pattern — the number of Hamiltonian paths is always ODD.**

This was proved by the Hungarian mathematician László Rédei in 1934. It's called **Rédei's Theorem**.

Think about what this means. There's no such thing as a tournament with 2, 4, 6, 8, or any even number of Hamiltonian paths. The possible counts are 1, 3, 5, 7, 9, 11, 13, ...

Why is this true? The proof uses an elegant symmetry: for any Hamiltonian path in a tournament, you can reverse it (flip all the arrows along the path) to get another path — unless the reversed path is the same as the original, which can only happen in special circumstances. This "path reversal" pairs up paths in twos, leaving exactly one unpaired. Hence the count is always 1 + 2k = odd.

---

## So What Are We Studying?

We write H(T) to mean "the number of Hamiltonian paths in tournament T."

By Rédei's theorem, H(T) is always odd. But:
- Which odd numbers actually occur?
- Can H(T) = 7? Can H(T) = 1,000,001?
- Given a tournament, how do we compute H(T) quickly?
- Which tournaments have the MOST Hamiltonian paths?

These questions are the heart of this research project.

---

## Small Examples

| Tournament size | Minimum H(T) | Maximum H(T) | All possible values |
|---|---|---|---|
| 3 players | 1 | 3 | 1, 3 |
| 4 players | 1 | 3 | 1, 3 |
| 5 players | 1 | 15 | 1, 3, 5, 7, 9, 11, 13, 15 (many) |
| 7 players | 1 | ~34 million | Exactly 77 distinct values |

Notice that at 5 players, the possible values already aren't just "all odd numbers" — some odd numbers are skipped. And at 7 players, there are gaps in which odd numbers are achievable.

---

## Why This Question Matters

Counting Hamiltonian paths is related to:

1. **Ranking problems**: How many consistent total orderings are implied by a set of pairwise preferences?
2. **Algorithm design**: Hamiltonian path problems are generally very hard to solve — tournaments are a special case where we can say more.
3. **Combinatorics**: The count H(T) encodes deep information about the structure of the tournament.

The number H(T) turns out to be a surprisingly faithful "fingerprint" of the tournament. Two very different-looking tournaments can have the same H, but two tournaments with the same H often share a lot of structure.

---

## Key Words

- **Hamiltonian path**: a path through every vertex of the tournament, always following arrows
- **H(T)**: the number of Hamiltonian paths in tournament T
- **Rédei's Theorem**: H(T) is always odd, for every tournament T
- **Transitive tournament**: the "ranking" tournament, where there's exactly one Hamiltonian path
