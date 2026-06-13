# Special Tournaments: Paley, Cyclic, and Self-Complementary

## Introduction

Among the enormous variety of possible tournaments, three families stand out as especially remarkable. Two of them — Paley tournaments and cyclic interval tournaments — are constructed using number theory. The third — self-complementary tournaments — are defined by a mirror symmetry. All three turn out to maximize or near-maximize the number of Hamiltonian paths.

---

## Paley Tournaments

### What They Are

Paley tournaments are built from **prime numbers** and a clever trick from number theory.

Take a prime number p where p = 4k + 3 for some integer k. So p could be 3, 7, 11, 19, 23, 31, ... (primes that leave remainder 3 when divided by 4).

The players in the tournament are labeled 0, 1, 2, ..., p-1 (the integers modulo p).

Now for the arrows: player A beats player B if and only if (A - B) is a **quadratic residue** mod p.

A number x is a quadratic residue mod p if there exists some y such that y² ≡ x (mod p) — in other words, x is a perfect square in the arithmetic of the clock with p hours.

For example, in mod 7 arithmetic:
- 1² = 1, 2² = 4, 3² = 2 (since 9 mod 7 = 2), 4² = 2, 5² = 4, 6² = 1
- The quadratic residues mod 7 are: {1, 2, 4}
- The non-residues are: {3, 5, 6}

So in the Paley tournament on 7 players: player A beats player B if (A-B) mod 7 is in {1, 2, 4}.

### Why They're Special

Paley tournaments have an extraordinary property called **rotational symmetry**: if you relabel all the players by shifting (player k becomes player k+1 mod p), the tournament looks exactly the same. Every player has an identical role — no one has any structural advantage over anyone else.

This perfect symmetry has a remarkable consequence: these tournaments are **self-complementary**. If you reverse all the arrows (turn every win into a loss and vice versa), you get a tournament that looks just like the original under some relabeling. The tournament contains its own mirror image.

### Their Path Counts

The Paley tournaments have some of the highest H(T) values known:

| Tournament | Players | H(T) |
|---|---|---|
| Paley T₃ | 3 | 3 |
| Paley T₇ | 7 | 189 |
| Paley T₁₁ | 11 | 95,095 |
| Paley T₁₉ | 19 | 1,172,695,746,915 |
| Paley T₂₃ | 23 | 15,760,206,976,379,349 |

For small primes (7 and 11), Paley tournaments achieve the **maximum possible** H(T) among all tournaments on that many players.

---

## Cyclic Interval Tournaments

### What They Are

Cyclic interval tournaments are also built mod a prime, but with a simpler rule.

Player A beats player B if (A - B) mod n falls in the range {1, 2, 3, ..., (n-1)/2}. In other words: player A beats the next ⌊(n-1)/2⌋ players "ahead" of them in the cyclic order.

Think of n players sitting around a clock face. Each player beats the half of the circle immediately clockwise from them.

### Why They're Special

Like Paley tournaments, these are rotationally symmetric — every player is in the same position relative to all others. But their construction is simpler and more regular.

For small primes, the cyclic tournament's path count is lower than Paley's. But at larger primes (p ≥ 13), the cyclic tournament's H(T) overtakes Paley's. The two families trade off which one wins, with a crossover point near n ≈ 13.

---

## Self-Complementary Tournaments

### What They Are

A tournament is **self-complementary (SC)** if reversing all its arrows gives a tournament that is isomorphic to the original — meaning they have the same structure up to relabeling.

Imagine a tournament where you can find a way to rename all the players such that the "new" tournament with all arrows flipped is actually identical to the "old" one.

Only tournaments on 4k or 4k+3 players (like 3, 4, 7, 8, 11, 12, ...) can be self-complementary.

### Why They're Special: The Symmetry Bonus

Self-complementary tournaments tend to have very high H(T). Here's the intuitive reason:

Reversing all arrows is an operation that maps some Hamiltonian path P to another path P'. For SC tournaments, this reversal is also a relabeling of the players. The map that simultaneously reverses arrows and relabels players creates a **pairing** among the odd cycles.

When cycles come in pairs like this, the way they interact — and the independence polynomial that computes H(T) — gets boosted in a predictable way. SC tournaments are essentially "optimized" for packing vertex-disjoint odd cycles.

### The Maximizer Conjecture

Computationally verified for all tournaments up to 8 players:

**Among all tournaments on n players, the maximum H(T) is always achieved by a self-complementary tournament.**

At n = 6: the unique maximizer is SC.
At n = 7: all global maximizers are SC.
At n = 8: all global maximizers are SC.

This has not been proved for all n yet, but the evidence is very strong.

---

## How These Families Connect

Both Paley and the cyclic interval tournament are SC (self-complementary) when the number of players is of the right form. This is no coincidence — the rotational symmetry of both families forces them into the SC world.

The three families form a nested structure:
- Cyclic interval tournaments ⊂ Rotationally symmetric tournaments ⊂ All tournaments
- Paley tournaments ⊂ Self-complementary tournaments ⊂ All tournaments
- Both Paley and cyclic are SC (when p ≡ 3 mod 4)

The fact that the highest H(T) values cluster in the SC corner of this picture is one of the deepest empirical observations of this project.

---

## Key Words

- **Paley tournament T_p**: built from quadratic residues mod a prime p ≡ 3 mod 4
- **Quadratic residue**: a number that is a perfect square in modular arithmetic
- **Rotational symmetry (circulant)**: shifting all player labels by 1 leaves the tournament unchanged
- **Self-complementary (SC)**: reversing all arrows gives an isomorphic tournament
- **H-maximizer**: a tournament achieving the highest possible path count for its size
