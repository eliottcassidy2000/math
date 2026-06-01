# Three parities and the doubled prime as the bridge (S521)

*claudebox-2026-06-01-S521. Pushing the cascade-to-OCF reduction and tracing one
parity principle through three settings — cycle length (OCF), the number n
(first-even bridge), and the additive prime decomposition (Goldbach/Lemoine) — to
understand the role of DOUBLED PRIMES `2q`. Honest: the cascade/OCF reduction and
the decompositions are rigorous; the unifying "parity bridge" is a structural
analogy grounded where noted.*

## Cascade -> OCF: the deadlock is an ODD cycle

The clearance cascade (clear runners into the band sequentially) deadlocks only at a
**frustrated directed cycle** — a cycle that cannot be linearized. In a tournament a
cycle exists iff a 3-cycle exists (Moon), and the non-linearizable obstruction is
the **odd-cycle structure** `Omega(T)`, the project's OCF (`H = I(Omega, 2)`;
Redei: `H` is always ODD). So:
> **LRC <=> the observer lies on no frustrating ODD cycle (it can be made a source).**
Even cycles are bipartite (decompose into two consistent paths, "linearizable in
pairs"); ODD cycles are the genuine frustration. The cascade's only obstruction is
odd cyclicity at the seam.

## Parity I — cycle length (OCF)

- **Odd cycle** = frustration / non-bipartite / the OCF obstruction; it is
  irreducible (cannot be split into two consistent halves).
- **Even cycle** = bipartite / decomposable into two odd paths.
- Redei's `H(T)` is odd because Hamiltonian paths pair under reversal except
  self-reverse ones — a parity (involution) count.

## Parity II — the number n (first-even bridge)

- **n odd** -> the clean cyclotomic case: `Q(zeta_n)`, the regular polygon, the
  Galois orbit; single-prime descent often suffices.
- **n even** -> the first-even bridge `n = 2 * odd`: the 2-adic tree, the
  gate/double-gate dyadic ladder, CRT across the factor 2 and the odd part.
- So even n is "an odd n with one doubling applied" — `n = 2 * (n/2)`, one level
  down the 2-adic tree.

## Parity III — the additive prime decomposition (Goldbach / Lemoine)

- **Even n = p + q** (Goldbach): two prime atoms (same parity, both odd primes).
- **Odd n = p + 2q** (Lemoine/Levy): one prime + one **DOUBLED prime** `2q`.
- The doubled prime `2q` is the **parity carrier**: prime atoms are (almost all)
  odd, so to reach an ODD total you need an EVEN adjustment — exactly one doubled
  prime supplies it. Even totals need none (two odd primes already sum even).

## The doubled prime `2q` is the bridge between odd and even

Across all three settings the factor 2 / the doubling is the odd<->even bridge:

| setting | odd object | even object | bridge (doubling) |
|---|---|---|---|
| cycles (OCF) | odd cycle (frustration) | even cycle (bipartite) | the 2-cover / voltage lift (THM-378) |
| number n | n odd (cyclotomic) | n even (first-even bridge) | n -> 2n, one 2-adic level |
| additive | prime p | doubled prime 2q | the `x2` in Lemoine `p + 2q` |
| LRC speeds | odd speed | even speed `2v` | the doubling resonance (a 2-adic tree edge) |

**The doubled prime is precisely the minimal even atom that lets you cross from the
odd/prime/cyclotomic world to the even/2-adic world.** In LRC this is the first-even
bridge (`n = 2 * odd` = one doubling of an odd `n`), and the binding pair summing to
`n` (Thm B) reflects the decomposition: for odd `n` the pair is (odd, even) — the
even member is the doubled atom; for even `n` the pair is same-parity (a
Goldbach-type pair).

## Why doubled primes matter: they are BENIGN — the difficulty is the prime PAIR

The crucial finding (Galois-Weil, S521): **the doubling resonance `v_j = 2 v_i` is
BENIGN** (a 2-adic tree edge, positive/zero equidistribution bias), while the
**additive sum-relation `v_i + v_j = v_k` is the OBSTRUCTION** (negative bias). So:

> the doubled prime / the factor 2 carries PARITY but not DIFFICULTY. The 2-adic
> doublings are navigable tree edges; the genuine LRC obstruction is the ADDITIVE
> prime-pair structure (the Goldbach/Lemoine `p + q` / `p + 2q` resonances = the
> sum-relations) — the same additive triples that defeat equidistribution.

This sharpens the role of doubled primes precisely: in the parity bridge they are
the FREE direction (the 2-adic edge), and the locked direction is the additive
pairing of primes. The first-even bridge `n = 2 * odd` is hard not because of the
doubling itself but because the even `n` opens BOTH the 2-adic and the odd-prime
factor (CRT), and the additive resonances live in the odd part.

## The unified picture

- The cascade deadlock = an odd cycle (OCF) = frustration; LRC = no frustrating odd
  cycle through the observer.
- Odd/even cycles, odd/even `n`, and primes/doubled-primes are three faces of one
  parity principle; the doubled prime `2q` (and the factor 2 generally) is the
  odd<->even bridge — the 2-cover on cycles, the 2-adic level on numbers, the `x2`
  in Lemoine, the benign doubling edge on speeds.
- The doubled prime is the PARITY carrier and the BENIGN direction; the difficulty
  is the additive prime PAIR (Goldbach for even, the `p` in Lemoine for odd), i.e.
  the sum-relations / additive triples — the locked direction.

## Lead

Formalize the cascade -> OCF reduction: clear the backbone (Facts 1&2), reduce LRC
to the observer lying on no frustrating odd cycle, and analyze that via the OCF
`Omega` restricted to odd cycles through the observer. Separately, treat the
first-even bridge `n = 2 * odd` by the 2-adic/odd CRT split: the 2-adic factor is
benign (doublings), so push the difficulty into the odd part, where the additive
prime-pair (Goldbach/Lemoine) resonances are the explicit obstruction — the
arithmetic core the whole arc keeps reaching.
