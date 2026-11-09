"""Re-check of the quotient-tournament parity lemma (KNOWN: THM-4922 (4), proof step 5 in
05-knowledge/results/first_creators_genealogy_20261107/audit/AUDIT.md section 4; re-derived here before that proof
was found, and used in THM-4951 (5)).

Let T be a tournament on {0, ..., 2c} with the involutive anti-automorphism iota(i) = 2c - i (u -> v iff iota(v) ->
iota(u)); iota fixes only c.  For each mirror pair p = {u, iota u} let p+ be the member with c -> p+ (exactly one
exists) and p- the other.  Define the quotient T_flat on the c pairs by
        p -> q in T_flat   iff   [p+ -> q+] + [p- -> q+] = 0 (mod 2)       (q+ beaten by both or neither member of p).
Lemma.  (a) T_flat is a tournament.  (b) For every sequence of distinct pairs (p_1, ..., p_h), the number of
mirror-free directed paths c -> v_1 -> ... -> v_h with v_j in p_j is congruent mod 2 to [p_1 -> p_2 -> ... -> p_h is a
path of T_flat].  (c) Hence the number of nonempty mirror-free paths from c is = sum over nonempty S of H(T_flat[S])
= 2^c - 1 = 1 (mod 2) (Redei).
This script checks (a), (b), (c) on random iota-tournaments, with and without a forced Hamiltonian path 0 -> 1 -> ... .
usage: python quotient_tournament_check.py
"""
import random
import itertools


def random_iota_tournament(n, rng, with_path):
    out = [0] * n
    decided = {}
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in decided:
                continue
            d = 1 if (with_path and j == i + 1) else rng.randrange(2)
            decided[(i, j)] = d
            decided[(n - 1 - j, n - 1 - i)] = d       # i -> j  iff  iota(j) -> iota(i)
    for (i, j), d in decided.items():
        if d:
            out[i] |= 1 << j
        else:
            out[j] |= 1 << i
    for i in range(n):
        for j in range(n):
            if i != j:
                assert (out[i] >> j & 1) == (out[n - 1 - j] >> (n - 1 - i) & 1)
    return out


def main():
    rng = random.Random(2026)
    for n in (3, 5, 7, 9):
        c = n // 2
        for with_path in (False, True):
            ok_a = ok_b = ok_c = 0
            trials = 150
            for _ in range(trials):
                out = random_iota_tournament(n, rng, with_path)
                beats = lambda x, y: out[x] >> y & 1
                plus = [i if beats(c, i) else n - 1 - i for i in range(c)]          # pair index i <-> {i, n-1-i}
                minus = [n - 1 - plus[i] if plus[i] == i else i for i in range(c)]
                D = [[0] * c for _ in range(c)]
                for p in range(c):
                    for q in range(c):
                        if p != q:
                            D[p][q] = 1 if (beats(plus[p], plus[q]) + beats(minus[p], plus[q])) % 2 == 0 else 0
                ok_a += all(D[p][q] + D[q][p] == 1 for p in range(c) for q in range(c) if p != q)
                good_b = True
                total = 0
                for h in range(1, c + 1):
                    for seq in itertools.permutations(range(c), h):
                        cnt = 0
                        for signs in itertools.product((0, 1), repeat=h):
                            vs = [plus[p] if s == 0 else minus[p] for p, s in zip(seq, signs)]
                            walk = [c] + vs
                            if all(beats(walk[k], walk[k + 1]) for k in range(h)):
                                cnt += 1
                        total += cnt
                        pathflat = all(D[seq[k]][seq[k + 1]] for k in range(h - 1))
                        if cnt % 2 != int(pathflat):
                            good_b = False
                ok_b += good_b
                ok_c += (total % 2 == 1)
            print(f"n = {n}, path {'forced' if with_path else 'free'}: (a) tournament {ok_a}/{trials}, "
                  f"(b) path parity {ok_b}/{trials}, (c) odd total {ok_c}/{trials}")


if __name__ == "__main__":
    main()
