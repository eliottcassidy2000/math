"""Independent replay of the bounded-order pure-blue atom classification.

Reads certificate data only; imports no producer or other repository code.
Every acceptance test raises explicitly and remains active under python -O.
"""
from collections import Counter
from functools import lru_cache
import gzip
import hashlib
import itertools
import json
from pathlib import Path


PACKET = Path(__file__).resolve().parent
SOURCE = PACKET.parent / "merged_metagraph_finite_pure_atoms_20261020"


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def load_certificates():
    manifest = {}
    for line in (SOURCE / "SHA256SUMS.txt").read_text().splitlines():
        if line.strip():
            digest, name = line.split(maxsplit=1)
            manifest[name.strip()] = digest

    def read(name):
        raw = (SOURCE / name).read_bytes()
        # The imported manifest hashes Git's LF text; Windows checks out CRLF.
        payload = raw if name.endswith(".gz") else raw.replace(b"\r\n", b"\n")
        need(hashlib.sha256(payload).hexdigest() == manifest[name],
             ("certificate hash", name))
        return json.loads(gzip.decompress(payload) if name.endswith(".gz") else payload)

    data = read("cycle_certificates.json")
    banks = {}
    for row in data["certificates"]:
        banks.setdefault(row["n"], []).append(row)
    need(sorted(banks) == list(range(3, 16, 2)), "base order universe")
    for n in range(17, 32, 2):
        bank = read("circulant_n%d_certificates.json.gz" % n)
        need(bank["n"] == n, ("extended order", n))
        need(bank["presentations"] == bank["excluded"] == 2 ** ((n - 1) // 2),
             ("extended header", n))
        banks[n] = bank["rows"]
    need(sorted(banks) == list(range(3, 32, 2)), "complete order universe")
    return banks


def check_circulants():
    total = negative = positive = 0
    for n, rows in sorted(load_certificates().items()):
        m = (n - 1) // 2
        need(sorted(row["mask"] for row in rows) == list(range(2 ** m)),
             ("mask universe", n))
        for row in rows:
            steps = {d if row["mask"] & (1 << (d - 1)) else n - d
                     for d in range(1, m + 1)}

            def arc(a, b):
                return a != b and (b - a) % n in steps

            need(len(steps) == m and not steps & {n - d for d in steps},
                 ("skew step set", n, row["mask"]))
            if n < 7:
                paths = 0
                for word in itertools.permutations(range(n)):
                    if not all(arc(a, b) for a, b in zip(word, word[1:])):
                        continue
                    paths += 1
                    reverse = dict(zip(word, word[::-1]))
                    need(all(arc(a, b) == arc(reverse[b], reverse[a])
                             for a in range(n) for b in range(n)),
                         ("positive purity", n, row["mask"], word))
                need(paths == (3 if n == 3 else 15), ("positive H", n, paths))
                positive += 1
            else:
                cert = row["certificate"]
                cycle, word = cert["cycle"], cert["bad_cut_path"]
                need(sorted(cycle) == list(range(n)) == sorted(word),
                     ("vertex words", n, row["mask"]))
                need(all(arc(cycle[i], cycle[(i + 1) % n]) for i in range(n)),
                     ("directed cycle", n, row["mask"]))
                need(any(word == cycle[i:] + cycle[:i] for i in range(n)),
                     ("cycle cut", n, row["mask"]))
                reverse = dict(zip(word, word[::-1]))
                a, b = cert["reversal_bad_pair"]
                need(a in reverse and b in reverse
                     and arc(a, b) != arc(reverse[b], reverse[a]),
                     ("nonblue witness", n, row["mask"]))
                successor = dict(zip(cycle, cycle[1:] + cycle[:1]))
                a, b = cert["successor_bad_pair"]
                need(a in successor and b in successor
                     and arc(a, b) != arc(successor[a], successor[b]),
                     ("nonautomorphic successor", n, row["mask"]))
                negative += 1
            total += 1
        print("n=%d exact masks=%d" % (n, len(rows)))
    need((total, negative, positive) == (65534, 65528, 6), "certificate totals")
    # All four order-five cyclic presentations are one ordinary isomorphism class.
    five = {frozenset(s) for s in itertools.product((1, 4), (2, 3))}
    need(five == {frozenset((a * d) % 5 for d in (1, 2)) for a in range(1, 5)},
         "order-five multiplier identification")
    print("PASS: 65528 negative reversal certificates; all paths in six positive presentations.")


def check_cutoff_and_grammar():
    d, s, r = 0, 2, 1
    while r < 4 * s * (d + 1):
        r *= 2
    need((s, r, 2 * s * r) == (2, 8, 32), "zero-defect bound specialization")
    print("Imported bound at d=0: n<=32; imported odd-order reduction gives n<=31.")
    # Symbolic class words retain component order; an optional center has size1,3,5.
    @lru_cache(None)
    def half_words(n):
        if n == 0:
            return ((),)
        return tuple((a,) + tail for a in (1, 3) if a <= n
                     for tail in half_words(n - a))

    limit = 32
    counts = Counter()
    signatures = set()
    for half_size in range(limit // 2 + 1):
        for left in half_words(half_size):
            for center in (None, 1, 3, 5):
                n = 2 * half_size + (center or 0)
                if n > limit:
                    continue
                word = left + (() if center is None else (center,)) + left[::-1]
                need(word not in signatures, ("duplicate SCC grammar word", word))
                signatures.add(word)
                counts[n] += 1
    c = [0] * (limit // 2 + 1)
    c[0] = 1
    for j in range(1, len(c)):
        c[j] = c[j - 1] + (c[j - 3] if j >= 3 else 0)
    for n in range(limit + 1):
        expected = (c[n // 2] if n % 2 == 0 else
                    sum(c[(n - a) // 2] for a in (1, 3, 5) if n >= a))
        need(counts[n] == expected, ("parity count", n))
        numerator = counts[n] - (counts[n - 2] if n >= 2 else 0)
        numerator -= counts[n - 6] if n >= 6 else 0
        need(numerator == int(n in (0, 1, 3, 5)), ("rational OGF", n))
    need([counts[n] for n in range(3, 14)] == [2, 1, 3, 2, 4, 3, 6, 4, 9, 6, 13],
         "inherited finite count control")
    print("Unique palindromic SCC grammar checked through order32.")
    print("Pure-blue class counts, orders0..16:", [counts[n] for n in range(17)])
    print("PASS: B(z)=z+z^3+z^5; P(z)=(1+z+z^3+z^5)/(1-z^2-z^6).")


if __name__ == "__main__":
    check_circulants()
    check_cutoff_and_grammar()
