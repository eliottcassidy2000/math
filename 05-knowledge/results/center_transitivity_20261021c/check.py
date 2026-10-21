"""Exact, search-free consumer for the complete Cayley-25 exclusion.

Run python -B check.py to save receipt.json; python -B -O check.py --verify
compares the recomputed deterministic receipt. Discovery is in produce.cpp.
No producer search counters are used to establish a mathematical claim.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
COUNT = 0

def require(value, message):
    global COUNT
    COUNT += 1
    if not value:
        raise ValueError(message)

def literal_int(x, lo, hi, name):
    if type(x) is not int or not lo <= x <= hi:
        raise ValueError(name)

def difference(x, y, kind):
    if kind == 0:
        return (y - x) % 25
    return (y % 5 - x % 5) % 5 + 5 * ((y // 5 - x // 5) % 5)

def source(kind, mask):
    literal_int(kind, 0, 1, "group kind")
    literal_int(mask, 0, 4095, "connection mask")
    reps = tuple(d for d in range(1, 25) if d < difference(d, 0, kind))
    connection = {difference(d, 0, kind) if mask >> i & 1 else d
                  for i, d in enumerate(reps)}
    return tuple(sum(1 << y for y in range(25)
                     if difference(x, y, kind) in connection) for x in range(25))

def validate_source(rows):
    if type(rows) is not tuple or not rows:
        raise ValueError("literal source tuple")
    n = len(rows)
    for i, row in enumerate(rows):
        if type(row) is not int or not 0 <= row < 1 << n or row >> i & 1:
            raise ValueError("source row")
        for j in range(i):
            if ((row >> j) ^ (rows[j] >> i)) & 1 != 1:
                raise ValueError("source tournament")

def validate_word(rows, word):
    n = len(rows)
    if type(word) is not tuple or len(word) != n or any(type(v) is not int for v in word):
        raise ValueError("literal word")
    if sorted(word) != list(range(n)):
        raise ValueError("word permutation")
    if not all(rows[a] >> b & 1 for a, b in zip(word, word[1:])):
        raise ValueError("word is not Hamilton")

def ordered_key(rows, word):
    """Full ordered adjacency, not a score word or a partial invariant."""
    validate_source(rows)
    validate_word(rows, word)
    key = 0
    for i, a in enumerate(word):
        for b in word[i + 1:]:
            key = 2 * key + ((rows[a] >> b) & 1)
    return key

def blue(rows, word):
    validate_source(rows)
    validate_word(rows, word)
    n = len(rows)
    return all(((rows[word[i]] >> word[j]) & 1) !=
               ((rows[word[n-1-i]] >> word[n-1-j]) & 1)
               for i in range(n) for j in range(i + 1, n))

def decode_word(text):
    if type(text) is not str or len(text) != 25 or set(text) != set("ABCDEFGHIJKLMNOPQRSTUVWXY"):
        raise ValueError("encoded source word")
    return tuple(ord(c) - ord('A') for c in text)

def audit_frame(kind, mask, encoded):
    rows = source(kind, mask)
    if type(encoded) is not tuple or len(encoded) != 7:
        raise ValueError("seven literal words")
    keys = []
    for text in encoded:
        word = decode_word(text)
        require(word[0] == 0, "declared translation normalization")
        keys.append(ordered_key(rows, word))
        require(not blue(rows, word), "black word required")
    require(len(set(keys)) == 7, "seven distinct full Aut orbits")
    return rows

def hamilton_words(rows):
    n = len(rows)
    return tuple(p for p in itertools.permutations(range(n))
                 if all(rows[a] >> b & 1 for a, b in zip(p, p[1:])))

def automorphisms(rows):
    n = len(rows)
    return tuple(p for p in itertools.permutations(range(n))
                 if all(((rows[i] >> j) & 1) == ((rows[p[i]] >> p[j]) & 1)
                        for i in range(n) for j in range(n)))

def small_controls():
    result = []
    for n, connection, expected in ((3, {1}, (3, 3, 1, 1)),
                                    (5, {1, 2}, (15, 5, 3, 3)),
                                    (7, {1, 2, 4}, (189, 21, 9, 3))):
        rows = tuple(sum(1 << j for j in range(n) if (j-i) % n in connection) for i in range(n))
        words = hamilton_words(rows)
        autos = automorphisms(rows)
        keys = {ordered_key(rows, w) for w in words}
        blue_keys = {ordered_key(rows, w) for w in words if blue(rows, w)}
        actual = (len(words), len(autos), len(keys), len(blue_keys))
        require(actual == expected, "independent small orbit census")
        # Independent orbit enumeration uses literal permutations, not adjacency keys.
        unseen = set(words)
        orbits = []
        while unseen:
            w = min(unseen)
            orbit = {tuple(p[v] for v in w) for p in autos}
            require(orbit <= unseen, "literal automorphism orbit partition")
            require(len({ordered_key(rows, v) for v in orbit}) == 1, "key on orbit")
            unseen -= orbit
            orbits.append(orbit)
        require(len(orbits) == len(keys), "keys iff full Aut orbits")
        result.append({"n": n, "H_A_h_b": actual})
    # A regular tournament need not be vertex-transitive. Its local 3-cycle
    # counts in outneighbourhoods are literal automorphism invariants.
    hostile = (14, 84, 56, 50, 97, 67, 13)
    validate_source(hostile)
    profile = []
    for v in range(7):
        neighbours = [j for j in range(7) if hostile[v] >> j & 1]
        profile.append(sum(all(sum((hostile[a] >> b) & 1 for b in s if b != a) == 1 for a in s)
                           for s in itertools.combinations(neighbours, 3)))
    require(all(r.bit_count() == 3 for r in hostile), "hostile regular")
    require(profile == [1, 1, 0, 1, 0, 1, 0], "regular is not transitive")
    return result, {"rows": hostile, "outneighbour_triangle_profile": profile}

def rejected(call):
    try:
        call()
    except (ValueError, TypeError):
        return
    raise ValueError("malformed certificate was accepted")

def run():
    global COUNT
    COUNT = 0
    data = (HERE / "certificates.txt").read_text(encoding="ascii").splitlines()
    require(len(data) == 8192, "complete two-group frame universe")
    source_hash = hashlib.sha256()
    seen_sources = set()
    total_nodes = max_nodes = max_words = 0
    first = None
    for index, line in enumerate(data):
        fields = line.split()
        require(len(fields) == 11, "literal frame schema")
        kind, mask = map(int, fields[:2])
        require((kind, mask) == divmod(index, 4096), "complete frame address")
        encoded = tuple(fields[2:9])
        if first is None:
            first = encoded
        rows = audit_frame(kind, mask, encoded)
        require(all(r.bit_count() == 12 for r in rows), "Cayley regularity")
        # Inversion is a literal anti-automorphism of each abelian Cayley source.
        inversion = tuple(difference(x, 0, kind) for x in range(25))
        require(all(((rows[x] >> y) & 1) != ((rows[inversion[x]] >> inversion[y]) & 1)
                    for x in range(25) for y in range(x+1, 25)), "inversion anti")
        raw = b"".join(r.to_bytes(4, "little") for r in rows)
        source_hash.update(raw)
        seen_sources.add(raw)
        nodes, words = map(int, fields[-2:])
        require(nodes > 0 and words >= 7, "nonproof discovery metadata")
        total_nodes += nodes
        max_nodes = max(max_nodes, nodes)
        max_words = max(max_words, words)
    controls, regular_hostile = small_controls()
    bad = [lambda: source(False, 0), lambda: source(0, 0.0),
           lambda: source(2, 0), lambda: source(0, 4096),
           lambda: audit_frame(0, 0, list(first)),
           lambda: audit_frame(0, 0, (first[0],) * 7),
           lambda: audit_frame(0, 0, (first[0][:-1] + first[0][0],) + first[1:]),
           lambda: audit_frame(0, 4095, first),
           lambda: ordered_key(source(0, 0), (False,) + tuple(range(1, 25))),
           lambda: ordered_key(source(0, 0), tuple(range(24)) + (-1,)),
           lambda: blue((2, 4, 1), (0.0, 1, 2)),
           lambda: ordered_key((True, 4, 1), (0, 1, 2))]
    for f in bad:
        rejected(f)
        require(True, "literal interface hostile")
    # This is the reduction from the frozen ten-order frontier, not a claim
    # about any concurrently produced order-twelve result.
    before = [12, 14, 15, 18, 20, 21, 24, 25, 28, 35]
    after = [n for n in before if n != 25]
    require(len(after) == 9 and after[-1] == 35, "own nine-order residual")
    return {"status": "FINITE-EXACT complete necessary Cayley-25 universe",
            "frames": 8192, "groups": ["Z/25", "(Z/5)^2"],
            "distinct_labelled_source_rows": len(seen_sources),
            "literal_black_words": 7 * 8192, "black_orbits_per_frame_at_least": 7,
            "source_adjacency_sha256": source_hash.hexdigest(),
            "certificate_lf_sha256": hashlib.sha256(("\n".join(data)+"\n").encode("ascii")).hexdigest(),
            "discovery_only": {"nodes": total_nodes, "max_nodes": max_nodes, "max_completed_words": max_words},
            "small_controls": controls, "regular_nontransitive_hostile": regular_hostile,
            "malformed_controls": len(bad), "residual_orders_after_this_packet_alone": after,
            "checks": COUNT}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--verify", action="store_true")
    args = p.parse_args()
    result = run()
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    path = HERE / "receipt.json"
    if args.verify:
        require(path.read_text(encoding="utf-8") == text, "saved receipt differs")
    else:
        path.write_text(text, encoding="utf-8", newline="\n")
    print(f"PASS {result['checks']} body checks; LF SHA256 {hashlib.sha256(text.encode()).hexdigest()}")

if __name__ == "__main__":
    main()
