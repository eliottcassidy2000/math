"""Complete H<=9 SCC-word family and finite local-operation dictionary.

Run from any directory: python -B classify_cap9.py [--output receipt.json].
No external inputs. All graphs are outgoing row masks; ordinary isomorphism
is used, not converse merging. Native lines use Hamilton-path frames.
"""
from __future__ import annotations
import argparse
import itertools as it
import json
from collections import Counter
from functools import lru_cache
from math import comb, prod
from pathlib import Path

CORES = {
    "x": (0,), "A": (2, 4, 1), "B": (2, 4, 9, 3),
    "C": (2, 4, 9, 19, 7), "D": (2, 4, 25, 3, 11),
}
WEIGHTS = {"x": 1, "A": 3, "B": 5, "C": 9, "D": 9}
KERNELS = ((), ("A",), ("B",), ("C",), ("D",), ("A", "A"))

def need(condition, message):
    if not condition:
        raise RuntimeError(message)

def falling(n, r):
    return prod(range(n-r+1, n+1)) if n >= r else 0

def cycles(p):
    unseen = set(range(len(p)))
    lengths = []
    while unseen:
        v = min(unseen)
        cycle = []
        while v not in cycle:
            cycle.append(v)
            v = p[v]
        unseen.difference_update(cycle)
        lengths.append(len(cycle))
    return lengths

def relabel(rows, order):
    return tuple(sum(1 << j for j, v in enumerate(order)
                     if (rows[u] >> v) & 1) for u in order)

def induced(rows, keep):
    return relabel(rows, tuple(keep))

def reverse(rows):
    mask = (1 << len(rows)) - 1
    return tuple(mask ^ (1 << i) ^ row for i, row in enumerate(rows))

@lru_cache(None)
def canonical(rows):
    return min(relabel(rows, p) for p in it.permutations(range(len(rows))))

CANON_TO_SYMBOL = {canonical(rows): c for c, rows in CORES.items()}

@lru_cache(None)
def automorphisms(rows):
    return tuple(p for p in it.permutations(range(len(rows)))
                 if relabel(rows, p) == rows)

def orbits(items, group, action):
    unseen = set(items)
    answer = []
    while unseen:
        seed = min(unseen)
        orbit = {action(p, seed) for p in group}
        unseen.difference_update(orbit)
        answer.append(sorted(orbit))
    return answer

def vertex_orbits(rows):
    return orbits(range(len(rows)), automorphisms(rows), lambda p, v: p[v])

def hp(rows):
    n = len(rows)
    if not n:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for v, count in enumerate(dp[mask]):
            if count:
                nxt = rows[v] & (full ^ mask)
                while nxt:
                    bit = nxt & -nxt
                    nxt -= bit
                    w = bit.bit_length() - 1
                    dp[mask | bit][w] += count
    return sum(dp[full])

def hp_permutations(rows):
    return sum(all((rows[p[i]] >> p[i+1]) & 1
                   for i in range(len(rows)-1))
               for p in it.permutations(range(len(rows))))

def sccs(rows):
    n = len(rows)
    reach = [row | (1 << i) for i, row in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if (reach[i] >> k) & 1:
                reach[i] |= reach[k]
    unused = set(range(n))
    parts = []
    while unused:
        v = min(unused)
        part = {u for u in unused if (reach[v] >> u) & 1
                and (reach[u] >> v) & 1}
        unused -= part
        parts.append(tuple(sorted(part)))
    # Between different SCCs all arcs have a common orientation.
    return sorted(parts, key=lambda a: -sum((rows[a[0]] >> b[0]) & 1
                                            for b in parts if b != a))

def word_of(rows):
    return tuple(CANON_TO_SYMBOL[canonical(induced(rows, part))]
                 for part in sccs(rows))

def from_word(word):
    n = sum(len(CORES[c]) for c in word)
    rows = []
    offset = 0
    for c in word:
        core = CORES[c]
        after = ((1 << n) - 1) ^ ((1 << (offset + len(core))) - 1)
        rows.extend((row << offset) | after for row in core)
        offset += len(core)
    return tuple(rows)

def compositions(n, r):
    if r == 1:
        yield (n,)
    else:
        for a in range(n+1):
            for tail in compositions(n-a, r-1):
                yield (a,) + tail

def family_words(n):
    for kernel in KERNELS:
        size = sum(len(CORES[c]) for c in kernel)
        if size > n:
            continue
        for gaps in compositions(n-size, len(kernel)+1):
            word = []
            for i, c in enumerate(kernel):
                word += ["x"] * gaps[i] + [c]
            word += ["x"] * gaps[-1]
            yield tuple(word)

def word_aut(word):
    groups, offset = [], 0
    for c in word:
        groups.append(tuple(tuple(offset+v for v in p)
                            for p in automorphisms(CORES[c])))
        offset += len(CORES[c])
    return tuple(tuple(v for part in choice for v in part)
                 for choice in it.product(*groups))

def birth_rules():
    rules = []
    for c in "ABCD":
        rows = CORES[c]
        group = automorphisms(rows)
        for orbit in vertex_orbits(rows):
            v = orbit[0]
            keep = tuple(u for u in range(len(rows)) if u != v)
            parent = induced(rows, keep)
            pword = word_of(parent)
            stabilizer = len(group) // len(orbit)
            rules.append(dict(source="".join(pword), target=c, root=v,
                              root_orbit=orbit, child_rows=rows,
                              parent_rows=parent,
                              parent_aut=len(word_aut(pword)),
                              child_root_stabilizer=stabilizer,
                              named_attachment_weight=len(word_aut(pword))//stabilizer,
                              named_deletion_weight=len(orbit)))
    return rules

def arc_rules():
    rules = []
    for n in range(2, 6):
        for word in family_words(n):
            rows = from_word(word)
            owner = [i for i, c in enumerate(word) for _ in CORES[c]]
            candidates = [(u, v) for u in range(n) for v in range(u+1, n)
                          if owner[u] == 0 and owner[v] == len(word)-1]
            group = word_aut(word)
            for orbit in orbits(candidates, group,
                                lambda p, e: tuple(sorted((p[e[0]], p[e[1]])))):
                u, v = orbit[0]
                target = list(rows)
                target[u] ^= 1 << v
                target[v] ^= 1 << u
                target = tuple(target)
                count = hp(target)
                if count <= 9:
                    rules.append(dict(source="".join(word),
                                      target="".join(word_of(target)),
                                      arc=[u, v], arc_orbit_size=len(orbit),
                                      source_h=hp(rows), target_h=count,
                                      source_rows=rows, target_rows=target))
    return rules

def native_lines():
    lines = []
    for n in range(3, 6):
        free = [(i, j) for i in range(n) for j in range(i+2, n)]
        full = (1 << len(free)) - 1
        def frame(code):
            rows = [0] * n
            for i in range(n-1):
                rows[i] |= 1 << (i+1)
            for k, (i, j) in enumerate(free):
                if (code >> k) & 1:
                    rows[i] |= 1 << j
                else:
                    rows[j] |= 1 << i
            return tuple(rows)
        pos = {edge: i for i, edge in enumerate(free)}
        for code in range((full+1)//2):
            rows, partner = frame(code), frame(full ^ code)
            h, h2 = hp(rows), hp(partner)
            if max(h, h2) > 9:
                continue
            reflected = sum(((code >> k) & 1) << pos[(n-1-j, n-1-i)]
                            for k, (i, j) in enumerate(free))
            lines.append(dict(n=n, frame=code, complement=full ^ code,
                              source="".join(word_of(rows)),
                              target="".join(word_of(partner)),
                              source_h=h, target_h=h2,
                              color="blue" if reflected == code else "black"))
    return lines

def all_labelled(n):
    pairs = list(it.combinations(range(n), 2))
    for code in range(1 << len(pairs)):
        rows = [0] * n
        for k, (u, v) in enumerate(pairs):
            if (code >> k) & 1:
                rows[u] |= 1 << v
            else:
                rows[v] |= 1 << u
        yield tuple(rows)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    core_data = {c: dict(rows=r, h=hp(r), aut=len(automorphisms(r)),
                         vertex_orbits=vertex_orbits(r), self_converse=canonical(r)==canonical(reverse(r)))
                 for c, r in CORES.items()}
    need([core_data[c]["h"] for c in CORES] == [1,3,5,9,9], "core path counts")
    need(all(v["self_converse"] for v in core_data.values()), "core converse types")
    universe = []
    for n in range(1, 7):
        seen = set()
        total = 0
        for rows in all_labelled(n):
            h = hp(rows)
            if n <= 5:
                need(h == hp_permutations(rows), "literal permutation count")
            if h <= 9:
                seen.add(word_of(rows))
                total += 1
        expected = set(family_words(n))
        need(seen == expected, f"class mismatch at {n}: {seen ^ expected}")
        universe.append(dict(n=n, labelled_tested=1 << comb(n,2),
                             labelled_under_cap=total, classes=len(seen)))
    rows_table = []
    for n in range(0, 51):
        words = list(family_words(n))
        aut_hist = Counter()
        rooted = frames = named_h_mass = self_converse = attachment_orbits = 0
        named_marked = [0,0,0,0]
        for word in words:
            group = word_aut(word)
            a = len(group)
            h = 1
            for c in word:
                h *= WEIGHTS[c]
            aut_hist[a] += 1
            rooted += sum(len(core_data[c]["vertex_orbits"]) for c in word)
            frames += h // a
            named_h_mass += h
            self_converse += word == word[::-1]
            attachment_orbits += sum(2**len(cycles(p)) for p in group)//a
            for r in range(4):
                named_marked[r] += sum(falling(sum(i==v for i,v in enumerate(p)),r)
                                       for p in group)//a
        if n >= 4:
            need(len(words) == (n*n-n-4)//2, "ordinary series")
            need(rooted == (n**3-5*n*n+28*n-72)//2, "rooted series")
            need(frames == (n*n+39*n-156)//2, "frame series")
        if n >= 6:
            need(attachment_orbits == 2**(n-3)*(n*n+19*n-68), "attachment orbit series")
            for r in range(4):
                formula = ((3*n-10)*falling(n,r)
                    +(n-2)*(falling(n,r)+2*falling(n-3,r))//3
                    +comb(n-4,2)*(falling(n,r)+4*falling(n-3,r)+4*falling(n-6,r))//9)
                need(named_marked[r] == formula, "ordered named marking series")
        rows_table.append(dict(n=n, ordinary=len(words), rooted=rooted,
                               frames=frames, named_h_sum=named_h_mass,
                               self_converse=self_converse,
                               converse_merged=(len(words)+self_converse)//2,
                               attachment_orbits=attachment_orbits,
                               ordered_named_markings_0_to_3=named_marked,
                               aut_hist=dict(aut_hist)))
    br, ar, native = birth_rules(), arc_rules(), native_lines()
    # Independent literal check of every one-vertex deletion through n=10.
    # The local birth dictionary must regenerate all strong-child deletions.
    rule_support = {(r["source"], r["target"]) for r in br}
    delete_checks = 0
    for n in range(1,11):
        for word in family_words(n):
            rows = from_word(word)
            for v in range(n):
                child = induced(rows, [u for u in range(n) if u != v])
                need(hp(child) <= hp(rows) <= 9, "deletion monotonicity")
                need(word_of(child) in set(family_words(n-1)), "deletion family")
                owner = 0
                while v >= sum(len(CORES[c]) for c in word[:owner+1]):
                    owner += 1
                symbol = word[owner]
                if symbol != "x":
                    expected_middle = word_of(induced(CORES[symbol],
                         [u for u in range(len(CORES[symbol]))
                          if u != v-sum(len(CORES[c]) for c in word[:owner])]))
                    need(("".join(expected_middle), symbol) in rule_support, "deletion rule")
                delete_checks += 1
    receipt = dict(status="PASS", cap=9, convention="ordinary classes; x=singleton; native literal lines n>=3",
                   cores=core_data, labelled_universe=universe,
                   literal_deletion_checks=delete_checks,
                   strong_birth_rules=br, arc_rules=ar, native_lines=native,
                   series_rows=rows_table,
                   boundaries={"all_order_strong_classification":"inherited audited five-core theorem",
                               "all_order_rewriting":"proved in FAMILIES.md; enumeration only audits finite dictionary",
                               "native_order_cutoff":5,
                               "unbounded_H_classification":"OPEN"})
    payload = json.dumps(receipt, indent=2) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
        columns = ("n","ordinary","rooted","frames","self_converse",
                   "converse_merged","attachment_orbits","two_named","three_named")
        lines = ["\t".join(columns)]
        for row in rows_table:
            out = dict(row,two_named=row["ordered_named_markings_0_to_3"][2],
                       three_named=row["ordered_named_markings_0_to_3"][3])
            lines.append("\t".join(str(out[c]) for c in columns))
        args.output.with_name("sequences.tsv").write_text("\n".join(lines)+"\n",encoding="utf-8")
    print(json.dumps({"status":"PASS", "labelled_graphs":sum(x["labelled_tested"] for x in universe),
                      "birth_rules":len(br), "arc_rules":len(ar),
                      "native_lines":len(native), "literal_deletion_checks":delete_checks,
                      "series_through":50}))

if __name__ == "__main__":
    main()
