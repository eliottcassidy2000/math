"""Independent exact controls for NATURAL_ROOTS.md; no external packages."""

from itertools import combinations, permutations
import json


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def labelled(n, code):
    rows = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        tail, head = (i, j) if code >> k & 1 else (j, i)
        rows[tail] |= 1 << head
    return tuple(rows)


def all_h(rows):
    n = len(rows)
    table = [0] * (1 << n)
    table[0] = 1
    ends = [[0] * n for _ in table]
    for v in range(n):
        ends[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            ways = ends[mask][v]
            if not ways:
                continue
            table[mask] += ways
            options = rows[v] & ~mask
            while options:
                bit = options & -options
                options -= bit
                ends[mask | bit][bit.bit_length() - 1] += ways
    return table


def literal_endpoints(rows):
    by_end = [0] * len(rows)
    for word in permutations(range(len(rows))):
        if all(rows[x] >> y & 1 for x, y in zip(word, word[1:])):
            by_end[word[-1]] += 1
    return by_end


def marked_g(rows, root, hs):
    full = (1 << len(rows)) - 1

    def walk(last, used):
        value = hs[full ^ used]
        options = rows[last] & ~used
        while options:
            bit = options & -options
            options -= bit
            value += walk(bit.bit_length() - 1, used | bit)
        return value

    return walk(root, 1 << root)


def double_birth(rows, root):
    n = len(rows)
    a, b = n, n + 1
    out = [row | (1 << a) for row in rows]
    for v in range(n):
        if v != root:
            out[v] |= 1 << b
    out.extend((1 << b, 1 << root))
    return tuple(out)


def is_left_padded_triangle(rows, root):
    n = len(rows)
    for triangle in combinations(range(n), 3):
        if root not in triangle:
            continue
        triangle_mask = sum(1 << x for x in triangle)
        if any((rows[x] & triangle_mask).bit_count() != 1 for x in triangle):
            continue
        prefix = tuple(v for v in range(n) if v not in triangle)
        if any(rows[v] & triangle_mask != triangle_mask for v in prefix):
            continue
        # A labelled tournament is transitive iff it has an all-forward order.
        if any(all(rows[x] >> y & 1 for x, y in combinations(word, 2))
               for word in permutations(prefix)):
            return True
    return False


def padded_triangle(left, right):
    n = left + 3 + right
    rows = [(1 << n) - (1 << (v + 1)) for v in range(n)]
    first, last = left, left + 2
    rows[first] ^= 1 << last
    rows[last] ^= 1 << first
    return tuple(rows)


def check_one(rows, root, base_hs=None, endpoints=None):
    n = len(rows)
    hs = all_h(rows) if base_hs is None else base_hs
    endpoints = literal_endpoints(rows) if endpoints is None else endpoints
    h, q, e = hs[-1], hs[((1 << n) - 1) ^ (1 << root)], endpoints[root]
    require(sum(endpoints) == h, ("base count", rows))
    g = marked_g(rows, root, hs)
    child = double_birth(rows, root)
    chs = all_h(child)
    child_endpoints = literal_endpoints(child)
    hp = chs[-1]
    require(sum(child_endpoints) == hp == h + 2 * g, ("whole count", rows, root))
    ga, gb, gu = (marked_g(child, v, chs) for v in (n, n + 1, root))
    require(ga == hp + h - e, ("new a", rows, root))
    require(gb == hp, ("new b", rows, root))
    require(gu == 4 * g - q, ("old u", rows, root))
    require(chs[((1 << (n + 2)) - 1) ^ (1 << root)] == q,
            ("old deletion", rows, root))
    classified = is_left_padded_triangle(rows, root)
    require((ga - hp == 2) == classified, ("classification", rows, root))
    if classified:
        require((h, g, e, hp) == (3, 3, 1, 9), ("exception", rows, root))
    if g - h == 2:
        require(ga - hp != 2 and gb - hp != 2 and gu - hp >= 4,
                ("prepared source", rows, root))
    return {"h": h, "g": g, "q": q, "endpaths": e, "child_h": hp,
            "defects_a_b_u": [ga - hp, gb - hp, gu - hp],
            "a_classified_ready": classified, "source_ready": g - h == 2}


def main():
    census = []
    for n in range(1, 5):
        checked = ready_sources = ready_a = 0
        for code in range(1 << (n * (n - 1) // 2)):
            rows = labelled(n, code)
            hs, endpoints = all_h(rows), literal_endpoints(rows)
            for root in range(n):
                result = check_one(rows, root, hs, endpoints)
                checked += 1
                ready_sources += result["source_ready"]
                ready_a += result["a_classified_ready"]
        census.append({"n": n, "rooted_labelled_cases": checked,
                       "ready_sources": ready_sources, "new_a_ready": ready_a})
    require(sum(row["rooted_labelled_cases"] for row in census) == 285,
            "incomplete universe")
    controls = []
    for left, right in ((0, 0), (1, 0), (3, 0), (0, 1), (1, 1), (2, 2)):
        rows = padded_triangle(left, right)
        hs, endpoints = all_h(rows), literal_endpoints(rows)
        for root in range(left, left + 3):
            result = check_one(rows, root, hs, endpoints)
            require(result["a_classified_ready"] == (right == 0),
                    ("padding boundary", left, right, root))
            controls.append({"left": left, "right": right, "root": root,
                             "rows": rows, **result})
    print(json.dumps({"status": "PASS", "scope": "N1-N3 finite controls",
                      "census": census, "padding_controls": controls}, indent=2))


if __name__ == "__main__":
    main()
