"""Independent exterior-word audit of the all-order two-run cone proof."""
from itertools import permutations, combinations
from pathlib import Path
import json
import check as base

HERE = Path(__file__).resolve().parent

def rows_of_code(n, code):
    rows = [0]*n
    for bit, (i, j) in enumerate(combinations(range(n), 2)):
        a, b = (i, j) if code >> bit & 1 else (j, i)
        rows[a] |= 1 << b
    return tuple(rows)

def hub_of(rows, root):
    incoming = [i for i in range(len(rows)) if rows[i] >> root & 1]
    outgoing = [o for o in range(len(rows)) if rows[root] >> o & 1]
    edges = [(o, i) for o in outgoing for i in incoming if rows[o] >> i & 1]
    if not edges or any(o != p and i != j for o, i in edges for p, j in edges):
        return None
    common_o = set(outgoing)
    for o, i in edges:
        common_o.intersection_update([o])
    if common_o:
        return rows, root, min(common_o), False
    n = len(rows)
    converse = tuple(((1 << n)-1) ^ row ^ (1 << i) for i, row in enumerate(rows))
    got = hub_of(converse, root)
    base.need(got is not None and not got[3], "star converse reduction")
    return got[0], root, got[2], True

def paths_on(rows, vertices):
    return [p for p in permutations(vertices) if all(rows[a] >> b & 1 for a, b in zip(p, p[1:]))]

def analyze(rows, root, full_skeleton=True):
    reduced = hub_of(rows, root)
    base.need(reduced is not None, "cyclic two-run guard")
    original = rows
    rows, root, hub, reversed_flag = reduced
    incoming = {i for i in range(len(rows)) if rows[i] >> root & 1}
    outgoing = {o for o in range(len(rows)) if rows[root] >> o & 1}
    exterior = sorted(incoming | outgoing)
    q0 = q1 = d0 = d1 = a1 = a2 = 0
    q1_words = set()
    for word in permutations(exterior):
        bad = [(a, b) for a, b in zip(word, word[1:]) if not (rows[a] >> b & 1)]
        backwards = sum(a in outgoing and b in incoming for a, b in zip(word, word[1:]))
        if not bad:
            base.need(backwards <= 1, "at most one legal back transition")
            if backwards:
                q1 += 1
                q1_words.add(word)
            else:
                q0 += 1
        elif len(bad) == 1 and bad[0][0] in incoming and bad[0][1] == hub:
            base.need(backwards <= 1, "at most one repaired back transition")
            if backwards:
                d1 += 1
            else:
                d0 += 1
        if full_skeleton:
            for j in (1, 2):
                for slots in combinations(range(len(word)+1), j):
                    expanded = []
                    for i in range(len(word)+1):
                        if i in slots:
                            expanded.append(root)
                        if i < len(word):
                            expanded.append(word[i])
                    if all(rows[a] >> b & 1 for a, b in zip(expanded, expanded[1:])):
                        if j == 1:
                            a1 += 1
                        else:
                            a2 += 1
    ipaths = paths_on(rows, sorted(incoming))
    opaths = paths_on(rows, sorted(outgoing))
    B = len(ipaths)*len(opaths)
    image = set()
    for pi in ipaths:
        split = next(i for i, v in enumerate(pi) if rows[hub] >> v & 1)
        for po in opaths:
            hub_index = po.index(hub)
            out = pi[:split]+po[:hub_index]+(hub,)+pi[split:]+po[hub_index+1:]
            base.need(out in q1_words, "literal injection lands in Q1")
            base.need(tuple(v for v in out if v in incoming) == pi, "I projection recovery")
            base.need(tuple(v for v in out if v in outgoing) == po, "O projection recovery")
            image.add(out)
    base.need(len(image) == B and q1 >= B > 0, "injective B-to-Q1 map")
    base.need(q0+d0 == B, "zero-transition partition")
    H = base.literal_h(original)
    base.need(H == q0+2*q1+d0+d1, "full Hamilton classification")
    if full_skeleton:
        base.need(a1 == H and a2 == q1+d1, "independent skeleton counts")
    else:
        a1, a2 = H, q1+d1
    base.need(a2+2 <= H <= 3*a2, "cone")
    base.need((H == 3*a2) == (q1 == B and d1 == 0), "upper equality boundary")
    return {"n": len(rows), "root": root, "hub_after_possible_converse": hub,
            "conversed": reversed_flag, "H": H, "a2": a2, "B": B,
            "classes": {"Q0": q0, "Q1": q1, "D0": d0, "D1": d1}}

def main():
    census = []
    for n in range(1, 6):
        count = upper = lower = 0
        for code in range(1 << (n*(n-1)//2)):
            rows = rows_of_code(n, code)
            for root in range(n):
                if hub_of(rows, root) is None:
                    continue
                result = analyze(rows, root)
                count += 1
                upper += result['H'] == 3*result['a2']
                lower += result['H'] == result['a2']+2
        census.append({"n": n, "rooted_cyclic_two_run": count, "upper_equality": upper, "lower_equality": lower})
    named = []
    for name, (rows, root) in [('S', base.SOURCE), ('B2', base.PLUS2), ('B4', base.PLUS4)]:
        result = analyze(rows, root)
        result['name'] = name
        named.append(result)
    report = {"status": "VERIFIED finite controls for independently audited TWO_RUN_CONE.md",
              "universe": "all rooted labelled tournaments n1..5 passing cyclic matching-star guard; three explicit incoming frames",
              "census": census, "named_frames": named}
    suffix = '' if __debug__ else '.optimized'
    (HERE / ('cone_receipt'+suffix+'.json')).write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({"status": "VERIFIED", "census": census, "named": [(r['name'], r['H'], r['a2']) for r in named]}))

if __name__ == '__main__':
    main()
