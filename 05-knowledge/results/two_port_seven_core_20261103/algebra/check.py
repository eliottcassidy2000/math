"""Exact full-response content and ambient-converse controls; stdlib only."""
from hashlib import sha256
from itertools import combinations, permutations, product
from math import factorial, gcd
from pathlib import Path
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def tournament(n, code):
    rows = [0] * n
    for bit, (u, v) in enumerate(combinations(range(n), 2)):
        if code >> bit & 1:
            rows[u] |= 1 << v
        else:
            rows[v] |= 1 << u
    return rows


def opposite(rows):
    full = (1 << len(rows)) - 1
    return [full ^ row ^ (1 << u) for u, row in enumerate(rows)]


def ends(rows):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    for u in range(n):
        dp[1 << u][u] = 1
    for mask in range(1, 1 << n):
        for u, count in enumerate(dp[mask]):
            if not count:
                continue
            targets = rows[u] & ~mask
            while targets:
                bit = targets & -targets
                targets -= bit
                dp[mask | bit][bit.bit_length() - 1] += count
    return dp


def h(rows):
    return sum(ends(rows)[-1]) if rows else 1


def literal_h(rows):
    return sum(all(rows[u] >> v & 1 for u, v in zip(p, p[1:]))
               for p in permutations(range(len(rows))))


def sidecar(rows):
    n, full = len(rows), (1 << len(rows)) - 1
    e, s = ends(rows), ends(opposite(rows))
    c = [[0] * n for _ in rows]
    for mask in range(1, full):
        for u, eu in enumerate(e[mask]):
            if not eu:
                continue
            for v, sv in enumerate(s[full ^ mask]):
                c[u][v] += eu * sv
    return sum(e[-1]) if rows else 1, s[-1], e[-1], c


def join(blocks):
    out, shift = [], 0
    n = sum(map(len, blocks))
    for block in blocks:
        after = ((1 << n) - 1) ^ ((1 << (shift + len(block))) - 1)
        out.extend((row << shift) | after for row in block)
        shift += len(block)
    return out


def inflate(q, module):
    r = len(module)
    out = module[:] + [0] * (len(q) - 1)
    for u in range(len(q)):
        us = range(r) if u == 0 else (r + u - 1,)
        for v in range(len(q)):
            if q[u] >> v & 1:
                vs = range(r) if v == 0 else (r + v - 1,)
                for a in us:
                    for b in vs:
                        out[a] |= 1 << b
    return out


def extend(rows, masks):
    n = len(rows)
    out = rows[:] + [0] * len(masks)
    for j, mask in enumerate(masks):
        for u in range(n):
            if mask >> u & 1:
                out[u] |= 1 << (n + j)
            else:
                out[n + j] |= 1 << u
        for k in range(j + 1, len(masks)):
            out[n + j] |= 1 << (n + k)
    return out


TRIANGLE = [2, 4, 1]


def profile(rows):
    """Ordered covers, with one fixed-slot endpoint; independent subset DP."""
    n, full = len(rows), (1 << len(rows)) - 1
    e, b = ends(rows), ends(opposite(rows))
    hs = [sum(x) for x in e]
    q = [[0] * (n + 1) for _ in range(1 << n)]
    q[0][0] = 1
    for mask in range(1, 1 << n):
        sub = mask
        while sub:
            for k in range(1, mask.bit_count() + 1):
                q[mask][k] += hs[sub] * q[mask ^ sub][k - 1]
            sub = (sub - 1) & mask
    es, bs = [[0] * n for _ in range(n + 1)], [[0] * n for _ in range(n + 1)]
    for mask in range(1, 1 << n):
        for k in range(1, n + 1):
            mult = q[full ^ mask][k - 1]
            for u in range(n):
                es[k][u] += e[mask][u] * mult
                bs[k][u] += b[mask][u] * mult
    for k in range(2, n + 1):
        for u in range(n):
            need((es[k][u] - bs[k][u]) % 2 == 0, "fixed-slot parity")
            need(es[k][u] % factorial(k - 1) == 0, "fixed-slot normalization")
    return es, bs


def channels(q):
    """Literal exterior run skeletons, end channel and each outgoing derivative."""
    m = len(q) - 1
    result = [[0] * (m + 2) for _ in q]
    for order in permutations(range(1, m + 1)):
        for gaps in range(1, 1 << (m + 1)):
            for channel in range(len(q)):
                if channel == 0:
                    if not gaps >> m & 1:
                        continue
                elif not gaps >> order.index(channel) & 1:
                    continue
                valid = True
                for j in range(m + 1):
                    if gaps >> j & 1:
                        if j and not q[order[j - 1]] & 1:
                            valid = False
                        if j < m and order[j] != channel and not q[0] >> order[j] & 1:
                            valid = False
                    elif 0 < j < m and not q[order[j - 1]] >> order[j] & 1:
                        valid = False
                if valid:
                    result[channel][gaps.bit_count()] += 1
    hq, sq, eq, cq = sidecar(q)
    need(result[0][1] == eq[0], "leading end channel")
    for w in range(1, len(q)):
        need(result[w][1] == cq[0][w], "leading cross channel")
    need(eq[0] + sum(cq[0][w] for w in range(1, len(q)) if q[0] >> w & 1) == hq,
         "quotient odd channel identity")
    return result


def ambient_controls():
    modules = [tournament(n, code) for n in range(1, 5)
               for code in range(1 << (n * (n - 1) // 2))]
    quotients = [tournament(n, code) for n in range(1, 4)
                 for code in range(1 << (n * (n - 1) // 2))]
    hostile = [4, 17, 42, 131, 45, 75, 31, 119]
    extras = [tournament(5, 152), tournament(6, 689), hostile, opposite(hostile)]
    # Alternating transitive exterior permits three or four nonempty module runs.
    deep = []
    for qn in (5, 7):
        q = [0] * qn
        for u in range(1, qn):
            for v in range(u + 1, qn):
                q[u] |= 1 << v
            if u % 2:
                q[0] |= 1 << u
            else:
                q[u] |= 1
        deep.append(q)
    cases = [(module, q) for module in modules + extras for q in quotients]
    cases += [(module, q) for module in modules[-8:] + extras[:2] for q in deep]
    profiles, skeletons = {}, {}
    channels_checked = skew_entries = 0
    largest_run = 0
    for module, q in cases:
        key = tuple(module)
        if key not in profiles:
            profiles[key] = profile(module), sidecar(module)
        (es, bs), (_, sm, em, cm) = profiles[key]
        if tuple(q) not in skeletons:
            skeletons[tuple(q)] = channels(q)
        ds = skeletons[tuple(q)]
        ht, st, et, ct = sidecar(inflate(q, module))
        hp, sp, ep, cp = sidecar(inflate(q, opposite(module)))
        need(ht == hp, "module reversal H")
        r = len(module)
        for ch, coeff in enumerate(ds):
            largest_run = max(largest_run, max((k for k, val in enumerate(coeff) if val), default=0))
            for u in range(r):
                predicted = sum(coeff[k] * es[k][u] for k in range(1, min(len(coeff), r + 1)))
                reverse = sum(coeff[k] * bs[k][u] for k in range(1, min(len(coeff), r + 1)))
                actual = et[u] if ch == 0 else ct[u][r + ch - 1]
                other = ep[u] if ch == 0 else cp[u][r + ch - 1]
                need(predicted == actual and reverse == other, "exact ambient fixed-slot channel")
                need((actual - other - coeff[1] * (em[u] - sm[u])) % 2 == 0,
                     "ambient endpoint descent")
                channels_checked += 1
        for u in range(r):
            for v in range(r):
                need(cp[u][v] == ct[v][u], "internal converse transpose")
                need((ct[u][v] - ct[v][u] - h(q) * (cm[u][v] - cm[v][u])) % 4 == 0,
                     "ambient ordered-cover mod4 descent")
                skew_entries += 1
    hh, sh, eh, ch = sidecar(hostile)
    need(all((x - y) % 2 == 0 for x, y in zip(sh, eh)), "hostile endpoint parity")
    need(ch[0][2] == 162 and ch[2][0] == 12, "hostile nonzero ordered-cover residue")
    return {"ambient_pairs": len(cases), "fixed_slot_channels": channels_checked,
            "internal_skew_entries": skew_entries, "largest_run_tested": largest_run,
            "complete_core_orders": [1, 4], "complete_quotient_orders": [1, 3],
            "additional_cores": extras, "deep_quotients": deep,
            "hostile": {"rows": hostile, "H": hh, "S": sh, "E": eh,
                        "C_02": ch[0][2], "C_20": ch[2][0]}}


def content_controls():
    records = []
    for s, m in ((0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2),
                 (2, 0), (2, 1), (2, 2), (3, 0), (3, 1)):
        rows = join([TRIANGLE] * s)
        n = 3 * s
        values = [h(extend(rows, [(code >> (j * n)) & ((1 << n) - 1) for j in range(m)]))
                  for code in range(1 << (m * n))]
        content = 0
        for value in values:
            content = gcd(content, value)
        need(content == 3 ** max(s - 2 * m, 0), "complete content table")
        coeff = values[:]
        for bit in range(m * n):
            for code in range(1 << (m * n)):
                if code >> bit & 1:
                    coeff[code] -= coeff[code ^ (1 << bit)]
        nonzero = 0
        for code, value in enumerate(coeff):
            hit = {((bit % n) // 3) for bit in range(m * n) if code >> bit & 1}
            need(value % (3 ** (s - len(hit))) == 0, "pointwise stabilizer divisibility")
            need(not value or code.bit_count() <= 2 * m, "full polynomial total degree")
            for j in range(m):
                need(not value or ((code >> (j * n)) & ((1 << n) - 1)).bit_count() <= 2,
                     "full polynomial port degree")
            nonzero += value != 0
        records.append({"s": s, "m": m, "queries": len(values), "gcd": content,
                        "nonzero_coefficients": nonzero,
                        "values_sha256": sha256(json.dumps(values).encode()).hexdigest()})
    # Literal permutation check of the local constants used in all-order sharpness.
    local = []
    for s, masks in ((1, [0, 1]), (2, [0, 1, 8, 9])):
        rows = join([TRIANGLE] * s)
        vals = [literal_h(extend(rows, [a])) for a in masks]
        need(vals == ([3, 5] if s == 1 else [9, 15, 45, 41]), "literal local coefficient")
        local.append({"triangles": s, "masks": masks, "values": vals})
    sharp = []
    for s, m in ((0, 3), (1, 3), (2, 3), (3, 2), (4, 2), (5, 2)):
        rows, n = join([TRIANGLE] * s), 3 * s
        t = min(s, 2 * m)
        masks = [(1 << (6 * j)) - 1 for j in range((t + 1) // 2)]
        masks.extend([(1 << n) - 1] * (m - len(masks)))
        actual = []
        coefficient = 0
        for chosen in range(1 << t):
            query = masks[:]
            for triangle in range(t):
                if chosen >> triangle & 1:
                    query[triangle // 2] |= 1 << (3 * triangle)
            value = h(extend(rows, query))
            actual.append({"masks": query, "value": value})
            coefficient += (-1) ** (t - chosen.bit_count()) * value
        prediction = 3 ** (s - t) * (-10) ** (t // 2) * 2 ** (t % 2)
        need(coefficient == prediction, "factorized full response sharpness coefficient")
        need(gcd(3 ** s, coefficient) == 3 ** max(s - 2 * m, 0), "content sharpness certificate")
        sharp.append({"s": s, "m": m, "coefficient": coefficient, "queries": actual})
    return {"complete_tables": records, "literal_local_constants": local, "sharpness": sharp}


def multiport_controls():
    count = 0
    for n in range(1, 4):
        for code in range(1 << (n * (n - 1) // 2)):
            b = tournament(n, code)
            hb, sb, eb, cb = sidecar(b)
            for mask in range(1, 1 << n):
                w = [v for v in range(n) if mask >> v & 1]
                j = sum(cb[u][v] for u in w for v in w)
                need(j >= len(w) * (len(w) - 1), "physical-port positive defect")
                need((j == 0) == (len(w) == 1), "exact singleton zero-ready boundary")
                fw = sum(eb[u] + sum(cb[u]) for u in w)
                for a in ([0], TRIANGLE):
                    ha, sa, ea, ca = sidecar(a)
                    pa2 = 2 * ha + sum(map(sum, ca))
                    extension = extend(join([a, b]), [mask << len(a)])
                    hm, sm, em, cm = sidecar(extension)
                    er2 = em[-1] + sum(cm[-1])
                    need(hm == ha * (hb - j) + pa2 * fw, "multi-in-port H identity")
                    need(er2 - hm == ha * j, "multi-in-port root defect identity")
                    count += 1
    return {"complete_suffix_orders": [1, 3], "prefixes": [[0], TRIANGLE], "cases": count}


def main():
    path = Path(__file__)
    source = sha256(path.read_bytes()).hexdigest()
    result = {"source_sha256": source, "status": "FINITE-EXACT controls of proved B1-B5; universal two-port OPEN",
              "content": content_controls(), "ambient": ambient_controls(),
              "multiport": multiport_controls()}
    need(source == sha256(path.read_bytes()).hexdigest(), "source changed while running")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
