"""Exact full-count controls for terminal scale/capacity statements."""
from itertools import permutations
from math import comb, gcd
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
Q = (2, 4, 9, 3)
SOURCE = ((2, 20, 9, 35, 45, 7), 0)
PLUS2 = ((2, 12, 9, 49, 135, 215, 159, 15), 4)

def need(test, message):
    if not test:
        raise RuntimeError(message)

def tt(n):
    return tuple(((1 << n)-1) ^ ((1 << (v+1))-1) for v in range(n))

def ordinal(*blocks):
    n = sum(map(len, blocks))
    rows = []
    offset = 0
    for block in blocks:
        later = ((1 << n)-1) ^ ((1 << (offset+len(block)))-1)
        rows.extend((row << offset) | later for row in block)
        offset += len(block)
    return tuple(rows)

def apex(word):
    n = len(word)
    rows = [tt(n)[i] | ((b == '0') << n) for i, b in enumerate(word)]
    rows.append(sum((b == '1') << i for i, b in enumerate(word)))
    return tuple(rows)

PLUS4 = (apex('10010110'), 4)

def inflate(frame, root, module):
    groups = []
    offset = 0
    for v in range(len(frame)):
        size = len(module) if v == root else 1
        groups.append(list(range(offset, offset+size)))
        offset += size
    rows = [0] * offset
    for v, group in enumerate(groups):
        if v == root:
            for i, actual in enumerate(group):
                rows[actual] |= sum(1 << group[j] for j in range(len(module)) if module[i] >> j & 1)
        for w, other in enumerate(groups):
            if frame[v] >> w & 1:
                mask = sum(1 << actual for actual in other)
                for actual in group:
                    rows[actual] |= mask
    return tuple(rows)

def valid(rows):
    n = len(rows)
    return all(not (row >> i & 1) for i, row in enumerate(rows)) and all(
        (rows[i] >> j & 1) + (rows[j] >> i & 1) == 1 for i in range(n) for j in range(i))

def strong(rows):
    reach = [row | (1 << i) for i, row in enumerate(rows)]
    for k in range(len(rows)):
        for i in range(len(rows)):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return all(x == (1 << len(rows))-1 for x in reach)

def all_h(rows):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    h = [0] * (1 << n)
    h[0] = 1
    for mask in range(1, 1 << n):
        h[mask] = sum(dp[mask])
        for v, count in enumerate(dp[mask]):
            if count:
                todo = rows[v] & ~mask
                while todo:
                    bit = todo & -todo
                    todo -= bit
                    dp[mask | bit][bit.bit_length()-1] += count
    return h

def full_cover_profile(rows):
    h = all_h(rows)
    size = len(h)
    previous = [0] * size
    previous[0] = 1
    answer = [0]
    for j in range(1, len(rows)+1):
        current = [0] * size
        for mask in range(1, size):
            sub = mask
            while sub:
                current[mask] += h[sub] * previous[mask ^ sub]
                sub = (sub-1) & mask
        answer.append(current[-1])
        previous = current
    return answer

def slot(profile, q):
    return sum(comb(q, j)*value for j, value in enumerate(profile) if j <= q)

def literal_h(rows):
    n = len(rows)
    return sum(all(rows[p[i]] >> p[i+1] & 1 for i in range(n-1)) for p in permutations(range(n)))

def main():
    profile_q = full_cover_profile(Q)
    need(profile_q == [0, 5, 22, 36, 24], "Q full ordered cover profile")
    need(literal_h(Q) == 5, "independent permutation control")
    need(all(slot(profile_q, q) == q**4+4*q*q for q in range(9)), "Q slot polynomial")
    module_checks = []
    for r in (1, 2, 3):
        module = ordinal(*([Q]*r))
        profile = full_cover_profile(module)
        need(all(slot(profile, q) == (q**4+4*q*q)**r for q in range(1, 10)), "ordinal full-profile product")
        need(gcd(*profile) == 1, "primitive full cover content")
        module_checks.append({"r": r, "rows": module, "profile": profile, "content": gcd(*profile),
                              "h": profile[1], "p": slot(profile, 2)})
        for j in range(1, len(module)+1):
            inverted = sum((-1)**(j-q)*comb(j, q)*(q**4+4*q*q)**r for q in range(1, j+1))
            need(inverted == profile[j], "independent spectral/binomial inversion")
    frame_checks = []
    for r, padding in ((1, 0), (2, 0), (1, 2)):
        module = ordinal(*([ordinal(Q, tt(padding))]*r))
        actual = []
        records = []
        for name, (frame, root) in [('S', SOURCE), ('B2', PLUS2), ('B4', PLUS4)]:
            inflated = inflate(frame, root, module)
            need(valid(inflated) and strong(inflated), "literal strong inflation")
            count = all_h(inflated)[-1]
            actual.append(count)
            records.append({"frame": name, "n": len(inflated), "H_DP": count,
                            "rows_sha256": hashlib.sha256(str(inflated).encode()).hexdigest()})
        p = 2**((5+padding)*r)
        expected = [15*p-5**r, 45*p-5**r, 45*p+5**r]
        need(actual == expected, "full approximate-transition counts")
        frame_checks.append({"r": r, "padding_per_atom": padding, "literal_frames": records, "residuals": [actual[1]-3*actual[0], actual[2]-3*actual[0]]})
    # Whole-profile order invariance, not merely equality of h or p.
    reorders = [ordinal(Q, tt(2), Q), ordinal(tt(1), Q, Q, tt(1)), ordinal(Q, Q, tt(2))]
    profiles = [full_cover_profile(m) for m in reorders]
    need(all(p == profiles[0] for p in profiles), "all covers survive atom reordering")
    clusters = []
    for r in range(1, 9):
        values = sorted({45*32**r + sign*5**j for sign in (-1, 1) for j in range(r+1)})
        need(len(values) == 2*(r+1), "cluster cardinality")
        need(max(b-a for a, b in zip(values, values[1:])) == 4*5**(r-1), "exact cluster gap")
        clusters.append({"r": r, "count": len(values), "min": values[0], "max": values[-1],
                         "largest_gap": 4*5**(r-1)})
    equal_slot = []
    for r in (1, 2):
        for j in range(r+1):
            module = ordinal(*([Q]*j), tt(5*(r-j)))
            h = all_h(module)
            p = sum(h[s]*h[(len(h)-1)^s] for s in range(len(h)))
            need(h[-1] == 5**j and p == 32**r, "literal equal-two-slot family")
            equal_slot.append({"r": r, "j": j, "n": len(module), "h_DP": h[-1], "p_DP": p})
    report = {"status": "VERIFIED finite controls; universal proofs in PROOF.md",
              "universe": "Q ordinal powers r1..3; all covers; three inherited strong frames with (r,padding)=(1,0),(2,0),(1,2); equal-p modules r1..2; clusters r1..8",
              "module_checks": module_checks, "frame_checks": frame_checks,
              "reordered_full_profile": profiles[0], "equal_slot_checks": equal_slot, "clusters": clusters}
    no_go_controls = []
    for r in range(3, 13):
        h, p = 5**r, 32**r
        k3 = 117**r-3*p+3*h
        target = 45*p-3*h
        need(h % 5 == 0 and target % 5 == 0 and all((target+d) % 5 != 0 for d in (2,4)),
             "correct acyclic-root obstruction is modulo five, not h")
        base_checks = (k3 > 45*p, p > 47*h+4, h > 4)
        need(all(base_checks), "all-frame no-go inequalities")
        no_go_controls.append({"r": r, "h": h, "p": p, "K3": k3, "targets": [target+2, target+4]})
    need(module_checks[-1]['profile'][3] == no_go_controls[0]['K3'], "literal K3 no-go control")
    report['all_frame_no_go_controls'] = no_go_controls
    small_cone = []
    for r in (1, 2):
        h, p = 5**r, 32**r
        for delta in (2, 4):
            target = 45*p-3*h+delta
            lo = (target+p+h-1)//(p+h)
            hi = (target-2*h)//(p-h)
            candidates = []
            for a in range(lo, hi+1):
                if (target-a*p) % h == 0:
                    b = (target-a*p)//h
                    if -a+2 <= b <= a:
                        candidates.append([a,b])
            if r == 2:
                need(not candidates, "r2 two-run-only obstruction")
            small_cone.append({"r": r, "delta": delta, "a_range": [lo,hi], "candidate_pairs": candidates})
    report['small_parameter_cone_boundaries'] = small_cone
    suffix = '' if __debug__ else '.optimized'
    (HERE / ('receipt'+suffix+'.json')).write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({"status": "VERIFIED", "frame_counts": [[x['H_DP'] for x in r['literal_frames']] for r in frame_checks],
                      "primitive_cover_contents": [r['content'] for r in module_checks]}))

if __name__ == '__main__':
    main()
