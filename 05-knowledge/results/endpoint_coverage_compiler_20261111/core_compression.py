"""Independent finite controls and all-order padding-transform coefficients."""
from collections import Counter,defaultdict
from itertools import combinations
from math import comb
from pathlib import Path
import json
import inflation_bank as ib
import construct as con

HERE = Path(__file__).resolve().parent


def components(rows):
    n = len(rows)
    reach = [rows[v] | (1 << v) for v in range(n)]
    for k in range(n):
        for v in range(n):
            if reach[v] >> k & 1:
                reach[v] |= reach[k]
    unused, out = set(range(n)),[]
    while unused:
        v = min(unused)
        group = sorted(w for w in unused if reach[v]>>w&1 and reach[w]>>v&1)
        out.append(group)
        unused.difference_update(group)
    return sorted(out,key=lambda group:-sum(bool(rows[group[0]] >> other[0] & 1) for other in out if other != group))


def main():
    records = json.loads((HERE/'inflation_bank_7.json').read_text())['records']
    core_k, all_k, core_p, all_p, core_q, all_q = (Counter() for _ in range(6))
    marks = defaultdict(list)
    for item in records:
        n = item['n']
        all_p[n] += 1
        strong = len(components(item['rows'])) == 1
        if strong:
            core_p[n] += 1
        marks[n,item['code'],item['selected']].append(item)
    for (n,code,selected),items in marks.items():
        pairs = comb(len(items),2)
        all_q[n] += pairs
        if len(components(items[0]['rows'])) == 1:
            core_q[n] += pairs
    bank, controls = [[]],0
    for n in range(1,8):
        codes = {ib.cc.canonical(ib.cc.child(r,m)) for r in bank for m in range(1 << (n-1))}
        bank = [ib.cc.decode(n,c) for c in sorted(codes)]
        for rows in bank:
            h,u,w,_,_,_ = ib.tc.path_tables(rows)
            parts = components(rows)
            for mask in range(1,1 << n):
                selected = [v for v in range(n) if mask >> v & 1]
                if h-sum(w[a][b] for a,b in combinations(selected,2)) != 1:
                    continue
                all_k[n] += 1
                if len(parts) == 1:
                    core_k[n] += 1
                active = [group for group in parts if any(v in selected for v in group)]
                ib.tc.require(len(active) == 1,'one selected SCC')
                core = active[0]
                ib.tc.require(all(group == core or len(group) == 1 for group in parts),'singleton exterior')
                subrows = [sum(1<<j for j,b in enumerate(core) if rows[a]>>b&1) for a in core]
                subselected = [i for i,v in enumerate(core) if v in selected]
                ch,cu,cw,_,_,_ = ib.tc.path_tables(subrows)
                ib.tc.require(ch-sum(cw[a][b] for a,b in combinations(subselected,2)) == 1,'core unit offset')
                trailing = sum(len(group) for group in parts[parts.index(core)+1:])
                ib.tc.require(sum(u[v] for v in selected) == (1<<trailing)*sum(cu[v] for v in subselected),'exact sink power')
                controls += 1
    for n in range(1,8):
        k = sum((n-m+1)*core_k[m] for m in range(1,n+1))
        p = sum((n-m+1)*(core_p[m]+(n-m)*core_k[m]) for m in range(1,n+1))
        q = sum((n-m+1)*(core_q[m]+(n-m)*core_p[m]+comb(n-m,2)*core_k[m]) for m in range(1,n+1))
        ib.tc.require((k,p,q) == (all_k[n],all_p[n],all_q[n]),'all three padding transforms')
    for family in con.FAMILIES:
        rows,_,_ = con.seed(family,1)
        ib.tc.require(len(components(rows)) == 1,'named base is itself strong')
    result = dict(ordinary_unit_seed_controls=controls,
                  all_seed_presentations=[all_k[n] for n in range(1,8)],
                  strong_core_seed_presentations=[core_k[n] for n in range(1,8)],
                  all_stable_ports=[all_p[n] for n in range(1,8)],
                  strong_core_stable_ports=[core_p[n] for n in range(1,8)],
                  all_stable_pairs=[all_q[n] for n in range(1,8)],
                  strong_core_stable_pairs=[core_q[n] for n in range(1,8)],
                  named_bases_strong=True)
    (HERE/'core_compression.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
