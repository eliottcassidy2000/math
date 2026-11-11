"""Probe whether separately stable ports are jointly stable."""
from pathlib import Path
import json
from itertools import combinations
import inflation_bank as ib

HERE = Path(__file__).resolve().parent


def defect(rows, selected):
    h, u, w, _, _, _ = ib.tc.path_tables(rows)
    return h - sum(w[a][b] for a, b in combinations(selected, 2)), sum(u[a] for a in selected)


def expand(rows, selected, sizes):
    labels = [v for v in range(len(rows)) for _ in range(sizes.get(v, 1))]
    out = [0] * len(labels)
    for i, j in combinations(range(len(labels)), 2):
        forward = labels[i] == labels[j] or bool(rows[labels[i]] >> labels[j] & 1)
        out[i if forward else j] |= 1 << (j if forward else i)
    return out, [i for i, v in enumerate(labels) if v in selected]


def main():
    records = json.loads((HERE / 'inflation_bank_7.json').read_text())['records']
    groups = {}
    for item in records:
        key = (item['n'], item['code'], item['selected'])
        groups.setdefault(key, []).append(item)
    count, failures = 0, []
    for key, items in groups.items():
        for left, right in combinations(items, 2):
            rows = left['rows']
            selected = [v for v in range(len(rows)) if key[2] >> v & 1]
            a, b = left['inflated_vertex'], right['inflated_vertex']
            expanded, selected_expanded = expand(rows, selected, {a: 2, b: 2})
            d, u = defect(expanded, selected_expanded)
            count += 1
            if d != 1:
                failures.append(dict(n=key[0], code=key[1], rows=rows, selected=key[2],
                                     marked=[a,b], doubled_defect=d, doubled_slope=u))
    output = {'separately_stable_port_pairs_checked': count, 'joint_defect_failures_at_2_2': len(failures),
              'failures': failures}
    (HERE / 'composition_probe.json').write_text(json.dumps(output, indent=2)+'\n', encoding='utf8')
    print(json.dumps({k:v for k,v in output.items() if k!='failures'}, indent=2))
    print(json.dumps(failures[:3], indent=2))


if __name__ == '__main__':
    main()
