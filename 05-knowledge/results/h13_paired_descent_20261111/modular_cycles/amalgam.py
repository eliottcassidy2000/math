"""Literal controls for rooted ordinal amalgamation and additive switches."""
from pathlib import Path
import importlib.util
import json

HERE = Path(__file__).resolve().parent
sp = importlib.util.spec_from_file_location("paired_probe", HERE / "probe.py")
p = importlib.util.module_from_spec(sp)
sp.loader.exec_module(p)


def converse(rows):
    full = (1 << len(rows)) - 1
    return [full ^ row ^ (1 << i) for i, row in enumerate(rows)]


def extended(rows, root):
    h, e, u, d = p.state(rows, root)
    hc, s, v, dc = p.state(converse(rows), root)
    assert (h, d) == (hc, dc)
    keep = [i for i in range(len(rows)) if i != root]
    deleted = [sum(((rows[i] >> j) & 1) << b for b, j in enumerate(keep)) for i in keep]
    hs = p.c.hp_subsets(deleted)
    full = (1 << len(deleted)) - 1
    k = sum(hs[m] * hs[full ^ m] for m in range(full + 1))
    return dict(H=h, E=e, S=s, U=u, V=v, d=d, K=k)


def amalgam(left, lr, right, rr):
    """Identify roots; every nonroot left vertex beats every nonroot right."""
    aa = [i for i in range(len(left)) if i != lr]
    bb = [i for i in range(len(right)) if i != rr]
    root = len(aa) + len(bb)
    rows = [0] * (root + 1)
    lm = {v: i for i, v in enumerate(aa)} | {lr: root}
    rm = {v: len(aa) + i for i, v in enumerate(bb)} | {rr: root}
    for source, mapping in ((left, lm), (right, rm)):
        for a in mapping:
            for b in mapping:
                if source[a] >> b & 1:
                    rows[mapping[a]] |= 1 << mapping[b]
    right_mask = sum(1 << rm[v] for v in bb)
    for a in aa:
        rows[lm[a]] |= right_mask
    return rows, root


def predicted(a, b):
    # Assumes both factors have a nonroot vertex; singleton is the identity.
    d, e, u = a['d'], a['E'], a['U']
    D, S, V = b['d'], b['S'], b['V']
    return dict(
        H=D * (a['H'] - e) + d * (b['H'] - S) + e * S + (u - d) * (V - D),
        E=d * b['E'], S=D * a['S'],
        U=(u - d) * b['K'] + d * b['U'],
        V=(V - D) * a['K'] + D * a['V'],
        d=d * D, K=a['K'] * b['K'])


def q(st):
    return [st['H'], st['U'], st['U'] - st['d'], st['H'] + st['U'] - st['E']]


def qmatrix(b):
    h, e, s, u, v, d, k = (b[t] for t in ('H', 'E', 'S', 'U', 'V', 'd', 'K'))
    first = [s, h - d, s + v - h - d, d - s]
    return [first, [0, u, k - u, 0], [0, u - d, k - u + d, 0],
            [first[0], first[1] + u - e, first[2] + k - u + e, first[3]]]


def main():
    Q = [18, 8, 3, 5, 14]
    R = [6, 12, 24, 17, 3]
    h13 = [p.c.construct(Q, cut, 0, 0) for cut in (6, 3)]
    c5 = [p.c.construct(R, cut, 0, 0) for cut in (3, 5)]
    c25 = [p.c.construct(p.add(Q, mask), cut, 0, 1) for mask, cut in ((14, 38), (8, 35))]
    # These controls include transitive roots at each endpoint, a cyclic root,
    # two distinct H13 boundary responses and the different H25 deletion cores.
    sources = [[2, 0], [0, 1], [2, 4, 1], *h13, *c5, *c25]
    pairs = [(i, j) for i in range(5) for j in range(5)] + [(5, 7), (6, 8)]
    controls = []
    for i, j in pairs:
        a, b = sources[i], sources[j]
        aa, bb = extended(a, len(a)-1), extended(b, len(b)-1)
        joined, root = amalgam(a, len(a)-1, b, len(b)-1)
        exact, formula = extended(joined, root), predicted(aa, bb)
        if exact != formula:
            raise RuntimeError(('amalgam', i, j, exact, formula))
        matrix = qmatrix(bb)
        if any(x < 0 for row in matrix for x in row):
            raise RuntimeError(('negative matrix', j, matrix))
        if q(exact) != [sum(x*y for x, y in zip(row, q(aa))) for row in matrix]:
            raise RuntimeError(('q response', i, j))
        controls.append({'left': i, 'right': j, 'state': exact})
    rectangles = []
    for label, right in [('C5', c5), ('H25', c25), ('converse-H25', list(map(converse, c25)))]:
        states = []
        for a in h13:
            row = []
            for b in right:
                t, r = amalgam(a, len(a)-1, b, len(b)-1)
                row.append(extended(t, r))
            states.append(row)
        cross = {key: states[1][1][key] - states[1][0][key] - states[0][1][key] + states[0][0][key]
                 for key in states[0][0]}
        if any(cross.values()):
            raise RuntimeError(('mixed response', label, cross))
        if label != 'converse-H25':
            sidecars = {(s['E'], s['U'], s['d'], s['K']) for row in states for s in row}
            if len(sidecars) != 1:
                raise RuntimeError(('pure H rectangle', label, sidecars))
        rectangles.append({'right_pair': label, 'states': states, 'mixed_difference': cross})
    old = [p.c.construct(Q, cut, 0, 4) for cut in (6, 3)]
    signed_order = []
    for branch, (before, after) in enumerate(zip(old, c25)):
        for i, gadget in enumerate(sources[:3]):
            for side in ('left', 'right'):
                states = []
                for source in (before, after):
                    args = (gadget, len(gadget)-1, source, len(source)-1) if side == 'left' else (
                        source, len(source)-1, gadget, len(gadget)-1)
                    graph, root = amalgam(*args)
                    states.append(extended(graph, root))
                diff = {key: states[0][key] - states[1][key] for key in states[0]}
                if any(diff[key] < 0 for key in ('H', 'U', 'V', 'K')) or any(
                    diff[key] > 0 for key in ('E', 'S', 'd')):
                    raise RuntimeError(('signed two-sided order', branch, i, side, diff))
                signed_order.append(dict(branch=branch, gadget=i, side=side, states=states, difference=diff))
    report = {'status': 'FINITE-EXACT literal controls; universal results in ADDITIVITY.md',
              'universe': '27 named ordered source pairs; three full 2x2 switch rectangles; twelve two-sided reset controls',
              'sources': [{'rows': t, 'root': len(t)-1, 'state': extended(t, len(t)-1)} for t in sources],
              'controls': controls, 'rectangles': rectangles,
              'reset_full_states': [[extended(t, len(t)-1) for t in group] for group in (old, c25)],
              'two_sided_reset_controls': signed_order}
    (HERE / 'amalgam_receipt.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'controls': len(controls), 'rectangles': [
        {'pair': r['right_pair'], 'H': [[s['H'] for s in row] for row in r['states']]}
        for r in rectangles]}, indent=2))


if __name__ == '__main__':
    main()
