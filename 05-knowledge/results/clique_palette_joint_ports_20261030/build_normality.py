"""Construct small exact decompositions for a frozen triangulation."""
from pathlib import Path
from collections import Counter
from itertools import product
import json
from explore import GOOD, palette, N

HERE = Path(__file__).resolve().parent

def det(a):
    a = [row[:] for row in a]
    sign, old, n = 1, 1, len(a)
    for k in range(n-1):
        pivot = next((j for j in range(k, n) if a[j][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign *= -1
        current = a[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                top = a[i][j]*current-a[i][k]*a[k][j]
                if top % old:
                    raise RuntimeError('inexact Bareiss division')
                a[i][j] = top//old
            a[i][k] = 0
        old = current
    return sign*a[-1][-1]

def decomposition(d, q):
    d, colors = tuple(d), []
    while any(d):
        active = sum(1 << i for i, x in enumerate(d) if x)
        for s in GOOD:
            if s & active != s:
                continue
            child = tuple(x-((s >> i) & 1) for i, x in enumerate(d))
            if palette(child)+1 == palette(d):
                colors.append(s)
                d = child
                break
        else:
            raise RuntimeError('no palette step')
    if len(colors) > q:
        raise RuntimeError('hole in the proposed normal semigroup')
    return colors+[0]*(q-len(colors))

if __name__ == '__main__':
    cert = json.loads((HERE/'triangulation.json').read_text())
    verts = cert['vertex_masks']
    generators = [[1]+[(s >> i) & 1 for i in range(7)] for s in verts]
    data, degrees = [], Counter()
    for index, simplex in enumerate(cert['simplices']):
        a = [[generators[j][i] for j in simplex] for i in range(8)]
        determinant = abs(det(a))
        if determinant == 1:
            continue
        if determinant not in (2, 3):
            raise RuntimeError('unsupported determinant in the small certificate')
        for residue in product(range(determinant), repeat=8):
            if not any(residue):
                continue
            numerator = [sum(x*y for x, y in zip(row, residue)) for row in a]
            if any(x % determinant for x in numerator):
                continue
            point = [x//determinant for x in numerator]
            colors = decomposition(point[1:], point[0])
            data.append({'simplex': index, 'residue': residue, 'point': point, 'colors': colors})
            degrees[point[0]] += 1
    cert['fundamental_points'] = data
    (HERE/'normality_certificate.json').write_text(json.dumps(cert, separators=(',', ':'))+'\n', encoding='utf-8')
    print(json.dumps({'nonzero_fundamental_points': len(data), 'heights': dict(degrees)}))
