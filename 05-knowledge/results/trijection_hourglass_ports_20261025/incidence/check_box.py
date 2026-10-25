"""Exact line-incidence controls for three-direction boxes. No dependencies."""
import json


def rank(rows, p):
    pivots = {}
    for raw in rows:
        row = [v % p for v in raw]
        for i in range(len(row)):
            v = row[i]
            if not v:
                continue
            if i in pivots:
                row = [(a - v*b) % p for a, b in zip(row, pivots[i])]
            else:
                iv = pow(v, -1, p)
                pivots[i] = [(a*iv) % p for a in row]
                break
    return len(pivots)


def rank_columns(rows, p):
    """Independent rectangular elimination: column scans and row swaps."""
    a = [[x % p for x in row] for row in rows]
    if not a:
        return 0
    r = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][col], -1, p)
        a[r] = [(x*inv) % p for x in a[r]]
        for i in range(r+1, len(a)):
            q = a[i][col]
            if q:
                a[i] = [(x-q*y) % p for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def check(p, w, length):
    directions = [(1, 0), (0, 1), w]
    points = sorted({((-a-b*w[0]) % p, (j-b*w[1]) % p)
                     for a in range(2) for b in range(2) for j in range(length)})
    if len(points) != 4*length:
        raise ValueError('box collided')
    m = len(points)
    matrices = []
    for d in range(3):
        rows = []
        for e in range(3):
            if d == e:
                continue
            x, y = directions[e]
            labels = [(x*v-y*u) % p for u, v in points]
            for label in sorted(set(labels)):
                rows.append([int(k == label) for k in labels])
        matrices.append(rows)
    sigma = [-w[0], -w[1], 1]
    whole = []
    for d, rows in enumerate(matrices):
        for row in rows:
            whole.append([0]*(d*m) + row + [0]*((2-d)*m))
    for i in range(m):
        row = [0]*(3*m)
        for d in range(3):
            row[d*m+i] = sigma[d]
        whole.append(row)
    dims = [m-rank(rows, p) for rows in matrices]
    nullity = 3*m-rank(whole, p)
    if [m-rank_columns(rows, p) for rows in matrices] != dims:
        raise ValueError('independent dimensions disagree')
    if 3*m-rank_columns(whole, p) != nullity:
        raise ValueError('independent coupled ranks disagree')
    # The adjacent-cube relations are actual vectors, not only dimensions.
    index = {point: i for i, point in enumerate(points)}
    basis = []
    for j in range(length-1):
        cube = [0]*m
        for a in range(2):
            for b in range(2):
                for c in range(2):
                    point = ((-a-b*w[0]) % p, (j+c-b*w[1]) % p)
                    cube[index[point]] += (-1)**(a+b+c)
        for d in (0, 1):
            vec = [0]*(3*m)
            for i, x in enumerate(cube):
                vec[d*m+i] = x
                vec[2*m+i] = -sigma[d]*x
            if any(sum(x*y for x,y in zip(row,vec)) % p for row in whole):
                raise ValueError('adjacent cube does not satisfy constraints')
            basis.append(vec)
    if rank_columns(basis, p) != 2*(length-1):
        raise ValueError('adjacent cubes dependent')
    return {'p': p, 'w': w, 'length': length, 'points': m,
            'dimensions': dims, 'coupled_rank': sum(dims)-nullity,
            'coupled_nullity': nullity,
            'adjacent_cube_rank': 2*(length-1),
            'extra_kernel_dimensions': nullity-2*(length-1)}


if __name__ == '__main__':
    rows = [check(29, (3, 9), length) for length in range(2, 13)]
    generic = [check(101, (17, 37), length) for length in range(2, 11)]
    for row in generic:
        length = row['length']
        if row['dimensions'] != [2*(length-1), length, 2*(length-1)]:
            raise ValueError(('generic dimensions', row))
        if row['coupled_nullity'] != 2*(length-1):
            raise ValueError(('generic nullity', row))
    if rows[2]['dimensions'] != [6, 4, 6] or rows[2]['coupled_nullity'] != 6:
        raise ValueError('inherited sixteen-point witness changed')
    if [r['extra_kernel_dimensions'] for r in rows] != [0]*9+[2,4]:
        raise ValueError('first coupled resonance boundary changed')
    print(json.dumps({'actual29': rows, 'generic101': generic}, indent=2))
