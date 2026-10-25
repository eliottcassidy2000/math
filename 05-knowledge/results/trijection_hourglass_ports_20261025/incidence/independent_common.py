"""Independent common-kernel reconstruction via physical parallelism."""
from collections import Counter
from itertools import product
import json


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def reduce_rows(rows, width, p):
    a = [[v % p for v in row] for row in rows]
    pivots = []
    for column in range(width):
        k = next((i for i in range(len(pivots), len(a)) if a[i][column]), None)
        if k is None:
            continue
        r = len(pivots)
        a[k], a[r] = a[r], a[k]
        inv = pow(a[r][column], p-2, p)
        a[r] = [(v*inv) % p for v in a[r]]
        for i in range(len(a)):
            if i != r and a[i][column]:
                factor = a[i][column]
                a[i] = [(x-factor*y) % p for x,y in zip(a[i],a[r])]
        pivots.append(column)
    return a[:len(pivots)], pivots


def nullspace(rows, width, p):
    a, pivots = reduce_rows(rows, width, p)
    basis = []
    for free in range(width):
        if free not in pivots:
            v = [0]*width
            v[free] = 1
            for row, pivot in zip(a,pivots):
                v[pivot] = -row[free] % p
            basis.append(v)
    return basis


def check(p, length):
    addresses = list(product(range(2), range(2), range(length)))
    points = [((-a-3*b) % p, (j-9*b) % p) for a,b,j in addresses]
    need(len(set(points)) == len(points), "literal points distinct")
    width = len(points)
    constraints = []
    line_sizes = []
    for direction in ((1,0),(0,1),(3,9)):
        remaining = set(range(width))
        while remaining:
            i = min(remaining)
            x0,y0 = points[i]
            line = {j for j,(x,y) in enumerate(points)
                    if ((x-x0)*direction[1]-(y-y0)*direction[0]) % p == 0}
            need(i in line and line <= remaining, "line partition")
            remaining -= line
            constraints.append([int(j in line) for j in range(width)])
            line_sizes.append(len(line))
    exact = nullspace(constraints, width, p)
    q = max(0,length-9)
    model = []

    def array(h, r):
        get = lambda seq,j: seq[j] if 0 <= j < len(seq) else 0
        values = []
        for a,b,j in addresses:
            P = -sum(get(r,j-s) for s in (3,6,9))
            R = get(r,j)
            S = -get(r,j-9)
            values.append(((h[j],P-h[j],R-h[j],S-P+h[j])[2*a+b]) % p)
        return values

    for j in range(length-1):
        h = [0]*length
        h[j],h[j+1] = 1,-1
        model.append(array(h,[0]*q))
    cubes = model[:]
    for j in range(max(0,q-1)):
        r = [0]*q
        r[j],r[j+1] = 1,-1
        model.append(array([0]*length,r))
    expected = length-1+max(0,length-10)
    need(len(exact) == len(model) == expected, "physical common dimension")
    need(len(reduce_rows(model,width,p)[1]) == expected, "model independent")
    need(all(sum(a*b for a,b in zip(row,v)) % p == 0
             for row in constraints for v in model), "literal line sums")
    need(len(reduce_rows(model+exact,width,p)[1]) == expected, "complete span")
    bridge_support = None
    if length >= 11:
        bridge = model[len(cubes)]
        bridge_support = sum(v != 0 for v in bridge)
        need(bridge_support == 12, "twelve-point bridge")
        need(len(reduce_rows(cubes+[bridge],width,p)[1]) == len(cubes)+1,
             "bridge outside cube span")
    return dict(p=p,L=length,J=len(exact),coupled=2*len(exact),
                cube_modes=len(cubes),bridge_modes=max(0,length-10),
                first_bridge_support=bridge_support)


records = [check(p,L) for p in (23,29,31) for L in range(2,p-8)]
need(len(records) == 53, "full pre-wrap universe")
print(json.dumps({"status":"PASS","universe":records,
                  "cases":len(records),"first_bridge_at_L":11,
                  "method":"parallelism-built lines plus full nullspace; no producer imports"},
                 indent=2,sort_keys=True))
