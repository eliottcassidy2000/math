"""Integral joint-character certificates for the three-star partial chart."""
from math import gcd
from random import Random
from three_stars import A, B, IDENTITY, PIVOT, SLOTS, make_chart, routes

PACKETS = (
    (24,51,20,25,12,29,44,6,33,41,53),
    (4,32,5,3,9,47,48,1,46,40,34),
    (49,31,14,17,13,22,50,16,38,15,21),
    (10,18,26,23,52,42,19,30,27,7,8),
    (35,36,2,54,11,45,55,39,37,43),
)


def relations(chart, packet):
    rows = []
    for h in SLOTS:
        if h not in packet:
            continue
        for x in range(1,56):
            row = [0]*165
            for j, vertex in enumerate((PIVOT,A,B)):
                f = chart[h,vertex]
                if f[x]:
                    row[55*j+f[x]-1] += 1
                if f[0]:
                    row[55*j+f[0]-1] -= 1
            rows.append(row)
    return rows


def square_minor(rows, seed):
    order = list(range(len(rows)))
    Random(seed).shuffle(order)
    pivots, chosen = {}, []
    for i in order:
        word = sum((v&1)<<j for j,v in enumerate(rows[i]))
        while word:
            j = word.bit_length()-1
            if j not in pivots:
                pivots[j] = word
                chosen.append(rows[i])
                break
            word ^= pivots[j]
        if len(chosen) == len(rows[0]):
            return chosen
    raise ValueError("no full-rank binary minor")


def determinant(matrix):
    # Small optional acceleration; the fallback is integer Bareiss elimination.
    try:
        from flint import fmpz_mat
    except ImportError:
        a = [r[:] for r in matrix]
        last, sign = 1, 1
        for i in range(len(a)-1):
            v = next(k for k in range(i,len(a)) if a[k][i])
            if i != v:
                a[i],a[v] = a[v],a[i]
                sign = -sign
            p = a[i][i]
            for k in range(i+1,len(a)):
                b = a[k][i]
                for j in range(i+1,len(a)):
                    a[k][j] = (p*a[k][j]-b*a[i][j])//last
            last = p
        return sign*a[-1][-1]
    return int(fmpz_mat(matrix).det())


def certificate(chart, packet):
    rows = relations(chart,packet)
    common = 0
    determinants = []
    for seed in range(3):
        d = abs(determinant(square_minor(rows,seed)))
        determinants.append(d)
        common = gcd(common,d)
        if common == 1:
            return determinants
    raise ValueError("chosen minors fail to certify saturation")


def old_control():
    chart = make_chart(())
    charts = (PIVOT,A,B)
    # Before trading, x_0-x_A and x_0-x_B are constant on each exterior branch.
    weights = ((1,-1,0),(1,0,-1))
    for prime in (2,7):
        constants = [[sum(t*chart[h,j][0] for t,j in zip(w,charts)) for h in SLOTS] for w in weights]
        assert all((sum(t*chart[h,j][x] for t,j in zip(w,charts))-c)%prime==0
                   for w,values in zip(weights,constants) for h,c in zip(SLOTS,values) for x in IDENTITY)
        failed = set()
        for b in SLOTS:
            omitted = [routes(chart,b,j) for j in charts]
            for x in IDENTITY:
                for w,values in zip(weights,constants):
                    available=sum(t*(sum(IDENTITY)-sum(r[x] for r in rr)) for t,rr in zip(w,omitted))
                    required=sum(values)-values[SLOTS.index(b)]
                    if (required-available)%prime:
                        failed.add((b,x))
        print("cyclic control",prime,"failing slots",len(failed))
        assert len(failed)=={2:2352,7:3024}[prime]


if __name__ == "__main__":
    assert sorted(x for packet in PACKETS for x in packet)==SLOTS
    chart = make_chart()
    for packet in (SLOTS, *PACKETS):
        print(len(packet),"branches, absolute maximal minors",certificate(chart,packet))
    old_control()
