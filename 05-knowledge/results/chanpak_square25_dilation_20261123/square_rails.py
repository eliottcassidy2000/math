"""Six-vertex path norm, boundary cases and the first physical arc transfers."""
from selector import readings, skeleton

SOURCE = (2,16,11,35,13,23)


def inflate(a,b):
    types = [i for i in range(6) for _ in range(a if i==0 else b if i==5 else 1)]
    rows = [0]*len(types)
    for j,v in enumerate(types):
        for i,u in enumerate(types[:j]):
            if u == v or SOURCE[u] >> v & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
    return rows


def unit_edges(rows, h):
    changed = []
    for v in range(len(rows)):
        for u in range(v):
            target = rows.copy()
            target[u] ^= 1<<v
            target[v] ^= 1<<u
            delta = readings(target)[0] - h
            if abs(delta) == 2:
                changed.append((u,v,delta))
    return changed


if __name__ == "__main__":
    table = skeleton(SOURCE, (0,5))
    assert table == {(1,1):23,(2,1):9,(1,2):11,(2,2):4}
    assert readings(inflate(0,0))[0]==5
    for a in range(1,4):
        for b in range(1,4):
            h = readings(inflate(a,b))[0]
            assert 4*h+7 == (4*2**a+3)*(4*2**b+1)
    for k,expected in enumerate(([(0,4,2)],[(0,2,2)],[]),1):
        root = 2**(k+1)+1
        edges = unit_edges(inflate(k,k), root**2-2)
        assert edges==expected
        print("root",root,"lower height",root**2-2,"all unit edges of this lower host",edges)
