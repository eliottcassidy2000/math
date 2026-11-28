"""Check literal shared edges, and feed forced choices back into the graph."""
from collections import Counter
from pathlib import Path
import argparse
import subprocess
from string import ascii_letters, digits

def read(s):
    lines = [list(map(int, x.split())) for x in s.splitlines()]
    n, = lines.pop(0)
    assert len(lines) == n
    A = []
    for i, (k, *row) in enumerate(lines):
        assert len(set(row)) == len(row) == k and i not in row
        A.append(set(row))
    assert all(u in A[v] for u, row in enumerate(A) for v in row)
    return A

def dump(A):
    return str(len(A))+'\n'+''.join(' '.join(map(str, [len(row), *sorted(row)]))+'\n' for row in A)

def distances(A):
    """Existence of exactly one positive-length route (length one or two) to each member."""
    R = []
    for u, row in enumerate(A):
        seen = set(row)
        for v in row:
            add = A[v] - {u}
            assert not (seen & add), (u, v, seen & add)
            seen |= add
        R.append(seen)
    return R

def audit(O, A):
    n = len(A)
    d = len(O[0])
    assert n == d*d+1 and all(O[u]<=A[u] for u in range(n))
    assert max(map(len, A)) == d
    R0, R = distances(O), distances(A)
    S = {u for u in range(n) if len(O[u])==d}
    assert S == {u for u in range(n) if len(A[u])==d}
    holes = [sum(n-1-len(Q[u]) for u in S) for Q in (R0, R)]
    assert all(S - {u} <= R0[u] <= R[u] for u in S)
    print("Uncovered pairs from saturated endpoints", *holes)
    U = {u for u in range(n) if len(O[u])<d}
    fixed_before = fixed_after = 0
    board_new = completed_new = completed_old = board_multi = 0
    participation = Counter()
    for u in sorted(S):
        for v in sorted(S - O[u]):
            if u>=v:
                continue
            mid, = O[u] & O[v]
            L, T = O[u]-{mid}, O[v]-{mid}
            pairs = [(x,y) for x in L for y in A[x]&T]
            assert len({x for x,y in pairs}) == len({y for x,y in pairs}) == len(pairs)
            old = sum(y in O[x] for x,y in pairs)
            gain = len(pairs)-old
            board_new += gain>0
            board_multi += gain>1
            completed_old += old==d-1
            completed_new += old<d-1 and len(pairs)==d-1
            fixed_before += old
            fixed_after += len(pairs)
            for x,y in pairs:
                if y not in O[x]:
                    participation[min(x,y),max(x,y)] += 1
    new = {(u,v) for u in U for v in A[u]-O[u] if u<v}
    assert set(participation) == new
    assert all(participation[u,v]==len(O[u]&S)*len(O[v]&S) for u,v in new)
    print('New edges',len(new),'total edges',sum(map(len,A))//2,'maximum degree',max(map(len,A)))
    print('Boards improved',board_new,'gained multiple hits',board_multi)
    print('Boards completed',completed_old,'->',completed_old+completed_new)
    print('Fixed incidences',fixed_before,'->',fixed_after)
    print('Physical-edge multiplicity in boards',sorted(Counter(participation.values()).items()))

def install_witness(O, data):
    A = [row.copy() for row in O]
    d = len(O[0])
    size = d - 1
    letters = ascii_letters + digits
    seen = set()
    for line in data.splitlines():
        b, c, word = line.split()
        b, c = int(b), int(c)
        assert 0 <= b < c < d and (b,c) not in seen
        seen.add((b,c))
        image = [letters.index(t) for t in word]
        assert sorted(image) == list(range(size))
        lo, ro = d+1+b*size, d+1+c*size
        for x,y in enumerate(image):
            assert all(z == ro+y for z in A[lo+x] if ro <= z < ro+size)
            assert all(z == lo+x for z in A[ro+y] if lo <= z < lo+size)
            A[lo+x].add(ro+y)
            A[ro+y].add(lo+x)
    distances(A)
    assert max(map(len, A)) == d
    print('Shared branch-pair matchings in certificate',len(seen),flush=True)
    return A

def main():
    p=argparse.ArgumentParser()
    p.add_argument('old',type=Path)
    p.add_argument('new',type=Path)
    p.add_argument('checker',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--matchings',action='store_true',help='read branch-pair witness instead of adjacency')
    a=p.parse_args()
    O = read(a.old.read_text())
    A = install_witness(O,a.new.read_text()) if a.matchings else read(a.new.read_text())
    for round in range(100):
        rr=subprocess.run([str(a.checker.resolve())],input=dump(A),text=True,capture_output=True)
        print(rr.stdout,end='',flush=True)
        assert rr.returncode==0, f'Wedge Hall failure at {round=}'
        force = [tuple(map(int,l.split()[1:])) for l in rr.stdout.splitlines() if l.startswith('force ')]
        if not force:
            break
        print('Install',len(force),'shared edges',flush=True)
        for u,v in force:
            assert v not in A[u]
            A[u].add(v); A[v].add(u)
        distances(A)
        assert max(map(len,A))==len(O[0])
    else:
        raise AssertionError('iteration cap')
    audit(O,A)
    a.output.write_text(dump(A))

if __name__=='__main__':
    main()
