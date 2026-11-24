"""Minimum physical dyadic square seams and the marked induced-containment obstruction."""
from itertools import combinations
from pathlib import Path
import subprocess
import sys
import tempfile

from bank import canonical, converse, height

# (root, outgoing vertex masks for the lower completion, flipped endpoints).
SEAMS = (
    (5, (30,60,8,48,36,5), (1,5)),
    (9, (14,124,120,48,33,65,25), (4,6)),
    (17, (110,252,24,240,33,196,148,21), (3,7)),
    (33, (414,508,120,304,480,449,137,268,68), (7,8)),
    (65, (382,892,1016,752,672,960,912,259,536,129), (5,9)),
    (129, (894,440,954,1904,1984,976,1926,1801,512,1026,295), (10,1)),
    (257, (1142,748,928,69,654,472,788,1865,1563,1065,126), (10,5)),
    (513, (1084,1149,992,1492,3300,3784,897,3843,1075,1307,68,1871), (11,1)),
)


def flip(rows, pair):
    i, j = pair
    new = list(rows)
    new[i] ^= 1 << j
    new[j] ^= 1 << i
    return tuple(new)


def check_tournament(rows):
    n = len(rows)
    assert all(0 <= row < 1 << n and not (row >> i & 1) for i, row in enumerate(rows))
    for i,j in combinations(range(n),2):
        assert ((rows[i] >> j) ^ (rows[j] >> i)) & 1


def strong(rows):
    reach = 1
    while True:
        new = reach
        for i,row in enumerate(rows):
            if (reach >> i) & 1:
                new |= row
        if new == reach:
            return reach == (1 << len(rows))-1
        reach = new


def check_seams():
    for root, rows, seam in SEAMS:
        upper = flip(rows, seam)
        check_tournament(rows)
        check_tournament(upper)
        for end in (rows, upper):
            assert strong(end) and strong(converse(end))
        assert (height(rows),height(upper)) == (root**2-2,root**2)
        print('root', root, 'minimum vertices', len(rows), 'counts', root**2-2, root**2)



def check_skip():
    rows = (322,488,451,341,487,461,256,73,128)
    upper = flip(rows,(4,5))
    for reading in (rows,upper):
        assert strong(reading) and strong(converse(reading))
    assert (height(rows),height(upper)) == (287,289)
    assert (height(tuple(x & 63 for x in rows[:6])),
            height(tuple(x & 63 for x in upper[:6]))) == (25,23)
    print('same arc (4,5) on first six / all nine', (25,23), (287,289))

def census():
    here = Path(__file__).resolve().parent
    inherited = here.parent / 'chanpak_arc_fiber_grammar_20261122' / 'birth_cut_cover.cpp'
    with tempfile.TemporaryDirectory() as tmp:
        obj, exe = Path(tmp)/'counts.o', Path(tmp)/'scan'
        subprocess.run(['c++','-O3','-std=c++17','-Dmain=birth_cut_cover_main','-c',
                        str(inherited),'-o',str(obj)],check=True)
        subprocess.run(['c++','-O3','-std=c++17',str(here/'bank_scan.cpp'),str(obj),
                        '-o',str(exe)],check=True)
        for r, next_cap, orders in ((5,9,(6,)),(9,11,(7,8,9))):
            records_for_lift = []
            for n in orders:
                with subprocess.Popen(['gentourng','-q',str(n-1)],stdout=subprocess.PIPE) as gen:
                    records = subprocess.check_output([str(exe),str(n-1),str(r*r-2)],
                                                      stdin=gen.stdout,text=True).splitlines()
                    assert gen.wait() == 0
                classes = set()
                for record in records:
                    *rows,i,j = (int(x) for x in record.split())
                    rows = tuple(rows)
                    if strong(rows) and strong(converse(rows)):
                        classes.add(min(canonical(rows,(i,j)),canonical(converse(rows),(i,j))))
                for mark, rows in sorted(classes):
                    pair = tuple(i for i,ch in enumerate(mark) if ch)
                    assert height(rows) == r*r-2 and height(flip(rows,pair)) == r*r
                    records_for_lift.append(' '.join(map(str, (*rows,*pair))))
                print('source root', r, 'order', n, 'raw', len(records),
                      'strong modulo converse', len(classes))
            assert len(records_for_lift) == {5:2,9:4}[r]
            lift = subprocess.check_output([str(exe),'--lift',str(next_cap),str((2*r-1)**2-2)],
                                           input='\n'.join(records_for_lift)+'\n',text=True)
            expected = {5: [(2,50,218),(2,49,208)], 9: [(2,88,437,2080),(2,86,430,2044)]}
            for line, visited in zip(lift.splitlines(),expected[r],strict=True):
                counts = tuple(map(int,line.split()))
                assert counts[:next_cap] == (0,)*(next_cap-len(visited)) + visited
                assert counts[next_cap:] == (0,)*(next_cap+1)
                print('expanded', visited, 'target hits', sum(counts[next_cap:]))
            if r == 5 and '--skip-census' in sys.argv:
                skip = subprocess.check_output([str(exe),'--lift','11','287'],
                                               input='\n'.join(records_for_lift)+'\n',text=True)
                expected_hits = [(0,)*12, (0,)*9+(6,48,360)]
                for line, expected_row in zip(skip.splitlines(),expected_hits,strict=True):
                    counts = tuple(map(int,line.split()))
                    assert counts[11:] == expected_row
                    print('two-step extension hits by order', counts[11:])


if __name__ == '__main__':
    check_seams()
    check_skip()
    if '--census' in sys.argv or '--skip-census' in sys.argv:
        census()
