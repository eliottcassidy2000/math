from pathlib import Path
from itertools import combinations
import sys
root=Path(sys.argv[1] if len(sys.argv)>1 else '../moore57_singer_repunit_six_partitions_20261126')
here=Path(__file__).parent
abc='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_'
pairs=list(combinations(range(6),2))
def verify(rows,boards):
    used=set()
    for row in rows:
        assert len(row)==6
        a=[abc.index(c) for c in row]
        assert max(a)<54
        for z,(i,j) in enumerate(pairs):
            p=(z,a[i],a[j]);assert p not in used and boards[z*54+a[i]]>>a[j]&1
            used.add(p)
for i in range(4):
    boards=[int(v,16) for v in (root/f'b{i}.txt').read_text().split()]
    text=(root/f'witness{i}.b64').read_text().strip()
    rows=[text[k:k+6] for k in range(0,len(text),6)]
    verify(rows,boards)
    for line in (here/'moves.txt').read_text().splitlines():
        t,*tokens=line.split()
        if int(t)!=i:continue
        j=tokens.index('>'); rm,add=tokens[:j],tokens[j+1:]
        assert len(add)==len(rm)+1
        for x in rm:rows.remove(x)
        rows+=add
        verify(rows,boards)
    (here/f'after{i}.b64').write_text(''.join(rows)+'\n')
    print(i,len(rows),len(rows)+(163 if i%2==0 else 162))
