"""Independent finite controls of bounded exceptions and TT clone monotonicity."""
from pathlib import Path
from itertools import combinations
import json
import retained_reset as audit


def graph(n, code):
    rows=[0]*n
    for bit,(i,j) in enumerate(combinations(range(n),2)):
        u,v=(i,j) if code>>bit&1 else (j,i)
        rows[u]|=1<<v
    return tuple(rows)


def defects(rows):
    hp,ends=audit.tables(rows)
    full=(1 << len(rows))-1
    return tuple(sum(ends[m][r]*hp[full ^ m] for m in range(1,full+1)) - hp[full]
                 for r in range(len(rows)))


def residual_sccs(rows, vertices):
    # Mutual reachability by independent bitset transitive closure.
    reach={v:(rows[v] & vertices)|(1 << v) for v in range(len(rows)) if vertices>>v&1}
    for mid in reach:
        for v in reach:
            if reach[v]>>mid&1:
                reach[v]|=reach[mid]
    left=vertices; out=[]
    while left:
        v=(left & -left).bit_length()-1
        cls=sum(1 << w for w in reach if left>>w&1 and reach[v]>>w&1 and reach[w]>>v&1)
        out.append(cls);left^=cls
    return out


def clone(rows,x):
    n=len(rows); out=list(rows)+[rows[x]|(1<<x)]
    # The new vertex is the smaller clone, and has identical external contacts.
    for v in range(n):
        if v != x and rows[v]>>x&1:
            out[v]|=1 << n
    return tuple(out)


def main():
    graphs=roots=clones=0
    for n in range(3,6):
        for code in range(1 << (n*(n-1)//2)):
            rows=graph(n,code)
            if not audit.strong(rows):continue
            graphs+=1
            bs=defects(rows)
            for r,b in enumerate(bs):
                ports=[v for v in range(n) if rows[v]>>r&1]
                d=len(ports)
                if d<2:continue
                residual=((1 << n)-1) ^ (1 << r) ^ sum(1 << x for x in ports)
                counts={}
                for cls in residual_sccs(rows,residual):
                    if cls & (cls-1):continue
                    v=cls.bit_length()-1
                    contact=sum(1 << j for j,a in enumerate(ports) if rows[v]>>a&1)
                    counts[contact]=counts.get(contact,0)+1
                if d==2:
                    c0,c1,c2=(counts.get(x,0) for x in (0,1,2))
                    lhs=c0+c1*(c1-1)//2+c2*(c2-1)//2
                else:
                    lhs=sum(m for c,m in counts.items() if c!=(1 << d)-1)
                audit.need(2*lhs<=b,'exceptional atom bound')
                roots+=1
                for x in range(n):
                    if not residual>>x&1:continue
                    bigger=clone(rows,x)
                    audit.need(audit.strong(bigger),'strongness under nonempty TT inflation')
                    audit.need(defects(bigger)[r]>=b,'positive defect clone injection')
                    clones+=1
    report=dict(status='PASS', universe='Every labelled strong tournament through five vertices, every branching root and every out-neighbour clone',
                graphs=graphs, branching_roots=roots, clones=clones,
                method='Independent endpoint subset recurrence, mutual reachability SCCs, literal new clone arcs')
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(report)


if __name__=='__main__':main()
