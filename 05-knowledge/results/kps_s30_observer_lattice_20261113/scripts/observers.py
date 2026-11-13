# separating power of tournament observers on merged classes (tournaments up to converse), orders 3-9
import collections, itertools
def decode(code,n):
    out=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1: out[i]|=1<<j
            else: out[j]|=1<<i
            k+=1
    return out
for n in range(3,10):
    info={}
    for line in open('cls_%d.txt'%n):
        p=line.split('|')[0].split(); info[int(p[0])]=dict(H=int(p[1]),c3=int(p[2]),strong=int(p[3]),sc=int(p[4]),aut=int(p[5]),cc=int(p[6]))
    for line in open('prof_%d.txt'%n):
        p=line.split(); c=int(p[1]); info[c]['C']=int(p[6]); info[c]['Q']=tuple(map(int,p[7:]))
    for line in open('arcms_%d.txt'%n):
        left,right=line.split('|'); a=left.split(); info[int(a[0])]['M']=(int(a[2]),tuple(map(int,right.split())))
    reps=sorted({min(c,v['cc']) for c,v in info.items()})
    for r in reps:
        out=decode(r,n); sc=tuple(sorted(bin(x).count('1') for x in out))
        # converse-invariant score: multiset of {s, n-1-s} pairs
        info[r]['S']=tuple(sorted(min(s,n-1-s) for s in sc)), tuple(sorted(sc))
        info[r]['Ssym']=tuple(sorted(tuple(sorted((s,n-1-s))) for s in sc))
    obs={'H':lambda v:v['H'],'scores(sym)':lambda v:v['Ssym'],'Q':lambda v:v['Q'],'arc multiset':lambda v:v['M'],
         'H,C,aut,c3':lambda v:(v['H'],v['C'],v['aut'],v['c3']),
         'scores+Q':lambda v:(v['Ssym'],v['Q']),'scores+arcs':lambda v:(v['Ssym'],v['M']),'Q+arcs':lambda v:(v['Q'],v['M']),
         'scores+arcs+aut':lambda v:(v['Ssym'],v['M'],v['aut'])}
    row=['n=%d merged %d'%(n,len(reps))]
    for name,f in obs.items():
        vals=collections.Counter(f(info[r]) for r in reps)
        row.append('%s %d'%(name,len(vals)))
    print(' | '.join(row),flush=True)
