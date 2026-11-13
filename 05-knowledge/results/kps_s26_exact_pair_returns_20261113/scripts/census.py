# full census of profile groups: mixed groups split by kind of decomposable member
import sys, collections
def load(n):
    rows=[]
    for line in open('prof_%d.txt'%n):
        p=list(map(int,line.split()))
        rows.append(dict(n=p[0],code=p[1],strong=p[2],nsink=p[3],nsrc=p[4],H=p[5],C=p[6],F=tuple(p[7:])))
    return rows
summary=[]
for n in map(int,sys.argv[1:]):
    rows=load(n); grp=collections.defaultdict(list)
    for r in rows: grp[r['F']].append(r)
    mixed=[g for g in grp.values() if any(r['strong'] for r in g) and any(not r['strong'] for r in g)]
    sing=[g for g in mixed if any((not r['strong']) and (r['nsink']>0 or r['nsrc']>0) for r in g)]
    block=[g for g in mixed if all(r['strong'] or (r['nsink']==0 and r['nsrc']==0) for r in g)]
    twins=[g for g in grp.values() if sum(r['strong'] for r in g)>=2]
    strongclasses=sum(r['strong'] for r in rows)
    strong_in_mixed=sum(sum(r['strong'] for r in g) for g in mixed)
    print('order %d: classes %d, strong %d, profiles %d | mixed groups %d (with singleton-end member %d, block-only %d), strong classes in mixed groups %d | strong-twin groups %d'%(
        n,len(rows),strongclasses,len(grp),len(mixed),len(sing),len(block),strong_in_mixed,len(twins)))
    for g in sorted(mixed,key=lambda g:g[0]['H']):
        st=[r for r in g if r['strong']]; ds=[r for r in g if not r['strong']]
        print('   H=%-5d strong %d (C=%s) | decomposable %d (sink-end %d, source-end %d, no singleton end %d)'%(
            g[0]['H'],len(st),sorted(set(r['C'] for r in st)),len(ds),sum(1 for r in ds if r['nsink']>0),sum(1 for r in ds if r['nsrc']>0),sum(1 for r in ds if r['nsink']==0 and r['nsrc']==0)))
