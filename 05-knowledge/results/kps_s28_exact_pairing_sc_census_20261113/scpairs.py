# exact pairings (equal full slot polynomial) between merged classes (tournaments up to converse), with self-converse status
import collections, sys
for n in range(3,10):
    info={}
    for line in open('cls_%d.txt'%n):
        p=line.split('|')[0].split()
        code=int(p[0]); H=int(p[1]); strong=int(p[3]); sc=int(p[4]); aut=int(p[5]); cc=int(p[6])
        info[code]=(H,strong,sc,aut,cc)
    prof={}
    for line in open('prof_%d.txt'%n):
        p=line.split(); prof[int(p[1])]=tuple(map(int,p[7:]))
    assert set(prof)==set(info)
    # merged classes: representative = min(code, ccode)
    merged={}
    for c,(H,st,sc,aut,cc) in info.items():
        r=min(c,cc); merged.setdefault(r,[]).append(c)
    for r,cs in merged.items():
        assert prof[cs[0]]==prof[cs[-1]]          # Q(T) = Q(T^op)
    groups=collections.defaultdict(list)
    for r,cs in merged.items(): groups[prof[r]].append(r)
    nSC=sum(1 for v in info.values() if v[2]); nM=len(merged)
    multi=[g for g in groups.values() if len(g)>=2]
    def issc(r): return info[r][2]==1
    def isst(r): return info[r][1]==1
    scsc=sum(1 for g in multi if sum(issc(r) for r in g)>=2)
    scn=sum(1 for g in multi if any(issc(r) for r in g) and any(not issc(r) for r in g))
    stmulti=[g for g in groups.values() if sum(isst(r) for r in g)>=2]
    stscmulti=[g for g in stmulti if sum(issc(r) and isst(r) for r in g)>=1]
    # SC classes whose Q is shared with nothing else (Q determines the SC class)
    scalone=sum(1 for g in groups.values() if len(g)==1 and issc(g[0]))
    print('n=%d classes %d SC %d merged %d | Q-classes %d, Q-classes with >=2 merged classes %d (SC-SC %d, SC-NSC %d); strong twin Q-classes (>=2 strong merged) %d, with an SC member %d; SC classes alone in their Q-class %d'%(
        n,len(info),nSC,nM,len(groups),len(multi),scsc,scn,len(stmulti),len(stscmulti),scalone),flush=True)
