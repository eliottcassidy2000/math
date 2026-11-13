import collections
for n in range(3,10):
    info={}
    for line in open('cls_%d.txt'%n):
        p=line.split('|')[0].split(); info[int(p[0])]=(int(p[3]),int(p[4]),int(p[6]))  # strong, sc, ccode
    prof={}
    for line in open('prof_%d.txt'%n):
        p=line.split(); prof[int(p[1])]=tuple(map(int,p[7:]))
    ms={}
    for line in open('arcms_%d.txt'%n):
        left,right=line.split('|'); a=left.split(); ms[int(a[0])]=(int(a[2]),tuple(map(int,right.split())))
    # converse has the same multiset (paths reverse)
    merged={}
    for c,(st,sc,cc) in info.items(): merged.setdefault(min(c,cc),[]).append(c)
    for r,cs in merged.items(): assert ms[cs[0]]==ms[cs[-1]]
    reps=list(merged)
    M=collections.defaultdict(list); J=collections.defaultdict(list); Qc=collections.defaultdict(list)
    for r in reps: M[ms[r]].append(r); J[(ms[r],prof[r])].append(r); Qc[prof[r]].append(r)
    # Q-classes with >=2 merged classes, separated by multiset?
    multiQ=[g for g in Qc.values() if len(g)>=2]
    sep=sum(1 for g in multiQ if len(set(ms[r] for r in g))==len(g))
    part=sum(1 for g in multiQ if 1<len(set(ms[r] for r in g))<len(g))
    none=sum(1 for g in multiQ if len(set(ms[r] for r in g))==1)
    # multiset determines Q?
    mdq=sum(1 for g in M.values() if len(set(prof[r] for r in g))==1)
    darkstrong=sum(1 for r in reps if info[r][0]==1 and ms[r][0]>0)
    strongreps=sum(1 for r in reps if info[r][0]==1)
    print('n=%d merged %d | distinct multisets %d, (multiset,Q) %d, Q-classes %d | multisets determining Q %d/%d | Q-classes with >=2 merged: %d (fully separated by multiset %d, partly %d, not at all %d) | strong merged %d with dark arcs %d'%(
        n,len(reps),len(M),len(J),len(Qc),mdq,len(M),len(multiQ),sep,part,none,strongreps,darkstrong),flush=True)
# Session 310 pair and S26 pairs
ms8={}
for line in open('arcms_8.txt'):
    left,right=line.split('|'); a=left.split(); ms8[int(a[0])]=(int(a[2]),tuple(map(int,right.split())))
for c in (33559553,100664336,142622722,167804930,168296452,193724545,168043267):
    print('order8',c,ms8[c][0],'zeros; top',ms8[c][1][:6],'... sum',sum(ms8[c][1]))
