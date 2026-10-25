#!/usr/bin/env python3
import hashlib,itertools,json

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def edges(n,kind):
    removed={(i,i+1) for i in range(n-1)}
    if kind=='cycle': removed.add((0,n-1))
    return [(i,j) for i in range(n) for j in range(i+1,n) if (i,j) not in removed]

def score(n,ee,bits):
    r=[0]*n
    for q,(a,b) in enumerate(ee): r[b if (bits>>q)&1 else a]+=1
    return tuple(r)

def cuts(n,ee,r):
    if any(x<0 for x in r) or sum(r)!=len(ee):return False
    for mask in range(1<<n):
        if sum(r[i] for i in range(n) if mask>>i&1)<sum((mask>>a&1) and (mask>>b&1) for a,b in ee): return False
    return True

def realize(n,ee,r):
    if min(r)<0 or sum(r)!=len(ee):return None
    slots=[]
    for i,x in enumerate(r):slots.extend([i]*x)
    match=[-1]*len(slots)
    adj=[[s for s,v in enumerate(slots) if v in e] for e in ee]
    def augment(e,seen):
        for s in adj[e]:
            if s in seen:continue
            seen.add(s)
            if match[s]<0 or augment(match[s],seen):match[s]=e;return True
        return False
    for e in range(len(ee)):
        if not augment(e,set()): return None
    out=[None]*len(ee)
    for s,e in enumerate(match):out[e]=slots[s]
    need(tuple(out.count(i) for i in range(n))==r,'matching wrong score')
    need(all(t in ee[i] for i,t in enumerate(out)),'matching wrong tail')
    return out

report={'exhaustive':[],'prescribed_sets':[],'hostiles':{},'internal':[]}
for n in range(3,7):
 for kind in ('path','cycle'):
    ee=edges(n,kind); actual={score(n,ee,b) for b in range(1<<len(ee))}
    deg=[sum(i in e for e in ee) for i in range(n)]
    candidates=0; feasible=0
    for r in itertools.product(*(range(d+1) for d in deg)):
        if sum(r)!=len(ee):continue
        candidates+=1
        c=cuts(n,ee,r)
        need(c==(r in actual),'cut criterion mismatch')
        need((realize(n,ee,r) is not None)==c,'matching criterion mismatch')
        feasible+=c
    report['exhaustive'].append(dict(n=n,kind=kind,orientations=1<<len(ee),candidate_scores=candidates,realized_scores=feasible,digest=hashlib.sha256(repr(sorted(actual)).encode()).hexdigest()))
for n in range(4,15,2):
 for kind in ('path','cycle'):
    p=n//2;ee=edges(n,kind); size=p+(kind=='path'); count=0; failed=[]
    for H in itertools.combinations(range(n),size):
        r=tuple(p-2+(i in H) for i in range(n))
        result=realize(n,ee,r)
        expect=not(n==4 and kind=='cycle') or (len(set(H)&{0,2})==1 and len(set(H)&{1,3})==1)
        need((result is not None)==expect,'high-set theorem mismatch')
        if result is None:failed.append(H)
        else:count+=1
    report['prescribed_sets'].append(dict(n=n,kind=kind,realizable=count,nonrealizable=failed))
for n in (4,6):
    p=n//2;ee=edges(n,'path');patterns={}; internal=0;balanced=0
    for bits in range(1<<len(ee)):
        r=score(n,ee,bits);d=tuple(r[i]+(i<n-1) for i in range(n))
        dc=tuple(n-1-d[i]+(i==0)-(i==n-1) for i in range(n))
        bal=all(x in (p-1,p) for x in d)
        if not bal:continue
        balanced+=1
        isinternal=all(x in (p-1,p) for x in dc)
        need(isinternal==(d[0]==p and d[-1]==p-1),'endpoint criterion')
        if not isinternal:continue
        internal+=1
        H=tuple(i for i,x in enumerate(r) if x==p-1)
        need(all(x in(p-2,p-1) for x in r) and len(H)==p+1 and 0 in H and n-1 in H,'internal residual')
        patterns[H]=patterns.get(H,0)+1
        comp=score(n,ee,((1<<len(ee))-1)^bits)
        target=tuple(i for i,x in enumerate(comp) if x==p-1)
        need(target==tuple(i for i in range(n) if i in (0,n-1) or i not in H),'high-set transport')
        # Build parent adjacency directly, check degree and complemented lift.
        a=[[0]*(n+1) for _ in range(n+1)]
        for i in range(n-1):a[i][i+1]=1
        for q,(i,j) in enumerate(ee):a[j][i]=1 if bits>>q&1 else 0;a[i][j]=1-a[j][i]
        for i in range(n):a[n][i]=int(d[i]==p);a[i][n]=1-a[n][i]
        need(all(sum(row)==p for row in a),'regular extension')
        cy={(n,0),(n-1,n)}|{(i,i+1) for i in range(n-1)}
        need(all(a[i][j] for i,j in cy),'cycle missing')
        residual=[[a[i][j]-(int((i,j) in cy)) for j in range(n+1)] for i in range(n+1)]
        need(all(sum(row)==p-1 for row in residual),'Euler residual')
        ac=[[a[i][j] if (i,j) in cy or (j,i) in cy or i==j else a[j][i] for j in range(n+1)] for i in range(n+1)]
        need(all(sum(row)==p for row in ac),'complement regularity')
        need(all(ac[n][i]==int(dc[i]==p) for i in range(n)),'complement extension')
    expected=set(itertools.combinations(range(1,n-1),p-1))
    need({tuple(i for i in h if i not in(0,n-1)) for h in patterns}==expected,'missing internal score pattern')
    report['internal'].append(dict(n=n,balanced_words=balanced,internal_words=internal,pattern_counts=[{'H':h,'words':v} for h,v in sorted(patterns.items())]))
report['odd_prescribed_sets']=[]
for n in range(3,16,2):
    m=(n-3)//2
    ee=edges(n,'path');good=[]
    for w in range(n):
        r=tuple(m+(i==w) for i in range(n))
        ok=realize(n,ee,r) is not None
        need(ok==(n>=5 or w in(0,n-1)),'odd high-set theorem')
        if ok:good.append(w)
    need(realize(n,edges(n,'cycle'),(m,)*n) is not None,'odd Eulerian prescription')
    report['odd_prescribed_sets'].append(dict(n=n,realizable_high_vertices=good))
report['hostiles']={'cycle4_bad_H':[0,2],'cycle4_bad_r':[1,0,1,0],'path4_r':[1,1,1,0],'path4_high_count':3,'cycle4_high_count':2,'odd_path3_total':1,'odd_path3_constant_target_total':0}
need(not cuts(4,edges(4,'cycle'),(1,0,1,0)),'hostile not rejected')
need(realize(4,edges(4,'path'),(1,1,1,0)) is not None,'path positive missing')
print(json.dumps(report,sort_keys=True,indent=2))
