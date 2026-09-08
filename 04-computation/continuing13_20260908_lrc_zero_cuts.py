"""Two-clock zero-cut consumer: complete native banks and all valid profile words.

No producer imports. All native ratios, quotient clocks, positional words and
complete zero-component cuts are retained. Analytic proof is in the report.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,combinations_with_replacement
from collections import Counter
from math import gcd
from hashlib import sha256
import json,sys
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve()
ROOT=HERE.parent.parent if HERE.parent.name=='04-computation' else Path.cwd()
DEST=ROOT/'05-knowledge/results' if HERE.parent.name=='04-computation' else HERE.parent
PROBE='continuing13_20260908_lrc_zero_clock_probe_certificate.json'
PROBE_SHA='a56f8d8bdee84f642d7f837a6c9e6398008eadbc0318b0778d6d0e435c0b0201'
gates=0
def need(ok,why):
    global gates
    gates+=1
    if not ok:raise ArithmeticError(why)
def canonical(v):return json.dumps(v,sort_keys=True,separators=(',',':')).encode()
def vp(v,p):
    k=0
    while v%p==0:v//=p;k+=1
    return k
def allowed(s):
    p=2
    while p*p<=s:
        e=0
        while s%p==0:s//=p;e+=1
        if e and (p%3!=2 or e>2):return False
        p+=1
    return s==1 or s%3==2
def geometry(p,q):
    L=14*p*q
    A=[(max(0,14*k*q-q),min(L,14*k*q+q)) for k in range(p+1)]
    B=[(max(0,14*k*p-p),min(L,14*k*p+p)) for k in range(q+1)]
    i=j=0;out=[]
    while i<len(A) and j<len(B):
        a,b=max(A[i][0],B[j][0]),min(A[i][1],B[j][1])
        if a<b:out.append((a,b))
        if A[i][1]<B[j][1]:i+=1
        elif A[i][1]>B[j][1]:j+=1
        else:i+=1;j+=1
    return L,[(out[-1][0]-L,out[0][1])]+out[1:-1]
def capacity(n,L,intervals):
    base=0;events={};cur=0
    for a,b in intervals:
        h,r=divmod(n*(b-a),L);base+=h
        A,B=n*a%L,n*b%L
        if not r:events.setdefault(A,[0,0,0])[2]+=1
        else:
            events.setdefault(A,[0,0,0])[0]+=1
            events.setdefault(B,[0,0,0])[1]+=1
            cur+=A>B
    best=10**9
    for wall,(starts,ends,holes) in sorted(events.items()):
        best=min(best,cur-ends-holes)
        cur+=starts-ends;best=min(best,cur)
    return base+best
def literal(n,L,intervals):
    walls=sorted({n*x%L for ab in intervals for x in ab})
    phases=[2*x for x in walls]+[(walls[i]+(walls[(i+1)%len(walls)]+(L if i+1==len(walls) else 0)))%(2*L) for i in range(len(walls))]
    counts=[sum(-((-(2*n*b-r))//(2*L))-(2*n*a-r)//(2*L)-1 for a,b in intervals) for r in phases]
    at=min(range(len(phases)),key=lambda i:counts[i])
    return counts[at],Q(phases[at],2*L)

raw=((ROOT/'05-knowledge/results'/PROBE) if HERE.parent.name=='04-computation' else HERE.with_name(PROBE)).read_bytes()
need(sha256(raw).hexdigest()==PROBE_SHA,'complete profile-word census pin')
J=json.loads(raw)
need(J['declared_clocks']==[9240,11088] and [r['t'] for r in J['clocks']]==[9240,11088],'exact declared clock universe')
AT=[(p,total-p,*geometry(p,total-p)) for total in range(3,357) if allowed(total) for p in range(1,(total+1)//2) if gcd(p,total-p)==1]
need(len(AT)==5855,'whole native coprime strict atlas')
expected_no_isolation={9240:[[6,8,8,8,12,15,15]],11088:[[6,8,8,8,9,9,36],[6,8,8,9,9,16,36],[6,8,9,9,16,16,36]]}
records=[]
for R in J['clocks']:
    t=R['t'];domain=R['domain'];pairkeys=list(combinations_with_replacement(domain,2))
    need(all(t%(7*d)==0 for d in domain),'entire divisor alphabet has exact zero budget')
    mins={ab:10**10 for ab in pairkeys};banks={ab:[] for ab in pairkeys};compatible=0
    for e in sorted({gcd(*ab) for ab in pairkeys}):
        n=t//e
        for p,q,L,intervals in AT:
            dp,dq=e*gcd(n,p),e*gcd(n,q);ab=tuple(sorted((dp,dq)))
            if ab not in mins:continue
            compatible+=1
            cap=capacity(n,L,intervals)
            need(0<=cap<=n,'complete unpruned native capacity bounds')
            mins[ab]=min(mins[ab],e*cap)
            if cap==0:
                # This evaluates every literal boundary and chamber for every zero ratio.
                direct,phase=literal(n,L,intervals)
                need(direct==0,'all-wall/chamber literal native zero witness')
                need(gcd(t,e*p)==dp and gcd(t,e*q)==dq,'each zero ratio has a separate exact native margin realization')
                banks[ab].append([p,q,dp,dq,phase.numerator,phase.denominator])
    inherited={tuple(ab):v[0] for ab,v in R['weights']}
    for ab in pairkeys:
        need(mins[ab]==inherited[ab],'all-ratio unpruned minimum agrees with full consumer')
        need((mins[ab]==0)==bool(banks[ab]),'zero adjacency iff full native bank nonempty')
    need(compatible==R['compatible_atlas_edges'],'complete compatible native universe agrees')
    shapes=Counter();isolated=0;exceptions=[];cuts=[]
    for word in R['words']:
        roots=list(range(7));degrees=[0]*7
        def find(i):
            while roots[i]!=i:i=roots[i]
            return i
        for i,j in combinations(range(7),2):
            if banks[tuple(sorted((word[i],word[j])))]:
                roots[find(i)]=find(j);degrees[i]+=1;degrees[j]+=1
        components={}
        for i in range(7):components.setdefault(find(i),[]).append(i)
        cc=sorted(components.values(),key=lambda I:(len(I),I))
        need(len(cc)>1,'every actual connected graph must cross a positive-capacity cut')
        shapes[tuple(sorted(map(len,cc)))]+=1
        isos=[i for i,d in enumerate(degrees) if not d]
        if isos:
            isolated+=1
            cut=[isos[0]]
        else:
            exceptions.append(word)
            cut=cc[0]
        credit=min(mins[tuple(sorted((word[i],word[j])))] for i in cut for j in range(7) if j not in cut)
        need(credit>0,'selected whole cut has strictly positive native capacity on every possible crossing')
        if not isos:
            cuts.append(dict(word=word,components=cc,component_margins=[[word[i] for i in C] for C in cc],cut=cut,minimum_crossing_credit=credit))
    need(exceptions==expected_no_isolation[t],'whole exceptional bank to the isolated-vertex shortcut')
    need(sum(shapes.values())==R['word_count'],'complete positional word coverage')
    zero_graph_components=[];parent={d:d for d in domain}
    def findval(d):
        while parent[d]!=d:d=parent[d]
        return d
    for ab,bank in banks.items():
        if bank:parent[findval(ab[0])]=findval(ab[1])
    groups={}
    for d in domain:groups.setdefault(findval(d),[]).append(d)
    zero_graph_components=sorted(groups.values())
    need(any(len(C)>1 and gcd(*C)==1 for C in zero_graph_components),'hostile: whole-alphabet zero connectivity cannot replace seven-label incidence')
    records.append(dict(clock=t,domain=domain,compatible_native_ratios=compatible,zero_ratio_count=sum(map(len,banks.values())),
        zero_banks=[[list(ab),bank] for ab,bank in banks.items() if bank],word_count=R['word_count'],isolated_word_count=isolated,
        shapes=[[list(shape),count] for shape,count in sorted(shapes.items())],no_isolation=cuts,alphabet_components=zero_graph_components))
    print('CLOCK',t,'NATIVE',compatible,'ZERO_RATIOS',sum(map(len,banks.values())),'WORDS',R['word_count'],'ISOLATED',isolated,'EXCEPTIONAL_CUTS',len(cuts),flush=True)
need(J['removed_scales']==[9240,11088] and len(J['new_scales'])==7618,'exact two-clock deletion only')
need(J['new_scales_sha256']=='89a3e7545cc77467c5f85fbe0bdaab1f71071cff65184bbc1212eb86c9f07177','final necessary array pin')
cert=dict(status='FINITE-EXACT complete native zero banks and all126-profile positional cuts; analytic connected-complement consumer in report',
    probe_sha256=PROBE_SHA,declared_clocks=[9240,11088],records=records,new_scales=J['new_scales'],new_scales_sha256=J['new_scales_sha256'],
    removed_scales=[9240,11088],gates=gates,scope='Actual connected strict graph on seven complementary labels; selected-six gcd clock; all inherited126profiles; weak clearance only.')
data=canonical(cert)+b'\n';(DEST/(HERE.stem+'_certificate.json')).write_bytes(data)
print('COMPLETE ZERO CUTS: all41,249 words split; no residual reaches a native-depth or endpoint obligation.')
print('NECESSARY CLOCKS',len(J['new_scales']),'MAXIMUM',max(J['new_scales']))
print('NEW_SCALES_SHA256',J['new_scales_sha256'])
print('CERTIFICATE_SHA256',sha256(data).hexdigest())
print('Always-active exact gates:',gates)
