"""Clock7560 closure by a forced zero-credit wedge with incompatible depths.

No producer imports. All native ratios, phases, margin roles and the complete
six-word residual are retained. Analytic proof is in the matching report.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
from math import gcd
from hashlib import sha256
import json,sys
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve()
ROOT=HERE.parent.parent if HERE.parent.name=='04-computation' else Path.cwd()
DEST=ROOT/'05-knowledge/results' if HERE.parent.name=='04-computation' else HERE.parent
PROBE='continuing12_20260908_lrc_zero_clock_probe_certificate.json'
PROBE_SHA='e8caa081f3d9e8bffb2c7e4a003256ae4b39e1116799fd97294ceac014febf2b'
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
need(sha256(raw).hexdigest()==PROBE_SHA,'complete first-residual producer pin')
J=json.loads(raw)
need([r['t'] for r in J['clocks']]==[6930,7560],'first residual stopping rule retained')
R=J['clocks'][1];weights={tuple(ab):row[0] for ab,row in R['weights']}
need(R['word_count']==37209 and len(R['survivors'])==6,'complete residual size')
expected=[(4,4,4,9,9,18,24),(4,4,4,9,18,24,27),(4,4,8,9,9,18,24),(4,4,8,9,18,24,27),(4,8,8,9,9,18,24),(4,8,8,9,18,24,27)]
need([tuple(x[0]) for x in R['survivors']]==expected,'whole six-word bank, not selected owners')
need(all(e==m==0 for _,e,m in R['survivors']),'all residual budgets exactly zero')
need(all(7560%(7*d)==0 for d in R['domain']),'zero budget for whole declared divisor alphabet')
topology=[]
for word in expected:
    mid=word.index(24);odd={i for i,d in enumerate(word) if d in (9,18,27)}
    cross=[];zero=[]
    for i,j in combinations(range(7),2):
        if weights[tuple(sorted((word[i],word[j])))]==0:
            zero.append([i,j])
            if (i in odd)!=(j in odd):cross.append([i,j])
    need(len(cross)==1 and {word[i] for i in cross[0]}=={18,24},'unique zero-edge bridge to odd side')
    four=[i for i,d in enumerate(word) if d==4]
    need(bool(four) and word.count(24)==word.count(18)==1,'actual unique middle and bridge label')
    for i in four:
        neighbors=[j for j in range(7) if j!=i and weights[tuple(sorted((word[i],word[j])))]==0]
        need(neighbors==[mid],'every4 label has only24 as a possible zero neighbor')
    topology.append(dict(word=word,zero_edges=zero,forced_first=[four[0],mid],forced_second=cross[0]))

AT=[(p,total-p) for total in range(3,357) if allowed(total) for p in range(1,(total+1)//2) if gcd(p,total-p)==1]
need(len(AT)==5855,'complete strict atlas')
bank_records=[]
for end in [4,18]:
    middle=24;e=gcd(end,middle);n=7560//e
    rows=[];compatible=0;positive=0
    for p,q in AT:
        dp,dq=e*gcd(n,p),e*gcd(n,q)
        if sorted([dp,dq])!=[end,middle]:continue
        compatible+=1
        L,intervals=geometry(p,q)
        count=capacity(n,L,intervals)
        need(count>=0,'actual strict capacity nonnegative')
        if count:positive+=1;continue
        direct,phase=literal(n,L,intervals)
        need(direct==count==0,'literal all-wall/chamber zero control')
        numerator,denominator=(p,q) if (dp,dq)==(end,middle) else (q,p)
        depth=vp(end,2)+vp(denominator,2)-vp(numerator,2)
        need(vp(end,2)<vp(7560,2),'endpoint dyadic valuation is unsaturated')
        need(gcd(7560,e*numerator)==end and gcd(7560,e*denominator)==middle,'each arm separately has a literal native realization')
        rows.append(dict(ratio=[numerator,denominator],primitive_pair=[p,q],middle_depth=depth,
            separate_realization=[e*numerator,e*denominator],zero_phase=[phase.numerator,phase.denominator]))
    bank_records.append(dict(end=end,middle=middle,e=e,n=n,compatible=compatible,positive=positive,zero=rows))
need([r['ratio'] for r in bank_records[0]['zero']]==[[307,48]],'complete first native zero bank')
need([r['ratio'] for r in bank_records[1]['zero']]==[[309,8],[303,16],[219,128],[201,152]],'complete second native zero bank')
depthA={r['middle_depth'] for r in bank_records[0]['zero']};depthB={r['middle_depth'] for r in bank_records[1]['zero']}
need(depthA=={6} and depthB=={4,5,8},'exact lifted dyadic depth requirements')
need(depthA.isdisjoint(depthB),'no actual zero-credit wedge exists')
# The smallest inherited clipped-depth hostile: both edges feasible, full row wrong.
need((gcd(2,1),gcd(2,4))==(1,2) and (gcd(2,2),gcd(2,3))==(2,1),'local arm controls in saturated-middle hostile')
need(tuple(gcd(2,x) for x in (1,4,6))!=(1,2,1),'joint depth is an additional requirement')
old=J['new_scales'];need(len(old)==7621 and 7560 in old,'only prior6930 deletion has been inherited')
new=[t for t in old if t!=7560]
need(len(new)==7620 and max(new)==11935,'exact final necessary-array update')
cert=dict(status='FINITE-EXACT native banks and topology; analytic valuation closure in companion report',
    probe_sha256=PROBE_SHA,clock=7560,profile_word_count=37209,residual_count=6,
    topology=topology,zero_banks=bank_records,middle_depth_sets=[sorted(depthA),sorted(depthB)],
    removed_scales=[7560],new_scales=new,new_scales_sha256=sha256(canonical(new)).hexdigest(),
    scope='Primitive thirteen distinct speeds, selected-six gcd clock, actual strict graph on complement connected; all inherited126profiles; weak clearance only.',gates=gates)
raw=canonical(cert)+b'\n';(DEST/(HERE.stem+'_certificate.json')).write_bytes(raw)
print('COMPLETE7560: all six full-profile residuals force actual zero-edge path4--24--18.')
print('NATIVE BANKS:',[(b['end'],b['middle'],b['compatible'],len(b['zero'])) for b in bank_records])
print('DYADIC MIDDLE DEPTHS:',sorted(depthA),sorted(depthB),'DISJOINT')
print('NECESSARY CLOCKS:',len(new),'MAXIMUM',max(new))
print('NEW_SCALES_SHA256',sha256(canonical(new)).hexdigest())
print('CERTIFICATE_SHA256',sha256(raw).hexdigest())
print('Always-active exact gates:',gates)
