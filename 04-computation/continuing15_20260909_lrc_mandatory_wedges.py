"""Complete remaining-zero-clock bridge catalogue and bounded actual-ratio decoder.

Run from the repository root. No historical producer is imported. The companion
enumerates complete seven-multisets and tests every one of the 126 profiles.
"""
from pathlib import Path
from math import gcd,lcm,comb
from fractions import Fraction as F
from collections import Counter,defaultdict
from functools import lru_cache
from itertools import combinations,permutations
from hashlib import sha256
import json,subprocess,sys,tempfile

sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve()
ROOT=HERE.parent.parent if HERE.parent.name=='04-computation' else Path.cwd()
DEST=ROOT/'05-knowledge/results' if HERE.parent.name=='04-computation' else HERE.parent
PROFILE_SHA='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f'
INPUT_SHA='863d33accf19a5515b3ecd2048aac35d77f30f4e004371c73f5760d345bbb240'
ARRAY_SHA='b603ca0186c1ebaa2318782348a7446a16f529543147eb53062f1ba360e8e57b'
gates=0
def need(ok,why):
    global gates
    gates+=1
    if not ok:raise ArithmeticError(why)
def canonical(v):return json.dumps(v,sort_keys=True,separators=(',',':')).encode()
def semantic(v):return sha256(canonical(v)).hexdigest()
def normalized(t,a,b,c):
    g=gcd(a,b,c)
    return (t//g,a//g,b//g,c//g)
def atlas_sum(s):
    p=2
    while p*p<=s:
        k=0
        while s%p==0:s//=p;k+=1
        if k and (p%3!=2 or k>2):return False
        p+=1
    return s==1 or s%3==2
def geometry(p,q):
    L=14*p*q
    A=[(max(0,(14*k-1)*q),min(L,(14*k+1)*q)) for k in range(p+1)]
    B=[(max(0,(14*k-1)*p),min(L,(14*k+1)*p)) for k in range(q+1)]
    out=[];i=j=0
    while i<len(A) and j<len(B):
        a,b=max(A[i][0],B[j][0]),min(A[i][1],B[j][1])
        if a<b:out.append((a,b))
        if A[i][1]<B[j][1]:i+=1
        elif A[i][1]>B[j][1]:j+=1
        else:i+=1;j+=1
    return L,[(out[-1][0]-L,out[0][1])]+out[1:-1]
def sweep(n,L,arcs):
    base=cur=0;events={}
    for a,b in arcs:
        h,r=divmod(n*(b-a),L);base+=h;A,B=n*a%L,n*b%L
        if not r:events.setdefault(A,[0,0,0])[2]+=1
        else:
            events.setdefault(A,[0,0,0])[0]+=1
            events.setdefault(B,[0,0,0])[1]+=1
            cur+=A>B
    closed=opened=10**30
    for _,(starts,ends,holes) in sorted(events.items()):
        closed=min(closed,cur-ends-holes)
        cur+=starts-ends;opened=min(opened,cur)
    out=(base+min(closed,opened),base+opened)
    need(0<=out[0]<=out[1]<=n,'strict walls and chambers are valid native counts')
    return out
@lru_cache(None)
def capacity(n,p,q):return sweep(n,*geometry(p,q))
def literal_counts(n,L,arcs):
    walls=sorted({n*x%L for ab in arcs for x in ab})
    def count(r):return sum(-((-(2*n*b-r))//(2*L))-(2*n*a-r)//(2*L)-1 for a,b in arcs)
    wall=min(count(2*w) for w in walls)
    chamber=min(count(2*w+1) for w in walls)
    phase=min([2*w for w in walls]+[2*w+1 for w in walls],key=count)
    return (min(wall,chamber),chamber),phase
def grid_sets(t,U,phase_num=1,phase_den=2):
    den=t*phase_den
    return [{j for j in range(t) if 14*min((u*(phase_den*j+phase_num))%den,den-(u*(phase_den*j+phase_num))%den)<den} for u in U]
def primitive(r,s):
    L=lcm(r.denominator,s.denominator)
    v=[int(L*r),L,int(L*s)];g=gcd(*v)
    return [x//g for x in v]
def bridge_wedges(word,R):
    adj=[[] for _ in word]
    for i,j in combinations(range(7),2):
        if R[word[i],word[j]]:adj[i].append(j);adj[j].append(i)
    seen=[-1]*7;low=[-1]*7;bridges=[];clock=0
    def visit(i,parent):
        nonlocal clock
        seen[i]=low[i]=clock;clock+=1
        for j in adj[i]:
            if j==parent:continue
            if seen[j]<0:
                visit(j,i);low[i]=min(low[i],low[j])
                if low[j]>seen[i]:bridges.append(tuple(sorted((i,j))))
            else:low[i]=min(low[i],seen[j])
    visit(0,-1)
    if -1 in seen:return None
    m=Counter(word)
    if any(not R[a,a] and m[a]>sum(m[b]*R[a,b] for b in m if b!=a) for a in m):return None
    star=defaultdict(list)
    for i,j in bridges:star[i].append(j);star[j].append(i)
    single={};bounded={}
    for middle,ends in star.items():
        for i,j in combinations(sorted(ends),2):
            if word[i]>word[j]:i,j=j,i
            a,b,c=word[i],word[middle],word[j];key=(a,b,c)
            if min(R[a,b],R[c,b])==1:single.setdefault(key,(i,middle,j))
            if R[a,b]*R[c,b]<=256:bounded.setdefault(key,(i,middle,j))
    return single,bounded,sorted(bridges)

def main():
    raw=(ROOT/'04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
    need(sha256(raw).hexdigest()==PROFILE_SHA,'full profile supplier raw pin')
    P=json.loads(raw);V=P['levels']['6']['gcds']
    raw=(ROOT/'05-knowledge/results/continuing14_20260908_lrc_zero_classification_certificate.json').read_bytes()
    need(sha256(raw).hexdigest()==INPUT_SHA,'complete inherited zero-classification certificate raw pin')
    inherited=json.loads(raw);baseline=inherited['new_scales']
    need(semantic(baseline)==ARRAY_SHA and len(baseline)==6889 and max(baseline)==11934,'complete inherited necessary array')
    clocks=[r for r in inherited['clocks'] if r[5]]
    need([r[0] for r in clocks]==inherited['remaining_E0']==[t for t in baseline if t%7==0],'entire remaining zero-budget clock slice')
    need(len(clocks)==356 and max(r[0] for r in clocks)==7056 and sum(r[5] for r in clocks)==47568,'declared clock/word universe')
    need(len(V)==42 and all(d%7 for d in V),'whole sheet alphabet, no factor seven')
    ids=sorted({r[1] for r in clocks})
    need(len(ids)==121,'all distinct full alphabets needed by remaining clocks')
    table={(n,a,b):o for n,a,b,bits,c,o in inherited['quotient_zero_pairs']}
    domains=inherited['domains']
    rows=[(int(k),c,w) for k,L in P['levels'].items() for c,w in L['profiles']]
    lines=[str(len(rows))]+[' '.join(map(str,[k,c,*w])) for k,c,w in rows]
    lines+=[str(len(ids))]+[' '.join(map(str,[d,len(domains[d]['domain']),*domains[d]['domain']])) for d in ids]
    wordbanks={};domain_cert=[]
    with tempfile.TemporaryDirectory(prefix='continuing15_lrc_') as directory:
        tmp=Path(directory).resolve()
        need(tmp.parent==Path(tempfile.gettempdir()).resolve(),'resolved private temporary output directory')
        inp=tmp/'input.txt';inp.write_bytes(('\n'.join(lines)+'\n').encode())
        exe=tmp/('native.exe' if sys.platform=='win32' else 'native')
        build=subprocess.run(['g++','-O2','-std=c++17',str(HERE.with_suffix('.cpp')),'-o',str(exe)],capture_output=True)
        need(build.returncode==0,'native compiler: '+build.stderr.decode(errors='replace'))
        run=subprocess.run([str(exe),str(inp),str(tmp)],capture_output=True)
        need(run.returncode==0,'complete profile consumer: '+run.stderr.decode(errors='replace'))
        native=run.stdout.decode().replace('\r\n','\n').strip()
        need('DOMAINS 121 RAW_MULTISETS 7105596 ACCEPTED_WORDS 291811' in native,'entire unpruned domain universe')
        for line in (tmp/'counts.txt').read_text().splitlines():
            d,raw_count,count=map(int,line.split());D=domains[d]['domain'];data=(tmp/f'words_{d}.json').read_bytes();words=json.loads(data)
            need(raw_count==comb(len(D)+6,7)==domains[d]['unpruned_multisets'],'no omitted seven-multiset')
            need(count==len(words)==domains[d]['word_count'],'full accepted profile bank size')
            need(sha256(data).hexdigest()==domains[d]['words_sha256'],'byte-exact independent regeneration of whole inherited word bank')
            wordbanks[d]=words;domain_cert.append(dict(index=d,domain=D,raw_multisets=raw_count,word_count=count,sha256=sha256(data).hexdigest()))
    need(set(wordbanks)==set(ids),'all requested alphabets regenerated')
    singles=Counter();bounded=Counter();survey=[];residual_banks={}
    for t,d,total,closed,chamber,expected,*_ in clocks:
        D=domains[d]['domain']
        need(D==[a for a in V if t%a==0] and all(t%(7*a)==0 for a in D),'complete divisor alphabet and zero marginal excess')
        R={(a,b):table.get((t//gcd(a,b),min(a,b)//gcd(a,b),max(a,b)//gcd(a,b)),0) for a in D for b in D}
        count=0;eligible=[];bank=[]
        for word in wordbanks[d]:
            data=bridge_wedges(word,R)
            if data is None:continue
            count+=1;bank.append(word);S,B,bridges=data
            for key in S:singles[(t,*key)]+=1
            for key in B:bounded[(t,*key)]+=1
            if S or B:eligible.append(dict(word=word,singleton=[[list(k),list(v)] for k,v in sorted(S.items())],bounded=[[list(k),list(v)] for k,v in sorted(B.items())],bridges=bridges))
        need(count==expected,'every complete clock word survives exactly the inherited graph/slot predicates')
        survey.append(dict(t=t,domain_index=d,profile_word_count=total,residual_count=count,residual_sha256=semantic(bank),eligible=eligible))
        residual_banks[t]=bank
    singleton_keys=sorted({normalized(*key) for key in singles})
    bounded_keys=sorted({normalized(*key) for key in bounded})
    eligible_count=sum(sum(bool(w['bounded']) for w in r['eligible']) for r in survey)
    eligible_clocks=sum(any(w['bounded'] for w in r['eligible']) for r in survey)
    need((len(singles),len(singleton_keys))==(94,53),'declared whole singleton-bank wedge catalogue')
    need((len(bounded),len(bounded_keys),eligible_count,eligible_clocks)==(561,304,5424,77),'complete cheap bounded-product topology survey')
    AT=[(p,s-p,*geometry(p,s-p)) for s in range(3,357) if atlas_sum(s) for p in range(1,(s+1)//2) if gcd(p,s-p)==1]
    need(len(AT)==5855,'entire strict coprime atlas')
    native_scans=0;bank_records={}
    @lru_cache(None)
    def bank(t,a,b):
        nonlocal native_scans
        e=gcd(a,b);n=t//e;out=[];typed=0
        for p,q,L,I in AT:
            dp,dq=e*gcd(n,p),e*gcd(n,q)
            if {dp,dq}!={a,b}:continue
            typed+=1;native_scans+=1
            if sweep(n,L,I)[1]:continue
            if (dp,dq)==(a,b):out.append(F(p,q))
            if (dq,dp)==(a,b):out.append(F(q,p))
        need(len(out)==table.get((n,min(a,b)//e,max(a,b)//e),0),'complete oriented chamber-zero bank equals inherited cardinality')
        need(len(set(out))==len(out),'every oriented ratio appears once')
        out=tuple(sorted(out));bank_records[t,a,b]=dict(key=[t,a,b],quotient_clock=n,typed_native_pairs=typed,ratios=list(map(str,out)))
        return out
    evaluated={}
    def evaluate(key):
        if key in evaluated:return evaluated[key]
        T,A,B,C=key;products=[];stats=Counter()
        need(gcd(A,B,C)==1 and all(T%x==0 for x in (A,B,C)),'normalized primitive wedge margins divide clock')
        for r in bank(T,A,B):
            for s in bank(T,C,B):
                v=primitive(r,s);actual=[gcd(T,x) for x in v];endpoint=r/s;cap=None;n=None
                need(gcd(*v)==1 and F(v[0],v[1])==r and F(v[2],v[1])==s,'primitive ratio propagation is exact')
                if len(set(v))<3:kind='duplicate'
                elif actual!=[A,B,C]:kind='depth'
                else:
                    p,q=sorted((endpoint.numerator,endpoint.denominator));n=T//gcd(A,C)
                    cap=capacity(n,p,q)
                    need(bool(cap[0])==bool(cap[1]),'no chamber-only composite positive occurs in the declared arithmetic universe')
                    kind='positive' if cap[0]>0 else 'residual'
                stats[kind]+=1
                products.append(dict(r=str(r),s=str(s),primitive_triple=v,actual_margins=actual,endpoint_ratio=str(endpoint),endpoint_quotient=n,kind=kind,capacity=cap))
        need(len(products)==len(bank(T,A,B))*len(bank(T,C,B)),'complete unpruned Cartesian product of both oriented banks')
        result=dict(key=list(key),banks=[list(map(str,bank(T,A,B))),list(map(str,bank(T,C,B)))],counts=dict(sorted(stats.items())),closed=not stats['residual'],products=products)
        evaluated[key]=result;return result
    for key in singleton_keys:evaluate(key)
    initial_closed={k for k,r in evaluated.items() if r['closed']}
    def killed(row,proved,kind='bounded'):
        return [w for w in row['eligible'] if any(normalized(row['t'],*k) in proved for k,_ in w[kind])]
    singleton_deletions=[r['t'] for r in survey if len(killed(r,initial_closed,'singleton'))==r['residual_count']]
    need(singleton_deletions==[5880,7056] and sum(len(r['products']) for r in evaluated.values())==1980 and len(initial_closed)==33,'complete singleton-wedge arithmetic stage')
    candidates=sorted([r['t'] for r in survey if sum(bool(w['bounded']) for w in r['eligible'])==r['residual_count'] and r['t'] not in singleton_deletions],reverse=True)
    need(candidates==[6804,5376,5124],'predeclared descending newly wholly eligible boundary clocks')
    boundary=[];stop=None
    for t in candidates:
        row=next(r for r in survey if r['t']==t)
        keys=sorted({normalized(t,*k) for w in row['eligible'] for k,_ in w['bounded']})
        for key in keys:evaluate(key)
        proved={k for k,r in evaluated.items() if r['closed']}
        survivors=[w for w in residual_banks[t] if w not in [v['word'] for v in killed(row,proved)]]
        boundary.append(dict(t=t,complete_keys=[list(k) for k in keys],remaining_words=survivors))
        if survivors:stop=t;break
    need([r['t'] for r in boundary]==[6804,5376] and stop==5376,'stop first complete arithmetic residual')
    need(boundary[-1]['remaining_words']==[[3,3,6,8,12,16,16]],'exact final stopping word')
    proved={k for k,r in evaluated.items() if r['closed']};stats=Counter()
    for row in evaluated.values():stats.update(row['counts'])
    need(len(evaluated)==64 and len(proved)==37 and dict(stats)==dict(positive=2158,depth=630,duplicate=9,residual=258),'entire evaluated primitive-key/product classification')
    need(sum(stats.values())==3055,'complete declared ratio-product universe')
    removals=[];deletions=[];clock_results=[]
    for row in survey:
        witnesses=[]
        for w in row['eligible']:
            options=[(k,v) for k,v in w['bounded'] if normalized(row['t'],*k) in proved]
            if options:
                k,v=min(options);witnesses.append(w['word'])
                removals.append(dict(t=row['t'],word=w['word'],positions=v,key=list(normalized(row['t'],*k))))
        surviving=[w for w in residual_banks[row['t']] if w not in witnesses]
        if not surviving:deletions.append(row['t'])
        clock_results.append(dict(t=row['t'],old_count=row['residual_count'],removed=len(witnesses),new_count=len(surviving),new_residual_sha256=semantic(surviving),first_remaining=surviving[0] if surviving else None))
    need(len(removals)==1119 and deletions==[5880,6804,7056],'only complete clock closures change the necessary array')
    new=[t for t in baseline if t not in deletions];zero=[t for t in new if t%7==0]
    need(len(new)==6886 and max(new)==11934 and len(zero)==353 and max(zero)==6552,'final necessary array and zero-budget frontier')
    example=evaluated[3528,2,12,9]
    need(len(example['products'])==7 and all(r['capacity']==(55,56) for r in example['products']),'7056 mandatory4--24--18 endpoint overlap is uniformly110, with chamber minimum112')
    need(evaluated[2268,4,6,9]['counts']=={'positive':171},'6804 complete positive composite bank')
    need(evaluated[5880,4,20,5]['counts']=={'depth':2},'5880 primitive global-depth obstruction')
    controls=[]
    for p,q in [(1,1),(1,13),(1,112),(1,355),(11,263),(23,323),(33,320),(48,307)]:
        for n in [1,2,7,14,29,651,1260,1890,1939,1953,2485]:
            L,I=geometry(p,q);literal,phase=literal_counts(n,L,I)
            need(sweep(n,L,I)==literal,'independent all-wall and all-chamber direct interval counts')
            den=2*L*n
            direct=sum(14*min(p*(2*L*j+phase)%den,den-p*(2*L*j+phase)%den)<den and 14*min(q*(2*L*j+phase)%den,den-q*(2*L*j+phase)%den)<den for j in range(n))
            need(direct==literal[0],'literal native grid attains the computed all-phase minimum')
            controls.append([p,q,n,*literal,phase])
    need(primitive(F(1,4),F(3,2))==[1,4,6] and [gcd(2,x) for x in [1,4,6]]==[1,2,2],'inherited local-ratio compatibility does not imply common clipped depths')
    v=primitive(F(289,58),F(295,22))
    need(v==[3179,638,8555] and [gcd(1792,x) for x in v]==[1,2,1],'stopping wedge has distinct depth-compatible endpoints')
    leaf=evaluated[1792,1,2,1]
    need(leaf['counts']=={'duplicate':2,'residual':2},'full stopping two-leaf ratio product, not one selected ratio')
    U=[157007631,422522895,31510182,129461800,36729660,8763568,5333680]
    t=5376;W=[gcd(t,u) for u in U]
    need(gcd(*U)==1 and len(set(U))==7 and W==[3,3,6,8,12,16,16],'physical stopping row is primitive, distinct and has every actual margin')
    dangers=grid_sets(t,U);edges=[];allpairs=[]
    for i,j in combinations(range(7),2):
        h=gcd(U[i],U[j]);p,q=sorted((U[i]//h,U[j]//h));strict=p+q<=356 and atlas_sum(p+q)
        intersection=len(dangers[i]&dangers[j]);allpairs.append([i,j,p,q,strict,intersection])
        if strict:
            need(intersection==0,'every actual strict edge vanishes at the same original phase')
            edges.append([i,j,p,q])
    need(edges==[[0,2,58,289],[1,2,22,295],[2,4,163,190],[3,5,22,325],[4,5,68,285],[4,6,44,303]],'the full actual strict atlas graph is exactly the declared six-edge tree')
    reach={0}
    for _ in range(7):reach|={j if i in reach else i for i,j,_,_ in edges if i in reach or j in reach}
    need(len(reach)==7 and len(edges)==6,'actual graph is connected and acyclic')
    sizes=[len(D) for D in dangers];safe=t-len(set.union(*dangers))
    need(sizes==[768]*7 and safe==1600,'common native-edge zero and exact marginal attainment leave1600 safe points')
    need(any(not strict and overlap for i,j,p,q,strict,overlap in allpairs),'the missing overlap is carried by actual non-atlas pairs')
    physical=dict(clock=t,speeds=U,margins=W,phase=[1,2],danger_sizes=sizes,actual_strict_edges=edges,all_actual_pairs=allpairs,safe_points=safe,danger_set_sha256=[semantic(sorted(D)) for D in dangers])
    need([r['t'] for r in boundary if r['remaining_words']]==[5376],'honest stopping state, no inference of an unsafe row')
    certificate=dict(status='FINITE-EXACT; analytic promotion requires independent audit',scope='necessary clocks for the inherited primitive13/selected-six actual connected-complement setting; general LRC14 remains open',
        input_sha256=INPUT_SHA,profile_sha256=PROFILE_SHA,old_array_sha256=ARRAY_SHA,old_scales=baseline,declared_zero_clocks=[r[0] for r in clocks],domains=domain_cert,native_profile_output=native,
        singleton_signatures=[[list(k),n] for k,n in sorted(singles.items())],singleton_normalized_keys=[list(k) for k in singleton_keys],
        bounded_signatures=[[list(k),n] for k,n in sorted(bounded.items())],bounded_normalized_keys=[list(k) for k in bounded_keys],bounded_eligible_words=eligible_count,bounded_eligible_clocks=eligible_clocks,
        full_clock_topology_survey=survey,native_banks=[bank_records[k] for k in sorted(bank_records)],typed_native_scans=native_scans,
        evaluated_signatures=[evaluated[k] for k in sorted(evaluated)],product_counts=dict(sorted(stats.items())),closed_primitive_signatures=[list(k) for k in sorted(proved)],
        singleton_clock_closures=singleton_deletions,boundary_declared_order=candidates,boundary_results=boundary,arithmetic_stop=stop,arithmetic_untested=[5124],
        removed_word_witnesses=removals,clock_results=clock_results,new_clock_closures=deletions,new_scales=new,new_array_sha256=semantic(new),remaining_E0=zero,
        interval_controls=controls,physical_stopping_hostile=physical,python_gates=gates)
    path=DEST/(HERE.stem+'_certificate.json');path.write_bytes(json.dumps(certificate,indent=2,sort_keys=True).encode()+b'\n')
    print('COMPLETE_PROFILE_DOMAINS 121 UNPRUNED_MULTISETS 7105596 ACCEPTED_WORDS 291811')
    print('ZERO_CLOCKS 356 RESIDUAL_WORDS 47568 SINGLETON_SIGNATURES 94 NORMALIZED 53')
    print('BOUNDED_TOPOLOGY_SIGNATURES 561 NORMALIZED 304 ELIGIBLE_WORDS 5424 CLOCKS 77')
    print('EVALUATED_KEYS 64 PRODUCTS 3055 POSITIVE 2158 DEPTH 630 DUPLICATE 9 RESIDUAL 258')
    print('CLOSED_KEYS 37 REMOVED_WORDS 1119 CLOSED_CLOCKS 5880 6804 7056')
    print('BOUNDARY_ORDER 6804 5376 5124 STOP 5376 UNTESTED 5124')
    print('REMAINING 6886 MAX 11934 E0 353 MAX_E0 6552 ARRAY_SHA256',semantic(new))
    print('PHYSICAL_5376 COMMON_ZERO_TREE_EDGES 6 MARGINALS 768 SAFE_POINTS 1600')
    print('TYPED_NATIVE_SCANS',native_scans,'PYTHON_GATES',gates)

if __name__=='__main__':main()
