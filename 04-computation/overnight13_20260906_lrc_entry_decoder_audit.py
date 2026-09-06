"""Independent complete-entry referee: complement boxes and literal ranks.

No producer imports. Negative membership is checked by bounded nonnegative
representations of Q(a+b)-x, with the opposite modular variable from producer.
"""
from itertools import combinations
from math import gcd
from functools import reduce
from pathlib import Path
import hashlib,json,sys

sys.stdout.reconfigure(newline='\n')
GATES=0
Q=91**6

def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)

def cgcd(xs):return reduce(gcd,xs,0)
def ceildiv(a,b):return -((-a)//b)

def complement_box(q,a,b,x):
    sign=-1 if x<0 else 1;x=abs(x)
    if x>q*(a+b):return False,None
    M=q*(a+b)-x
    v0=0 if a==1 else (pow(b,-1,a)*M)%a
    u0=(M-b*v0)//a
    lo=max(0,ceildiv(u0-2*q,b))
    hi=min((2*q-v0)//a,u0//b)
    if lo>hi:return False,None
    u=u0-b*lo;v=v0+a*lo
    return True,(sign*(q-u),sign*(q-v))

def support(q,row,i,j,k):
    if row[i]>row[j]:i,j=j,i
    D=gcd(row[i],row[j]);a=row[i]//D;b=row[j]//D
    delta=gcd(D,row[k]);c=D//delta;x=row[k]//delta
    if c>q:return False,None
    yes,coeff=complement_box(q,a,b,x)
    if not yes:return False,None
    r,s=coeff;w=[0]*len(row);w[i]=-r;w[j]=-s;w[k]=c
    return True,w

def literal_rows(q,row):
    """All nonzero support<=3 relations, unique up to overall sign.

    The first occupied coefficient is positive. Enumerate two coefficients
    and solve the last by divisibility; no box or minimal-coefficient lemma.
    """
    out=[];n=len(row)
    for i,j in combinations(range(n),2):
        for r in range(1,q+1):
            value=-r*row[i]
            if value%row[j]:continue
            s=value//row[j]
            if 0<abs(s)<=q:
                w=[0]*n;w[i]=r;w[j]=s;out.append(w)
    for i,j,k in combinations(range(n),3):
        for r in range(1,q+1):
            for s in range(-q,q+1):
                if not s:continue
                value=-r*row[i]-s*row[j]
                if value%row[k]:continue
                t=value//row[k]
                if 0<abs(t)<=q:
                    w=[0]*n;w[i]=r;w[j]=s;w[k]=t;out.append(w)
    return out

def rational_rank(rows):
    # Fraction-free elimination followed by integer content removal.
    basis={}
    for w in rows:
        w=list(w)
        for p in sorted(basis):
            if w[p]:
                a=w[p];b=basis[p][p]
                w=[b*x-a*y for x,y in zip(w,basis[p])]
                d=cgcd(abs(x) for x in w)
                if d:w=[x//d for x in w]
        p=next((i for i,x in enumerate(w) if x),None)
        if p is not None:
            d=cgcd(abs(x) for x in w)
            if w[p]<0:d=-d
            basis[p]=[x//d for x in w]
    return len(basis)

def allowed_sums():
    # Sieve primes, then multiplicatively build sums; no per-pair factorization.
    prime=[True]*357;prime[0]=prime[1]=False
    for p in range(2,19):
        if prime[p]:
            for v in range(p*p,357,p):prime[v]=False
    allowed={1}
    for p in range(2,357):
        if not prime[p] or p%3!=2:continue
        old=list(allowed)
        allowed|={v*p**e for v in old for e in (1,2) if v*p**e<=356}
    return allowed

SUMS=allowed_sums()

def graph(row):
    n=len(row);parent=list(range(n));edges=[];vectors=[]
    def find(i):
        while i!=parent[i]:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i,j in combinations(range(n),2):
        d=gcd(row[i],row[j]);a=row[i]//d;b=row[j]//d
        if a+b in SUMS:
            parent[find(i)]=find(j)
            edges.append((i,j,min(a,b),max(a,b)))
            w=[0]*n;w[i]=b;w[j]=-a;vectors.append(w)
    components=sorted([sorted(i for i in range(n) if find(i)==p)
                       for p in {find(i) for i in range(n)}],key=lambda x:(len(x),x))
    return components,edges,vectors

def orientations(I,J):
    return [(i,j,k) for i,j in combinations(I,2) for k in J]+[(i,j,k) for i,j in combinations(J,2) for k in I]

def check_rel(q,row,w,I,label):
    need(len(w)==len(row) and sum(x*y for x,y in zip(w,row))==0,label+' literal kernel')
    need(0<sum(x!=0 for x in w)<=3 and max(map(abs,w))<=q,label+' literal support and height')
    need(sum(row[i]*w[i] for i in I)!=0,label+' actual nonzero component partial sum')

def main():
    here=Path(__file__).parent
    names={
      'overnight13_20260906_lrc_entry_decoder.py':'9b8fe6804a037a0f93a29396940e510691ee3b2765be1eeca1271ed28e5b7b6c',
      'overnight13_20260906_lrc_entry_decoder.out':'8932a324412f5c0b7ce8cfc243146af374454a2239aa6059030f00af261bd57d',
      'overnight13_20260906_lrc_entry_decoder_certificates.json':'af177ea9703390a8169d9466d04a5e4d8bb53f66833cfbfd223c94a0006cb847'}
    for name,h in names.items():
        p=here/name
        if not p.exists() and name.endswith('.out'):p=here.parent/'05-knowledge'/'results'/name
        need(hashlib.sha256(p.read_bytes()).hexdigest()==h,'frozen producer '+name)
    data=json.loads((here/'overnight13_20260906_lrc_entry_decoder_certificates.json').read_text())
    need(data['Q']==Q,'actual coefficient height')
    boxcases=0
    for q in range(2,8):
        for b in range(2,q+1):
            for a in range(1,b):
                if gcd(a,b)>1:continue
                raw={a*r+b*s for r in range(-q,q+1) for s in range(-q,q+1)}
                L=q*(a+b)
                for x in range(-L-1,L+2):
                    yes,w=complement_box(q,a,b,x)
                    need(yes==(x in raw),'opposite-congruence complement membership versus full coefficient set')
                    if yes:
                        r,s=w
                        need(abs(r)<=q and abs(s)<=q and a*r+b*s==x,'complement witness')
                    boxcases+=1
    toyrows=splits=literalcount=equalitycount=0
    for n,ceiling,qs in ((4,13,range(4,9)),(6,10,range(5,10))):
        for q in qs:
            for row in combinations(range(1,ceiling+1),n):
                if cgcd(row)!=1 or sum(row)>q*q:continue
                rels=literal_rows(q,row);rank=rational_rank(rels)
                need(rank>=n-2,'literal finite-box rank bound')
                literalcount+=len(rels);toyrows+=1
                # n=4 uses each unoriented 2+2 partition once; n=6 uses all 4+2.
                for J in combinations(range(n),2):
                    if n==4 and 0 not in J:continue
                    I=tuple(i for i in range(n) if i not in J)
                    inside=all(sum(row[i]*w[i] for i in I)==0 for w in rels)
                    actual_equal=inside and rank==n-2
                    heights=all(max(row[i],row[j])//gcd(row[i],row[j])<=q
                                for part in (I,J) for i,j in combinations(part,2))
                    predicted=False
                    if heights:
                        positives=[support(q,row,i,j,k)[0] for i,j,k in orientations(I,J)]
                        predicted=not any(positives)
                    need(predicted==actual_equal,'complete arithmetic criterion versus full literal rational span')
                    equalitycount+=actual_equal;splits+=1
    controls=[]
    for case in data['cases']:
        row=case['speeds'];n=len(row)
        need(n==13 and len(set(row))==13 and min(row)>0 and cgcd(row)==1,'actual primitive row type')
        comps,edges,vecs=graph(row)
        need(sorted(map(len,comps))==[2,11],'independent all-scale actual components')
        need(rational_rank(vecs)==11,'literal weighted-edge rational rank')
        need(all(max(map(abs,w))<=355 for w in vecs),'decoder span lies inside coefficient budget')
        domain=sum(row)<=Q*Q
        need(case['domain']==domain,'physical-box domain status')
        J,I=comps
        if not domain:
            need(case['reason']=='outside_physical_box' and 'entry' not in case,'outside domain is not a false equality rejection')
            K=max(row[i] for i in I);g=min(row[i] for i in J)
            need(g>2*Q*K,'outside-box all-orientation crossing exclusion by exact dominance')
            need(max(row[i]//gcd(row[i],row[j]) for i in I for j in I)>Q,'outside-box internal-height hostile')
            controls.append(dict(name=case['name'],domain=False,decoder_rank=11,actual_equal=True))
            continue
        need(case['components']==comps,'full component certificate')
        need(case['decoder_edges']==[list(e) for e in edges],'all 78 edge decisions, retained positive edges')
        expected_heights={(i,j):max(row[i],row[j])//gcd(row[i],row[j])
                          for part in (I,J) for i,j in combinations(part,2)}
        need({tuple(r['indices']):r['height'] for r in case['internal_heights']}==expected_heights,'complete internal-height certificate')
        need(max(expected_heights.values())<=Q,'all minimal-coefficient hypotheses hold')
        required=set(orientations(I,J));seen=set();positive=[];counts=[0,0]
        for rec in case['supports']:
            i,j=rec['pair_indices'];k=rec['distinguished_index'];key=(i,j,k)
            need(key in required and key not in seen,'one complete mixed support, no duplicate or omitted orientation')
            seen.add(key);D=gcd(row[i],row[j]);a=row[i]//D;b=row[j]//D
            delta=gcd(D,row[k]);c=D//delta;x=row[k]//delta
            need([rec[z] for z in ('a','b','pair_gcd','delta','c','x')]==[a,b,D,delta,c,x],'physical clearing data')
            r0=rec['r0'];s0=rec['s0']
            need(0<=r0<b and a*r0+b*s0==x,'producer base solution is exact')
            need(rec['lower']==max(ceildiv(-Q-r0,b),ceildiv(s0-Q,a)) and
                 rec['upper']==min((Q-r0)//b,(s0+Q)//a),'claimed integer certificate interval endpoints')
            yes,w=support(Q,row,i,j,k)
            need(yes==rec['crossing'],'complement engine independently replays every support outcome')
            orientation=0 if k in J else 1;counts[orientation]+=yes
            need(rec['orientation']==('two_core_one_pair' if orientation==0 else 'one_core_two_pair'),'distinguished-label typing')
            if yes:
                check_rel(Q,row,w,I,'independent positive witness')
                check_rel(Q,row,rec['witness'],I,'frozen positive witness')
                positive.append(w)
            else:
                need('witness' not in rec,'negative record has no false witness')
        need(seen==required and len(seen)==121,'all 121 supports exhaust both orientations')
        need(case['entry']==(not positive),'full entry decision')
        rank=rational_rank(vecs+positive)
        need(rank==(11 if not positive else 12),'literal decoder-plus-crossing rational rank')
        controls.append(dict(name=case['name'],domain=True,rank=rank,positive_counts=counts))
    # Separate small exact versions of the two theorem hostiles.
    q=41;t=3*q+1;row=[t,4*t,6*t,1,3];I=(0,1,2);J=(3,4)
    need(sum(row)<=q*q and cgcd(row)==1,'small opposite-orientation hostile lies in finite box')
    comps,edges,vecs=graph(row)
    need(comps==[list(J),list(I)],'small hostile uses actual decoder components')
    rels=literal_rows(q,row)
    need(rational_rank(rels)==4 and rational_rank(vecs)==3,'small exact crossing increases literal relation rank')
    need(not any(support(q,row,i,j,k)[0] for i,j in combinations(I,2) for k in J),'small analogue of all 110 negative')
    need(sum(support(q,row,*J,k)[0] for k in I)==1,'exactly one omitted orientation is positive')
    check_rel(q,row,[1,0,0,-1,-q],I,'small omitted-orientation witness')
    q=7;row=[1,3,9,127,381];I=(0,1,2);J=(3,4)
    rels=literal_rows(q,row);comps,edges,vecs=graph(row)
    need(sum(row)>q*q and 9>q,'small physical-box boundary hostile')
    need(comps==[list(J),list(I)] and rational_rank(vecs)==3,'small outside-box actual decoder graph and rank')
    need(rational_rank(rels)==3 and all(sum(row[i]*w[i] for i in I)==0 for w in rels),
         'outside box literal equality survives a failed internal-height gate')
    positive_toys=[]
    for q,core in ((53,[1,3]),(101,[1,3,9])):
        g=2*q*max(core)+1;row=core+[g,3*g];n=len(row)
        I=tuple(range(n-2));J=(n-2,n-1)
        rels=literal_rows(q,row);comps,edges,vecs=graph(row)
        need(sum(row)<=q*q and cgcd(row)==1,'positive toy is a physical finite-box row')
        need({tuple(c) for c in comps}=={I,J},'positive toy has the actual decoder graph')
        need(rational_rank(vecs)==n-2 and rational_rank(rels)==n-2,'positive toy exact equality of rational ranks')
        need(all(sum(row[i]*w[i] for i in I)==0 for w in rels),'every literal positive-toy relation remains internal')
        need(not any(support(q,row,i,j,k)[0] for i,j,k in orientations(I,J)),
             'positive toy passes every arithmetic orientation')
        positive_toys.append(dict(Q=q,row=row,rank=n-2,literal_relations=len(rels)))
    # An additional in-box actual 13-label height-gate branch. This is supplied
    # separately from the producer JSON and checked without a box membership
    # call on the inadmissible large primitive pair.
    core=[355**j for j in range(6)]+[3,9,27,81,243]
    K=max(core);g=1000*K+1;row=core+[g,3*g];I=tuple(range(11))
    comps,edges,vecs=graph(row)
    need(len(set(row))==13 and cgcd(row)==1 and sum(row)<=Q*Q,'actual finite-box height-gate control')
    need(comps==[[11,12],list(I)] and rational_rank(vecs)==11,'actual height-gate graph and decoder rank')
    need(K>Q,'actual height-gate obstruction is triggered')
    w=[0]*13;w[0]=-1;w[5]=-1000;w[11]=1
    check_rel(Q,row,w,I,'actual height-gate explicit crossing')
    need(rational_rank(vecs+[w])==12,'actual height-gate crossing adds the missing rank')
    print('STATUS: INDEPENDENT ANALYTIC AND EXACT PASS; complete actual 11+2 equality decoder')
    print('UNIVERSES:',boxcases,'complement-box integers;',toyrows,'literal physical toy rows;',splits,'complete partitions;',literalcount,'literal relation rows;',equalitycount,'toy equality partitions')
    print('FULL CERTIFICATES:',json.dumps(controls,sort_keys=True))
    print('POSITIVE LITERAL TOYS:',json.dumps(positive_toys,sort_keys=True))
    print('SMALL HOSTILES: Q41 opposite orientation raises rank 3 to 4; Q7 outside-box equality rank 3 survives internal height 9')
    print('ACTUAL HEIGHT-GATE CONTROL: K=355^5>Q; g=1000K+1; decoder rank 11; explicit crossing raises rank to 12 inside Q^2 box')
    print('SCOPE: actual-entry equality only; no safe-phase conclusion; finite-box necessity is retained')
    print('SEMANTIC SHA256:',hashlib.sha256(json.dumps(controls,sort_keys=True).encode()).hexdigest())
    print('ACTIVE GATES:',GATES)

if __name__=='__main__':main()
