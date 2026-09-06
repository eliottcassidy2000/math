"""Bounded one-clock stopping certificate, not a complete new clock sieve.
All988 original candidates retained; only first107 located profiles evaluated.
The first located survivor is then removed by exact forced-word enumeration.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement
from collections import Counter
from math import gcd,comb
import hashlib,json,sys
sys.stdout.reconfigure(newline='\n')
T=16704;CAPS=(90,30,9,4,2,1,1);G=0

def need(ok,label):
    global G
    G+=1
    if not ok:raise ArithmeticError(label)

def digest(x):return hashlib.sha256(json.dumps(x,separators=(',',':')).encode()).hexdigest()
def atlas(n):
    p=2
    while p*p<=n:
        e=0
        while n%p==0:n//=p;e+=1
        if e and(p%3!=2 or e>2):return False
        p+=1
    return n==1 or n%3==2

def geometry(p,q):
 L=14*p*q
 a=[(max(0,14*k*q-q),min(L,14*k*q+q)) for k in range(p+1)]
 b=[(max(0,14*k*p-p),min(L,14*k*p+p)) for k in range(q+1)]
 i=j=0;I=[]
 while i<len(a) and j<len(b):
  lo=max(a[i][0],b[j][0]);hi=min(a[i][1],b[j][1])
  if lo<hi:I.append((lo,hi))
  if a[i][1]<b[j][1]:i+=1
  elif b[j][1]<a[i][1]:j+=1
  else:i+=1;j+=1
 return L,[(I[-1][0]-L,I[0][1])]+I[1:-1]

def located(n,p,q):
 L,I=geometry(p,q);w=sorted({n*x%L for arc in I for x in arc})
 pts=[(x,1) for x in w]+[(a+w[(i+1)%len(w)]+(L if i==len(w)-1 else 0),2) for i,a in enumerate(w)]
 def count(x,s):return sum(-((x-s*n*b)//(s*L))-(s*n*a-x)//(s*L)-1 for a,b in I)
 vals=[(count(x,s),F(x,s*L)%1) for x,s in pts]
 return min(vals),max(vals),len(w)

def main():
    here=Path(__file__).resolve().parent
    inherited='overnight12_20260906_lrc_decoder_descent_inherited_profiles.json'
    path=here/inherited
    if not path.is_file():path=Path('C:/w/s0905/04-computation')/inherited
    raw=path.read_bytes()
    need(hashlib.sha256(raw).hexdigest()=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','complete inherited profile pin')
    data=json.loads(raw)
    profiles={int(k):{(c,tuple(w)) for c,w in v['profiles']} for k,v in data['levels'].items()}
    allowed=data['levels']['6']['gcds']
    need(len(allowed)==42 and allowed==sorted({c for c,w in profiles[6]}),'complete allowed-state domain')
    tablepath=here.parent/'05-knowledge/results/third_20260906_grid_full_words.out' if here.name=='04-computation' else Path('C:/w/s0905/05-knowledge/results/third_20260906_grid_full_words.out')
    rows=[json.loads(s[8:]) for s in tablepath.read_text().splitlines() if s.startswith('MAXIMUM ')]
    need(digest([[t,E] for t,E,*_ in rows])=='ca6b6f562db1fc3632f8b7570b89a16020a981ae8aa130be200dc1bdcb4264ca','independently audited complete maximum table semantic pin')
    match=[row for row in rows if row[0]==T]
    need(len(match)==1 and match[0][1]==188,'exact inherited M16704')
    E=188
    avail={d:sum(d<=c for c in CAPS) for d in allowed if T%d==0}
    cost={d:d*((-(T//d))%7) for d in avail}
    entries=sorted(((cost[d],d) for d in avail for _ in range(avail[d])),reverse=True)
    need(len(avail)==18,'all eighteen clock divisors in allowed-state domain')
    for d,c in cost.items():need(c==7*d*((T+7*d-1)//(7*d))-T,'literal ceiling cost')
    def paired(a,b):
        if a not in avail or b not in avail or(a==b and avail[a]<2):return None
        u=Counter((a,b));rest=[]
        for c,d in entries:
            if u[d]:u[d]-=1
            else:rest.append((c,d))
        value=cost[a]+cost[b]+sum(c for c,d in rest[:5])
        need(len(rest)>=5 and value%7==0,'correct five-slot paired marginal maximum')
        return value//7
    candidates=[];pairs=components=0;filtered={'component_cost':0,'forced_margin':0,'paired_cost':0}
    es=[e for e in range(1,7) if T%e==0]
    need(es==[1,2,3,4,6],'every supplied e<=6 divisor retained')
    for p in range(1,178):
        for q in range(p+1,357-p):
            if gcd(p,q)>1 or not atlas(p+q):continue
            pairs+=1
            D=14*p*q
            numerators=[2*p]+[min(2*p,p+q-14*k) for k in range(1,(p+q-1)//14+1) for _ in range(2)]
            need(all(0<a<=2*p for a in numerators),'strict complete clipped length numerators')
            for e in es:
                n=T//e;sep=e*sum((n*a+D-1)//D-1 for a in numerators)
                if sep>E:filtered['component_cost']+=1;continue
                components+=1
                dp,dq=e*gcd(n,p),e*gcd(n,q)
                need(gcd(dp,dq)==e,'forced margins retain common sheet gcd')
                ep=paired(dp,dq)
                if ep is None:filtered['forced_margin']+=1;continue
                if sep>min(E,ep):filtered['paired_cost']+=1;continue
                candidates.append([e,p,q,dp,dq,sep,ep])
    candidates.sort()
    need(pairs==5855 and len(candidates)==988,'complete exact current candidate universe')
    need(sum(filtered.values())+len(candidates)==pairs*len(es),'every base candidate accounted for exactly once')
    prefix=[];survivor=None
    for index,row in enumerate(candidates):
        e,p,q,dp,dq,sep,ep=row
        low,high,walls=located(T//e,p,q)
        credit=e*low[0]
        need(credit>=sep,'minimum of sum dominates sum of minima')
        prefix.append(row+[credit,low[0],str(low[1]),high[0],str(high[1]),walls])
        if credit<=min(E,ep):survivor=(index,row,low,high,walls);break
    need(survivor is not None and len(prefix)==107,'explicit stopping at first located survivor')
    index,row,low,high,walls=survivor
    need(index==106 and row==[4,23,323,4,4,180,207],'first lexicographic survivor and exact old costs')
    need(low==(45,F(0)) and high==(92,F(75651,104006)) and walls==98,'exact first-survivor profile')
    e,p,q,dp,dq,sep,ep=row;n=T//e;L,I=geometry(p,q)
    need(len(I)==49 and sum(b-a for a,b in I)==2122,'complete49 intervals and raw length sum')
    need(Counter(b-a for a,b in I)==Counter({46:43,38:2,24:2,10:2}),'exact clipped interval distribution')
    # Separate literal danger tests at every critical phase and every cell.
    w=sorted({n*x%L for arc in I for x in arc})
    phases=[F(a,L) for a in w]+[F(a+w[(i+1)%len(w)]+(L if i==len(w)-1 else 0),2*L)%1 for i,a in enumerate(w)]
    for alpha in phases:
        den=n*alpha.denominator
        direct=sum(all(14*min((v*(alpha.numerator+j*alpha.denominator))%den,(-v*(alpha.numerator+j*alpha.denominator))%den)<den for v in (p,q)) for j in range(n))
        count=sum(-((alpha.numerator*L-n*b*alpha.denominator)//(L*alpha.denominator))-(n*a*alpha.denominator-alpha.numerator*L)//(L*alpha.denominator)-1 for a,b in I)
        need(count==direct,'complete native196-phase first-survivor check')
    # A sharply bounded full-word test, with no branch pruning or producer oracle.
    tests=[(k,J,tuple(j for j in range(7) if j not in J)) for k in range(1,7) for J in combinations(range(7),k)]
    def valid(S):
        if gcd(*S)!=1:return False
        for k,J,K in tests:
            c=gcd(*(S[i] for i in J));word=tuple(sorted(gcd(c,S[j]) for j in K))
            if(c,word) not in profiles[7-k]:return False
        return True
    best=-1;owner=None;total=valid_count=0
    for tail in combinations_with_replacement(sorted(avail),5):
        S=(4,4)+tail;total+=1
        if not valid(S):continue
        valid_count+=1;value=sum(d*((T+7*d-1)//(7*d)) for d in S)-T
        need(7*value==sum(cost[d] for d in S),'ceiling and residue costs agree for every full-valid completion')
        if value>best:best=value;owner=S
    need(total==comb(22,5)==26334 and valid_count==2422,'complete unpruned forced-word universe')
    need(best==134 and owner==(4,4,3,9,36,58,64),'exact forced4,4 maximum and attaining owner')
    for k,J,K in tests:
        c=gcd(*(owner[i] for i in J));word=tuple(sorted(gcd(c,owner[j]) for j in K))
        need((c,word) in profiles[7-k],'all126 words of the maximizing owner')
    need(not valid((4,4,72,58,64,36,29)),'independent bag maximum owner is not a full word')
    need(180>best and 180-best==46,'exact owner conditioning repairs the first located survivor')
    print('STATUS FINITE-EXACT bounded stopping result; no whole-clock elimination claim')
    print('UNIVERSE atlas5855 divisors',es,'old_candidates988 filters',json.dumps(filtered,sort_keys=True))
    print('FIRST index106 prefix107 edge4,23,323 margins4,4 old_sep180 located180 min45 at0 max92 Epair207 M188')
    print('CONDITIONED forced4,4 complete26334 valid2422 maximum134 owner',owner,'conditional_safe_surplus46')
    print('CANDIDATE_SCHEMA [e,p,q,dp,dq,Csep,Epair]')
    print('CANDIDATES',json.dumps(candidates,separators=(',',':')))
    print('CANDIDATE_SHA256',digest(candidates))
    print('PREFIX_SCHEMA candidate+[Cloc,min_count,min_phase,max_count,max_phase,wall_count]')
    print('LOCATED_PREFIX',json.dumps(prefix,separators=(',',':')))
    print('PREFIX_SHA256',digest(prefix))
    print('SCOPE candidates107..987 untested; first located survivor is removed by exact conditioning, not a physical hostile')
    print('PASS',G,'always-active exact gates; raw LF')

if __name__=='__main__':main()
