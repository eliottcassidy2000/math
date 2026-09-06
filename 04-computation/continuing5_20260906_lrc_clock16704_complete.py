"""Complete exclusion of one clock by compatible words and located pair counts.
No other clock is newly tested. Every full word and every eligible edge retained.
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
    name='overnight12_20260906_lrc_decoder_descent_inherited_profiles.json'
    path=here/name
    if not path.is_file():path=Path('C:/w/s0905/04-computation')/name
    raw=path.read_bytes();profile_pin='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f'
    need(hashlib.sha256(raw).hexdigest()==profile_pin,'complete inherited profile pin')
    data=json.loads(raw)
    profiles={int(k):{(c,tuple(w)) for c,w in v['profiles']} for k,v in data['levels'].items()}
    allowed=data['levels']['6']['gcds'];values=[d for d in allowed if T%d==0]
    need(len(allowed)==42 and allowed==sorted({c for c,w in profiles[6]}),'all inherited state values')
    need(values==[1,2,3,4,6,8,9,12,16,18,24,29,32,36,48,58,64,72],'complete18-divisor universe')
    tests=[(k,J,tuple(j for j in range(7) if j not in J)) for k in range(1,7) for J in combinations(range(7),k)]
    def valid(S):
        if gcd(*S)!=1:return False
        for k,J,K in tests:
            c=gcd(*(S[i] for i in J));word=tuple(sorted(gcd(c,S[j]) for j in K))
            if(c,word) not in profiles[7-k]:return False
        return True
    # One complete, unpruned multiset enumeration updates all pair maxima at once.
    conditional={};owners={};valid_words=[];total=0;best=-1
    for S in combinations_with_replacement(values,7):
        total+=1
        if not valid(S):continue
        valid_words.append(S)
        E=sum(d*((T+7*d-1)//(7*d)) for d in S)-T
        need(7*E==sum(d*((-(T//d))%7) for d in S),'native ceiling and residue cost for every valid word')
        best=max(best,E)
        for pair in set(combinations(S,2)):
            if E>conditional.get(pair,-1):conditional[pair]=E;owners[pair]=S
    need(total==comb(24,7)==346104 and len(valid_words)==19073,'complete full-word enumeration')
    need(best==188 and len(conditional)==161,'global and conditioned maxima reconstructed without prior optimizer')
    allpairs=list(combinations_with_replacement(values,2))
    need(len(allpairs)==171,'all unordered margin types, including repeated values')
    for pair,S in owners.items():
        need(valid(S) and all(S.count(d)>=pair.count(d) for d in set(pair)),'attaining owner contains the forced margins with multiplicity')
        need(sum(d*((T+7*d-1)//(7*d)) for d in S)-T==conditional[pair],'literal conditional owner cost')
    need(conditional[4,4]==134 and owners[4,4]==(3,4,4,9,36,58,64),'crossed boundary full-word result')
    need(conditional[12,16]==188,'other crossed boundary remains at full global cost')
    # Reconstruct the old988 candidates, rather than importing a consumer.
    avail={d:sum(d<=c for c in CAPS) for d in values}
    costs={d:d*((-(T//d))%7) for d in values}
    entries=sorted(((costs[d],d) for d in values for _ in range(avail[d])),reverse=True)
    def paired(a,b):
        if a not in avail or b not in avail or(a==b and avail[a]<2):return None
        used=Counter((a,b));rest=[]
        for c,d in entries:
            if used[d]:used[d]-=1
            else:rest.append(c)
        v=costs[a]+costs[b]+sum(rest[:5])
        need(v%7==0 and len(rest)>=5,'exact old forced marginal-cost relaxation')
        return v//7
    candidates=[];atlas_size=0;filters={'component_cost':0,'forced_margin':0,'paired_cost':0};es=[e for e in range(1,7) if T%e==0]
    need(es==[1,2,3,4,6],'all supplied small-edge gcd divisors')
    for p in range(1,178):
        for q in range(p+1,357-p):
            if gcd(p,q)>1 or not atlas(p+q):continue
            atlas_size+=1;D=14*p*q
            nums=[2*p]+[min(2*p,p+q-14*k) for k in range(1,(p+q-1)//14+1) for _ in range(2)]
            need(all(0<a<=2*p for a in nums),'full strict clipped length list')
            for e in es:
                n=T//e;sep=e*sum((n*a+D-1)//D-1 for a in nums)
                if sep>best:filters['component_cost']+=1;continue
                dp,dq=e*gcd(n,p),e*gcd(n,q);pair=tuple(sorted((dp,dq)))
                need(gcd(dp,dq)==e,'exact common sheet gcd')
                ep=paired(dp,dq)
                if ep is None:
                    need(pair not in conditional,'invalid old margins are absent from every full word')
                    filters['forced_margin']+=1;continue
                if sep>min(best,ep):
                    need(pair not in conditional or conditional[pair]<=ep,'old paired bound remains a necessary relaxation')
                    filters['paired_cost']+=1;continue
                candidates.append([e,p,q,dp,dq,sep,ep])
    candidates.sort()
    need(atlas_size==5855 and len(candidates)==988 and sum(filters.values())+988==5855*5,'complete inherited candidate accounting')
    need(digest(candidates)=='8610d17932b51afbaaeec4699f5a21e2362441f58b9c28235c48182106d45bc1','same complete988-candidate universe as crossed checkpoint')
    decisions=[];located_profiles=[];separate=0
    for row in candidates:
        e,p,q,dp,dq,sep,ep=row;pair=tuple(sorted((dp,dq)))
        need(pair in conditional and conditional[pair]<=min(best,ep),'actual full-word upper bound fits old relaxation')
        bound=conditional[pair]
        if sep>bound:
            credit=sep;mode='separate';separate+=1
        else:
            n=T//e;low,high,walls=located(n,p,q);credit=e*low[0];mode='located'
            need(credit>=sep,'same-phase minimum dominates separate minima')
            need(credit>bound,'every remaining candidate closes by its full located profile')
            located_profiles.append(dict(candidate=row,conditional=bound,minimum=low[0],min_phase=str(low[1]),maximum=high[0],max_phase=str(high[1]),walls=walls,credit=credit,surplus=credit-bound))
        need(credit>bound,'positive weak-safe grid surplus for every candidate')
        decisions.append(row+[bound,mode,credit,credit-bound])
    need(separate==835 and len(located_profiles)==153,'all988 candidates exhausted by two coupled tests')
    worst=min(located_profiles,key=lambda r:(r['surplus'],r['candidate']))
    need(worst['surplus']==18 and worst['candidate']==[6,132,221,72,6,174,226] and worst['conditional']==180 and worst['minimum']==33,'sharp observed located surplus control')
    # Literal grid controls for the least-surplus pair and both crossed boundaries.
    native_controls=[]
    for e,p,q in [(6,132,221),(4,3,308),(4,23,323)]:
        n=T//e;L,I=geometry(p,q)
        need(all(a<b for a,b in I),'all literal components nonempty')
        w=sorted({n*x%L for arc in I for x in arc})
        phases=[F(a,L) for a in w]+[F(a+w[(i+1)%len(w)]+(L if i==len(w)-1 else 0),2*L)%1 for i,a in enumerate(w)]
        vals=[]
        for alpha in phases:
            den=n*alpha.denominator
            direct=sum(all(14*min((v*(alpha.numerator+j*alpha.denominator))%den,(-v*(alpha.numerator+j*alpha.denominator))%den)<den for v in (p,q)) for j in range(n))
            count=sum(-((alpha.numerator*L-n*b*alpha.denominator)//(L*alpha.denominator))-(n*a*alpha.denominator-alpha.numerator*L)//(L*alpha.denominator)-1 for a,b in I)
            need(count==direct,'all critical walls/cells by independent native modular predicates');vals.append(count)
        low,high,nw=located(n,p,q)
        need(min(vals)==low[0] and max(vals)==high[0] and nw==len(w),'native control full extrema')
        native_controls.append([e,p,q,len(I),len(w),low[0],high[0]])
    # Change exactly one element of the already audited necessary clock set.
    path=here.parent/'05-knowledge/results/third_20260906_grid_refined.out' if here.name=='04-computation' else Path('C:/w/s0905/05-knowledge/results/third_20260906_grid_refined.out')
    matches=[json.loads(s[10:]) for s in path.read_text().splitlines() if s.startswith('SCALE_SET ') and json.loads(s[10:])['name']=='full_words_coupled']
    need(len(matches)==1,'unique inherited complete refined scale set')
    oldset=matches[0]['survivors'];oldpin='4f29481d984ead40d0144556ce1c45dce210e30b964bb65835a7904ca6353e59'
    need(digest(oldset)==oldpin and len(oldset)==8202 and max(oldset)==16704,'complete audited input scale set pin')
    newset=[t for t in oldset if t!=T]
    need(oldset==sorted(set(oldset)) and len(newset)==8201 and max(newset)==14904,'delete only proven clock16704, new maximum14904')
    record={'status':'PROOF CANDIDATE pending independent full-clock audit','clock':T,'profile_sha256':profile_pin,'allowed_divisors':values,'word_universe':total,'valid_word_count':len(valid_words),'valid_word_sha256':digest(valid_words),'global_maximum':best,'conditional_maxima':[[a,b,conditional.get((a,b)),owners.get((a,b))] for a,b in allpairs],'atlas_count':atlas_size,'edge_gcds':es,'old_filters':filters,'candidates':candidates,'candidate_sha256':digest(candidates),'decisions':decisions,'located_profiles':located_profiles,'native_controls':native_controls,'old_scale_sha256':oldpin,'new_scales':newset,'new_scale_sha256':digest(newset)}
    dest=here.parent/'05-knowledge/results' if here.name=='04-computation' else here
    output=dest/(Path(__file__).stem+'_certificate.json')
    output.write_bytes((json.dumps(record,sort_keys=True,indent=2)+'\n').encode())
    print('CANDIDATE PROOF: every retained connected-seven-complement row atclock16704 is weakly safe; independent audit pending')
    print('WORDS complete346104 valid19073 conditional_margin_types161/171 globalmax188 valid_sha',digest(valid_words))
    print('EDGES atlas5855 times5 divisors; filters',json.dumps(filters,sort_keys=True),'remaining988; Csep835 located153; no survivors')
    print('MINIMUM_LOCATED_SURPLUS',json.dumps(worst,sort_keys=True))
    print('MINIMUM_ALL_SURPLUS',min(row[-1] for row in decisions),'NATIVE_CONTROLS',native_controls)
    print('SCALE_SET deletion16704 count8201 maximum14904 sha256',digest(newset))
    print('CERTIFICATE_SHA256',hashlib.sha256(output.read_bytes()).hexdigest())
    print('SCOPE one new clock only; general chosen-six connected-seven complement, weak safety, no LRC14 closure')
    print('PASS',G,'always-active exact gates; raw LF')

if __name__=='__main__':main()
