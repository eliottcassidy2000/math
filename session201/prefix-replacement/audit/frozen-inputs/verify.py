"""Constructive maximal-prefix signed repair. Python standard library only."""
from pathlib import Path
from itertools import product
from functools import lru_cache
import argparse,hashlib,json
T=(1,0,2,1)
def need(ok,info):
    if not ok:raise RuntimeError(info)
def step(s,c):
    D,E,U,d=s
    if c=='O':return D,0,2*U,d
    if c=='Z':return D-E+U,d,2*U-d,d
    if c=='J':return D-E+U,D-E,D+U,D
    raise ValueError(c)
def state(word,s=T):
    for c in word:s=step(s,c)
    return s
def f(s):return s[2]-s[0]+s[1]-s[3]
def signed(s):return 6*s[2]-7*s[0]
def section(s):return s[1]==2*(s[2]-s[0]) and s[3]==3*(s[2]-s[0])
def bank(seed,cap,Dcap=None):
    if seed[2]>cap or Dcap is not None and seed[0]>Dcap:return {}
    out={seed:''};queue=[seed]
    for s in queue:
        for c in 'OZJ':
            q=step(s,c)
            need(q[2]>s[2] and q[0]>=s[0],'physical growth')
            if q[2]<=cap and (Dcap is None or q[0]<=Dcap) and q not in out:
                out[q]=out[s]+c;queue.append(q)
    return out

class TargetTable:
    def __init__(self,n):
        need(n>0 and n%2==1,'positive odd target required')
        self.n=n;self.cap=2*n-3
        self.nodes=bank(T,self.cap,n) if n>1 else {}
        self.good={}
        for s in sorted(self.nodes,key=lambda x:x[2],reverse=True):
            if s[0]==n:self.good[s]='';continue
            choices=[c+self.good[q] for c in 'OZJ' if (q:=step(s,c)) in self.good]
            if choices:self.good[s]=min(choices,key=lambda w:(len(w),w))
    def D_suffix(self,s):
        if self.n<s[0]:return None
        if self.n==s[0]:return ''
        if s[2]>self.cap:return None
        return self.good.get(s)
    def completion(self,s):
        options=[];tail=self.D_suffix(s)
        if tail is not None:options.append((tail+'OJJJJ','descendant D'))
        if f(s)==0 and s[3]==self.n:options.append(('ZJJJ','deletion'))
        for j in range(3):
            word='J'*j;q=state(word,s)
            if section(q) and signed(q)==self.n:options.append((word,'short completion'))
        return min(options,key=lambda x:(len(x[0]),x[0])) if options else None

@lru_cache(maxsize=None)
def table(n):return TargetTable(n)

def repair(word,n=None):
    need(all(c in 'OZJ' for c in word),'physical alphabet')
    prefixes=[T]
    for c in word:prefixes.append(step(prefixes[-1],c))
    if n is None:n=signed(prefixes[-1])
    tab=table(n);queries=1
    if tab.completion(T) is None:
        return dict(status='NO_C1_SECTION',target=n,original_word=word,target_table_profiles=len(tab.nodes),queries=queries)
    lo,hi=0,len(word)+1
    while hi-lo>1:
        mid=(lo+hi)//2;queries+=1
        if tab.completion(prefixes[mid]) is None:hi=mid
        else:lo=mid
    ext,route=tab.completion(prefixes[lo]);newword=word[:lo]+ext;q=state(newword)
    need(section(q) and signed(q)==n,'constructed section witness')
    unsigned=tab.D_suffix(T)
    need(unsigned is not None and state(unsigned)[0]==n,'unsigned source witness')
    c1word=None
    if n>1:
        core=unsigned.rstrip('O');need(core and core[-1] in 'ZJ','positive-count birth')
        c1word=core[:-1]+'J';need(state(c1word)[0]==n,'C1 reconstruction count')
    return dict(status='REPAIRED',target=n,original_word=word,retained_prefix_length=lo,
        removed_letters=len(word)-lo,retained_prefix=word[:lo],extension=ext,compiler_route=route,
        repaired_word=newword,repaired_profile=q,unsigned_source_word=unsigned,C1_word=c1word,
        target_table_profiles=len(tab.nodes),queries=queries)

def graph(word):
    rows=[2,0];root=0
    for c in word:
        n=len(rows);rows.append(0)
        for v in range(n):
            if c!='O' and v==root:rows[n]|=1<<v
            else:rows[v]|=1<<n
        if c=='J':root=n
    return rows,root

def literal(rows,root):
    n=len(rows);full=(1<<n)-1
    def dp(starts):
        out=[{} for _ in range(full+1)]
        for v in starts:out[1<<v][v]=1
        for mask in range(1,full+1):
            for v,c in out[mask].items():
                bits=rows[v]&~mask&full
                while bits:
                    bit=bits&-bits;bits-=bit;j=bit.bit_length()-1
                    dest=out[mask|bit];dest[j]=dest.get(j,0)+c
        return out
    A=dp(range(n));P=dp([root]);H=[sum(x.values())for x in A];H[0]=1
    starts=[sum(x.values())for x in P]
    U=sum(starts[m]*H[full^m]for m in range(1,full+1)if m>>root&1)
    prof=(H[full],A[full].get(root,0),U,H[full^(1<<root)])
    cycles=sum(c for v,c in P[full].items()if rows[v]>>root&1)
    return dict(profile=prof,root_start=starts[full],cycles=cycles,vertices=n)

def tests(max_length,target_cap):
    need(max_length>=6 and target_cap>=301,'controls included')
    full=bank(T,2*target_cap-3);Dvalues={s[0]for s in full};Dvalues.add(1)
    checked=success=failed=0;sample=[]
    for length in range(max_length+1):
        for letters in product('OZJ',repeat=length):
            word=''.join(letters);n=signed(state(word))
            if not 0<n<=target_cap:continue
            result=repair(word);tab=table(n)
            prefixes=[T]
            for c in word:prefixes.append(step(prefixes[-1],c))
            flags=[tab.completion(s) is not None for s in prefixes]
            need(all(not flags[i+1] or flags[i] for i in range(length)),'initial interval')
            need((result['status']=='REPAIRED')==(n in Dvalues),'complete root membership')
            if result['status']=='REPAIRED':
                expected=max(i for i,x in enumerate(flags)if x)
                need(result['retained_prefix_length']==expected,'binary-search maximal cut')
                need(result['queries']<=1+(length+1).bit_length(),'logarithmic queries')
                success+=1
            else:
                need(not any(flags),'no later prefix can repair');failed+=1
            checked+=1
            if word in ['ZOZ','ZOOZ','OZZ','OZJJ','OZ']:sample.append(result)
    b=state('Z');base=step(b,'O');inventory=bank(base,63)
    values=sorted({s[0]for s in inventory if s[0]<=33})
    need(len(inventory)==55 and values==[3,9,15,19,27,31],'complete 33 tail exclusion')
    short_words=0
    for length in range(7):
        for letters in product('OZJ',repeat=length):
            q=state(''.join(letters),b)
            need(not(section(q) and signed(q)==33),'minimum seven-letter repair')
            short_words+=1
    need(short_words==1093,'complete short-word count')
    for r in range(41):
        scale=3**r;vword='Z'+'OJ'*r;v=state(vword)
        need(v==tuple(scale*x for x in b),'physical scaled boundary')
        old1=state('OZ',v);q1=state('OJJ',v)
        old2=state('OOZ',v);q2=state('ZZZZOJJ',v)
        need(q1==tuple(scale*x for x in (15,6,18,9)),'first repair identity')
        need(q2==tuple(scale*x for x in (165,66,198,99)),'second repair identity')
        need(section(q1) and signed(q1)==signed(old1)==3*scale,'first target preserved')
        need(section(q2) and signed(q2)==signed(old2)==33*scale,'second target preserved')
    literal_results=[]
    for word in ['ZOJJ','ZZZZZOJJ','Z','ZZZZZ']:
        data=literal(*graph(word));need(data['profile']==state(word),'literal profile equality')
        need(data['root_start']==data['cycles']==1,'literal start-one and C1')
        if word=='Z':need(data['profile'][0]==3,'ordinary target3')
        if word=='ZZZZZ':need(data['profile'][0]==33,'ordinary target33')
        literal_results.append(dict(word=word,**data))
    return dict(status='PASS',word_length_bound=max_length,target_cap=target_cap,
        complete_global_membership_bank=len(full),tested_signed_words=checked,successful_repairs=success,
        exact_no_C1_controls=failed,distinct_target_tables=table.cache_info().currsize,
        total_target_table_profiles=sum(table(n).nodes.__len__()for n in range(1,target_cap+1,2)),
        prefix33_inventory=len(inventory),prefix33_D_values=values,minimum_length_word_checks=short_words,
        infinite_repair_formula_r=[0,40],literal_graphs=literal_results,sample_repairs=sample,
        verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        dependencies={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted((Path(__file__).parent/'dependencies').glob('*.md'))},
        scope='Constructive maximal-prefix replacement conditional on exact C1 membership; no universal occupancy claim.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--max-length',type=int,default=6);p.add_argument('--target-cap',type=int,default=501)
    p.add_argument('--repair-word');p.add_argument('--target',type=int);a=p.parse_args()
    out=repair(a.repair_word,a.target)if a.repair_word is not None else tests(a.max_length,a.target_cap)
    print(json.dumps(out,indent=2,sort_keys=True))
