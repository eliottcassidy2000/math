"""Independent forward-path and boundary-cover audit of finite source tables."""
import itertools, json, math, pathlib, random, sys
from functools import lru_cache
ROOT=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'source'))
import verify as source

def induced(rows, vs):
    return tuple(sum(1<<j for j,y in enumerate(vs) if rows[x]&(1<<y)) for x in vs)

def endpoints(rows):
    n=len(rows); full=(1<<n)-1
    dp=[[0]*n for _ in range(full+1)]
    for v in range(n):dp[1<<v][v]=1
    for S in range(1,full+1):
        for v,x in enumerate(dp[S]):
            if not x:continue
            left=rows[v]&~S
            while left:
                bit=left&-left;left-=bit
                dp[S|bit][bit.bit_length()-1]+=x
    return dp

def data(rows):
    n=len(rows); full=(1<<n)-1
    end=endpoints(rows)
    converse=tuple(sum(1<<u for u in range(n) if rows[u]&(1<<v)) for v in range(n))
    start=endpoints(converse)
    H=[1]+[sum(end[S]) for S in range(1,full+1)]
    R=[0]*n;B=[[0]*n for _ in range(n)]
    for Q in range(1,full+1):
        P=full^Q
        for x,a in enumerate(start[Q]):
            if a:
                R[x]+=H[P]*a
                for y,b in enumerate(end[P]):
                    if b:B[y][x]+=b*a
    profiles=[]
    for r in range(n):
        profiles.append((H[full],end[full][r],R[r],H[full^(1<<r)]))
    return profiles,H,B,start,end

def components(rows):
    n=len(rows);reach=[]
    for s in range(n):
        seen={s};stack=[s]
        while stack:
            u=stack.pop()
            for v in range(n):
                if rows[u]&(1<<v) and v not in seen:seen.add(v);stack.append(v)
        reach.append(seen)
    todo=set(range(n));cs=[]
    while todo:
        u=min(todo);c={v for v in todo if v in reach[u] and u in reach[v]}
        todo-=c;cs.append(c)
    # Independent topological sort rather than component score sorting.
    order=[]
    while cs:
        c=next(c for c in cs if all(c is d or rows[min(c)]&(1<<min(d)) for d in cs))
        cs.remove(c);order.append(sorted(c))
    return order

def literal(rows,r):
    n=len(rows);h=e=R=0
    for p in itertools.permutations(range(n)):
        bad=[i for i in range(n-1) if not rows[p[i]]&(1<<p[i+1])]
        if not bad:h+=1;e+=p[-1]==r
        cut=p.index(r)
        R+=not bad or (len(bad)==1 and bad[0]==cut-1)
    S=induced(rows,[v for v in range(n) if v!=r])
    d=sum(all(S[p[i]]&(1<<p[i+1]) for i in range(len(S)-1)) for p in itertools.permutations(range(len(S))))
    return h,e,R,d

def normalize(rows,r):
    vs=[v for v in range(len(rows)) if v!=r];S=induced(rows,vs)
    O={v for v in range(len(S)) if rows[r]&(1<<vs[v])};cs=components(S)
    indices=[i for i,c in enumerate(cs) if O.intersection(c)]
    lo=min(indices) if indices else -1;hi=max(indices) if indices else -1
    keep=[];tail=0;middle=0
    for i,c in enumerate(cs):
        if len(c)>1 or O.intersection(c) or lo<i<hi:
            keep+=c;middle+=len(c)==1 and not O.intersection(c)
        elif indices and i>hi:tail+=1
    C=induced(rows,[vs[v] for v in keep]+[r])
    return C,len(C)-1,tail,middle

def check(rows,literal_check=False):
    pp,H,B,starts,ends=data(rows);n=len(rows)
    assert pp==source.profiles(list(rows))[0]
    for x in range(n):
        for y in range(x):assert B[x][y]+B[y][x]>=H[((1<<n)-1)^(1<<x)^(1<<y)]
    for r,z in enumerate(pp):
        if literal_check:assert z==literal(rows,r)
        h,e,R,d=z;f=R+e-h-d;U=R-d
        C,c,t,u=normalize(rows,r);cp=data(C)[0][c]
        ch,ce,cR,cd=cp;cf=cR+ce-ch-cd;cU=cR-cd
        predicted=cp if not t else ((cU<<t)+cd-cf,cd,cd+(cU<<t),cd)
        assert predicted==z
        assert d==cd and f==cf
        assert len(C)<=1+d+(1+math.isqrt(1+8*f))//2+(max(1,f).bit_length()-1)
        assert not u or (1<<u)<=f
        if U:assert e<=h-e and h<=2*(U-f)-1
        S=induced(rows,[v for v in range(n) if v!=r]);sp,sH,sB,ss,se=data(S)
        O=[i for i,v in enumerate(v for v in range(n) if v!=r) if rows[r]&(1<<v)]
        assert U==sum(sp[x][2] for x in O)
        assert f==sum(sB[y][x] for x in O for y in O if x!=y)
        assert f>=len(O)*(len(O)-1)//2
    return n

def deletion(core,before,after):
    n=before+len(core)+after
    rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
    for i in range(len(core)):
        for j in range(len(core)):
            rows[before+i]&=~(1<<(before+j))
            if core[i]&(1<<j):rows[before+i]|=1<<(before+j)
    return tuple(rows)

def table(core,d,beta):
    f=d-beta;k=(1+math.isqrt(1+8*f))//2;limit=k+max(1,f).bit_length()-1
    values={};examined=0
    for singles in range(limit+1):
        for before in (range(singles+1) if core else [0]):
            S=deletion(core,before,singles-before);n=len(S)
            profiles,H,B,starts,ends=data(S)
            assert H[-1]==d
            cs=components(S);belongs={v:i for i,c in enumerate(cs) for v in c}
            singleton=[c[0] for c in cs if len(c)==1]
            for out in range(1,1<<n):
                O=[v for v in range(n) if out&(1<<v)]
                if len(O)>k:continue
                first=min(belongs[v] for v in O);last=max(belongs[v] for v in O)
                if any(v not in O and not(first<belongs[v]<last) for v in singleton):continue
                examined+=1
                actual_f=sum(B[y][x] for x in O for y in O if x!=y)
                if actual_f!=f:continue
                U=sum(profiles[x][2] for x in O)
                e=d-sum(ends[-1][x] for x in O)
                T=tuple(S[v]|((1<<n) if v not in O else 0) for v in range(n))+(out,)
                z=(U+e-f,e,d+U,d)
                if U not in values:
                    assert data(T)[0][-1]==z
                    values[U]={'rows':T,'root':n,'profile':z}
    return {'examined':examined,'bases':sorted(values),'witnesses':values}

def main():
    checks=0;strong={3:set(),5:set()}
    for n in range(1,6):
        for rows in source.labelled(n):
            rows=tuple(rows);checks+=check(rows,True)
            if n>1 and len(components(rows))==1:
                h=data(rows)[1][-1]
                if h in strong:strong[h].add(min(induced(rows,p) for p in itertools.permutations(range(n))))
    report={'all_labelled_root_profile_and_normalization_checks_through_n5':checks,'strong_H3_H5_isomorphism_types':{str(k):sorted(v) for k,v in strong.items()}}
    print('Small profiles passed',report,flush=True)
    rng=random.Random(197555);random_checks=0
    for n,count in ((6,100),(7,50),(8,25),(9,10)):
        for _ in range(count):
            rows=[0]*n
            for i,j in itertools.combinations(range(n),2):
                a,b=(i,j) if rng.randrange(2) else (j,i);rows[a]|=1<<b
            random_checks+=check(tuple(rows))
    report['random_larger_root_checks']=random_checks
    expected={(1,1):[1],(1,-1):[3,5],(3,1):[6],(3,-1):[9,15],(5,1):[11,12,16],(5,-1):[15,39]}
    tables={}
    for core,d in (((),1),((4,1,2),3),((8,1,3,6),5)):
        for beta in (1,-1):
            r=table(core,d,beta);assert r['bases']==expected[d,beta]
            tables[f'{d}:{beta}']=r
            print('Independent boundary-cover table',d,beta,r['bases'],r['examined'],flush=True)
    report['independent_unit_tables']=tables
    E=(2,12,9,17,7);h13=[]
    for m in range(1,7):
        # New clone construction with clones first, unlike the source.
        labels=[0]*m+list(range(1,5));rows=[0]*len(labels)
        for i,j in itertools.combinations(range(len(labels)),2):
            if labels[i]==labels[j] or E[labels[i]]&(1<<labels[j]):rows[i]|=1<<j
            else:rows[j]|=1<<i
        h=data(tuple(rows))[1][-1];assert h==6*(1<<m)+1;h13.append(h)
        if m<=4:assert literal(tuple(rows),0)[0]==h
    report['independent_H13_clone_counts']=h13
    (ROOT/'AUDIT.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
