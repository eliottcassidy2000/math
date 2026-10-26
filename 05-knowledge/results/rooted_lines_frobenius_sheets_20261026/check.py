"""Exact marked-line composition; literal enumeration versus rooted convolution."""
from collections import Counter, defaultdict
from functools import lru_cache
import itertools as it
import json

def need(b,s):
    if not b: raise RuntimeError(s)

def edge(t,i,j): return t[i]>>j&1
def opposite(t): return tuple(((1<<len(t))-1)^x^(1<<i) for i,x in enumerate(t))
def reindex(t,p): return tuple(sum(edge(t,p[i],p[j])<<j for j in range(len(t))) for i in range(len(t)))
def sigma(t): return reindex(opposite(t),tuple(reversed(range(len(t)))))
def tau(t):
    a=list(opposite(t))
    for i in range(len(t)-1):a[i]^=1<<(i+1);a[i+1]^=1<<i
    return tuple(a)
def presented(n,mask):
    rows=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            forward=j==i+1 or mask>>k&1
            if j>i+1:k+=1
            if forward: rows[i]|=1<<j
            else: rows[j]|=1<<i
    return tuple(rows)

@lru_cache(None)
def canonical(t,roots=()):
    fixed=list(dict.fromkeys(roots)); pat=tuple(fixed.index(x) for x in roots)
    groups=defaultdict(list)
    for i in range(len(t)):
        if i not in fixed: groups[t[i].bit_count()].append(i)
    perms=[list(it.permutations(groups[k])) for k in sorted(groups)]
    rows=min(reindex(t,tuple(fixed)+tuple(it.chain.from_iterable(parts))) for parts in it.product(*perms))
    return rows,pat

def root_sig(t): return canonical(t,(0,len(t)-1))
def reflect_sig(h):
    t,(s,e)=h
    return canonical(opposite(t),(e,s))
def dual_class(a):return canonical(opposite(a[0]))

@lru_cache(None)
def glue(h,k):
    a,(s,t)=h;b,(u,v)=k;n=len(a);m=len(b)
    rows=list(a)+[(x<<n)|((1<<n)-1) for x in b]
    rows[t]|=1<<(n+u);rows[n+u]^=1<<t
    return canonical(tuple(rows),(s,n+v))

def fold(seq):
    h=seq[0]
    for k in seq[1:]:h=glue(h,k)
    return h

def components(t):
    n=len(t)
    cuts=[0]+[k for k in range(1,n) if all(edge(t,i,j) for i in range(k) for j in range(k,n))]+[n]
    return [tuple((t[i]>>a)&((1<<(b-a))-1) for i in range(a,b)) for a,b in zip(cuts,cuts[1:])]

def paths(t):
    def visit(p,seen):
        if len(p)==len(t):yield tuple(p);return
        for v in range(len(t)):
            if not seen>>v&1 and (not p or edge(t,p[-1],v)):
                yield from visit(p+[v],seen|(1<<v))
    yield from visit([],0)

def hamilton(t):
    @lru_cache(None)
    def f(mask,last):
        if mask==(1<<len(t))-1:return 1
        return sum(f(mask|(1<<j),j) for j in range(len(t)) if not mask>>j&1 and edge(t,last,j))
    return sum(f(1<<i,i) for i in range(len(t)))

def convolution_check():
    F=defaultdict(Counter);B=defaultdict(Counter);observed=defaultdict(Counter);blue=defaultdict(Counter)
    candidates=defaultdict(list);totals=[]
    for n in range(1,7):
        counts=Counter()
        for mask in range(1<<((n-1)*(n-2)//2)):
            x=presented(n,mask);y=tau(x);h=root_sig(y);atoms=components(x)
            word=tuple(canonical(a) for a in atoms);isblue=(sigma(x)==x)
            observed[word][h]+=1
            if isblue:blue[word][h]+=1
            if len(atoms)==1:
                F[word[0]][h]+=1
                if isblue:B[word[0]][h]+=1
                key=(word[0],canonical(y),canonical(y,(0,)),canonical(y,(n-1,)),isblue)
                candidates[key].append((h,mask,x,y))
            counts['presentations']+=1;counts['blue_presentations']+=isblue
            counts['join_presentations']+=len(atoms)>1
        totals.append({'n':n,**dict(counts)})
    checked=0;blue_words=0
    for word,actual in observed.items():
        if len(word)<2:continue
        pred=Counter()
        for terms in it.product(*(list(F[a].items()) for a in word)):
            hs,weights=zip(*terms);prod=1
            for weight in weights:prod*=weight
            pred[fold(hs)]+=prod
        need(pred==actual,'two-port total convolution')
        bp=Counter()
        if word==tuple(dual_class(a) for a in reversed(word)):
            half=len(word)//2
            choices=[list(F[a].items()) for a in word[:half]]
            if len(word)%2:choices.append(list(B[word[half]].items()))
            for terms in it.product(*choices):
                hs=[h for h,_ in terms];prod=1
                for _,weight in terms:prod*=weight
                left=hs[:half];middle=hs[half:]
                sequence=left+middle+[reflect_sig(h) for h in reversed(left)]
                bp[fold(sequence)]+=prod
        need(bp==blue[word],'reflection-fold blue convolution')
        need(all(bp[h]<=pred[h] for h in bp),'black = total minus blue')
        checked+=1;blue_words+=bool(bp)
    # Same source, target, both separate rooted types and colour; different joint ports.
    hostile=None;unit=(((0,), (0,0)))
    for key,items in candidates.items():
        if hostile:break
        for a,b in it.combinations(items,2):
            if a[0]==b[0]:continue
            ga=fold([unit,a[0],unit]);gb=fold([unit,b[0],unit])
            ha,hb=hamilton(ga[0]),hamilton(gb[0])
            if ha!=hb:
                hostile={'factor_n':len(a[2]),'masks':[a[1],b[1]],
                         'source_rows':[list(a[2]),list(b[2])],
                         'target_rows':[list(a[3]),list(b[3])],
                         'same_separate_root_types':True,'same_source_target_and_colour':True,
                         'singleton_sandwich_path_counts':[ha,hb],
                         'joined_target_rows':[list(ga[0]),list(gb[0])]}
                break
    need(hostile is not None,'find a joint-port hostile')
    # Product of opposite reflection-paired factors: blue uses diagonal weights.
    diagonal=None
    for a,hist in F.items():
        if any(v>1 for v in hist.values()) and 2*len(a[0])<=12:
            diagonal={'atom_n':len(a[0]),'source_rows':list(a[0]),
                      'correct_blue_total':sum(hist.values()),
                      'incorrect_independent_matching':sum(v*v for v in hist.values())}
            break
    return {'universe':totals,'join_class_words_checked':checked,'words_with_blue':blue_words,
            'joint_port_hostile':hostile,'reflection_diagonal_control':diagonal}

def substitute(q,bs):
    starts=[0]
    for b in bs:starts.append(starts[-1]+len(b))
    rows=[];owner=[]
    for i,b in enumerate(bs):
        ext=sum(((1<<len(bs[j]))-1)<<starts[j] for j in range(len(bs)) if edge(q,i,j))
        rows.extend((r<<starts[i])|ext for r in b);owner.extend([i]*len(b))
    return tuple(rows),owner

def forest_check():
    c3=(2,4,1);tt2=(2,0);ap=(10,12,9,0)
    cases=[(c3,[tt2,tt2,tt2]),(c3,[ap,(0,),(0,)]),(c3,[c3,c3,(0,)])]
    out=[]
    for q,bs in cases:
        t,owner=substitute(q,bs);count=0;maxruns=0;minimal=None
        for p in paths(t):
            runs=[]
            for v in p:
                if not runs or owner[runs[-1][-1]]!=owner[v]:runs.append([v])
                else:runs[-1].append(v)
            rr=Counter(owner[r[0]] for r in runs);maxruns=max(maxruns,max(rr.values()))
            direct=list(opposite(t))
            for a,b in zip(p,p[1:]):direct[a]^=1<<b;direct[b]^=1<<a
            rebuilt=[0]*len(t)
            # Build complemented child forests and reversed quotient cross arcs.
            internal={tuple(sorted((a,b))) for run in runs for a,b in zip(run,run[1:])}
            hinges={tuple(sorted((a[-1],b[0]))) for a,b in zip(runs,runs[1:])}
            for i,j in it.combinations(range(len(t)),2):
                old=edge(t,i,j)
                if owner[i]==owner[j]: forward=(1-old)^int((i,j) in internal)
                else: forward=int(edge(q,owner[j],owner[i]))^int((i,j) in hinges)
                if forward:rebuilt[i]|=1<<j
                else:rebuilt[j]|=1<<i
            need(tuple(direct)==tuple(rebuilt),'run forest reconstructs literal partner')
            if max(rr.values())>1 and minimal is None:
                minimal={'path':list(p),'run_word':[owner[r[0]] for r in runs],
                         'runs':runs,'partner_rows':direct}
            count+=1
        out.append({'sizes':[len(b) for b in bs],'paths_checked':count,'max_runs_in_one_block':maxruns,
                    'reentry_witness':minimal})
    return out

def nonclosure():
    rows=[]
    for m in [3,5,7,9,11,13]:
        q=tuple(sum(1<<j for j in range(m) if 0<(j-i)%m<=m//2) for i in range(m))
        t,_=substitute(q,[(10,12,9,0)]*m)
        partner=tau(t);hist=Counter(r.bit_count() for r in partner)
        need(hist=={2*m-1:3*m-1,2*m:2,2*m+1:m-1},'unlabelled carousel nonclosure')
        need(all(hist.get(2*m-2,0)!=a or hist.get(2*m,0)!=3*a for a in range(m+1)), 'no child recolouring recovers partner')
        rows.append({'m':m,'target_degree_multiplicities':dict(sorted(hist.items()))})
    return rows

if __name__=='__main__':
    print(json.dumps({'composition':convolution_check(),'forest':forest_check(),'unlabelled_nonclosure':nonclosure()},sort_keys=True,indent=2))
