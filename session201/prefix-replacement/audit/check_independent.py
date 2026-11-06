"""Independent session 201 audit: positive coordinates and direct graph DP.

No author arithmetic, inventory, or graph-counting functions are reused.
Author module is imported only as the system under test after independent
reference results have been produced.
"""
from collections import deque
from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib
import importlib.util
import json
import sys

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
PACKET = HERE / 'frozen-inputs'
EXPECTED = {
    'PROOF.md': '3ac0eefd2815134456101cd6ba7841b756bc13667ab6e149417b529b7c2281d9',
    'verify.py': '2868922f90f5ba2cc8d6122408d87df692bfd64d1afa5747b814ca4f070af8dc',
    'VALIDATION-normal.json': 'bb10aecd4f968a2868737cd88525d0c5b374ff1119e7e0e21e3fe4197dec7da1',
    'VALIDATION-optimized.json': 'bb10aecd4f968a2868737cd88525d0c5b374ff1119e7e0e21e3fe4197dec7da1',
}

def require(ok, label):
    if not ok:
        raise RuntimeError(label)

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

# All independent arithmetic uses y=(a,f,e,p), unlike the author's raw states.
T = (0,0,0,1)
def raw(y):
    a,f,e,p = y
    return (a+2*e+p, e, a+f+2*e+2*p, e+p)

def move(y,c):
    a,f,e,p = y
    if c == 'O': return (a+e, a+2*f+e+2*p, 0, e+p)
    if c == 'Z': return (2*a+f+e+p, f, e+p, 0)
    if c == 'J': return (f+p, 0, a+e+p, e)
    raise ValueError(c)

def apply(w,y=T):
    for c in w: y=move(y,c)
    return y

def height(y): return raw(y)[2]
def target(y):
    D,E,U,d=raw(y)
    return 6*U-7*D

def section(y):
    # Slice written in positive coordinates; unlike author's raw equations.
    a,f,e,p=y
    return f==0 and e==2*p

def bank(y,cap,d_cap=None):
    if height(y)>cap or d_cap is not None and raw(y)[0]>d_cap: return {}
    states={y:''}; q=deque([y])
    while q:
        s=q.popleft()
        for c in 'JZO':
            v=move(s,c)
            require(height(v)>height(s) and raw(v)[0]>=raw(s)[0], 'monotone physical edges')
            if height(v)<=cap and (d_cap is None or raw(v)[0]<=d_cap) and v not in states:
                states[v]=states[s]+c;q.append(v)
    return states

@lru_cache(None)
def descendant_targets(y,cap):
    return {raw(v)[0] for v in bank(y,2*cap-3) if raw(v)[0]<=cap} | {raw(y)[0]}

def feasible(y,n,cap):
    if n in descendant_targets(y,cap): return True
    if y[1]==0 and raw(y)[3]==n: return True
    return any(section(apply('J'*j,y)) and target(apply('J'*j,y))==n for j in range(3))

def word_list(maxlen):
    for n in range(maxlen+1):
        for w in product('OZJ',repeat=n): yield ''.join(w)

# Independent literal graph realization, using an adjacency matrix, not masks.
def make_graph(w):
    g=[[False,True],[False,False]]; root=0
    for c in w:
        k=len(g)
        new_beats_old=[c!='O' and v==root for v in range(k)]
        for v in range(k):g[v].append(not new_beats_old[v])
        g.append(new_beats_old+[False])
        if c=='J':root=k
    return g,root

def count_graph(w):
    g,r=make_graph(w); n=len(g);full=(1<<n)-1
    @lru_cache(None)
    def ending(S,v):
        smaller=S^(1<<v)
        if not smaller:return 1
        return sum(ending(smaller,u) for u in range(n) if smaller>>u&1 and g[u][v])
    @lru_cache(None)
    def from_root(S,v):
        if not S>>r&1:return 0
        smaller=S^(1<<v)
        if not smaller:return int(v==r)
        return sum(from_root(smaller,u) for u in range(n) if smaller>>u&1 and g[u][v])
    @lru_cache(None)
    def paths(S):
        return sum(ending(S,v) for v in range(n) if S>>v&1) if S else 1
    D=paths(full);E=ending(full,r);d=paths(full^(1<<r))
    U=sum(paths(full^S)*sum(from_root(S,v) for v in range(n) if S>>v&1)
          for S in range(1,full+1) if S>>r&1)
    starts=sum(from_root(full,v) for v in range(n))
    cycles=sum(from_root(full,v) for v in range(n) if g[v][r])
    profile=(D,E,U,d)
    require(profile==raw(apply(w)),('literal profile',w,profile,raw(apply(w))))
    require(starts==1,('literal root-start',w))
    return dict(word=w,order=n,profile=profile,root_start=starts,cycles=cycles)

def main():
    for name,h in EXPECTED.items():require(digest(PACKET/name)==h,('frozen input',name))
    b=apply('Z');bo=apply('ZO')
    inventory=bank(bo,63)
    require(len(inventory)==55, '55 U-only states')
    filtered={s:w for s,w in inventory.items() if raw(s)[0]<=33}
    require(len(filtered)==34,'34 U,D-filtered states')
    vals=sorted({raw(s)[0] for s in filtered})
    require(vals==[3,9,15,19,27,31], 'exact target33 exclusion')
    short=list(word_list(6))
    require(len(short)==1093, '1093 short words')
    require(not any(section(apply(w,b)) and target(apply(w,b))==33 for w in short), 'minimum7')
    length7=[w for w in word_list(7) if len(w)==7 and section(apply(w,b)) and target(apply(w,b))==33]
    require('ZZZZOJJ' in length7, 'claimed length7 witness')
    require(not any(section(apply(w,bo)) for w in word_list(1)), 'minimum2')
    for r in range(101):
        v=apply('Z'+'OJ'*r);scale=3**r
        require(raw(v)==tuple(scale*x for x in (3,1,3,1)), 'scaled physical prefix')
        for old,new,n,profile in [('OZ','OJJ',3,(15,6,18,9)),('OOZ','ZZZZOJJ',33,(165,66,198,99))]:
            q=apply(new,v)
            require(raw(q)==tuple(scale*x for x in profile), 'repair scaled identity')
            require(section(q) and target(q)==target(apply(old,v))==n*scale,'signed repair')
    spec=importlib.util.spec_from_file_location('author_sut',PACKET/'verify.py')
    author=importlib.util.module_from_spec(spec);spec.loader.exec_module(author)
    query_sources=bank(T,127);checks=0;positive=0
    for y in query_sources:
        expected=descendant_targets(y,63)
        for n in range(1,64,2):
            tail=author.table(n).D_suffix(raw(y))
            require((tail is not None)==(n in expected),('Dtail query',raw(y),n,tail))
            if tail is not None:
                require(raw(apply(tail,y))[0]==n,'Dtail reconstruction');positive+=1
            checks+=1
    # Equality remains valid arbitrarily far beyond the target table cutoff.
    for base in ['','Z','ZZZZZ','ZOJJ']:
        y=apply(base+'O'*80);n=raw(y)[0]
        require(author.table(n).D_suffix(raw(y))=='','large-height equality branch')
    # Direct forward section enumeration under the inherited complete bound,
    # not a table/classification comparison with itself.
    forward_sources=bank(T,31);forward_checks=0;forward_states=0
    for y in forward_sources:
        cap=max(34*31-42,14*height(y)+6*31)
        reachable=bank(y,cap);forward_states+=len(reachable)
        targets={target(v) for v in reachable if section(v)}
        for n in range(1,32,2):
            result=author.table(n).completion(raw(y))
            require((result is not None)==(n in targets),('direct section comparison',raw(y),n))
            require(feasible(y,n,31)==(n in targets),'independent exact classification')
            if result is not None:
                out=apply(result[0],y)
                require(section(out) and target(out)==n,'section completion witness')
            forward_checks+=1
    # Exhaust every original word through length7, using independent flags.
    repairs=success=failure=0;total_prefixes=0
    for w in word_list(7):
        y=apply(w);n=target(y)
        if not 0<n<=301:continue
        flags=[feasible(apply(w[:i]),n,301) for i in range(len(w)+1)]
        require(flags==sorted(flags,reverse=True),'initial interval')
        answer=author.repair(w);repairs+=1;total_prefixes+=len(flags)
        if any(flags):
            success+=1;cut=max(i for i,x in enumerate(flags) if x)
            require(answer['status']=='REPAIRED' and answer['retained_prefix_length']==cut,'maximal cut')
            q=apply(answer['repaired_word'])
            require(section(q) and target(q)==n,'repair witness')
            require(raw(apply(answer['unsigned_source_word']))[0]==n,'unsigned witness')
            if n>1:
                require(answer['C1_word'].endswith('J') and raw(apply(answer['C1_word']))[0]==n,'C1 recipe')
        else:
            failure+=1;require(answer['status']=='NO_C1_SECTION','no C1 exact control')
    require(author.repair('OZ')['status']=='NO_C1_SECTION','signed7 hostile')
    # Explicit requested targets, including empty words, already-section
    # prefixes, no-op repairs, and root-only repairs.
    explicit=[]
    for w in ['', 'OZ', 'ZOZ', 'ZOOZ', 'ZOJJ', 'ZZZZZOJJ', 'OOOOOOOOOO']:
        for n in [1,3,7,9,21,33,63]:
            flags=[feasible(apply(w[:i]),n,63) for i in range(len(w)+1)]
            a=author.repair(w,n)
            require((a['status']=='REPAIRED')==any(flags),'explicit target feasibility')
            if any(flags):require(a['retained_prefix_length']==max(i for i,x in enumerate(flags) if x),'explicit target cut')
            explicit.append({'word':w,'target':n,'status':a['status'],'cut':a.get('retained_prefix_length')})
    graph_words=['ZOZ','ZOOZ','ZOJJ','ZZZZZOJJ','Z','ZZZZZ','J','ZZZZJ',
                 'ZOJOJJ','ZOJZZZZOJJ','ZOJ','ZOJZZZZ', 'JZOZJ']
    graphs=[count_graph(w) for w in graph_words]
    for x in graphs:
        if x['word'] in ['ZOJJ','ZZZZZOJJ','Z','ZZZZZ','J','ZZZZJ','JZOZJ']:
            require(x['cycles']==1,('literal C1',x['word']))
    for name,h in EXPECTED.items():require(digest(PACKET/name)==h,('unchanged frozen input',name))
    return dict(status='PASS',proof_sha256=EXPECTED['PROOF.md'],author_sha256=EXPECTED['verify.py'],
        verifier_sha256=digest(Path(__file__)),
        dependencies={p.name:digest(p) for p in sorted((PACKET/'dependencies').glob('*.md'))},
        inventory_U63_profiles=len(inventory),inventory_U63_D33_profiles=len(filtered),inventory_D_values=vals,
        inventory=[{'raw':raw(s),'suffix_from_bO':w} for s,w in sorted(inventory.items(),key=lambda item:raw(item[0]))],
        short_words_checked=len(short),minimal_length7_witnesses=length7,infinite_identity_checked_r=[0,100],
        Dtail_query_sources=len(query_sources),Dtail_target_cap=63,Dtail_checks=checks,Dtail_successes=positive,
        forward_section_sources=len(forward_sources),forward_section_target_cap=31,
        forward_section_profile_visits=forward_states,forward_section_queries=forward_checks,
        exhaustive_original_word_length=7,exhaustive_target_cap=301,repair_cases=repairs,
        successful_repairs=success,no_C1_cases=failure,checked_prefixes=total_prefixes,
        explicit_target_cases=explicit,literal_graphs=graphs,
        algorithm_sample=author.repair('ZOOZ'),
        scope='Exact C1-conditional maximal-prefix repair, not universal C1/C2/C3 occupancy; family suffix length is a separate optimization.')

if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True))
