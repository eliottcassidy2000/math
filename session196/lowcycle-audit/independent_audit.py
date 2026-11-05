"""Independent audit: no imports from the proposed packet except comparison calls.
Uses sparse slot records with terminal slots restricted to E/G, independently
built partition-refinement, adjacency subset counts, and regex full matches.
"""
import itertools, collections, json, hashlib, random, sys, time, re
from pathlib import Path
BASE=Path(__file__).resolve().parent
SRC=BASE/'source'
sys.path.insert(0,str(SRC))
import apex_transfer as proposed

# E=0; G=1; B=2, but terminal slots only have E/G.
BOUNDS={'c':[(0,1),(1,0)],'ab':[(None,0),(0,1),(1,None)],'ba':[(None,1),(1,0),(0,None)]}
EMPTY={'c':(True,False),'ab':(True,True,True),'ba':(True,False,True)}
def start(kind):return (((0,)*len(BOUNDS[kind]),1),)
def advance(v,kind,letter,cap=None):
    accum=collections.defaultdict(int)
    for old,count in v:
        for slot,(left,right) in enumerate(BOUNDS[kind]):
            if old[slot]==0 and left is not None and not (letter&(1<<left)):continue
            new=list(old);new[slot]=1 if right is None or not (letter&(1<<right)) else 2
            accum[tuple(new)]+=count
    return tuple(sorted((s,min(cap,n) if cap is not None else n) for s,n in accum.items()))
def end(v,kind,cap=None):
    val=sum(n for s,n in v if all(t==1 or (t==0 and EMPTY[kind][i]) for i,t in enumerate(s)))
    return min(val,cap) if cap is not None else val

def slot_counts(word):
    vs={k:start(k) for k in BOUNDS}
    for a in word:vs={k:advance(v,k,a) for k,v in vs.items()}
    return end(vs['ab'],'ab')+end(vs['ba'],'ba'),end(vs['c'],'c')

def graph_counts(word):
    # Vertex 0=a, 1=b, 2+i=x_i; independently construct a matrix.
    n=len(word)+2
    arcs=[[False]*n for _ in range(n)];arcs[0][1]=True
    for j,a in enumerate(word,2):
        for root in (0,1):
            if (a>>root)&1:arcs[root][j]=True
            else:arcs[j][root]=True
        for i in range(2,j):arcs[i][j]=True
    # Recurrence strips the last vertex; first component all starts, second a only.
    D={1<<v:{v:(1,int(v==0))} for v in range(n)}
    for mask in range(1,1<<n):
        if mask.bit_count()==1:continue
        row={}
        for v in range(n):
            if mask&(1<<v):
                prev=D.get(mask^(1<<v),{})
                row[v]=(sum(x[0] for u,x in prev.items() if arcs[u][v]),sum(x[1] for u,x in prev.items() if arcs[u][v]))
        D[mask]=row
    row=D[(1<<n)-1]
    return sum(x[0] for x in row.values()),sum(x[1] for v,x in row.items() if arcs[v][0])

def automaton():
    vs=[start('c')];ids={vs[0]:0};T=[]
    for v in vs:
        row=[]
        for a in range(4):
            w=advance(v,'c',a,3)
            if w not in ids:ids[w]=len(vs);vs.append(w)
            row.append(ids[w])
        T.append(row)
    O=[end(v,'c',3) for v in vs]
    # Refine sets explicitly; do not assume proposed numbering.
    blocks=[{i for i,o in enumerate(O) if o==c} for c in sorted(set(O))]
    while True:
        label={s:b for b,B in enumerate(blocks) for s in B};new=[]
        for B in blocks:
            buckets=collections.defaultdict(set)
            for s in B:buckets[tuple(label[t] for t in T[s])].add(s)
            new.extend(buckets.values())
        if len(new)==len(blocks):break
        blocks=new
    label={s:b for b,B in enumerate(blocks) for s in B}
    t=[[label[w] for w in T[min(B)]] for B in blocks];o=[O[min(B)] for B in blocks]
    return vs,T,O,label[0],t,o

def bounded(B,init,T,O):
    live={s for s,o in enumerate(O) if o in (1,2)}
    while True:
        larger=live|{s for s,ts in enumerate(T) if any(t in live for t in ts)}
        if larger==live:break
        live=larger
    states=[(init,start('ab'),start('ba'))];index={states[0]:0};parents=[None];depth=[0];edges=0;hits={};h_by_c=collections.defaultdict(set)
    max_coordinate=0
    for i,(s,ab,ba) in enumerate(states):
        h=end(ab,'ab')+end(ba,'ba')
        assert h<=B
        if O[s] in (1,2):
            hits.setdefault(h,i);h_by_c[O[s]].add(h)
        for a,t in enumerate(T[s]):
            if t not in live:continue
            nab=advance(ab,'ab',a,B+1);nba=advance(ba,'ba',a,B+1)
            max_coordinate=max(max_coordinate,*(n for _,n in nab),*(n for _,n in nba))
            if end(nab,'ab')+end(nba,'ba')>B:continue
            state=(t,nab,nba);edges+=1
            if state not in index:
                index[state]=len(states);states.append(state);parents.append((i,a));depth.append(depth[i]+1)
    witnesses={}
    for h,i in hits.items():
        word=[]
        while parents[i] is not None:i,a=parents[i];word.append(a)
        word=word[::-1]
        assert slot_counts(word)[0]==h and slot_counts(word)[1] in (1,2)
        witnesses[h]=word
    return {'bound':B,'states':len(states),'transitions':edges,'max_discovery_depth':max(depth),'max_coordinate_before_pruning':max_coordinate,'missing_strong_odd':[h for h in range(3,B+1,2) if h not in hits],'H_by_C':{str(c):sorted(hs) for c,hs in h_by_c.items()},'witnesses':witnesses}

def main():
    result={}
    start_time=time.time();tested=0
    for m in range(7):
        for word in itertools.product(range(4),repeat=m):
            a=graph_counts(word);b=slot_counts(word);c=proposed.counts(word)
            assert a==b==c,(word,a,b,c)
            tested+=1
        print('independent adjacency check through base',m,'total',tested,'elapsed',time.time()-start_time,flush=True)
    result['graph_cross_checks']=tested
    rng=random.Random(196)
    for m in (7,8,9,10):
        for _ in range(25):
            word=tuple(rng.randrange(4) for _ in range(m));assert graph_counts(word)==slot_counts(word)==proposed.counts(word)
    result['additional_random_graph_checks']=100
    vs,T,O,init,t,o=automaton();assert (len(vs),len(t))==(251,32)
    result['raw_cycle_states']=len(vs);result['minimal_moore_states']=len(t)
    # Stored machine equivalence, which proposed verify_language does not read.
    saved=json.loads((SRC/'cycle_dfa.json').read_text());todo=[(init,saved['initial_state'])];seen=set(todo)
    for u,v in todo:
        assert o[u]==saved['outputs'][v]
        for a in range(4):
            pair=(t[u][a],saved['transitions'][v][a])
            if pair not in seen:seen.add(pair);todo.append(pair)
    assert len({u for u,v in seen})==32 and len({v for u,v in seen})==32
    result['saved_DFA_equivalence_states']=len(seen)
    for B in (67,501):
        b=bounded(B,init,t,o)
        saved=json.loads((SRC/f'bounded_image_{B}.json').read_text())
        for key in ('bound','states','transitions','max_discovery_depth','missing_strong_odd'):assert b[key]==saved[key],(key,b[key],saved[key])
        assert set(map(int,saved['witnesses']))==set(b['witnesses'])
        for h,w in saved['witnesses'].items():assert slot_counts(w)[0]==int(h) and slot_counts(w)[1] in (1,2)
        result[f'bounded_{B}']=b
        print('independent bounded result',B,b['states'],b['transitions'],b['missing_strong_odd'],flush=True)
    # Exact and saturated transfers agree for hostile large multiplicities.
    for kind,pkind in [('c','cycle'),('ab','path_ab'),('ba','path_ba')]:
        for cap in (1,2,3,68):
            for _ in range(100):
                word=[rng.randrange(4) for _ in range(40)];v=start(kind);u=start(kind)
                for a in word:v=advance(v,kind,a);u=advance(u,kind,a,cap)
                assert u==tuple((s,min(n,cap)) for s,n in v)
                assert end(u,kind,cap)==min(cap,end(v,kind))
    result['saturation_random_word_checks']=1200
    result['elapsed_seconds']=time.time()-start_time
    (BASE/'audit_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',result['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
