"""Exact one-birth response completion and holonomy on the actual 81/27 fibre."""
from itertools import combinations
from collections import deque
import json,importlib.util
from support import HERE,PREV,ROOT,need,save,junction
from hostile_embedding import tournament,group,power

def push(T,p):
    inv=[p.index(i) for i in range(len(p))];return junction.sub(T,inv)
def polynomial(T):
    n=len(T);full=(1<<n)-1;E=junction.endpoints(T);R=tuple(full^x^(1<<i) for i,x in enumerate(T));S=junction.endpoints(R)
    H=sum(E[-1]);K=[[0]*n for _ in T]
    for mask in range(1,full):
        for v,x in enumerate(E[mask]):
            if x:
                for u,y in enumerate(S[full^mask]):K[v][u]+=x*y
    linear=[E[-1][v]-S[-1][v]+sum(K[v]) for v in range(n)]
    quadratic={(v,u):K[v][u]+K[u][v] for v,u in combinations(range(n),2)}
    return H,linear,quadratic
def responses(T):
    H,L,Q=polynomial(T);n=len(T);vals=[H]*(1<<n)
    for m in range(1,1<<n):
        bit=m&-m;v=bit.bit_length()-1;rest=m^bit
        vals[m]=vals[rest]+L[v]-sum(Q[min(v,u),max(v,u)] for u in range(n) if rest>>u&1)
    return vals
def cover(values,initial):
    pairs=[(a,b) for a,b in combinations(range(len(initial)),2) if initial[a]==initial[b]]
    full=(1<<len(pairs))-1;dictionary={}
    for name,col in values:
        mask=sum((col[a]!=col[b])<<i for i,(a,b) in enumerate(pairs))
        if mask:dictionary.setdefault(mask,name)
    todo=deque([0]);known={0:()}
    while todo:
        m=todo.popleft()
        if m==full:return known[m]
        for add,name in dictionary.items():
            out=m|add
            if out not in known:known[out]=known[m]+(name,);todo.append(out)
    return None
def main():
    d=json.loads((PREV/'hostile_embedding.json').read_text());G=sorted(group([tuple(d['translation_permutation']),tuple(d['doubling_permutation'])]))
    A,B=tournament(False),tournament(True);orbits=[sorted({push(T,p) for p in G}) for T in (A,B)];states=sum(orbits,[])
    need(list(map(len,orbits))==[2,6] and len(set(states))==8,'actual disjoint physical orbit fibre')
    path=ROOT/'05-knowledge/results/regular_parent_symmetry_20261029/verify.py';spec=importlib.util.spec_from_file_location('fresh81',path);h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    need((len(h.iso_maps(A,A)),len(h.iso_maps(B,B)))==(81,27),'fresh full hostile automorphisms')
    rows=[responses(T) for T in states];baseline=[row[0] for row in rows];cuts=[(m,[row[m] for row in rows]) for m in range(1<<12)]
    chosen=cover(cuts,baseline);need(chosen is not None,'one-birth contexts jointly complete on fibre')
    frames=[]
    for a,b in combinations(range(12),2):
        menu=[(m,col) for m,col in cuts if not(m&~((1<<a)|(1<<b)))];picked=cover(menu,baseline)
        frames.append({'ports':[a,b],'completion':None if picked is None else list(picked)})
    action=[[states.index(push(T,g)) for T in states] for g in G]
    levels=[]
    for label,menu in [('baseline',[]),('all_singletons',[1<<v for v in range(12)]),('all_one_birth',list(range(1<<12)))]:
        invisible=[i for i,perm in enumerate(action) if all(baseline[j]==baseline[perm[j]] and all(rows[j][m]==rows[perm[j]][m] for m in menu) for j in range(8))]
        classes={}
        for j in range(8):classes.setdefault((baseline[j],)+tuple(rows[j][m] for m in menu),[]).append(j)
        levels.append({'menu':label,'response_classes':list(classes.values()),'invisible_group_order':len(invisible),'elements':[G[i] for i in invisible]})
    widths=[]
    for width in range(1,5):
        answer=cover([(m,col) for m,col in cuts if m.bit_count()<=width],baseline)
        widths.append({'incoming_width':width,'optimal_contexts':answer})
    kernel=group([power(tuple(d['translation_permutation']),3),power(tuple(d['doubling_permutation']),2)])
    need(len(kernel)==9 and all(set(map(tuple,row['elements']))==kernel for row in levels[1:]),'actual affine reduction kernel')
    for j,T in enumerate(states):
        for cut in set([0,1,3,85,4095]+list(chosen)):
            need(junction.H(junction.ear(T,cut))==rows[j][cut],'independent literal one-birth DP')
        for cut in (1,288,291):need(junction.strong(junction.ear(T,cut)),'chosen read heads preserve strongness')
        junction.endpoints.cache_clear()
    save('completion_hostile.json',{'status':'FINITE_EXACT_FULL_ONE_BIRTH_MENU_ON_EIGHT_ACTUAL_SOURCES','states':states,
      'baseline_H':baseline,'orbit_sizes':[2,6],'group_order':54,'optimal_contexts':list(chosen),'optimal_responses':[[row[m] for m in chosen] for row in rows],
      'singleton_completion':cover([(m,col) for m,col in cuts if m.bit_count()==1],baseline),
      'completion_by_incoming_width':widths,
      'two_old_port_frames':frames,'holonomy_filtration':levels,
      'distinct_one_birth_response_columns':len({tuple(c) for m,c in cuts}),
      'scope':'fixed old labels; arbitrary incoming cuts for one named newborn; no strong-output filter; group transformations are not free native commands'})
    print('PASS hostile completion',chosen,'holonomy orders',[x['invisible_group_order'] for x in levels],flush=True)
if __name__=='__main__':main()
