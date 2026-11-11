"""Actual named arc-flip transport of singleton gaps past C3 components."""
from pathlib import Path
from itertools import product, combinations
import json
from extension_content import C3,tt,join,h,coefficients,need,path_table

def profile(adj,r):
    n=len(adj);full=(1<<n)-1;end=path_table(adj)
    reverse=tuple(sum(1<<j for j in range(n) if adj[j]>>i&1) for i in range(n))
    start=path_table(reverse);hm=[sum(v) for v in end];hm[0]=1
    U=sum(start[m][r]*hm[full^m] for m in range(full+1) if m>>r&1)
    V=sum(end[m][r]*hm[full^m] for m in range(full+1) if m>>r&1)
    deleted=full^(1<<r)
    K=sum(hm[m]*hm[deleted^m] for m in range(full+1) if not(m>>r&1))
    return [hm[full],end[full][r],start[full][r],U,V,hm[deleted],K]

def flip(adj,u,v):
    out=list(adj);out[u]^=1<<v;out[v]^=1<<u
    return tuple(out)

def start(gaps):
    r=len(gaps)-1; components=[]; next_label=0
    for i,g in enumerate(gaps):
        for _ in range(g):components.append([next_label]);next_label+=1
        if i<r:components.append(list(range(next_label,next_label+3)));next_label+=3
    blocks=[C3 if len(c)==3 else (0,) for c in components]
    return join(blocks),components

def normalize(gaps):
    adj,components=start(gaps);r=len(gaps)-1;original=adj
    records=[]
    # Bubble every singleton left, retaining all labels and all other arcs.
    while True:
        at=next((i for i in range(1,len(components))
                 if len(components[i])==1 and len(components[i-1])==3),None)
        if at is None:break
        tri=components[at-1];v=components[at][0]
        before=adj
        trace=[h(adj)]
        for u in tri:
            need(adj[u]>>v&1,'old triangle beats adjacent singleton')
            adj=flip(adj,u,v);trace.append(h(adj))
        need(trace==[3**r,5*3**(r-1),5*3**(r-1),3**r],
             'exact macro count trace')
        low=flip(before,tri[2],tri[0]);low_trace=[h(before),h(low)]
        for u in reversed(tri):low=flip(low,u,v);low_trace.append(h(low))
        low=flip(low,tri[2],tri[0]);low_trace.append(h(low))
        need(low==adj,'same endpoint for upper and lower excursions')
        need(low_trace==[3**r]+[3**(r-1)]*4+[3**r],'exact lower excursion')
        records.append({'triangle':tri,'singleton':v,'trace':trace,'lower_trace':low_trace})
        components[at-1],components[at]=components[at],components[at-1]
    order=[v for c in components for v in c]
    canonical=tuple(sum((adj[u]>>v&1)<<j for j,v in enumerate(order)) for u in order)
    expected=join([tt(sum(gaps))]+[C3]*r)
    need(canonical==expected,'full named endpoint graph')
    need(len(records)==sum(i*g for i,g in enumerate(gaps)),'exact crossing count')
    return original,adj,records,order

def main():
    controls=[]
    for r in range(1,4):
        for gaps in product(range(3),repeat=r+1):
            if sum(gaps)>3:continue
            old,new,records,order=normalize(gaps)
            expected=3**(r-2 if r>=2 else 0)
            need(coefficients(new)[3]==expected,'normalized first content')
            controls.append({'r':r,'gaps':gaps,'macros':len(records),
                             'peak':max([3**r]+[max(x['trace']) for x in records])})
    old,new,records,order=normalize((0,1,0))
    first_steps=[{'arc':[i,j],'H':h(flip(old,i,j))}
                 for i,j in combinations(range(len(old)),2)]
    need(min(x['H'] for x in first_steps if x['H']>=9)==15,
         'minimal seven-vertex H-nondecreasing first-step barrier')
    need(coefficients(old)[3]==3 and coefficients(new)[3]==1,'content repair')
    rooted_old=join((old,(0,)));rooted_new=join((new,(0,)))
    old_profile=profile(rooted_old,len(old));new_profile=profile(rooted_new,len(new))
    need(old_profile==new_profile==[9,9,0,9,288,9,288],'equal full seven-coordinate rooted states')
    need(coefficients(rooted_old)[3]==3 and coefficients(rooted_new)[3]==1,
         'same full rooted state does not fix arbitrary extension content')
    out={'status':'PROVED macro controls / FINITE-EXACT sharp minimal-source peak',
         'source':old,'target':new,'named_target_order':order,'macro_records':records,
         'minimal_source_all21_first_steps':first_steps,'family_controls':controls,
         'full_rooted_state_collision':{'first_rows':rooted_old,'second_rows':rooted_new,
             'root':len(old),'profile':old_profile,'first_C1':3,'second_C1':1},
         'control_count':len(controls)}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'result':'PASS','family_controls':len(controls),
                      'minimal_peak':15,'source_content':3,'target_content':1}))

if __name__=='__main__':main()
