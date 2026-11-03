"""What contraction retains: scalar counts, named ranks and completion budgets."""
from itertools import permutations,product,combinations
import json
from support import HERE,PREV,need,save,junction
from completion_hostile import cover

def hierarchy_orders(labels):
    if len(labels)==1:return [tuple(labels)]
    m=len(labels)//2;out=[]
    for A in hierarchy_orders(labels[:m]):
        for B in hierarchy_orders(labels[m:]):out.extend((A+B,B+A))
    return out
def main():
    trials=[]
    for m in (2,4,8):
        orders=list(permutations(range(m))) if m<=4 else None
        for label,bank in [('bare_fibre',orders),('unordered_binary_hierarchy',hierarchy_orders(list(range(m))))]:
            if bank is None:continue
            # TT_m: the one-root incoming cut has H=1+2^rank-[rank=0].
            values=[(v,[1+2**P.index(v)-(P.index(v)==0) for P in bank]) for v in range(m)]
            query=cover(values,[1]*len(bank));expected=m-1 if label=='bare_fibre' else m//2
            need(len(query)==expected,'exact finite minimum single-mark completion')
            if m<=4:
                for P in bank:
                    T=tuple(sum(1<<P[j] for j in range(i+1,m)) for v in range(m) for i in [P.index(v)])
                    for v,col in values:need(junction.H(junction.ear(T,1<<v))==col[bank.index(P)],'literal transitive-fibre readout')
            trials.append({'m':m,'retained':label,'source_count':len(bank),'minimum_queries':expected,'query_labels':query})
    packet=json.loads((PREV/'halving_normalization.json').read_text());controls=[]
    for row in packet['programs']:
        if row['k']>2:continue
        T=tuple(row['final_rows']);U=tuple(row['first_quotient']);W=tuple(row['odd_core']);fibres=row['odd_core_original_fibres']
        profiles=[junction.profile(T,r) for r in range(len(T))];read=[h+F-ss for h,ss,q,F in profiles]
        ordered=[]
        for fibre in fibres:
            order=sorted(fibre,key=lambda r:sum(T[r]>>v&1 for v in fibre),reverse=True);ordered.append(order)
            need(all(read[a]<read[b] for a,b in zip(order,order[1:])),'strict position readout within actual strong fibre')
            for a,b in zip(order,order[1:]):
                keep=[v for v in range(len(T)) if v!=b];C=junction.sub(T,keep);rr=keep.index(a);h,ss,q,F=junction.profile(C,rr)
                need(read[b]-read[a]==F+ss,'exact adjacent-twin response difference')
        for r in range(len(T)):need(junction.H(junction.ear(T,1<<r))==read[r],'literal newborn versus rooted observer')
        controls.append({'k':row['k'],'flipped':row['flipped'],'H_before_normalization':junction.H(tuple(row['source_rows'])),
          'H_after_normalization':junction.H(T),'H_first_halving':junction.H(U),'H_odd_core':junction.H(W),
          'ordered_original_fibres':ordered,'single_port_response':read,'bare_fibre_completion':3*len(fibres),'hierarchy_completion':2*len(fibres)})
        junction.endpoints.cache_clear()
    save('halving_response.json',{'status':'EXACT_FINITE_CONTROLS_FOR_ALL_SIZE_TWIN_COMPLETION_THEOREM','finite_atlases':trials,'actual_normalization_controls':controls,
      'scope':'single-old-vertex incoming cuts; quotient and named fibre membership known; all-size minima proved separately, not extrapolated'})
    print('PASS halving completion: bare fibres m-1; binary hierarchy m/2; ordered fibres0',flush=True)
if __name__=='__main__':main()
