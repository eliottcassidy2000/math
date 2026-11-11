"""All-height two-triangle-core inverse and a named no-deletion+4 escape."""
from pathlib import Path
import importlib.util,json,sys
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parent.parent
sys.path.insert(0,str(RESULTS/'uniform_return_metagraph_20261111/arithmetic'))
import chart_inverse as inv
spec=importlib.util.spec_from_file_location('literal_compiler',RESULTS/'session251_creative_completion_20261111/ternary/compiler.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
sys.path.insert(0,str(RESULTS/'uniform_return_metagraph_20261111/reset'))
from classify_modules import H,strong,need,inherited as a,two_birth as b
from gap_return_probe import BASE,flip
from finite_bridge import replace
from funnel_six_state import a as funnel
CORE=(58,60,57,16,32,8)

def main():
    rows=[];hits={5993:[],5995:[]}
    for cut in range(64):
        d,kp,km,ip,im,op,om,jp,jm=c.signature(CORE,cut)
        need(d==9 and d+km==0,'exact two-triangle signature type')
        coeff=(d+jp+jm,d,d+kp,ip,op,-im,-om)
        record=dict(cut=cut,coefficients=coeff,solutions={})
        for N in hits:
            pairs=inv.chart_inverse(N,coeff)
            record['solutions'][str(N)]=pairs
            for A,B in pairs:
                T=c.construct(CORE,cut,A,B)
                need(H(T)==N and strong(T),'independent all normalized inverse positive witnesses')
                hits[N].append(dict(cut=cut,A=A,B=B,rows=T))
        rows.append(record)
    need(not hits[5993] and len(hits[5995])==18,'complete all-height normalized core inverse')
    standard=a.birth_beaten_only_by(a.birth_beaten_only_by(BASE,0),len(BASE))
    T=standard;trace=[H(T)];macros=[]
    # Named source residual order: newtriangle, gap1, oldtriangle, six right singletons.
    blocks=[[0,12,13],[1],[2,3,4],[5],[6],[7],[8],[9],[10]]
    for singleton in (1,5,6):
        while True:
            i=next(j for j,B in enumerate(blocks) if B==[singleton])
            if i==0 or len(blocks[i-1])==1:break
            tri=blocks[i-1];local=[H(T)]
            for v in tri:T=flip(T,v,singleton);local.append(H(T));trace.append(H(T))
            blocks[i-1],blocks[i]=blocks[i],blocks[i-1]
            macros.append(dict(singleton=singleton,triangle=tri,full_graph_trace=local))
    order=[v for B in blocks for v in B]
    need(order==[1,5,6,0,12,13,2,3,4,7,8,9,10] and len(macros)==5,'exact named gap redistribution')
    apex=11;cut=3273;changed=[];star_trace=[H(T)]
    for j,v in enumerate(order):
        wanted=(cut>>j)&1
        if (T[apex]>>v&1)!=wanted:
            T=flip(T,apex,v);changed.append(v);star_trace.append(H(T))
    need(H(T)==5995==3*H(BASE)+4 and strong(T),'actual all-name-retained+4 target')
    need(len(T)==len(standard)==14,'no birth/deletion in repair')
    for tri in ((0,12,13),(2,3,4)):
        need(a.induced(T,tri)==a.induced(standard,tri),'both originaltriangle orientations retained')
    need(T[11]&1,'original in-neighbour still beats original root')
    profile=a.root_profiles(T)[0]
    ready=funnel.unrank(13,578);ready=a.birth_beaten_only_by(ready,0)
    kept=[0,2,3,7,8,9,10,11,12];mapping=[13,4,5,8,9,10,11,0,3]
    ready,names,return_ready=replace(T,[f'old{i}' for i in range(14)],ready,kept,mapping,
                                   ['new1','new2','new6','new7','new12'])
    need(names[13]=='old0' and names[0]=='old11','same root and original in-neighbour survive ready return')
    need(H(ready)==5995 and a.root_profiles(ready)[13]['gap']==1 and strong(ready),'literal1-ready5995 return')
    following=a.birth_beaten_only_by(a.birth_beaten_only_by(ready,13),14)
    need(H(following)==17987,'one further lawful+2 exit')
    # Singleton-module slope obstruction in both normalized flank patterns.
    slopes=[]
    for A,B in ((4,9),(9,4)):
        witness=next(x for x in hits[5995] if x['A']==A and x['B']==B)
        U=tuple(witness['rows']);left=A.bit_length();base=len(U)-1
        entry=[]
        for v in list(range(left))+list(range(left+6,base)):
            runs=b.module_runs(U,v)
            need(not any(runs[3:]),'all singleton flank responses are dyadic')
            need(runs[2]!=3402,'no preserved single-module uniform+4 coefficient')
            entry.append(dict(vertex=v,slope=runs[2]))
        slopes.append(dict(A=A,B=B,coefficients=entry))
    result=dict(status='PASS',scope='All64two-trianglecore cuts, exact unbounded binaryflank inverse; named14vertex+4repair',
                core=CORE,inverse_rows=rows,normalized_hits={str(N):h for N,h in hits.items()},
                standard_rows=standard,transport_macros=macros,residual_target_order=order,
                apex=apex,target_cut=cut,changed_apex_neighbors=changed,apex_recut_trace=star_trace,
                target_rows=T,old_root_profile=profile,singleton_module_slopes=slopes,
                same_root_ready_return=return_ready,ready5995_rows=ready,next_target=17987)
    (HERE/'two_triangle_inverse.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: all-height5993two-triangleexclusion;18normalized5995charts;namedzero-loss+4repair')
    print('MACRO_TRACES',[x['full_graph_trace'] for x in macros]);print('STAR_TRACE',star_trace)
    print('OLD_ROOT_PROFILE',profile)
    print('READY_RETURN_TRACE',return_ready['trace'])

if __name__=='__main__':main()
