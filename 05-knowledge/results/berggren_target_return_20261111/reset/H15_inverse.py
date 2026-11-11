"""Complete H15-menu target inverse and an actual S4 defect-five escape."""
from pathlib import Path
import importlib.util,json,sys
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parent.parent
sys.path.insert(0,str(RESULTS/'uniform_return_metagraph_20261111/arithmetic'))
import chart_inverse as inv
spec=importlib.util.spec_from_file_location('H15_literal',RESULTS/'session251_creative_completion_20261111/ternary/compiler.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
sys.path.insert(0,str(RESULTS/'uniform_return_metagraph_20261111/reset'))
from classify_modules import H,strong,need,inherited as a,two_birth as b
from gap_return_probe import BASE
from finite_bridge import replace
from target_ready_return import expand_map

def join(A,B):
    return tuple(x|(((1<<len(B))-1)<<len(A)) for x in A)+tuple(x<<len(A) for x in B)

def main():
    S4=(2,4,9,3);C3=(2,4,1);report=[];small=[]
    for name,core in [('S4_C3',join(S4,C3)),('C3_S4',join(C3,S4)),('S4',S4)]:
        rows=[];hits=[]
        for cut in range(1<<len(core)):
            d,kp,km,ip,im,op,om,jp,jm=c.signature(core,cut)
            need(d+km==0,'positive observer parity-cross cancellation')
            coeff=(d+jp+jm,d,d+kp,ip,op,-im,-om)
            pairs=inv.chart_inverse(5993,coeff)
            rows.append(dict(cut=cut,coefficients=coeff,solutions=pairs))
            for A,B in pairs:
                T=tuple(c.construct(core,cut,A,B))
                need(H(T)==5993 and strong(T),'independent positive normalized witness')
                hit=dict(cut=cut,A=A,B=B,rows=T);hits.append(hit)
                if name=='S4':small.append(hit)
        need(len(hits)==(6 if name=='S4' else 0),'complete all-height finite-core inverse answers')
        report.append(dict(name=name,core=core,cut_universe=len(rows),inverse=rows,hits=hits))
        print('CORE',name,'cuts',len(rows),'hits',len(hits),flush=True)
    module_checks=[];total=0
    for hit in small:
        T=hit['rows'];data=[]
        for v in range(len(T)):
            co=b.module_runs(T,v);total+=1
            need(not(co[2]==3402 and not any(co[3:])),
                 'no single-placeholder lift of any normalizedS4witness has required ternary slope')
            data.append(dict(vertex=v,runs=co))
        module_checks.append(dict(cut=hit['cut'],A=hit['A'],B=hit['B'],vertices=data))
        print('MODULE',hit['cut'],hit['A'],hit['B'],'order',len(T),flush=True)
    need(total==84,'complete six-witness all-vertex placeholder universe')
    target=next(tuple(h['rows']) for h in small if (h['cut'],h['A'],h['B'])==(1,3,57))
    kept=[0,1,2,4,5,6,7,8,10,11];mapping=[0,1,2,5,3,6,9,10,11,12]
    target,names,rewrite=replace(BASE,[f'old{i}' for i in range(12)],target,kept,mapping,['new4','new7','new8'])
    need(names[0]=='old0' and names[12]=='old11','same physical root and sole in-neighbour')
    need(a.root_profiles(target)[0]['gap']==5,'actual defect-five target')
    tr=b.module_runs(target,1)
    need(tr[1:3]==[5993,3684] and not any(tr[3:]),'complete uniform target family count')
    controls=[]
    for m in range(1,5):
        S=b.inflate(BASE,1,m);T=b.inflate(target,1,m)
        ids,p=expand_map(kept,mapping,1,1,m)
        need(a.induced(S,ids)==a.induced(T,p),'all-k retained module and old-arc map')
        D=1<<m;sp=a.root_profiles(S)[0];tp=a.root_profiles(T)[0]
        need(sp['gap']==3 and tp['gap']==5 and strong(S) and strong(T),'actual uniform interface types')
        need((H(S),H(T))==(1134*D-271,3684*D-1375),'independent uniform count controls')
        need(H(T)-(3*H(S)+2)==282*(D-2),'exact drift and unique first-size target equality')
        controls.append(dict(m=m,source_H=H(S),target_H=H(T),surplus=282*(D-2),kept=ids,map=p))
    result=dict(status='PASS',scope='All128cuts ofeach H15core and16S4cuts withunbounded exact flanks;84complete module expansions',
                cores=report,module_checks=module_checks,module_universe=total,
                actual_direct_plus2=rewrite,root0_target_profile=a.root_profiles(target)[0],
                uniform_target_runs=tr,uniform_controls=controls)
    (HERE/'H15_inverse.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: H15menuexcluded; actual1997->5993same-rootdefect3->5; uniformsurplus282(2^m-2)')
    print('TRACE',rewrite['trace'])

if __name__=='__main__':main()
