"""Full joint certificate for the first new five-port strong birth."""
from common import *
from joint_nine import kernel,expo,cones
from synchrony_probe import expand
from sparse_covers import kernel as cover_kernel
import time

def run():
    r=read('rank_candidates.json.gz')['records'][0];P=r['stable_ports'];start=time.time()
    K=kernel(r['rows'],r['selected'],P);print('signed kernel',len(K),'seconds',time.time()-start,flush=True)
    D,U=cover_kernel(r['rows'],P,[v for v in range(len(r['rows'])) if r['selected']>>v&1]);need(K==D,'independent complete joint kernel');print('independent kernel',len(D),len(U),'seconds',time.time()-start,flush=True)
    need(K=={(1,)*len(P):1},'joint stability of all five original ports')
    controls=[]
    for lengths in [[2,1,2,1,2],[1,2,1,2,1],[2]*5]:
        rr,O,_=expand(r['rows'],r['selected'],dict(zip(P,lengths)));h,u,w,_=tables(np.array(rr,dtype=np.int64));d=defect(h,w,[v for v in range(len(rr)) if O>>v&1]);need(d==1,'independent nonuniform literal joint control');controls.append({'lengths':lengths,'order':len(rr),'H':int(h),'defect':d})
    out={'scope':'All positive independent lengths of five named unselected ports in one explicit strong order13 birth. Full signed tensor and independent path-cover compiler, plus literal controls. No minimum-order claim.','rows':r['rows'],'selected':r['selected'],'ports':P,'defect_Stirling':[[*k,c] for k,c in sorted(K.items())],'U_Stirling':[[*k,c] for k,c in sorted(U.items())],'literal_controls':controls}
    save('rank_kernel.json.gz',out);print('complete',time.time()-start,flush=True)
if __name__=='__main__':run()
