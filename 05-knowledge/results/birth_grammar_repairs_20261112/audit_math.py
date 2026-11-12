"""Independent small-universe laws and positive/hostile controls."""
from common import *
from birth_extensions import quadratic
from extension_heads import add,end_counts,attachments
from hamilton_clocks import graphs,skeleton,trim
from new_axis_bank import slope_skeleton
from sparse_covers import kernel as covers
from joint_nine import kernel,cones,expo,evaluate
from collections import Counter
from math import factorial

def run():
    stats=Counter()
    for rows in graphs(5):
        n=len(rows);a=np.array(rows,dtype=np.int64);h,u,w,_=tables(a)
        H=attachments(h,u,w,end_counts(a));literal=[]
        for cut in range(1<<n):
            rr=add(rows,cut);hh,uu,ww,_=tables(np.array(rr,dtype=np.int64));need(int(hh)==int(H[cut]),'every small attachment quadratic');stats['attachment_controls']+=1
            literal.append((int(hh),ww))
        for O in range(1<<n):
            cube=quadratic(rows,O)
            for cut,(hh,ww) in enumerate(literal):
                for mark in [0,1]:
                    selected=[v for v in range(n+1) if (O|(mark<<n))>>v&1]
                    need(int(cube[mark][cut])==defect(hh,ww,selected),'entire small defect cubes');stats['defect_cube_controls']+=1
        for p in range(n):
            hc=[int(x) for x in skeleton(a,p)];ih,_,_=ib.coefficients(rows,p)
            need(hc==[int(x) for x in ih[1:len(hc)+1]] and not any(ih[len(hc)+1:]),'independent small whole H skeleton');stats['Hamilton_skeleton_controls']+=1
            O=((1<<n)-1)^(1<<p);D,U=covers(rows,[p],[v for v in range(n) if O>>v&1]);uc=[int(x) for x in slope_skeleton(a,O,p)]
            need({(j+1,):c for j,c in enumerate(uc) if c}==U,'positive U skeleton versus endpoint convolution');stats['U_skeleton_controls']+=1
        if n>=2:
            P=[0,1];O=list(range(2,n));D,U=covers(rows,P,O);oldD,oldU=cones.kernel(rows,P,O)
            need((D,U)==(oldD,oldU),'sparse versus dense endpoint convolution');need(D==kernel(rows,sum(1<<v for v in O),P),'signed versus positive endpoint kernels');stats['joint_kernel_controls']+=1
    # A real mixed obstruction is mandatory, not only unit-positive examples.
    rows=[16,1,67,7,142,31,187,303,127];D,U=covers(rows,[5,8],[1,3]);need(evaluate(expo(D),[2,2])==-3 and D=={(1,1):1,(2,2):-1},'J1 hostile joint gate')
    sync=read('synchrony_kernels.json.gz')['records'][0];D,U=covers(sync['rows'],sync['ports'],[v for v in range(len(sync['rows'])) if sync['selected']>>v&1])
    need([[*b,c] for b,c in sorted(expo(D).items())]==sync['defect_exponential'] and [[*b,c] for b,c in sorted(expo(U).items())]==sync['U_exponential'],'complete sparse synchrony replay')
    # Positive-length exponential continuation is not vertex deletion.
    h=[int(x) for x in skeleton(np.array([2,4,1],dtype=np.int64),0)];need(expcoeff(h)==[1,1] and sum(expcoeff(h))==2,'C3 zero-clock hostile');need(int(tables(np.array([2,0],dtype=np.int64))[0])==1,'actual deleted C3 is transitive')
    # Fixed-cut two-cover upper bound; every literal producer here uses n<=21.
    bounds={n:factorial(n)*(2*n-1)//2**(n//2) for n in range(2,22)}
    need(max(bounds.values())<2**63,'all n<=21 endpoint counts fit signed64 under disjoint-pair bound')
    need(read('hamilton_clocks_9.json.gz')['source_counts']=={'1':1,'2':1,'3':2,'4':4,'5':12,'6':56,'7':456,'8':6880,'9':191536},'complete ordinary profile-bank source counts')
    out={'scope':'Every ordinary tournament through5, every new orientation and every selected subset for the quadratic defect law; every small pointed H/U kernel; separate sparse/dense and signed/endpoint methods; physical mixed hostile and zero-clock hostile.','statistics':dict(stats),'hostiles':['J1 double inflation has defect-3','C3 exponential value at length0 is2, actual deletion has H1'],'largest_signed64_endpoint_bound':max(bounds.values()),'maximum_literal_order_supported':21}
    save('audit_math.json.gz',out);print(json.dumps(out),flush=True)
if __name__=='__main__':run()
