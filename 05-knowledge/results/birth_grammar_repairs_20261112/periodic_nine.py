"""A concrete formal repair interface versus every admissible65-bank period."""
from common import *
from repair_obstructions import det
from collections import Counter
import time
sys.path.insert(0,str(HERE.parent/'inverse65_zeckendorf_20261103'))
from bank import MATS

def run():
    stats=Counter();hits=[];one=0;start=time.time()
    u=[6,9,5]
    need(all(5*sum(B[i][j]*u[j] for j in range(3))>=9*u[i] for B in MATS.values() for i in range(3)),'complete inherited positive norm-growth certificate')
    def visit(A,word):
        nonlocal one
        if word:
            stats[len(word)]+=1
            if det([[int(i==j)-A[i][j] for j in range(3)] for i in range(3)])==0:
                one+=1
                if det([[9*int(i==j)-A[i][j] for j in range(3)] for i in range(3)])==0:hits.append({'word':word,'matrix':A})
        if len(word)==3:return
        for label,B in MATS.items():visit([[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)],word+[label])
    visit([[int(i==j) for j in range(3)] for i in range(3)],[])
    need(sum(stats.values())==278915,'entire period universe through3')
    from construct import tail
    band={r['c']:r for r in read('coefficient_band.json.gz',True)['records']};witnesses=[{'m':0,'H':1,'rows':[0]}];small={}
    for r in read('hamilton_clocks_9.json.gz')['records']:
        H=value(r['exponential_coefficients'],1)
        if H not in small:small[H]=r['rows']
    for m in range(1,7):
        N=(9**m+1)//2;s=v2(N-1);c=(N-1)//2**s;need(s==2+v2(m),'exact LTE invoice for formal family');r=band[c]
        rows=tail(r['rows'],[v for v in range(len(r['rows'])) if r['selected']>>v&1],s-r['v2U']) if s>r['v2U'] else small[N]
        need(len(rows)<=21 and int(tables(np.array(rows,dtype=np.int64))[0])==N,'literal physical finite control outside the excluded periodic strategy')
        witnesses.append({'m':m,'H':N,'c':c,'s':s,'rows':rows})
    out={'scope':'All65-bank periods of lengths1..3, exact eigenvalue1/9 screen. The inherited positive norm growth9/5 per letter excludes longer physical periods for any output asymptotic to9^m, so this test addresses the formal candidate(9^m+1)/2 under fixed prefix/period/suffix. It does not exclude adaptive words or arbitrary tournaments. Actual graphs are supplied only for the finite controls m0..6.','period_counts':dict(stats),'eigenvalue_one_periods':one,'joint_eigenvalue_one_nine_periods':hits,'formal_candidate_first_values':[(9**m+1)//2 for m in range(6)],'physical_finite_controls':witnesses}
    save('periodic_nine.json.gz',out);print(json.dumps(out),flush=True);print('seconds',time.time()-start,flush=True)
if __name__=='__main__':run()
