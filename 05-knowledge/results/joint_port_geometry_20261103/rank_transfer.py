"""Resolve the physical-source rank boundary before enlarging a census."""
import json,random
import sympy as s
from return_probe import HERE,junction,need,observer
import effects as e

def main():
    p=1009;basis={};witnesses=[];examined=0
    def reduce(row):
        row=[int(x)%p for x in row]
        for pivot,b in sorted(basis.items()):
            c=row[pivot]
            if c:row=[(x-c*y)%p for x,y in zip(row,b)]
        return row
    def add(U,ports,provenance):
        nonlocal examined
        V=observer(U,ports);v=[int(V.get(k,0)) for k in e.STATES];row=[x%p for x in v];examined+=1
        row=reduce(row)
        pivot=next((j for j,x in enumerate(row) if x),None)
        if pivot is not None:
            inv=pow(row[pivot],-1,p);basis[pivot]=[(x*inv)%p for x in row]
            witnesses.append({'rows':U,'ports':ports,'cap_vector':v,'provenance':provenance})
            print('rank',len(basis),'examined',examined,flush=True)
        return pivot is not None
    for n in (3,4,5):
        for number,T in enumerate(junction.all_tournaments(n)):
            if not junction.strong(T):continue
            for r in range(n):
                U=junction.ear(junction.ear(T,1<<r),1<<n)
                for other in list(range(n))+[n+1]:add(U,(n,other),['seal',n,number,r,other])
            if len(basis)>=31:break
        if len(basis)>=31:break
    N=e.macro(0);ops={'3':e.macro(3),'5':e.macro(5),'N':N,**{'J'+str(m):e.dense((0,m,2))*N for m in range(4)}}
    W=s.Matrix([x['cap_vector'] for x in witnesses])
    print('actual initial',len(basis),flush=True)
    for source in list(witnesses):
        v=s.Matrix([source['cap_vector']])
        for m in range(4):
            row=v*ops['J'+str(m)]
            if any(reduce(list(row))):
                T=tuple(source['rows']);r,b=source['ports'];n=len(T);inc=((1<<r) if m&1 else 0)|((1<<b) if m&2 else 0)
                A=junction.ear(T,inc);A=junction.ear(A,1<<r);U=junction.ear(A,1<<(n+1))
                need(add(U,(n+1,b),['J',m,source['provenance']]),'actual attachment realizes new source direction')
                need(witnesses[-1]['cap_vector']==list(row),'literal cap vector matches formal joint operator')
                W=s.Matrix([x['cap_vector'] for x in witnesses]);print('expanded actual rank',W.rows,flush=True)
        if W.rows==32:break
    cases=json.loads((HERE/'joint_probe.json').read_text())['cases'];tests=[]
    for c in cases:
        from algebra import mod_rank
        rank=mod_rank(W*s.Matrix(c['basis']),p)
        need(rank==c['sealed_rank'],'actual-source lower bound meets universal upper bound')
        tests.append({'masks':c['masks'],'actual_hankel_rank':rank})
    (HERE/'rank_transfer.json').write_text(json.dumps({'status':'FINITE_EXACT_ACTUAL_RANK','source_rank':W.rows,'examined':examined,
        'witnesses':witnesses,'observable_tests':tests},indent=2)+'\n',encoding='utf8',newline='\n')
    print(tests,flush=True)
if __name__=='__main__':main()
