"""Independent literal and complete-axis audit of the1997/5995 escape."""
from pathlib import Path
import json,hashlib
from check import need,H,induced,birth,root_profile,strong,state,path_tables
from audit_reset import coefficient_polynomial,inflate
from wide_return import converse

HERE=Path(__file__).resolve().parent;RESET=HERE.parent/'reset'
BASE=(2046,2044,2024,2032,4068,4032,3968,3840,3584,1024,2048,527)
CORE=(58,60,57,16,32,8)


def literal(cut,A,B):
    left=list(map(int,bin(A)[2:])) if A else []
    right=[1-int(x) for x in reversed(bin(B)[2:])] if B else []
    l=len(left);m=l+6+len(right);T=[0]*(m+1)
    for i in range(m):
        for j in range(i+1,m):
            forward=(CORE[i-l]>>(j-l)&1) if l<=i<j<l+6 else 1
            if forward:T[i]|=1<<j
            else:T[j]|=1<<i
    for i,bit in enumerate(left+[(cut>>i)&1 for i in range(6)]+right):
        if bit:T[m]|=1<<i
        else:T[i]|=1<<m
    return tuple(T)


def g(n):
    return state(''.join('O' if b=='1' else'I' for b in bin(n)[2:]))[0] if n else 1


def recover(cut):
    probes={(a,b):H(literal(cut,a,b)) for a,b in((0,0),(1,0),(2,0),(3,0),(0,1),(0,2),(0,3),(1,1),(2,2))}
    h=probes[0,0];i=(probes[3,0]-probes[1,0])//2;o=(probes[0,3]-probes[0,1])//2
    s=probes[1,0]-h-i;e=probes[0,1]-h-o
    d=(probes[2,0]-h-2*i)//2
    need(probes[0,2]==h+2*d+2*o,'independent transverse consistency')
    L=(probes[2,2]-probes[2,0]-probes[0,2]+h)//4
    need(probes[1,1]-probes[1,0]-probes[0,1]+h==L,'no cross-parity coordinate')
    return(h,d,L,i,o,s,e)


def complete_small_inverse(N,co):
    """Direct exhaustive bounded rectangle plus axes; no recursive inverse."""
    h,d,L,i,o,s,e=co
    need(d==9 and L==144 and min(i,o,s,e)>=0 and max(s,e)<=d,'calibrated positive chart')
    found=set()
    def value(A,B):return h+d*(g(A)+g(B)-2)+L*A*B+i*A+o*B+s*(A&1)+e*(B&1)
    # Both positive: nonnegative remaining terms imply AB <=(N-h)/L.
    if N>=h:
        cap=(N-h)//L
        for A in range(1,cap+1):
            for B in range(1,cap//A+1):
                if value(A,B)==N:found.add((A,B))
    for transpose in(False,True):
        slope=i if transpose else o;parity=s if transpose else e
        if slope:
            for u in range(max(0,(N-h)//slope)+1):
                pair=(u,0) if transpose else(0,u)
                if value(*pair)==N:found.add(pair)
        else:
            need(h==d and parity==0,'uniform-axis structure')
            need(N%d!=0,'entire uniform axis excluded by divisibility, including all trailing-one representatives')
    return sorted(found)


def replacement_trace(source,target,keep,mapping,missing):
    need(induced(source,keep)==induced(target,mapping),'all retained literal arcs')
    labels=list(mapping);trace=[H(source),H(induced(target,labels))]
    for v in missing:labels.append(v);trace.append(H(induced(target,labels)))
    need(set(labels)==set(range(len(target))),'all target vertices inserted exactly once')
    return trace


def all_defects(T):
    end,totals=path_tables(T);full=len(totals)-1
    return [sum(end[mask][r]*totals[full^mask] for mask in range(full+1) if mask>>r&1)-totals[-1] for r in range(len(T))]


def main():
    raw=json.loads((RESET/'two_triangle_inverse.json').read_text())
    calibrated=[];solutions={5993:[],5995:[]}
    for row in raw['inverse_rows']:
        cut=row['cut'];co=recover(cut)
        need(list(co)==row['coefficients'],'independent literal reconstruction of every signature')
        item={'cut':cut,'coefficients':co,'solutions':{}}
        for N in solutions:
            ans=complete_small_inverse(N,co)
            need([list(p) for p in ans]==row['solutions'][str(N)],'independent complete rectangle-plus-axis inverse')
            item['solutions'][str(N)]=ans
            solutions[N].extend((cut,A,B) for A,B in ans)
        calibrated.append(item)
    need(not solutions[5993] and len(solutions[5995])==18,'all-height normalized outcomes')
    p=root_profile(BASE,0)
    need(H(BASE)==1997 and p['F']==2000 and p['indegree']==1 and BASE[11]&1,'literal defect3 root and port')
    co=coefficient_polynomial(BASE,1)
    need(co[:2]==[1997,1134] and not any(co[2:]),'independent complete infinite-source polynomial')
    for m in range(1,4):
        T=inflate(BASE,1,m);p=root_profile(T,0)
        need(strong(T) and p['H']==1134*(1<<m)-271 and p['F']-p['H']==3,'literal source family')
    readydata=json.loads((RESET/'target_ready_return.json').read_text())
    Q=tuple(readydata['controls'][0]['ready665']);source=birth(birth(Q,len(Q)-1),len(Q))
    keep=[0,1,2,3,4,6,7,9,10,11];mapping=[11,2,5,6,7,8,9,10,0,4]
    first_trace=replacement_trace(source,BASE,keep,mapping,[1,3])
    need(first_trace==[1997,241,619,1997],'independent1997 rebase trace')
    T=birth(birth(BASE,0),len(BASE));need(list(T)==raw['standard_rows'],'actual5997 source')
    macro_traces=[]
    for rec in raw['transport_macros']:
        trace=[H(T)]
        for v in rec['triangle']:
            R=list(T);u=rec['singleton'];R[v]^=1<<u;R[u]^=1<<v;T=tuple(R);trace.append(H(T))
        need(trace==rec['full_graph_trace'],'independent full-tournament macro prices')
        macro_traces.append(trace)
    need(H(T)==7237,'pre-recut target')
    need(raw['changed_apex_neighbors']==[8],'one exact apex incidence')
    R=list(T);R[11]^=1<<8;R[8]^=1<<11;T=tuple(R)
    need(list(T)==raw['target_rows'] and H(T)==5995 and strong(T),'all14-name retained+4 endpoint')
    p=root_profile(T,0)
    need((p['H'],p['S'],p['d'],p['F'],p['F']-p['H'])==(5995,339,1897,19759,13764),'old-root interface loss')
    need(all(b not in(1,2) for b in all_defects(T)+all_defects(converse(T))),'neither gauge ready at this exact endpoint')
    ready=tuple(raw['ready5995_rows']);rec=raw['same_root_ready_return']
    second_trace=replacement_trace(T,ready,rec['kept'],rec['map'],rec['insertion_target_labels'])
    need(second_trace==[5995,107,253,545,1085,2165,5995],'independent5995 ready return')
    need(rec['map'][rec['kept'].index(0)]==13 and rec['map'][rec['kept'].index(11)]==0,'same physical root and original port')
    p=root_profile(ready,13)
    need(strong(ready) and p['H']==5995 and p['F']==5996 and p['indegree']==1 and ready[0]>>13&1,'recovered ready root')
    following=birth(birth(ready,13),14);need(H(following)==17987,'lawful subsequent ternary exit')
    paths=[RESET/name for name in('DEFECT_THREE_ESCAPE.md','two_triangle_inverse.py','two_triangle_inverse.json','low_defect_escape.py')]
    out={'status':'PASS','scope':'576 fresh literal signature probes; independent complete bounded positive-product and axis enumeration for all64cuts; actual graphs, root profiles, full prices and both retained-source return traces; independent full run polynomial. Inherited strong classification and ready census not rerun here.',
         'inverse_rows':calibrated,'normalized_solutions':solutions,'source_runs':co,
         'rebase1997_trace':first_trace,'macro_traces':macro_traces,'ready5995_trace':second_trace,'next_H':17987,
         'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (HERE/'defect_escape_audit.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print('PASS: complete independent5993/5995 rectangle+axis inverse; literal1997/5995/17987 constructions and same-root return')


if __name__=='__main__':main()
