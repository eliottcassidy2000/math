"""Two physically composable star repairs, with a retained arbitrary bulk module."""
from pathlib import Path
import importlib.util,json
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('owned_one_stage_family',HERE/'star_family.py');f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f);e=f.e
C1=((183,1,19,1509),(90,0,10,1484),(0,0,0,468))
C2=((205,1,65,3107),(102,0,34,2950),(0,0,0,890))
def induced(T,keep):return tuple(sum(1<<j for j,v in enumerate(keep) if T[u]>>v&1) for u in keep)
def construct(B):
    a=len(B);module=tuple(row|(1<<a) for row in B)+(0,);T=f.source(module);P=list(T)
    patch1=[(0,1),(0,5),(0,7),(1,2),(1,4)]+[(1,8+j) for j in range(a+1)]
    for u,v in patch1:P[u]^=1<<v;P[v]^=1<<u
    P=tuple(P);Z,r=e.birth(P,0,True);Q=list(Z)
    patch2=[(0,2),(0,8+a),(0,5),(0,6),(0,r),(2,4),(3,4),(1,4),(4,r)]+[(8+j,4) for j in range(a+1)]
    for u,v in patch2:Q[u]^=1<<v;Q[v]^=1<<u
    Q=tuple(Q);e.need(Q[6]==0,'second parent has literal sink6');keep=[v for v in range(len(Q)) if v!=6];R=induced(Q,keep);W,wr=e.birth(R,0,True)
    return dict(source=T,parent1=P,ready1=Z,ready1_root=r,parent2_with_sink=Q,parent2=R,ready2=W,ready2_root=wr,patch1=patch1,patch2=patch2)
def skeleton(Q):
    n=len(Q)-1;M=n;full=(1<<n)-1;cap=sum(bool(Q[v]>>M&1) for v in range(n))+2
    def table(rows):
        end=[[[0]*(n+1) for _ in range(1<<n)] for _ in range(cap+1)]
        for v in range(n):end[0][1<<v][v]=1
        end[1][0][M]=1
        for j in range(cap+1):
            for mask in range(1<<n):
                for last,value in enumerate(end[j][mask]):
                    if not value:continue
                    for v in range(n):
                        if not(mask>>v&1) and rows[last]>>v&1:end[j][mask|(1<<v)][v]+=value
                    if j<cap and last!=M and rows[last]>>M&1:end[j+1][mask][M]+=value
        h=[[sum(end[j][mask]) for mask in range(1<<n)] for j in range(cap+1)];h[0][0]=1
        return h,end
    H,E=table(Q);_,S=table(e.opposite(Q));out=[]
    for j in range(1,cap+1):
        F=sum(E[a][mask][0]*H[j-a][full^mask] for a in range(j+1) for mask in range(1<<n) if mask&1)
        out.append((H[j][full],S[j][full][0],H[j][full^1],F))
    return out
def main():
    seed=construct((0,));q1=induced(seed['parent1'],[0,1,2,3,4,5,6,7,9,8]);q2=induced(seed['parent2'],[0,1,2,3,4,5,6,8,9,7])
    c1=skeleton(q1);c2=skeleton(q2)
    e.need(tuple(c1[:3])==C1 and all(x==(0,0,0,0) for x in c1[3:]),'first exact skeleton coefficients')
    e.need(tuple(c2[:3])==C2 and all(x==(0,0,0,0) for x in c2[3:]),'second exact skeleton coefficients')
    controls=0
    for n in range(1,4):
        for B in e.labelled(n):
            x=construct(B);k=f.K(B);h,v,w=k
            e.need(e.state(x['parent1'],0)==f.linear(C1,k),'first arbitrary bulk profile')
            e.need(e.state(x['parent2'],0)==f.linear(C2,k),'second arbitrary bulk profile')
            e.need(e.state(x['parent2_with_sink'],0)==e.state(x['parent2'],0),'neutral sink deletion')
            a=e.state(x['ready1'],x['ready1_root']);b=e.state(x['ready2'],x['ready2_root'])
            e.need(a[0]==1691*h+1574*v+468*w and b[0]==3311*h+3052*v+890*w,'composed output counts')
            e.need(a[2]-a[1]==b[2]-b[1]==h,'same exact readiness defect')
            e.need(e.strong(x['parent1']) and e.strong(x['parent2']) and e.strong(x['ready1']) and e.strong(x['ready2']),'strong physical endpoints')
            controls+=1;e.tables.cache_clear()
    family=[]
    for a in range(1,6):
        x=construct(e.transitive(a));z=e.state(x['ready1'],x['ready1_root']);zz=e.state(x['ready2'],x['ready2_root'])
        e.need(z==(468*3**a+170*2**a-53,90*2**a+2,90*2**a+3,468*3**a+170*2**a-52),'first ready TT family')
        e.need(zz==(890*3**a+382*2**a-123,102*2**a,102*2**a+1,890*3**a+382*2**a-122),'second ready TT family')
        family.append(dict(bulk_size=a,source_H=e.state(x['source'],0)[0],first_ready=z,second_ready=zz,offset=2*z[0]-zz[0]));e.tables.cache_clear()
    result=dict(status='PASS',first_parent_quotient=q1,second_parent_quotient=q2,first_coefficients=c1,second_coefficients=c2,arbitrary_bulk_controls=controls,family=family,scope='Two physically composable repairs for every nonempty bulk B; defectH(B), unit-ready iffBtransitive; no third-stage closure asserted')
    (HERE/('composed_family.json' if __debug__ else 'composed_family.optimized.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(status='PASS',controls=controls,first_coefficients=c1,second_coefficients=c2,first_two=[(x['source_H'],x['first_ready'][0],x['second_ready'][0]) for x in family[:2]])))
if __name__=='__main__':main()
