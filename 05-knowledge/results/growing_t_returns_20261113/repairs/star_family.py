"""Three-arc repair family: exact run skeletons plus independent literal graph counts."""
from pathlib import Path
from itertools import product
import importlib.util,json
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('independent_old_counter',HERE.parent.parent/'changing_core_returns_20261113/preparation/check.py');e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
A=(518,388,456,1011,487,455,387,257,1,502)
OUTSIDE=(0,1,2,3,6,7,8,9)
PATCH=((0,2),(0,5),(0,7))
def source(B):
    n=8+len(B);T=[0]*n
    for i,u in enumerate(OUTSIDE):
        for j,v in enumerate(OUTSIDE):
            if A[u]>>v&1:T[i]|=1<<j
        for j in range(len(B)):
            if A[u]>>4&1:T[i]|=1<<(8+j)
            else:T[8+j]|=1<<i
    for i,row in enumerate(B):T[8+i]|=row<<8
    return tuple(T)
def repair(T):
    U=list(T)
    for a,b in PATCH:U[a]^=1<<b;U[b]^=1<<a
    return tuple(U)
def skeleton_table(Q,cap=4):
    # Outside vertices0..7 are distinct; repeated symbol8 is a nonempty module run.
    n=8;full=(1<<n)-1;end=[[[0]*9 for _ in range(1<<n)] for _ in range(cap+1)]
    for v in range(n):end[0][1<<v][v]=1
    end[1][0][8]=1
    for j in range(cap+1):
        for mask in range(1<<n):
            for last,value in enumerate(end[j][mask]):
                if not value:continue
                for v in range(n):
                    if not(mask>>v&1) and Q[last]>>v&1:end[j][mask|(1<<v)][v]+=value
                if j<cap and last!=8 and Q[last]>>8&1:end[j+1][mask][8]+=value
    H=[[sum(end[j][mask]) for mask in range(1<<n)] for j in range(cap+1)];H[0][0]=1
    return H,end
def coefficients(Q):
    H,E=skeleton_table(Q);_,S=skeleton_table(e.opposite(Q));full=255;out=[]
    for j in range(1,5):
        F=sum(E[a][mask][0]*H[j-a][full^mask] for a in range(j+1) for mask in range(256) if mask&1)
        out.append((H[j][full],S[j][full][0],H[j][full^1],F))
    return out
def K(B):
    h,p,t=e.profile(B);return h,p-2*h,t-3*p+3*h
C_SOURCE=((369,8,15,766),(256,3,6,630),(52,0,0,156))
C_PARENT=((85,1,15,591),(34,0,6,474),(0,0,0,114))
def linear(C,profile):return tuple(sum(C[j][i]*profile[j] for j in range(3)) for i in range(4))
def main():
    cs=coefficients(source((0,)));cp=coefficients(repair(source((0,))))
    e.need(tuple(cs[:3])==C_SOURCE and cs[3]==(0,0,0,0),'source skeleton certificate')
    e.need(tuple(cp[:3])==C_PARENT and cp[3]==(0,0,0,0),'repaired skeleton certificate')
    controls=0
    for n in range(1,4):
        for B in e.labelled(n):
            T=source(B);U=repair(T);k=K(B)
            e.need(e.state(T,0)==linear(C_SOURCE,k),'arbitrary module source')
            e.need(e.state(U,0)==linear(C_PARENT,k),'arbitrary module repaired parent')
            e.need(e.strong(T) and e.strong(U),'strong substitutions')
            W,r=e.birth(U,0,True);z=e.state(W,r)
            e.need(z[0]==675*k[0]+508*k[1]+114*k[2] and z[2]-z[1]==k[0],'newborn exact count and defect')
            controls+=1;e.tables.cache_clear()
    family=[]
    for m in range(1,7):
        T=source(e.transitive(m));U=repair(T);W,r=e.birth(U,0,True)
        a=e.state(T,0);b=e.state(U,0);c=e.state(W,r)
        expected_a=(52*3**m+100*2**m+13,3*2**m+2,6*2**m+3,156*3**m+162*2**m-26)
        expected_b=(34*2**m+17,1,6*2**m+3,114*3**m+132*2**m-15)
        e.need(a==expected_a and b==expected_b,'transitive family profiles')
        e.need(c[0]==114*3**m+166*2**m+1 and c[2]-c[1]==1 and e.strong(W),'literal unit-ready family')
        e.need(a[0]<c[0]<52*3**(m+1)+100*2**(m+1)+13,'output strictly between source clocks')
        family.append(dict(m=m,source_state=a,parent_state=b,output_state=c,ray_error=c[0]-2*a[0]+71))
        e.tables.cache_clear()
    result=dict(status='PASS',source_quotient=source((0,)),repaired_quotient=repair(source((0,))),outside_old_names=OUTSIDE,module_old_names=[4,5],patch_in_family_indices=PATCH,source_coefficients=cs,parent_coefficients=cp,arbitrary_module_controls=controls,family=family,scope='All nonempty module B proved by finite run-skeleton coefficients; unit-ready iff H(B)=1 iff B transitive')
    (HERE/('star_family.json' if __debug__ else 'star_family.optimized.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(status='PASS',arbitrary_module_controls=controls,family_targets=[x['output_state'][0] for x in family],source_coefficients=cs,parent_coefficients=cp)))
if __name__=='__main__':main()
