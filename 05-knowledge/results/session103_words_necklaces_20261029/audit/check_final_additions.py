"""Independent recomputation using recursive ordered-path growth and rational elimination."""
from check_independent import covers,need
from itertools import permutations
from fractions import Fraction
from pathlib import Path
import json

def det(A):
    a=[[Fraction(v) for v in row] for row in A];v=Fraction(1)
    for j in range(len(a)):
        p=next((i for i in range(j,len(a)) if a[i][j]),None)
        if p is None:return 0
        if p!=j:a[p],a[j]=a[j],a[p];v=-v
        q=a[j][j];v*=q
        for i in range(j+1,len(a)):
            r=a[i][j]/q
            a[i]=[x-r*y for x,y in zip(a[i],a[j])]
    need(v.denominator==1,'noninteger determinant');return int(v)

def observer(A):
    ports=list(permutations(range(len(A)),2));D=covers(A,2)
    return ports,[[D[s,t] for t in ports] for s in ports]

def hc(A):
    n=len(A);return sum(all(A[p[i]][p[(i+1)%n]] for i in range(n)) for p in ((0,)+q for q in permutations(range(1,n))))

def main():
    root=Path(__file__).resolve().parents[1]
    h=json.loads((root/'necklaces/results.json').read_text())['hostile'];out={}
    for name in ('singular_parent','unimodular_parent','singular_child','unimodular_child','repair_child'):
        A=h[name];ports,M=observer(A);out[name]={'det':det(M),'hc':hc(A)}
        need(out[name]['hc']==1,('cycle',name))
        if name=='singular_child':
            idx={p:i for i,p in enumerate(ports)};vs=[]
            for terms in (((0,2,1),(1,2,-1)),((0,3,1),(1,3,-1)),((0,2,1),(2,5,1)),((0,3,1),(2,3,1),(3,5,1))):
                v=[0]*len(ports)
                for i,j,z in terms:v[idx[i,j]]+=z;v[idx[j,i]]-=z
                need(all(sum(a*b for a,b in zip(row,v))==0 for row in M),'kernel');vs.append(v)
            minor=[[v[idx[p]] for v in vs] for p in ((2,1),(3,1),(5,2),(5,3))]
            need(abs(det(minor))==1,'kernel saturation minor');out['kernel_minor_det']=det(minor)
    need([out[k]['det'] for k in ('singular_parent','unimodular_parent','singular_child','unimodular_child','repair_child')]==[1,1,0,1,-1],'det table')
    R5=[[int(i!=j and (j-i)%5 in (1,2)) for j in range(5)] for i in range(5)]
    need(hc(R5)==2,'exact R5 cycle count');out['R5_exact_hc']=2
    T=((1,1,1),(1,1,0),(1,0,0));v=(1,2,0);seq=[]
    for m in range(15):
        seq.append(v[0]);v=tuple(sum(v[i]*T[i][j] for i in range(3)) for j in range(3))
    need(seq[:9]==[1,3,7,16,36,81,182,409,919],'boundary coefficients')
    for m in range(3,len(seq)):need(seq[m]==2*seq[m-1]+seq[m-2]-seq[m-3],'boundary recurrence')
    out['marked_walk_coefficients']=seq
    (root/'audit/final_additions_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
