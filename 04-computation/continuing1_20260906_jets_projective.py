"""Projective complete higher-jet precision and an intrinsic p31 residue packet.

No mathematical producer imports. Literal homogeneous observations are checked
against cardinal forms, integer Smith forms, and modular unit-pivot peeling.
"""
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import comb,gcd,lcm,prod
import hashlib,json,sys

sys.stdout.reconfigure(newline='\n')
GATES=0

def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise ArithmeticError(label)

def det(v,w):return v[0]*w[1]-v[1]*w[0]
def cgcd(xs):return reduce(gcd,xs,0)

def tangent(v):
    a,b=v
    if b==0:return (0,a)
    y=0 if abs(b)==1 else pow(a,-1,abs(b))
    return ((a*y-1)//b,y)

def mul(a,b,limit=None):
    size=len(a)+len(b)-1
    if limit is not None:size=min(size,limit)
    out=[0]*size
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<size:out[i+j]+=x*y
    return out

def power(a,n):
    out=[1]
    for _ in range(n):out=mul(out,a)
    return out

def vp(x,p):
    if x==0:return 10**9
    if isinstance(x,F):return vp(x.numerator,p)-vp(x.denominator,p)
    x=abs(x);n=0
    while x%p==0:n+=1;x//=p
    return n

def observer(vectors,mults,tangents=None):
    D=sum(mults)-1
    if tangents is None:tangents=[tangent(v) for v in vectors]
    rows=[]
    for (a,b),(u,v),m in zip(vectors,tangents,mults):
        ap=[a**j for j in range(D+1)];bp=[b**j for j in range(D+1)]
        up=[u**j for j in range(D+1)];vpowers=[v**j for j in range(D+1)]
        for r in range(m):
            row=[]
            for k in range(D+1):
                value=0
                for s in range(max(0,r-k),min(r,D-k)+1):
                    value+=comb(D-k,s)*comb(k,r-s)*ap[D-k-s]*up[s]*bp[k-r+s]*vpowers[r-s]
                row.append(value)
            rows.append(row)
    return rows

def reciprocal_packet(vectors,mults,tangents=None):
    if tangents is None:tangents=[tangent(v) for v in vectors]
    packets=[];global_den=1
    for i,(vi,wi,mi) in enumerate(zip(vectors,tangents,mults)):
        q=[1]
        for j,(vj,mj) in enumerate(zip(vectors,mults)):
            if i==j:continue
            q=mul(q,power([det(vi,vj),det(wi,vj)],mj),mi)
        q+= [0]*(mi-len(q))
        h=[F(1,q[0])]
        for r in range(1,mi):h.append(-sum(q[j]*h[r-j] for j in range(1,r+1))/q[0])
        den=lcm(*(x.denominator for x in h))
        # Independent integral numerator recurrence and common denominator.
        f=[1]
        for r in range(1,mi):
            f.append(-sum(q[j]*q[0]**(j-1)*f[r-j] for j in range(1,r+1)))
        numerators=[q[0]**(mi-1-r)*f[r] for r in range(mi)]
        exact=abs(q[0])**mi//gcd(abs(q[0])**mi,cgcd(abs(x) for x in numerators))
        need(den==exact,'integral reciprocal numerator recurrence')
        global_den=lcm(global_den,den)
        packets.append(dict(q=q,h=h,denominator=den))
    return packets,global_den

def cardinal_columns(vectors,mults,tangents,packets):
    M=sum(mults);out=[]
    for i,(vi,wi,mi) in enumerate(zip(vectors,tangents,mults)):
        Qpoly=[1]
        for j,(vj,mj) in enumerate(zip(vectors,mults)):
            if i!=j:Qpoly=mul(Qpoly,power([vj[1],-vj[0]],mj))
        need(cgcd(abs(x) for x in Qpoly)==1,'primitive homogeneous vanishing factor')
        s=[wi[1],-wi[0]];t=[-vi[1],vi[0]]
        for k in range(mi):
            H=[F(0)]*mi
            for r in range(mi-k):
                term=mul(power(s,mi-1-k-r),power(t,k+r))
                for j,c in enumerate(term):H[j]+=packets[i]['h'][r]*c
            phi=mul(Qpoly,H)
            need(len(phi)==M,'cardinal homogeneous degree')
            need(lcm(*(x.denominator for x in phi))==lcm(*(x.denominator for x in packets[i]['h'][:mi-k])),
                 'primitive multiplication attains the entire truncated denominator')
            out.append(phi)
    return out

def smith_p(A,p,bound):
    mod=p**(bound+1)
    B=[[(x.numerator*pow(x.denominator,-1,mod))%mod if isinstance(x,F) else int(x)%mod
        for x in row] for row in A]
    out=[];level=0
    while B:
        n=len(B);pivot=next(((i,j) for i in range(n) for j in range(n) if B[i][j]%p),None)
        if pivot is None:
            need(mod>p,'known precision bound prevents invisible residual block')
            B=[[x//p for x in row] for row in B];mod//=p;level+=1
            continue
        i,j=pivot;B[0],B[i]=B[i],B[0]
        for row in B:row[0],row[j]=row[j],row[0]
        inv=pow(B[0][0],-1,mod);top=B[0]
        B=[[(row[k]-row[0]*inv*top[k])%mod for k in range(1,n)] for row in B[1:]]
        out.append(level)
    return out

def precision_bound(vectors,mults,p):
    bounds=[]
    for i,mi in enumerate(mults):
        depths=[vp(det(vectors[i],vectors[j]),p) for j in range(len(vectors)) if i!=j]
        D=sum(mults[j]*vp(det(vectors[i],vectors[j]),p) for j in range(len(vectors)) if i!=j)
        bounds.append(D+(mi-1)*max(depths,default=0))
    return max(bounds,default=0)

def transformed(g,v):return (g[0][0]*v[0]+g[0][1]*v[1],g[1][0]*v[0]+g[1][1]*v[1])

def residue_class(v,p):
    return ('infinity',) if v[1]%p==0 else (v[0]*pow(v[1],-1,p)%p,)

def sidecar(vectors,z,p):
    v0,v1,v2=vectors
    e=vp(det(v0,v1),p)
    need(e>=1 and all(vp(det(v,w),p)==e for v,w in combinations(vectors,2)),
         'equilateral positive-depth projective triple')
    need(det(v0,z)%p!=0,'transverse frame is a p-adic unit chart')
    U=det(v0,v2)*det(v1,z)//p**e
    V=det(v0,v1)*det(v2,z)//p**e
    need(U%p and V%p and (U-V)%p,'bracket packet lies in the admissible torus')
    return U*pow(V,-1,p)%p

def main():
    from sympy import Matrix,ZZ
    from sympy.matrices.normalforms import smith_normal_form
    configurations=[]
    for m in range(1,6):configurations.append(([(1,0)],[m]))
    configurations += [
      ([(1,0),(0,1)],[3,2]),
      ([(1,0),(1,2)],[3,4]),
      ([(1,0),(1,4),(0,1)],[2,3,1]),
      ([(0,1),(4,1),(8,1)],[3,3,3]),
      ([(0,1),(4,1),(12,1)],[3,3,3]),
      ([(0,1),(8,1),(4,1)],[2,2,1]),
      ([(4,1),(12,1),(0,1)],[2,2,1]),
      ([(1,0),(0,1),(1,1),(1,2)],[3,1,2,1]),
      ([(-1,2),(3,1),(2,-3)],[2,3,2])]
    small_records=[];smith_checks=0;symmetry_checks=0;value_hostiles=0
    for vectors,mults in configurations:
        need(all(cgcd(v)==1 for v in vectors),'global primitive directions')
        need(all(det(v,w)!=0 for v,w in combinations(vectors,2)),'distinct rational projective directions')
        tangents=[tangent(v) for v in vectors]
        need(all(det(v,w)==1 for v,w in zip(vectors,tangents)),'integer unimodular tangents')
        A=observer(vectors,mults,tangents);packets,den=reciprocal_packet(vectors,mults,tangents)
        inverse=cardinal_columns(vectors,mults,tangents,packets);M=len(A)
        for i in range(M):
            for j in range(M):need(sum(x*y for x,y in zip(A[i],inverse[j]))==int(i==j),'literal observer times cardinal inverse')
        matrix=Matrix(A);sf=smith_normal_form(matrix,domain=ZZ)
        need(abs(int(sf[-1,-1]))==den,'independent global integer Smith largest factor')
        determinant=prod(abs(det(vectors[i],vectors[j]))**(mults[i]*mults[j]) for i,j in combinations(range(len(vectors)),2))
        need(abs(int(matrix.det()))==determinant,'weighted projective confluent determinant')
        for p in (2,3,5):
            exps=smith_p(A,p,precision_bound(vectors,mults,p))
            need(exps==[vp(int(sf[i,i]),p) for i in range(M)],'modular peeling versus full integer Smith')
            need(exps[-1]==vp(den,p),'reciprocal precision is sharp at the actual matrix')
            smith_checks+=1
            if exps[-1]>0:
                index=max(range(len(mults)),key=lambda i:vp(packets[i]['denominator'],p))
                offset=sum(mults[:index]);Dcolumn=packets[index]['denominator']
                scaled=[p*Dcolumn*x for x in inverse[offset]]
                need(all(x.denominator==1 for x in scaled),'value-only hostile is an actual integral coefficient vector')
                polynomial=[int(x) for x in scaled]
                need(min(vp(x,p) for x in polynomial)==1,'value-only integral hostile fails coefficient recovery modulo p^2')
                values=[sum(x*y for x,y in zip(row,polynomial)) for row in A]
                need(values[offset]==p*Dcolumn and all(x==0 for j,x in enumerate(values) if j!=offset),
                     'all derivative data are unchanged in the sharp precision hostile')
                need(vp(values[offset],p)==exps[-1]+1,'value perturbation is invisible at one-less observation digit')
                value_hostiles+=1
        # Each value cardinal alone contains its node's maximal denominator.
        offset=0
        for i,mi in enumerate(mults):
            need(lcm(*(x.denominator for x in inverse[offset]))==packets[i]['denominator'],
                 'value-only column attains complete local inverse precision')
            offset+=mi
        target=[abs(int(sf[i,i])) for i in range(M)]
        moved_tangents=[(w[0]+(i-2)*v[0],w[1]+(i-2)*v[1]) for i,(v,w) in enumerate(zip(vectors,tangents))]
        alternatives=[observer(vectors,mults,moved_tangents)]
        for g in (((1,2),(1,3)),((0,1),(1,0))):
            alternatives.append(observer([transformed(g,v) for v in vectors],mults))
        for alternative in alternatives:
            nf=smith_normal_form(Matrix(alternative),domain=ZZ)
            need([abs(int(nf[i,i])) for i in range(M)]==target,'full integer Smith covariance under tangent changes and GL2Z')
            symmetry_checks+=1
        small_records.append(dict(vectors=vectors,mults=mults,largest_factor=den))
    # Higher-multiplicity projective complete-residue controls.
    residue_records=[]
    for p in (2,3,5):
        for e in (1,2):
            vectors=[(1,0),(1,p**e)]+[(a,1) for a in range(p)]
            mults=[3,4]+[1]*p;A=observer(vectors,mults)
            full=smith_p(A,p,precision_bound(vectors,mults,p))
            pair=smith_p(observer(vectors[:2],[3,4]),p,6*e)
            need(len({residue_class(v,p) for v in vectors})==p+1,'all projective residue classes obstruct a common integral affine chart')
            need(full==sorted([0]*p+pair),'complete higher-jet projective residue splitting')
            packets,den=reciprocal_packet(vectors,mults)
            need(full[-1]==vp(den,p),'higher-jet largest denominator with all residue classes occupied')
            residue_records.append(dict(p=p,e=e,exponents=full))
    # Literal higher-jet hostiles inherited from THM-4443, transported to infinity.
    for xs,expected in [([0,4,8],18),([0,4,12],19)]:
        vectors=[(1,-x) for x in xs]
        A=observer(vectors,[3]*3)
        need(smith_p(A,2,24)[-1]==expected,'projective higher-jet metric-only precision hostile')
    A=observer([(1,0),(0,1)],[3,2],[(0,1),(-1,0)])
    bad=observer([(2,0),(0,1)],[3,2],[(0,1),(-1,0)])
    need(smith_p(A,2,5)==[0]*5 and smith_p(bad,2,5)==[0,0,2,3,4],
         'nonprimitive representative without its content changes higher-jet observer')
    partial=observer([(1,0),(0,1)],[2,2]);partial[1]=[0,0,1,0]
    need(Matrix(partial).det()==0,'same number of incomplete local jets need not be an invertible observer')
    local_vectors=[(1,0),(1,3),(0,1)];local_mults=[3,2,1]
    local_tangents=[tangent(v) for v in local_vectors]
    expected=smith_p(observer(local_vectors,local_mults),3,6)
    moved_vectors=[(2*v[0],v[1]) for v in local_vectors]
    moved_tangents=[(w[0],F(w[1],2)) for w in local_tangents]
    need(all(det(v,w)==1 for v,w in zip(moved_vectors,moved_tangents)),'local GL2 unit chart retains unimodular rational tangents')
    need(smith_p(observer(moved_vectors,local_mults,moved_tangents),3,6)==expected,
         'GL2Zp covariance includes locally primitive representatives and unit denominators')
    # p31 complete-class observer: new projective transport of inherited ideals.
    exceptional={3,11,15,17,21,29};p=31;p31_records=[]
    raw=[0]*15
    for j in range(15):
        c=comb(14+j,j)*comb(28-j,14-j)
        for k in range(j+1):raw[14-j+k]+=c*comb(j,k)*(-1)**(j-k)
    content=cgcd(abs(x) for x in raw);q=[x//content for x in raw]
    need(content==19380,'inherited odd packet reconstructed from its universal binomial identity')
    def binary(U,V):return sum(c*U**r*V**(14-r) for r,c in enumerate(q))
    need({a for a in range(2,31) if binary(a,1)%31==binary(1,a)%31==0}==exceptional,
         'homogeneous bracket-packet zero locus is exactly the inherited orbit')
    for a in (3,4,879,-20):
        e=1 if a in (3,4) else 2
        triple=[(1,0),(1,p**e),(1,p**e*a)]
        for z in ((0,1),(1,1),(-2,1),(31,1)):
            need(sidecar(triple,z,p)==a%p,'intrinsic residue independent of transverse choice')
            U=det(triple[0],triple[2])*det(triple[1],z)//p**e
            V=det(triple[0],triple[1])*det(triple[2],z)//p**e
            need(min(vp(binary(U,V),p),vp(binary(V,U),p))==int(a%p in exceptional),
                 'intrinsic homogeneous packet retains the common all-lift cancellation, not an individual value')
        for g in (((1,2),(1,3)),((0,1),(-1,0))):
            moved=[transformed(g,v) for v in triple]
            need(sidecar(moved,transformed(g,(0,1)),p)==a%p,'projective GL2 covariance of bracket residue')
        vectors=triple+[(r,1) for r in range(p)]
        mults=[16]*3+[1]*p
        full=smith_p(observer(vectors,mults),p,47*e)
        residual=full[p+16:];kappa=int(a%p in exceptional)
        need(full[:p+16]==[0]*(p+16),'full projective class splitting supplies the additional unit factors')
        need(sum(residual[:28])==588*e+2 and sum(residual[:29])==631*e+1+kappa and
             sum(residual[:30])==675*e+1,'intrinsic p31 ideal block at its correct full-module indices')
        need(residual[28:30]==[43*e-1+kappa,44*e-kappa],'two identified projective Smith factors')
        need(sum(full)==768*e,'full bracket determinant valuation')
        need(len({residue_class(v,p) for v in vectors})==p+1,'p31 control has no common integral affine chart')
        p31_records.append(dict(a=a,e=e,kappa=kappa,D75=sum(full[:75]),D76=sum(full[:76]),D77=sum(full[:77]),
                                pair=full[75:77],kernel43=sum(min(43,x) for x in full)))
    print('STATUS: PASS; projective arbitrary complete-jet precision and intrinsic p31 ideal sidecar')
    print('SMALL GLOBAL SMITH CONTROLS:',len(small_records),'; prime comparisons',smith_checks,
          '; full Smith symmetry checks',symmetry_checks,'; value-only sharp hostiles',value_hostiles)
    print('COMPLETE-RESIDUE HIGHER-JET CONTROLS:',json.dumps(residue_records,sort_keys=True))
    print('P31 FULL 79-DIMENSION CONTROLS:',json.dumps(p31_records,sort_keys=True))
    print('SCOPE: reciprocal jets retained; arbitrary higher-jet precision is not metric-only; p31 block transported, not a full partition law')
    print('SEMANTIC SHA256:',hashlib.sha256(json.dumps([small_records,residue_records,p31_records],sort_keys=True).encode()).hexdigest())
    print('ACTIVE GATES:',GATES)

if __name__=='__main__':main()
