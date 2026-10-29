"""Independent native-observer certificate for the explicitly marked order-five hostile."""
from native_kernel import *
from check_insertion_update import insertion_update
from pathlib import Path
import hashlib,json

def sparse(ports,v):
    return [{'port':list(p),'coefficient':int(z)} for p,z in zip(ports,v) if z]

def unit_pivots(M):
    a=[row[:] for row in M];r=0;trace=[]
    while r<len(a):
        pivot=next(((i,j) for i in range(r,len(a)) for j in range(r,len(a[0])) if abs(a[i][j])==1),None)
        if pivot is None:break
        i,j=pivot;trace.append([r,i,j,a[i][j]])
        a[r],a[i]=a[i],a[r]
        for row in a:row[r],row[j]=row[j],row[r]
        if a[r][r]<0:a[r]=[-v for v in a[r]]
        for i in range(r+1,len(a)):
            z=a[i][r]
            if z:a[i]=[x-z*y for x,y in zip(a[i],a[r])]
        r+=1
    need(all(v==0 for row in a[r:] for v in row[r:]),'nonzero residual after unit pivots')
    return r,trace

def main():
    path=Path(__file__).resolve().parents[1]/'necklaces/results.json'
    data=path.read_bytes();h=json.loads(data)['hostile'];inc=[0,1,1,1,1]
    records={}
    for name in ('singular','unimodular'):
        A=h[name+'_parent'];B=h[name+'_child']
        need(B==[row+[inc[u]] for u,row in enumerate(A)]+[[1-v for v in inc]+[0]],'different marked extension')
        for suffix,T in (('parent',A),('child',B)):
            ports,M=observer(T)
            need((ports,M)==direct_observer(T),'subset versus permutation mismatch')
            H=hp_table(T)[(1<<len(T))-1]
            hc=sum(H[0][t]*T[t][0] for t in range(len(T)))
            need(hc==1,'not unique HC')
            records[name+'_'+suffix]={'determinant':det(M),'hamilton_cycles':hc,'scores':sorted(map(sum,T))}
        need(insertion_update(A,inc)==observer(B),'universal update mismatch')
    A=h['singular_parent'];B=h['singular_child'];ports,M=observer(B);idx={p:i for i,p in enumerate(ports)}
    other=h['unimodular_parent']
    need(not any(all(A[i][j]==other[p[i]][p[j]] for i in range(5) for j in range(5)) for p in permutations(range(5))), 'parents unexpectedly isomorphic')
    need(records['singular_parent']['scores']==records['unimodular_parent']['scores']==[1,1,2,3,3],'score hostile lost')
    need([records[x]['determinant'] for x in ('singular_parent','unimodular_parent','singular_child','unimodular_child')]==[1,1,0,1],'determinant hostile lost')
    def qsum(terms):
        v=[0]*len(ports)
        for i,j,c in terms:v[idx[i,j]]+=c;v[idx[j,i]]-=c
        return v
    basis=[qsum(t) for t in [[(0,2,1),(1,2,-1)],[(0,3,1),(1,3,-1)],[(0,2,1),(2,5,1)],[(0,3,1),(2,3,1),(3,5,1)]]]
    for v in basis:need(all(sum(a*b for a,b in zip(row,v))==0 for row in M),'native kernel residual')
    need(len(rref_kernel(M))==4,'nullity mismatch')
    # Their projection onto four coordinates is unimodular. Together with
    # nullity four this proves they generate the full integer kernel lattice.
    free=[(2,1),(3,1),(5,2),(5,3)]
    need(abs(det([[v[idx[p]] for p in free] for v in basis]))==1,'kernel lattice not primitive')
    P,S=sectors(ports,M)
    need(det(P)==16 and len(rref_kernel(P))==0 and len(rref_kernel(S))==4,'sector location mismatch')
    old,new,N,R,Z=insertion_schur(A,inc);residual=rref_kernel(R);lifts=[]
    need(len(residual)==4,'Schur nullity mismatch')
    for v in residual:
        w=[-sum(a*b for a,b in zip(row,v)) for row in Z]
        full=[0]*len(ports)
        for p,z in zip(old+new,v+w):full[idx[p]]=z
        need(all(sum(a*b for a,b in zip(row,full))==0 for row in M),'Schur lift residual')
        need(all(z.denominator==1 for z in v+w),'nonintegral Schur basis')
        lifts.append({'residual':sparse(old,v),'new_port_correction':sparse(new,w)})
    rank,trace=unit_pivots(M);need(rank==26,'Smith unit rank')
    out={'status':'FINITE-EXACT independent certificate','input_sha256':hashlib.sha256(data).hexdigest(),'scope':'Explicit common marked insertion, not an identification with missing b/e grammar','records':records,'native_dimension':30,'rank':26,'nullity':4,'kernel_basis': [sparse(ports,v) for v in basis],'kernel_coordinate_convention':'right kernel, ordered END-pair columns; q_ij=e_(i,j)-e_(j,i)','sector':{'plus_determinant':16,'plus_nullity':0,'minus_nullity':4},'schur_lifts':lifts,'native_smith_form':{'1':26,'0':4},'unit_pivot_trace':trace}
    Path(__file__).with_name('marked_hostile_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('kernel_basis','schur_lifts','unit_pivot_trace')},indent=2))
if __name__=='__main__':main()
