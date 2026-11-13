"""Symmetry, budget topology, and exact endpoint contents on Paley controls."""
from common import *
from math import gcd,comb
from content_checks import child,hp

def factors(n):
    out=[];p=2
    while p*p<=n:
        e=0
        while n%p==0:n//=p;e+=1
        if e:out.append([p,e])
        p+=1
    if n>1:out.append([n,1])
    return out

def run():
    rows=[]
    for q in [3,7,11,19]:
        squares={i*i%q for i in range(1,q)}
        rr=[sum(1<<j for j in range(q) if (j-i)%q in squares) for i in range(q)]
        H,U,W,masks=tables(np.array(rr,dtype=np.int64));h=int(H);w=int(W[0,1]);u=int(U[0])
        need(all(int(W[i,j])==w for i in range(q) for j in range(i+1,q)),'pair-transitive constant weight')
        need(all(int(x)==u for x in U) and q*u==h+comb(q,2)*w,'vertex-transitive slope identity')
        r=next(k for k in range(2,q+1) if comb(k,2)*w>h-1);d=r-2
        sizes=[k for k in range(1,q+1) if h-comb(k,2)*w==1]
        need(len(masks)==sum(comb(q,k) for k in sizes),'unit sets depend only on cardinality')
        chi=sum((-1)**(k-1)*comb(q,k) for k in range(1,d+2))-1
        need(chi==(-1)**d*comb(q-1,d+1),'complete-skeleton Euler rank')
        extension=[h+(u-h//q)*k-w*comb(k,2) for k in range(q+1)]
        need(all(v==h+w*k*(q-k)//2 for k,v in enumerate(extension)),'endpoint subtraction in the extension formula')
        content=0
        for v in extension:content=gcd(content,v)
        need(content==gcd(h,w),'actual one-extension content')
        controls=0
        if q<=11:
            for k,want in enumerate(extension):
                need(hp(child(rr,(1<<k)-1))==want,'independent literal child Hamilton count')
                controls+=1
        rows.append({'q':q,'H':h,'pair_weight':w,'per_vertex_U':u,'defect_content':gcd(h,w),'one_extension_content':content,'extension_counts_by_new_outdegree':extension,'literal_child_controls':controls,'budget_dimension':d,'reduced_homology_rank':comb(q-1,d+1),'unit_selection_sizes':sizes,'H_factors':factors(h),'W_factors':factors(w),'U_factors':factors(u)})
    need(rows[2]['H']==5*7*11*13*19 and rows[2]['pair_weight']==2*3*11*19*23,'exact Paley11 prime identities')
    save('paley_resources',{'scope':'All-order structural reduction for any tournament with automorphism group transitive on unordered vertex pairs is proved in TOURNAMENTS.md. Literal controls are Paley tournaments3,7,11,19 only. Homology uses the complete-skeleton proof, not a homology enumeration. Extension content uses the inherited one-vertex formula.','records':rows})
    print(rows,flush=True)
if __name__=='__main__':run()
