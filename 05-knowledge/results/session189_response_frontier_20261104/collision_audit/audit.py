#!/usr/bin/env python3
import hashlib,itertools,json,math,pathlib,sys
HERE=pathlib.Path(__file__).resolve().parent
SOURCE=pathlib.Path(sys.argv[1]).resolve() if len(sys.argv)>1 else HERE.parent/'session189-collision'
def need(b,msg):
    if not b:raise RuntimeError(msg)
def hc_end(rows):
    n=len(rows);N=1<<n;E=[[0]*n for _ in range(N)]
    for mask in range(1,N):
        for v in range(n):
            if not(mask>>v&1):continue
            rest=mask^(1<<v)
            E[mask][v]=sum(E[rest][u]for u in range(n)if rest>>u&1 and rows[u]>>v&1)if rest else 1
    return E
def joint(rows):
    n=len(rows);full=(1<<n)-1;E=hc_end(rows);S=hc_end([sum(1<<j for j in range(n)if rows[j]>>i&1)for i in range(n)])
    return[[sum(E[M][j]*S[full^M][k]for M in range(1,full))for k in range(n)]for j in range(n)]
def delv(rows,v):
    keep=[i for i in range(len(rows))if i!=v]
    return[sum(1<<j for j,b in enumerate(keep)if rows[a]>>b&1)for a in keep],{x:i for i,x in enumerate(keep)}
def shares(rows,m):
    ans=[]
    for i in range(m):
        q,ix=delv(rows,i);J=joint(q)
        for j in range(m):
            for k in range(j+1,m):
                if i not in(j,k):ans.append([i,j,k,J[ix[j]][ix[k]]+J[ix[k]][ix[j]]])
    return ans
def ambient(core,R,into):
    n=len(core);r=core[:]+[x<<n for x in R]
    for e in range(len(R)):
        for v in range(n):
            if into>>e&1:r[n+e]|=1<<v
            else:r[v]|=1<<(n+e)
    return r
def transitive(order):
    n=len(order);r=[0]*n
    for i,a in enumerate(order):
        for b in order[i+1:]:r[a]|=1<<b
    return r
def compositions(n,k):
    if k==1:yield(n,);return
    for a in range(n+1):
        for rest in compositions(n-a,k-1):yield(a,)+rest
def theta_literal(R,into,rmax):
    # Enumerate ordered path bins directly: prefix,suffix,excursion_1,... .
    e=len(R);ans=[0]*(rmax+1)
    def path(p):return all(R[a]>>b&1 for a,b in zip(p,p[1:]))
    for p in itertools.permutations(range(e)):
        for r in range(rmax+1):
            for lens in compositions(e,r+2):
                if any(z==0 for z in lens[2:]):continue
                bins=[];pos=0
                for z in lens:bins.append(p[pos:pos+z]);pos+=z
                if not all(path(q)for q in bins):continue
                pre,suf=bins[:2]
                if pre and not(into>>pre[-1]&1):continue
                if suf and into>>suf[0]&1:continue
                if any((into>>q[0]&1)or not(into>>q[-1]&1)for q in bins[2:]):continue
                ans[r]+=1
    return ans
def f(n,d,r):return(r+1)*sum((-1)**(r-s)*math.comb(r,s)*(s+1)**(n-d-1)*(s**(d-1)+(s+2)**(d-1))for s in range(r+1))
def converse(r):
    full=(1<<len(r))-1
    return[full^row^(1<<i)for i,row in enumerate(r)]
def inflate(Q,blocks):
    off=[];n=0
    for b in blocks:off.append(n);n+=len(b)
    r=[0]*n
    for i,b in enumerate(blocks):
        for v,row in enumerate(b):
            r[off[i]+v]|=row<<off[i]
            for j,c in enumerate(blocks):
                if Q[i]>>j&1:r[off[i]+v]|=((1<<len(c))-1)<<off[j]
    return r
def main():
    manifest=json.loads((HERE/'frozen_source_manifest.json').read_text())
    for name,z in manifest.items():need(hashlib.sha256((SOURCE/name).read_bytes()).hexdigest()==z['sha256'],'Frozen source changed: '+name)
    bank=json.loads((SOURCE/'bank.json').read_text());lines=(HERE/'literal_coefficients.txt').read_text().splitlines();need(len(lines)==len(bank['cases']),'Bank length')
    degree_counts={i:0 for i in range(5)}
    for case,line in zip(bank['cases'],lines):
        a=list(map(int,line.split()));L=len(case['masks'])
        need(a[:2]==[case['id'],L]and a[2:]==case['coefficients'],'Independent complete literal polynomial')
        for i,u in enumerate(case['masks']):
            for j,v in enumerate(case['masks']):degree_counts[u.bit_count()+v.bit_count()]+=1
        # Reconstruct every saved actual value from the independent coefficients.
        for i,u in enumerate(case['masks']):
            for j,v in enumerate(case['masks']):
                val=sum(a[2+ii*L+jj]for ii,uu in enumerate(case['masks'])if uu&u==uu for jj,vv in enumerate(case['masks'])if vv&v==vv)
                need(val==case['values'][i*L+j],'Independent full value reconstruction')
    need((HERE/'normal_receipt.json').read_bytes()==(HERE/'optimized_receipt.json').read_bytes()==(SOURCE/'receipt.json').read_bytes(),'Normal/optimized/frozen receipts mismatch')
    extra=[]
    for m,e in((5,4),(6,6)):
        R=transitive(list(range(e)));into=sum(1<<i for i in range(e//2,e));theta=theta_literal(R,into,m-3)
        need(theta[m-3]>0,'Higher excursion control failed to exercise top term')
        for order in(list(range(m)),list(reversed(range(m))),[0,2,1]+list(range(3,m))):
            rows=ambient(transitive(order),R,into);obs=shares(rows,m);pos={v:i for i,v in enumerate(order)}
            for i,j,k,value in obs:
                lo,hi=sorted((pos[j],pos[k]));d=hi-lo-int(lo<pos[i]<hi)
                need(value==sum(theta[r]*f(m-1,d,r)for r in range(m-2)),'Higher-excursion ambient transfer mismatch')
            extra.append({'m':m,'exterior_rows':R,'into':into,'theta':theta,'order':order,'rows':rows,'shared_quartics':obs})
    # Quotient-converse restoration on a genuinely nontransitive core;
    # target A has a nontrivial physical reorder, other children reverse.
    A=transitive(list(range(5)));Ap=transitive([0,2,1,3,4]);B=[2,4,1];C=[2,0];Q=[2,4,1]
    U=inflate(converse(Q),[Ap,converse(B),converse(C)])
    V=inflate(Q,[converse(Ap),B,C])
    restoration=[]
    for R,into in (([],0),([2,0],1),([2,4,1],3)):
        u=ambient(U,R,into);v=ambient(V,R,into);us=shares(u,5);vs=shares(v,5)
        need(us==vs,'All-order collision-restoration control')
        restoration.append({'exterior_rows':R,'into':into,'target_rows':u,'restored_rows':v,'common_shared_quartics':us})
    result={'status':'PASS_INDEPENDENT_REVIEW','frozen_source_manifest':manifest,'literal_bank_cases':len(lines),'complete_coefficients':sum(len(z['coefficients'])for z in bank['cases']),'coefficient_counts_by_total_degree':degree_counts,'normal_optimized_frozen_equal':True,'higher_excursion_controls':extra,'quotient_restoration_controls':restoration,'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in(HERE/'audit.py',HERE/'literal_bank.cpp')}}
    (HERE/'audit_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: corrected r+1; all 38888 literal coefficients; normal/-O; nonzero theta_2/theta_3 ambients; quotient-restoration controls.')
if __name__=='__main__':main()
