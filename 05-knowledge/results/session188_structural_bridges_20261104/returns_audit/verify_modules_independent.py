"""Independent arbitrary-module audit using literal covers and path enumeration."""
from itertools import combinations, permutations
from pathlib import Path
import hashlib,json
from verify_independent import BASE,OUTER,path_response,signature,require,SOURCE


def substitute(A):
    m=len(A);names=[None]*m+list(OUTER);rows=[0]*len(names)
    for i,j in combinations(range(len(names)),2):
        if i<m and j<m:direction=bool(A[i]>>j&1)
        else:
            oi=3 if i<m else names[i];oj=3 if j<m else names[j]
            direction=bool(BASE[oi]>>oj&1)
        if direction:rows[i]|=1<<j
        else:rows[j]|=1<<i
    return rows,{m+OUTER.index(u) for u in (1,5)}


def all_covers(A,k):
    n=len(A);out=[]
    for p in permutations(range(n)):
        for cuts in combinations(range(1,n),k-1):
            boundaries=(0,)+cuts+(n,)
            paths=tuple(p[boundaries[i]:boundaries[i+1]] for i in range(k))
            if all(all(A[x]>>y&1 for x,y in zip(P,P[1:])) for P in paths):out.append(paths)
    require(len(out)==len(set(out)),'duplicate literal path covers')
    return out


def run():
    prior=json.loads((SOURCE/'module_extension.receipt.json').read_text())
    for f,h in prior['source_sha256'].items():
        require(hashlib.sha256((SOURCE/f).read_bytes()).hexdigest()==h,'source hash '+f)
    require(hashlib.sha256((SOURCE/'PROOF.md').read_bytes()).hexdigest()==prior['proof_sha256'],
            'proof changed since producer module receipt')
    require((SOURCE/'module_extension.normal.stdout').read_bytes()==
            (SOURCE/'module_extension.optimized.stdout').read_bytes(),'module run mismatch')
    modules=masks=vertices=injections=0;triangle=None
    for n in (1,2,3,4):
        edges=list(combinations(range(n),2))
        for mask in range(1<<len(edges)):
            A=[0]*n
            for k,(u,v) in enumerate(edges):
                if mask>>k&1:A[u]|=1<<v
                else:A[v]|=1<<u
            covers={k:all_covers(A,k) if k<=n else [] for k in (1,2,3)}
            h,q2=len(covers[1]),len(covers[2]);B,W=substitute(A);R,_=path_response(B)
            require(R[0]==17*h+9*q2,'module Hamilton formula')
            require(signature(R,W)[4]==8*h+6*q2,'module external joint formula')
            one={}
            for u in range(n):
                E=[sum(C[0][-1]==u for C in covers[k]) for k in (1,2,3)]
                S=[sum(C[0][0]==u for C in covers[k]) for k in (1,2)]
                # Test actual injectivity, not just a cardinality inequality.
                images=set()
                for (P,) in covers[1]:
                    j=P.index(u)+1
                    target=(P,) if j==n else (P[:j],P[j:])
                    require(target in covers[len(target)],'Hamilton injection target invalid')
                    require(target[0][-1]==u and target not in images,'Hamilton injection collision')
                    images.add(target);injections+=1
                images=set()
                for C in covers[2]:
                    bit=0 if u in C[0] else 1;P,Q=C[bit],C[1-bit];j=P.index(u)+1
                    target=(P,Q) if j==len(P) else (P[:j],P[j:],Q)
                    key=(bit,target)
                    require(target in covers[len(target)] and target[0][-1]==u,'two-cover target invalid')
                    require(key not in images,'two-cover injection collision')
                    images.add(key);injections+=1
                E1,E2,E3=E;S1,S2=S
                c=9*h+3*q2-25*E1-44*E2-18*E3-S1-5*S2
                sig=signature(R,W|{u});require(sig[0]-sig[4]==c,'decorated coefficient formula')
                one[u]=c;vertices+=1
                if n>=2:
                    require(E2>=1 and S2>=1 and E3>=2**(n-1)-2,'universal lower bounds')
                    require(c<=-10-12*2**(n-1),'universal constant bound')
                if A==[2,4,1]:
                    require((h,q2,E1,E2,E3,S1,S2)==(3,6,1,2,2,1,2),'C3 decorated control')
                    require(c==-115,'C3 constant');triangle=c
            for ports in range(1,1<<n):
                U={u for u in range(n) if ports>>u&1};sig=signature(R,W|U);c=sig[0]-sig[4]
                require(all(c<=one[u] for u in U),'multiport monotonicity')
                if n>=4:require(c<=-106,'order-four threshold')
                if n>=2 and c>=-103:
                    for target in (-103,-101):
                        require(3*c<target<c and c-min(sig[2:4])<target,'small context reduction')
                masks+=1
            modules+=1
    out={'status':'PASS','labelled_modules':modules,'nonempty_masks':masks,
         'singleton_vertices':vertices,'explicit_injection_images_checked':injections,
         'triangle_singleton_constant':triangle,'producer_proof_sha256':prior['proof_sha256'],
         'method':'Literal module permutations with nonempty cuts; explicit injection collision checks; literal paths in substituted graph. No producer imports.'}
    print(json.dumps(out,indent=2,sort_keys=True))


if __name__=='__main__':run()
