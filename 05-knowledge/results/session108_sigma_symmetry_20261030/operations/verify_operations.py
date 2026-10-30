"""Exact observer/update symmetry and aliquot reconstruction controls."""
from itertools import permutations,product
from math import gcd,isqrt,prod
from pathlib import Path
import hashlib,json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def sigma(n):
    total=0
    for d in range(1,isqrt(n)+1):
        if n%d==0:
            total+=d
            if d*d!=n:total+=n//d
    return total
def sf(n):
    result=1;p=2
    while p*p<=n:
        parity=0
        while n%p==0:n//=p;parity^=1
        if parity:result*=p
        p+=1
    return result*n

def main():
    primes=(3,5,19)
    ns=[prod(p for j,p in enumerate(primes) if x>>j&1) for x in range(8)]
    observed=[sf(sigma(n)) for n in ns]
    # Quotient coordinates are the two nonzero-prime selection bits.
    need(all((observed[x]==observed[y])==((x>>1)==(y>>1)) for x in range(8) for y in range(8)),'four fibers of size two')
    preservation=[];central=[];affine=[]
    for perm in permutations(range(8)):
        if not all(observed[perm[x]]==observed[x] for x in range(8)):continue
        preservation.append(perm)
        if all(perm[x^(1<<j)]==perm[x]^(1<<j) for x in range(8) for j in range(3)):central.append(perm)
        c=perm[0]
        if all(perm[x^y]^c==perm[x]^perm[y] for x in range(8) for y in range(8)):affine.append(perm)
    need(len(preservation)==16 and len(central)==2 and len(affine)==8,'different declared symmetry groups')
    need(set(central)=={tuple(x^k for x in range(8)) for k in (0,1)},'only uniform zero-bit toggle')
    hostile=lambda x:x^(((x>>1&1)*(x>>2&1)))
    need(all(observed[hostile(x)]==observed[x] for x in range(8)),'conditional toggle preserves observer')
    need(any(hostile(x^2)!=hostile(x)^2 for x in range(8)),'conditional toggle does not commute with named5')
    need(any(hostile(x^4)!=hostile(x)^4 for x in range(8)),'conditional toggle does not commute with named19')
    # Eight linear-syndrome states, exact sparse-factor realization.
    for x in range(8):
        need(observed[x]==(6 if x&2 else 1)*(5 if x&4 else 1),'sigma squareclass formula')
    # Neutral exponent grid: arbitrary observer symmetries versus named-port ones.
    grid=list(product((0,1,4),(0,3))); typed=0
    for perm in permutations(range(6)):
        if all((grid[i][0]==grid[j][0])==(grid[perm[i]][0]==grid[perm[j]][0]) and (grid[i][1]==grid[j][1])==(grid[perm[i]][1]==grid[perm[j]][1]) for i in range(6) for j in range(6)):
            typed+=1
    need(typed==12,'typed prime-port grid S3 x C2')
    probes=0;coprime_checks=0
    for a in range(2,40):
        sa=sigma(a)-a
        for n in range(1,500):
            if gcd(a,n)!=1:continue
            sn=sigma(n)-n;san=sigma(a*n)-a*n
            need(san==sa*n+sigma(a)*sn,'coprime aliquot affine identity')
            need((san-sigma(a)*sn)//sa==n and (san-sigma(a)*sn)%sa==0,'exact source recovery')
            probes+=1
            if a==3:
                need(sf(sigma(3*n))==sf(sigma(n)),'guarded neutral sigma attachment')
                coprime_checks+=1
    need(sigma(3*3)!=sigma(3)*sigma(3),'shared-prime guard hostile')
    need([sigma(n)-n for n in (3,5,19)]==[1,1,1],'aliquot prime collapse')
    need([sigma(n)-n for n in (15,57)]==[9,23],'same3 attachment reveals5 versus19')
    # Equal rank fails closure even under intrinsic source deletion.
    gl_a={0,1};gl_b={0,3}
    delete_first=lambda W:{x>>1 for x in W}
    need(len(gl_a)==len(gl_b)==2,'GL equal-rank inputs')
    need(len(delete_first(gl_a))==1 and len(delete_first(gl_b))==2,'GL source-deletion rank split')
    affine_a={0,6,7,1};affine_b={0,4,7,3}
    need(len(affine_a)==len(affine_b)==4,'affine equal-rank inputs')
    need(len(delete_first(affine_a))==2 and len(delete_first(affine_b))==4,'affine source-deletion rank split')
    return {'status':'PASS','prime_bank':primes,'integers':ns,'sigma_squareclasses':observed,'observer_symmetry_order':len(preservation),'affine_observer_symmetry_order':len(affine),'named_update_symmetry_order':len(central),'neutral_grid_observer_symmetry_order':720,'neutral_grid_prime_port_symmetry_order':typed,'coprime_probe_checks':probes,'free3_controls':coprime_checks,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Exact finite controls; general statements proved in OPERATIONS_AND_OBSERVERS.md'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
