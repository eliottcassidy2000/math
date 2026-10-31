"""Independent arithmetic audit of the necessary native-core grammar."""
from fractions import Fraction as Q
import json

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def floors(n,m,a,b):
    # Recursive lattice reflection, distinct from producer's iterative routine.
    if n<=0:return 0
    qa,a=divmod(a,m);qb,b=divmod(b,m)
    base=qa*n*(n-1)//2+qb*n
    if a==0:return base
    h=(a*(n-1)+b)//m
    if not h:return base
    # Count levels y=1,...,h: x>=ceil((m*y-b)/a).
    return base+h*n-floors(h,a,m,m+a-1-b)

def count_clock(N,M,d,b,r):
    a0=r if r else 8
    n=max(0,(N-a0)//8+1)
    return floors(n,M,8*d,d*a0+b-1)-floors(n,M,8*(d-1),(d-1)*a0+b)

def comp(F,G,q):
    # F after G.
    return (F[0]*G[0]%q,(F[0]*G[1]+F[1])%q)

def main():
    eta9=Q(2832031,781250);g=eta9-1;D=8*5**10;T=2*5**9
    need(D*(Q(21,8)-g)==25,'exceptional clock cap')
    small=0
    for A in range(1,26):
        for L in range(1,A):
            need(pow(3,L,5**10)!=pow(2,A,5**10),'small return survived');small+=1
    need(g>Q(20,8) and 2*Q(22,8)-1>3,'unique-core separation')
    need(Q(5,8)*T-sum(T//(2*5**(e-1)) for e in range(1,11))==Q(1,4),'all-one block charge')
    tests=0
    for n in range(0,15):
        for m in range(1,13):
            for a in range(-3,16):
                for b in [-20,-1,0,1,17]:
                    need(floors(n,m,a,b)==sum((a*i+b)//m for i in range(n)),'signed lattice reflection')
                    tests+=1
    q=5**14;M=4*5**13;d=2883827743
    need(pow(3,d,q)==2 and pow(3,M//2,q)!=1 and pow(3,M//5,q)!=1,'discrete log/order')
    h=pow(3*pow(2,-1,q)%q,T,q)
    need(pow(h,5,5**11)==1 and pow(h,1,5**11)!=1,'five-block first return')
    rows=[]
    for j,cap0 in [(23,48828467),(24,97656592)]:
        for a in range(6):
            jc=j+8-2*a
            if jc<22:continue
            cap=min(int(D*(Q(jc,8)-g)),cap0-a*T)
            residue=5*jc%8;b=a*T*(d-1)%M;target=pow(h,-a,q)
            count=count_clock(cap,M,d,b,residue)
            placements=[]
            for p in range(a+1):
                if p>=5 or a-p>=5:continue
                hp=pow(h,p,q);hs=pow(h,a-p,q)
                F=(target,(target-1)%q)
                need(comp((hs,(hs-1)%q),comp(F,(hp,(hp-1)%q),q),q)==(1,0),'affine inversion')
                placements.append(p)
            need(cap<M,'count range must fit one length modulus')
            def select(rank):
                lo,hi=1,cap
                while lo<hi:
                    mid=(lo+hi)//2
                    if count_clock(mid,M,d,b,residue)>=rank:hi=mid
                    else:lo=mid+1
                L=(d*lo+b)%M
                need(0<L<lo and lo%8==residue,'actual clock pair')
                need(pow(3,L,q)*pow(pow(2,lo,q),-1,q)%q==target,'actual multiplier')
                return [lo,L]
            rows.append(dict(j11=j,all_one10_blocks=a,j10=jc,core_clock_cap=cap,
                 core_clock_mod8=residue,prefix_block_counts=placements,required_core_map=[target,(target-1)%q],
                 pairs=count,first_pair=select(1) if count else None,last_pair=select(count) if count else None))
    need(sum(len(r['prefix_block_counts']) for r in rows)==34,'grammar count')
    need(sum(r['pairs'] for r in rows)==423868,'independent pair count')
    eta10=Q(183105426,48828125);lam12=Q(5,8)-Q(1,8*5**11);lam14=Q(5,8)-Q(1,8*5**13)
    need(Q(24,25)*Q(5,8)+Q(1,25)*lam12==lam14,'convex slopes')
    fee33=Q(24,25)*Q(33,8)+eta10/25
    fee32=Q(24,25)*4+eta10/25
    need(fee33>4 and fee32<4,'conditional certificate threshold')
    return dict(status='PASS',scope='conditional necessary grammar; heavy certificate not assumed',
        small_clock_pairs=small,signed_floor_controls=tests,grammar_forms=34,
        multiplier_pairs_without_placement=423868,patterns=rows,
        conditional_terminal_fee_if_33=str(fee33),insufficient_fee_if_32=str(fee32),
        inherited_certificates_rerun=False)

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
