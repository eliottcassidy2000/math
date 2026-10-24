import json
from fractions import Fraction
q=2**59;a=3**37;b=23*q-27*a;d=q-a

def need(x,s):
    if not x:raise RuntimeError(s)
def T(n):return (3*n+1)//2 if n%2 else n//2
def v2(n):
    n=abs(n)
    if not n:raise ValueError('infinite valuation')
    return (n&-n).bit_length()-1
def bits(n,length):
    w=[]
    for j in range(length):w.append(n%2);n=T(n)
    return w,n
def source_from_bits(w):
    # Incremental parity lift, independent of the closed affine inversion.
    r=0;cur=0;odd=0
    for j,e in enumerate(w):
        if cur%2!=e:r+=2**j;cur+=3**odd
        need(cur%2==e,'lift parity')
        cur=T(cur);odd+=e
    return r,cur

def run():
    w,end=bits(27,59)
    need(sum(w)==37 and end==23,'base word')
    need(8*d<b<9*d and Fraction(b,d).denominator==d,'noninteger fixedpoint')
    rows=[];checked=0;prior=0
    for k in range(1,33):
        Q=q**k;A=a**k;S=(Q-A)//d;B=b*S
        need((Q-A)%d==0,'geometric sum')
        r=(-B*pow(A,-1,Q))%Q;e=(A*r+B)//Q
        need(r>0 and r%q==27,'source class')
        need(r==b*pow(d,-1,Q)%Q,'fixed point cylinder')
        need(source_from_bits(w*k)==(r,e),'independent binary lift')
        need(k==1 or r%(q**(k-1))==prior,'nesting');prior=r
        need(A*pow(Q,-1,7)%7==(1 if k%2==0 else 6),'quotient slope')
        need(B*pow(Q,-1,7)%7==(0 if k%2==0 else 1),'quotient carry')
        for h in list(range(7))+[2**80+1]:
            n=r+Q*h;cur=n
            need(d*n-b>0,'real side')
            K=v2(d*n-b)//59
            need(K>=k,'fuel lower bound')
            for block in range(k):
                initial=cur;wc=[];below=[];ordinary=0
                for j in range(59):
                    bit=cur%2;wc.append(bit)
                    # Direct ordinary-step expansion, independent of T.
                    if bit:cur=3*cur+1;ordinary+=1
                    cur//=2;ordinary+=1
                    if cur<initial:below.append(j+1)
                need(wc==w and below==[59] and ordinary==96,'fresh block descent')
            need(cur==e+A*h,'affine endpoint')
            need(cur%7==((n if k%2==0 else 1-n)%7),'C2 residue')
            ww,_=bits(n,59*(K+1))
            need(ww[:59*K]==w*K and ww!=w*(K+1),'exact fuel')
            checked+=1
        if k<=4:rows.append({'k':k,'source':r,'endpoint':e,'max_blocks_at_source':v2(d*r-b)//59})
    # Negative controls: the original 27 is not itself in the twice repeated cylinder.
    need(v2(d*27-b)//59==1 and bits(27,118)[0]!=w*2,'27 fails repeat two')
    # A mod-7 return after two copies is not an integer return.
    r=rows[1]['source'];e=rows[1]['endpoint']
    need(r%7==e%7 and e<r and e!=r,'quotient not native involution')
    # Confirm fuel formula also where no full copy is available.
    for n in range(1,2049):
        K=v2(d*n-b)//59
        ww,_=bits(n,59*(K+1))
        need(ww[:59*K]==w*K and ww!=w*(K+1),'small-source exact fuel')
    # Exact finite interval counts: nested arithmetic progressions, with CRT sidecar.
    def apcount(X,r,M):return (X-r)//M-(-r)//M
    def residue(k):
        M=q**k
        return b*pow(d,-1,M)%M if k else 0
    for X in (0,26,27,28,100,2048):
        for k in range(4):
            rk=residue(k);rn=residue(k+1)
            predicted=apcount(X,rk,q**k)-apcount(X,rn,q**(k+1))
            need(predicted==sum(v2(d*n-b)//59==k for n in range(1,X+1)),'exact distribution count')
            for z in range(7):
                def joint(level):
                    M=q**level;r=residue(level)
                    R=r+M*((z-r)*pow(M,-1,7)%7)
                    return apcount(X,R,7*M)
                need(joint(k)-joint(k+1)==sum(v2(d*n-b)//59==k and n%7==z for n in range(1,X+1)),'joint CRT count')
    return {'parameters':{'q':q,'a':a,'b':b,'d':d},'fixedpoint':str(Fraction(b,d)),'base_word':''.join(map(str,w)),'base_odd_count':sum(w),'examples':rows,'checks':{'depths':32,'lifts_per_depth':8,'lifted_sources':checked,'small_sources':2048,'distribution_interval_bounds':[0,26,27,28,100,2048],'distribution_levels':[0,1,2,3],'joint_residues':7},'hostiles':{'27_max_full_blocks':1,'two_block_source':r,'two_block_endpoint':e},'status':'ALL CHECKS PASSED'}
if __name__=='__main__':print(json.dumps(run(),sort_keys=True,indent=2))
