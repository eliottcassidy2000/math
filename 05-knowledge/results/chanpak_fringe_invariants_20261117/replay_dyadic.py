from itertools import product
from replay_joint_phase import direct

def coefficient(L,t,mod):
    return (pow(2,t+1,mod)*(pow(3,L-2-t,mod)+1))%mod

def compiler(a,H,P):
    """A physical one-apex word of odd length, profile (H,P) mod 2^a."""
    if a<1 or H%2==0 or P%gcd8(a)!=0:raise ValueError((a,H,P))
    mod=1<<a; P%=mod
    # Room for nonoverlapping fringes, H gadgets, and a zero guard after each.
    L=4*a*a+6*a+9
    digit=(P//8)%4
    target_length_mod4=(1 if digit in (0,3) else 3)
    L+= (target_length_mod4-L)%4
    bits=[0]*L; b=[0]*a; c=[0]*a; b[0]=c[0]=1
    if a>=2 and digit in (0,2) and a>3:b[1]=1
    for k in range(1,(a-2)//2):
        i=2*k
        if i+3>=a:break
        old=sum(x*y*coefficient(L,r+s,mod) for r,x in enumerate(b)
                for s,y in enumerate(c))%mod
        step=1<<(i+3)
        need=((P-old)//step)%min(4,mod//step)
        e=b if c[1]==0 else c
        correct= ((pow(3,L-2-i,16)+1)//4)%4
        vals=next((e0,x,y) for e0,x,y in product(range(2),repeat=3)
                  if (correct*e0+x+y-need)%min(4,mod//step)==0)
        e[i],b[i+1],c[i+1]=vals
    for i in range(a):
        bits[i]=b[i];bits[-i-1]=1-c[i]
    word=''.join(map(str,bits)); h,p=direct(word,mod)
    if p!=P:raise ValueError(('cover phase',a,H,P,h,p))
    delta=(h-H)%mod
    at=2*a+1
    for j in range(1,a):
        if delta>>j&1:
            bits[at]=bits[at+j+1]=1
        at+=2*a+2
    word=''.join(map(str,bits))
    if direct(word,mod)!=(H%mod,P) or not (word[0]=='1' and word[-1]=='0' and L%2):
        raise ValueError(('joint phase',a,H,P))
    return word

def gcd8(a):return 1<<min(a,3)

if __name__=='__main__':
    for a in range(1,8):
        k=0
        for h in range(1,1<<a,2):
            for p in range(0,1<<a,gcd8(a)):
                compiler(a,h,p); k+=1
        print(a,k)
    for a in (8,9,13,20):
        word=compiler(a,1+(1<<a)//2,8*((1<<a)//8-1))
        print("deeper",a,len(word))
