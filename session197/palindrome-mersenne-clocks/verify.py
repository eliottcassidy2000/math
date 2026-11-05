import argparse
import hashlib
import json
import math
from pathlib import Path

def need(condition,message):
    if not condition:raise RuntimeError(message)

def f(n):
    value=1
    while n:
        x=n//2
        if n%2==0:value+=x+x%2
        n=x
    return value

def word_f(word):
    n=0;value=1
    for bit in word:
        if bit=='0':value+=n+n%2
        n=2*n+int(bit)
    return n,value

def word(c,p,q,t):
    return '1'+'0'*(p+t)+bin(c)[2:]+'0'*(q+t)

def parameters(c,p,q):
    ell=c.bit_length();reverse=int(bin(c)[2:][::-1],2)
    return 1<<(p+q+ell+1),(c<<(q+1))-(reverse<<p),c-f(c)-1

def count_formula(c,p,q,t):
    A,B,delta=parameters(c,p,q)
    return A*(4**t)+B*(2**t)-delta

def graph_from_word(w):
    n=int(w,2);s=len(w);b=s+1;T=[0]*(s+2)
    for i in range(s+1):
        for j in range(i+1,s+1):T[i]|=1<<j
        if i>0 and n>>(i-1)&1:T[i]|=1<<b
        else:T[b]|=1<<i
    return T

def graph(c,p,q,t):return graph_from_word(word(c,p,q,t))

def fixed_births(T):
    n=len(T);a=0;b=n-1;v=n-2
    X=[row|((1<<n) if i in(a,b) else 0) for i,row in enumerate(T)]
    X.append(((1<<n)-1)^(1<<a)^(1<<b))
    Y=[row|((1<<(n+1)) if i!=v else 0) for i,row in enumerate(X)]
    Y.append(1<<v)
    mapping=[0]+[i+1 for i in range(1,n-2)]+[n,n+1,1,n-1]
    return Y,mapping

def adjacency_map(A,B,p):
    n=len(A)
    return len(B)==n and sorted(p)==list(range(n)) and all(((A[i]>>j)&1)==((B[p[i]]>>p[j])&1) for i in range(n) for j in range(n))

def hp_cycle(T):
    n=len(T);size=1<<n;full=size-1
    def run(cycle):
        dp=[[0]*n for _ in range(size)]
        for v in ([n-1] if cycle else range(n)):dp[1<<v][v]=1
        for S in range(1,size):
            for v in range(n):
                if not dp[S][v]:continue
                out=T[v]&~S
                while out:
                    bit=out&-out;out-=bit;j=bit.bit_length()-1
                    dp[S|bit][j]+=dp[S][v]
        return sum(dp[full][v] for v in range(n) if not cycle or T[v]>>(n-1)&1)
    return run(False),run(True)

def prime_certificate(n):
    need(n>=2,'positive prime candidate')
    remainders=[[d,n%d] for d in range(2,math.isqrt(n)+1)]
    need(all(r for d,r in remainders),'complete trial division')
    return {'prime':n,'tested_every_integer_from':2,'through':math.isqrt(n),'remainders':remainders}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--witness');args=parser.parse_args()
    formulas=ports=arcs=literals=offsets=mersennes=order_checks=0;digest=hashlib.sha256()
    for c in range(1,128,2):
        ell=c.bit_length();reverse=int(bin(c)[2:][::-1],2)
        need(f(c)<= (c+1)//2,'odd-input bound')
        for p in range(1,5):
            for q in range(4):
                A,B,delta=parameters(c,p,q)
                need((B==0)==(p==q+1 and c==reverse),'exact cancellation criterion')
                for t in range(7):
                    w=word(c,p,q,t);n,g=word_f(w)
                    H=n+f(n)+n%2
                    need(g==f(n),'independent binary recurrence')
                    need(H==count_formula(c,p,q,t),'all two-block formula terms')
                    need(count_formula(c,p,q,t+1)-4*H==-2*B*(2**t)+3*delta,'affine-return boundary')
                    formulas+=1;digest.update(f'{c},{p},{q},{t},{H}\n'.encode())
                    if c<=31 and p<=3 and q<=2 and t<=3:
                        T=graph(c,p,q,t);born,mapping=fixed_births(T);target=graph(c,p,q,t+1)
                        need(adjacency_map(born,target,mapping),'literal fixed three-port return')
                        ports+=1;arcs+=len(born)*(len(born)-1)
                    if c<=15 and p<=2 and q<=1 and t<=1:
                        T=graph(c,p,q,t);hh,C=hp_cycle(T)
                        need(hh==H and C==1,'literal Hamilton paths and unique cycle')
                        literals+=1
    for delta in range(1,102,2):
        exact=[c for c in range(delta+2,2*delta+4,2) if bin(c)[2:]==bin(c)[2:][::-1] and c-f(c)-1==delta]
        larger=[c for c in range(1,10001,2) if bin(c)[2:]==bin(c)[2:][::-1] and c-f(c)-1==delta]
        need(exact==larger,'complete finite offset inverse control');offsets+=1
    need([c for c in range(1,10001,2) if bin(c)[2:]==bin(c)[2:][::-1] and c-f(c)-1==-1]==[1],'positive unit constant classification')
    need([c for c in range(1,10001,2) if bin(c)[2:]==bin(c)[2:][::-1] and c-f(c)-1==1]==[3,5],'negative unit constant classification')
    for k in list(range(4,257))+[511,512,1023,1024]:
        c=3 if k%2==0 else 5;t=(k-4)//2 if c==3 else(k-5)//2
        w=word(c,1,0,t);n,g=word_f(w)
        need(len(w)==k and n+g+n%2==(1<<k)-1,'exact Mersenne family')
        need(n+1==(1<<(k-1))+c*(1<<t)+1,'explicit L formula');mersennes+=1
    for L in range(1,1<<14):
        g=f(L-1)
        need(g<=L and (L%2 or g<=L//2),'all-order-bound numerical control')
        order_checks+=1
    for k in range(4,15):
        need(all(f(L-1)+L-L%2<(1<<k)-1 for L in range(2,(1<<(k-1))+1)),'complete smaller feedback-one orders')
    witnesses=[]
    for k in (4,5,7,9,13):
        c=3 if k%2==0 else 5;t=(k-4)//2 if c==3 else(k-5)//2
        w=word(c,1,0,t);T=graph_from_word(w);H,C=hp_cycle(T)
        need(H==(1<<k)-1 and C==1,'literal named Mersenne witness')
        witnesses.append({'exponent':k,'H':H,'C':C,'L':int(w,2)+1,'word':w,'rows':T,'apex':len(T)-1,'fixed_ports':{'a':0,'b':len(T)-1,'v':len(T)-2},'order':len(T)})
    # The source proof's global-order hostile: R5 is smaller at H15.
    need(hp_cycle([6,12,24,17,3])==(15,2),'unrestricted-order hostile')
    certificates=[prime_certificate(p) for p in (31,127,73,8191)]
    seven=[]
    for k in range(4,65):
        H=(1<<k)-1
        if H%7:continue
        j=k//3;a=(1<<j)-1;b=(1<<(2*j))+(1<<j)+1
        need(a*b==H and a>1 and b>1,'seven-prime factor split')
        if a==7 or b==7:seven.append((k,H))
    need(seven==[(9,511)],'unique 7p member control')
    artifact={'family':'T(L) from one-apex integer compiler','witnesses':witnesses,'prime_certificates':certificates,'scope':'Mersenne exact counts and prime examples; no infinite-primality or unrestricted order-minimality assertion.'}
    witness_text=json.dumps(artifact,sort_keys=True,indent=2)+'\n'
    if args.witness:Path(args.witness).write_text(witness_text)
    print(json.dumps({'status':'PASS','two_block_formula_controls':formulas,'fixed_port_returns':ports,'fixed_port_directed_incidences':arcs,'literal_small_H_C_pairs':literals,'finite_offset_inverse_controls':offsets,'exact_Mersenne_exponents':mersennes,'gadget_order_bound_controls':order_checks,'named_literal_Mersenne_graphs':len(witnesses),'prime_certificates':len(certificates),'seven_prime_member':511,'formula_digest':digest.hexdigest(),'witness_sha256':hashlib.sha256(witness_text.encode()).hexdigest(),'scope':'All-order proof with literal three-port construction; minimum order only among strong feedback-one tournaments.'},sort_keys=True,indent=2))

if __name__=='__main__':main()
