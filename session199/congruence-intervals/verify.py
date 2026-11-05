"""Standalone finite controls for fixed-length residue forcing. No dependencies."""
from functools import lru_cache
from hashlib import sha256
from itertools import product
from math import gcd, isqrt
from pathlib import Path
import json
import random

ROOT=Path(__file__).resolve().parent

def check(ok, message):
    if not ok:
        raise RuntimeError(message)

# Independent evaluation by stripping least significant bits, not word reading.
def f(n):
    value=1
    while n:
        x=n//2
        if n%2==0:
            value+=x+x%2
        n=x
    return value


def step(n,F,bit,M):
    return (2*n+bit)%M, (F+(n+n%2 if bit==0 else 0))%M


def read(n,F,word,M):
    for b in word:
        n,F=step(n,F,int(b),M)
    return n,F


def params(M,dyadic=True):
    check(M>=2 and M%2==0,'even modulus')
    k=0;q=M
    while q%2==0:
        k+=1;q//=2
    b=(M-1).bit_length()
    a=2
    while pow(2,a,q)!=1%q:
        a+=1
    T=(M//2-1)*(a+k)
    compressed=(q==1 and dyadic)
    if compressed:
        T=0 if k==1 else (5*k*k-3*k-6)//2
    return {'M':M,'k':k,'q':q,'b':b,'a':a,'T':T,'L':2*b+T,
            'compressed':compressed}


def suffix(n,F,u,z,M,dyadic=True):
    par=params(M,dyadic)
    k,b,a,T=par['k'],par['b'],par['a'],par['T']
    check(F%2==1 and z%2==1,'odd f states')
    reset=format((-pow(2,b,M)*n)%M,f'0{b}b')
    n0,F0=read(n%M,F%M,reset,M)
    check(n0==0,'reset')
    r=(u-1)%M
    finish=format(r,f'0{b}b')
    need=(z-F0-f(r)+1)%M
    check(need%2==0,'even required increment')
    if par['compressed']:
        j=(-need//2)%(M//2)
        loops=[]
        for i in range(max(0,k-1)):
            if (j>>i)&1:
                loops.append('1'*(2 if i==0 else i+1)+'0'*k)
                if i>=1:
                    loops.append('11'+'0'*k)
    else:
        delta=(2**a-1)*(2**k-1)+1
        check(gcd(delta,M)==2,'loop gcd')
        j=0 if M==2 else ((need//2)*pow(delta//2,-1,M//2))%(M//2)
        loops=['1'*a+'0'*k]*j
    middle=''.join(loops)
    check(len(middle)<=T,'loop budget')
    word=reset+'0'*(T-len(middle))+middle+finish
    check(len(word)==par['L'],'fixed suffix length')
    check(read(n%M,F%M,word,M)==((u-1)%M,z%M),'suffix endpoint')
    return word


def realize(P,u,z,M,dyadic=True):
    check(P>=1,'positive prefix')
    word=suffix(P,f(P),u,z,M,dyadic)
    n=(P<<len(word))+int(word,2)
    check(n>>len(word)==P,'prefix retained')
    check((n+1)%M==u%M and f(n)%M==z%M,'integer endpoint')
    return n+1


def source(word):
    D,E,U,d=1,0,2,1
    for c in word:
        if c=='O':D,E,U,d=D,0,2*U,d
        elif c=='Z':D,E,U,d=D-E+U,d,2*U-d,d
        elif c=='J':D,E,U,d=D-E+U,D-E,D+U,D
        else:raise ValueError(c)
    return D,E,U,d


def word_for(U):
    check(U>=2,'nontrivial source')
    return ''.join('Z' if b=='0' else 'O' for b in bin(U-1)[3:])


def approximate(N,m):
    M=2*m
    par=params(M)
    B=1<<par['L']
    check(N%2==1 and N>=9*B,'approximation domain')
    Q=1
    while B*(2*Q)**2<=N:
        Q*=2
    P=N//(B*(Q+1))-1
    U=realize(P,1,(N-Q)%M,M)
    D=f(U-1)
    H=D+Q*U
    word=word_for(U)+'O'*(Q.bit_length()-2)+'OJ'
    actual=source(word)
    check(actual[0]==actual[2]==H,'balanced compiler output')
    check(0<=N-H and (N-H)%M==0,'lower bound and exact congruence')
    check((N-H)**2<9*B*N,'strict squared error bound')
    check(len(word)+2<=N.bit_length()+2,'vertex bound')
    check(word.count('J')==1,'one root change')
    check(D<=U,'physical D upper bound')
    check((Q+1)*U<=N and (Q+1)*(U+2*B)>N,'interval inequalities')
    return {'N':N,'m':m,'L':par['L'],'Q':Q,'U':U,'D':D,
            'H':H,'deficit':N-H,'vertices':len(word)+2,'word':word}


def graph(word):
    rows=[2,0];root=0
    for c in word:
        n=len(rows);rows.append(0)
        for i in range(n):
            if c!='O' and i==root:rows[n]|=1<<i
            else:rows[i]|=1<<n
        if c=='J':root=n
    return rows,root


def literal_counts(rows,root):
    n=len(rows);full=(1<<n)-1
    # Fix the cycle start at vertex 0 so each directed cycle is counted once.
    @lru_cache(None)
    def path(mask,last,cycle):
        if mask==full:
            return int(bool(rows[last]&1)) if cycle else 1
        total=0
        todo=rows[last]&~mask
        while todo:
            bit=todo&-todo;todo-=bit
            total+=path(mask|bit,bit.bit_length()-1,cycle)
        return total
    H=sum(path(1<<v,v,False) for v in range(n))
    C=path(1,0,True)
    starts=path(1<<root,root,False)
    return H,C,starts


def main():
    report={}
    count=0
    for M in range(2,15,2):
        for n in range(M):
            for F in range(1,M,2):
                for u in range(M):
                    for z in range(1,M,2):
                        suffix(n,F,u,z,M)
                        count+=1
    report['arbitrary_modular_start_and_target_checks']=count
    integers=0
    for M in range(2,21,2):
        for P in range(1,9):
            for u in range(M):
                for z in range(1,M,2):
                    U=realize(P,u,z,M)
                    check(source(word_for(U))==(f(U-1),U%2,U,1),'actual source integer')
                    integers+=1
    rng=random.Random(20261105)
    for k in range(1,15):
        M=1<<k
        for _ in range(20):
            P=rng.getrandbits(160)+1
            realize(P,rng.randrange(M),2*rng.randrange(M//2)+1,M)
            integers+=1
    for _ in range(200):
        M=2*rng.randrange(1,61)
        P=rng.getrandbits(80)+1
        realize(P,rng.randrange(M),2*rng.randrange(M//2)+1,M)
        integers+=1
    report['actual_integer_prefix_checks']=integers
    report['dyadic_suffix_lengths']={str(k):params(1<<k)['L'] for k in range(1,15)}
    for U in range(2,10001):
        check(source(word_for(U))==(f(U-1),U%2,U,1),'one apex dictionary')
        check(1<=f(U-1)<=U,'one apex cone')
    report['one_apex_sources_through_U']=10000
    # Both general and compressed constructors are tested at dyadic moduli.
    general=0
    for M in (2,4,8,16,32,64):
        for P in (1,2,3,17,100):
            for u,z in ((0,1),(1,3),(M-1,M-1)):
                realize(P,u,z,M,False)
                general+=1
    report['general_dyadic_constructor_checks']=general
    outputs=[]
    for N in range(37,10002,2):
        outputs.append(approximate(N,1))
    for m in range(1,31):
        B=1<<params(2*m)['L']
        for ratio in (9,10,15,16,17,25,64,100,1000,10**6,10**12):
            N=ratio*B+1
            outputs.append(approximate(N,m))
    report['approximation_checks']=len(outputs)
    report['sample_outputs']=[{k:v for k,v in row.items() if k!='word'} for row in outputs[:4]]
    report['sample_word']=outputs[100]
    literal=0
    for row in outputs:
        if row['vertices']<=11 and literal<30:
            H,C,starts=literal_counts(*graph(row['word']))
            check((H,C,starts)==(row['H'],1,1),'literal balanced graph counts')
            literal+=1
    report['literal_balanced_graphs']=literal
    # The old loop alone does not reset an arbitrary prefix to n=0.
    check(read(1,1,'110',6)[0]!=0,'hostile nonzero prefix control')
    report['hostile_unreset_prefix']={'M':6,'initial':[1,1],
        'word':'110','final':read(1,1,'110',6)}
    # A positive-density omission is compatible with the approximation bound.
    values={f(U-1)+2*U for U in range(2,501)}
    report['small_zero_J_missing_odds']=[n for n in range(3,100,2) if n not in values]
    report['proof_sha256']=sha256((ROOT/'PROOF.md').read_bytes()).hexdigest()
    report['verifier_sha256']=sha256(Path(__file__).read_bytes()).hexdigest()
    report['status']='PASS: all-order statements proved in text; finite controls only'
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
