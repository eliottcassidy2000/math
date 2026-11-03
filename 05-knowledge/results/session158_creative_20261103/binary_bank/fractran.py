#!/usr/bin/env python3
"""Compile 1...0 words to ordered FRACTRAN programs; literal + path-DP tests."""
import itertools
import sys

def primes():
    n=7
    while True:
        if all(n%d for d in range(2,int(n**0.5)+1)): yield n
        n+=1

def compile_word(w):
    if len(w)<2 or w[0]!='1' or w[-1]!='0' or any(c not in '01' for c in w):
        raise ValueError('word must start 1 and end 0')
    P=primes(); states=[tuple(next(P) for _ in range(4)) for i in range(1,len(w)-1)]
    final_p, final_q=next(P),next(P)
    program=[]
    for k,(p,q,r,s) in enumerate(states):
        i=k+1; bit=int(w[i]); prev=int(w[i-1]); delta=int(prev==0 and bit==1)
        pnext=states[k+1][0] if k+1<len(states) else final_p
        # Ordered tests: transfer B→2T plus optional H, then restore T→B.
        program.extend([(q*3**2*5**(1-bit),p*2), (p,q), (r,p),
                        (s*2,r*3), (r,s), (pnext*2**bit*5**delta,r)])
    program.extend([(final_q*5,final_p*2), (final_p,final_q)])
    initial=(states[0][0] if states else final_p)*2*5**2
    return program,initial,final_p

def run(program,n,limit=1000000):
    for step in range(limit):
        for a,b in program:
            if (n*a)%b==0:
                n=n*a//b
                break
        else: return n,step
    raise RuntimeError('step limit exceeded')

def paths(w):
    n=len(w)+1; out=[0]*n
    for i,c in enumerate(w,1):
        if c=='1': out[0]|=1<<i
        else: out[i]|=1
        for j in range(i+1,n): out[i]|=1<<j
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for i in range(n):
            if not dp[mask][i]: continue
            bits=out[i]&~mask
            while bits:
                bit=bits&-bits; bits-=bit
                dp[mask|bit][bit.bit_length()-1]+=dp[mask][i]
    return sum(dp[-1])

def exponent5(result, control):
    if result%control: raise ValueError('missing final control')
    result//=control; e=0
    while result%5==0 and result:
        result//=5; e+=1
    if result!=1: raise ValueError('nonempty garbage registers')
    return e

def main():
    if len(sys.argv)>1:
        w=sys.argv[1]; prog,start,final=compile_word(w); result,steps=run(prog,start)
        print('word:',w,'initial:',start,'fractions:',', '.join(f'{a}/{b}' for a,b in prog))
        print('final integer:',result,'steps:',steps,'v_5:',exponent5(result,final),'path-DP:',paths(w))
        return
    total=0
    for L in range(2,9):
        for bits in itertools.product('01',repeat=L-2):
            w='1'+''.join(bits)+'0'; prog,start,final=compile_word(w)
            if len(prog)!=6*(L-2)+2: raise ValueError(('length',w))
            result,steps=run(prog,start)
            if exponent5(result,final)!=paths(w): raise ValueError(('mismatch',w))
            total+=1
    print('PASS:',total,'literal ordered FRACTRAN runs matched independent subset path DP')
    main_example='1010'; prog,start,final=compile_word(main_example)
    res,steps=run(prog,start)
    print('example:',main_example,'initial',start,'fractions',len(prog),'v5',exponent5(res,final),'steps',steps)
    print(' '.join(f'{a}/{b}' for a,b in prog))

if __name__=='__main__': main()
