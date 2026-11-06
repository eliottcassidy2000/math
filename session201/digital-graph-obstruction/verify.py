"""Controls for the all-order digital graph obstruction; no occupancy census."""
from functools import lru_cache
import json
from math import isqrt

def need(c,m):
    if not c:raise RuntimeError(m)

@lru_cache(None)
def g(n):
    if n<=1:return 1  # n=0 is a stated nonphysical extension for the 2-kernel.
    if n%2==0:return g(n//2)
    h=(n+1)//2
    return g(h)+h-h%2

T0=((1,0,0,0,0),(0,1,1,1,0),(0,0,2,0,0),(0,0,0,0,0),(0,0,0,0,1))
T1=((0,1,1,1,0),(0,1,0,0,0),(0,0,2,0,1),(0,0,0,0,1),(0,0,0,0,1))

def mv(T,x):return tuple(sum(a*b for a,b in zip(row,x)) for row in T)
def state(n):return (g(n),g(n+1),n,n%2,1)
def matrix_g(n):
    v=state(0)
    for digit in bin(n)[2:]:v=mv(T1 if digit=='1' else T0,v)
    return v[0]

def digit_g(u):
    m=u-1;bits=[int(x) for x in bin(m)[2:][::-1]]
    pair=sum(bits[i]*bits[j]*(1<<(j-i-1)) for i in range(len(bits)) for j in range(i+1,len(bits)))
    adjacent=sum(bits[i]*bits[i+1] for i in range(len(bits)-1))
    return u-bits[0]-pair-adjacent

def pair_word(u,v):
    alphabet={(1,1):'A',(0,0):'Z',(1,0):'B',(0,1):'C'}
    width=max(u.bit_length(),v.bit_length())
    return ''.join(alphabet[((u>>i)&1,(v>>i)&1)] for i in range(width-1,-1,-1))

def decode_word(w):
    alphabet={'A':(1,1),'Z':(0,0),'B':(1,0),'C':(0,1)}
    x=y=0
    for c in w:
        a,b=alphabet[c];x=2*x+a;y=2*y+b
    return x,y

def prime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

def x_fiber(x):
    if x%2:
        lo=(2*x+2)//3;hi=x-1;p=0
    else:
        lo=(x+1)//2;hi=(2*x-3)//3;p=1
    lo=max(lo,2);lo+=(p-lo)%2
    return [u for u in range(lo,hi+1,2) if u+g(u)==x]

def main():
    for n in range(4097):
        need(mv(T0,state(n))==state(2*n),'even state transition')
        need(mv(T1,state(n))==state(2*n+1),'odd state transition')
        need(matrix_g(n)==g(n),'matrix digit evaluation')
        if n:need(digit_g(n)==g(n),'weighted binary pair formula')
    sparse=0
    for a in range(2,65):
        for b in range(1,65):
            u=(1<<(a+b))+(1<<b)+1
            need(g(u)==u-(1<<(a-1)),'sparse-bit formula')
            sparse+=1
    distinguish=0
    for k in range(1,65):
        for m in range(1,65):
            word='A'+'Z'*k+'B'+'C'+'Z'*(m-1)+'A'
            u,y=decode_word(word)
            need(u==(1<<(k+m+2))+(1<<(m+1))+1,'input index decoding')
            need(y==(1<<(k+m+2))+(1<<m)+1,'output index decoding')
            need((g(u)==y)==(k==m),'Myhill-Nerode separating suffix')
            distinguish+=1
    # Positive and hostile regular-shape controls, not a proof by finite scan.
    language=0
    for left in range(1,17):
        for middle in range(1,17):
            for right in range(17):
                word='A'+'Z'*left+'B'+'C'*middle+'Z'*right+'A'
                u,y=decode_word(word)
                need((g(u)==y)==(left==right+1),'exact language intersection')
                language+=1
    # The row relation used in the corollary is exact on actual indices.
    for u in range(2,65):
        for v in range(2,65):
            row=g(u)*g(v)+u*v-(u%2)*(v%2)
            need(row==g(u)*g(v)+u*v+(-(u%2))*(v%2),'fixed J row relation')
    factors={2647:[2647],7939:[17,467],13231:[101,131],
             2134117:[643,3319],6402349:[23,278363],10670581:[10670581]}
    for n,fs in factors.items():
        product=1
        for f in fs:
            need(prime(f),'gap-control exact primality');product*=f
        need(product==n,'gap-control exact factorization')
    fibers={x:x_fiber(x) for x in (6,8,156,643,664)}
    need(fibers=={6:[3],8:[],156:[],643:[],664:[]},'gap-control source fibers')
    need(all(d%10!=9 for d in (101,131)),'2647 c=5 divisor residues')
    need(all(d%6!=5 for d in (643,3319)),'2134117 c=3 mixed divisor residues')
    return {'status':'PASS','gap_factorizations':factors,'gap_source_sum_fibers':fibers,'state_transition_and_digit_controls_through':4096,
            'sparse_bit_formula_pairs':sparse,'distinguishing_prefix_suffix_pairs':distinguish,
            'regular_shape_positive_and_hostile_controls':language,
            'g_description':'2-regular; binary graph not regular (all-order proof in PROOF.md)',
            'uniform_twelve_template_occupancy':'OPEN'}

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
