"""Exact periodic-shadow first-descent compiler. Standard library only."""
import json
from itertools import product
from pathlib import Path

def ck(ok,msg):
    if not ok: raise RuntimeError(msg)
def v2(x):return (x & -x).bit_length()-1
def affine(w):
    B=S=0
    for a in w:B=3*B+(1<<S);S+=a
    return 3**len(w),B,1<<S

def rotate_good(w):
    for i in range(len(w)):
        v=w[i:]+w[:i]
        if all(3**j>=(1<<sum(v[:j])) for j in range(1,len(v)+1)): return v
    raise RuntimeError('no growth rotation')
def compile_word(w):
    P,B,Q=affine(w)
    ck(Q>P,'contractive word required')
    M=2*Q;r=((Q-B)*pow(P,-1,M))%M
    # n > B/(Q-P), and n positive, ensures strict descent.
    threshold=B//(Q-P)+1
    t=max(0,(threshold-r+M-1)//M)
    return r,M,r+t*M

def first_descent(n,w):
    x=n;first=None
    for j,a in enumerate(w,1):
        actual=v2(3*x+1);ck(actual==a,'cylinder valuation')
        x=(3*x+1)>>a
        if first is None and x<n:first=j
    ck(first==len(w),'wrong first-descent index')
    return x

def main():
    cases=checks=0;examples=[]
    for L in range(1,6):
      for w0 in product(range(1,4),repeat=L):
        if 3**L<=1<<sum(w0):continue
        w=rotate_good(list(w0))
        for q in [1,2,5,20,100]:
            prefix=w*q;P,B,Q=affine(prefix)
            # First number of 2-steps making the multiplicative slope < 1.
            r=0;Pn=P;Qn=Q
            while Pn>=Qn:r+=1;Pn*=3;Qn*=4
            full=prefix+[2]*r
            residue,mod,n0=compile_word(full)
            for t in range(3):
                n=n0+t*mod;y=first_descent(n,full);checks+=len(full)
                ck(y<n,'descent');cases+=1
            # Single-reset barrier is exact at the affine slope level.
            b=1
            while (1<<(sum(prefix)+b))<=3**(len(prefix)+1):b+=1
            one=prefix+[b];_,_,n1=compile_word(one)
            first_descent(n1,one);checks+=len(one);cases+=1
            if w==[1,2] and q in [1,5,100]:
                examples.append({'word':w,'repetitions':q,'two_reset_length':r,'first_descent_index':len(full),'single_reset_minimum':b,'source_bits':n0.bit_length(),'source':str(n0),'target':str(first_descent(n0,full))})
    # Direct bounded-valuation special family and long exact stress controls.
    long=[]
    for q in [1,2,10,100,1000]:
        w=[1,2]*q;r=0
        while 3**(2*q+r)>=2**(3*q+2*r):r+=1
        full=w+[2]*r;res,mod,n=compile_word(full)
        P,B,Q=affine(full)
        ck(B==5*9**q*3**r+8**q*4**r-6*8**q*3**r,'closed carry')
        first_descent(n,full)
        ck(v2(n+5)==3*q+1,'S179 first mismatch precision')
        long.append({'q':q,'r':r,'first_descent':2*q+r,'source_bits':n.bit_length()})
    print(json.dumps({'status':'PASS','checked_sources':cases,'checked_steps':checks,'examples':examples,'long_controls':long},indent=2))
if __name__=='__main__':main()
