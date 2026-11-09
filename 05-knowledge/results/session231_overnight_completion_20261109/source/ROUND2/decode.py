"""Exact membership and literal graph constructor for the inherited mirrored fan subfamily."""
from math import isqrt

def f(q):
    total=0
    while q:
        u=q//2
        if q%2: total+=u+u%2
        q=u
    return total

def decode(N):
    if not isinstance(N,int) or N<5 or N%2==0:return None
    m=(N-2).bit_length()
    d=N-(1<<m)
    if d>1:return None
    if m%2==0:
        T=1-d;q=isqrt(T);l=q.bit_length()
        if q*q+2*f(q)+q%2!=T or m<2*l+2:return None
        beta=None
    else:
        T=(1-d)//2;q=isqrt(T);l=q.bit_length()
        if m<2*l+3:return None
        choices=[b for b in [0,1] if q*q+f(q)+(1-b)*(q+q%2)==T]
        if not choices:return None
        beta=choices[0]
    left=format((1<<l)-1-q,f'0{l}b') if l else ''
    middle=left+('' if beta is None else str(beta))+left[::-1]
    s=(m-len(middle))//2
    word='1'*s+middle+'1'*s+'0'
    return {'N':N,'exponent':m,'q':q,'center_bit':beta,'clock':s,'chain_word':word}

def graph(certificate):
    word=certificate['chain_word'];n=len(word)+1;out=[0]*n
    for i in range(1,n):
        for j in range(i+1,n):out[i]|=1<<j
        if word[i-1]=='1':out[0]|=1<<i
        else:out[i]|=1
    return out
