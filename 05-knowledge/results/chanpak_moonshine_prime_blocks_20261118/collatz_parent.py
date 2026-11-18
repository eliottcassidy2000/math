"""Replay the decreasing physical carry parent and the actual C7 splice."""
from pathlib import Path
from math import isqrt
from fractions import Fraction
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent/'moonshine_prime_blocks'),str(HERE.parent/'dihedral_height')]
from rooted_pair import rooted,append,TT
from frame_replay import join,add
from sole_exit_reset import src,reset
from replay import profile

def parent(N,carry):
    """Unique penultimate-tail arithmetic inverse on its primitive interior."""
    h=carry
    if not 0<h<N:return None
    z=(N-h)//h
    if z<4:return None
    k=z.bit_length()-1
    u,r=divmod(N-h,1<<k)
    eta=u-h
    if r or not 0<eta<h:return None
    return h,eta,k

def rho(N,sgn):
    x=3*N+sgn
    return x//(x&-x)

def run(fixtures=False):
    assert [rho(n,-1) for n in (1,3,5,7,21,31,23)]==[1,1,7,5,31,23,17]
    assert [rho(n,-1) for n in (17,25,37,55,41,61,91)]==[25,37,55,41,61,91,17]
    assert [rho(n,1) for n in (1,7,11,17,13,5,21)]==[1,11,17,13,5,1,1]

    def step_trace(x,y):
        sg=1 if x%4==1 else -1
        n=3*x+sg
        a=(n&-n).bit_length()-1
        d=n>>a
        if d!=y:return d
        if a<6:return None
        w=rho(x,-sg); z=rho(w,-sg); v=rho(z,-sg)
        assert (w,z,v)==(y*2**(a-1)-sg,3*y*2**(a-3)-sg,9*y*2**(a-5)-sg)
        return v

    for y,forbidden,seeds in [(7,(7,21),(1,9,19,37,75)),(1,(1,3),(5,11))]:
        for x in range(1,10000,2):
            if x in forbidden or x in seeds:continue
            parent_x=step_trace(x,y)
            assert parent_x is not None and 0<parent_x<x and parent_x not in forbidden
    assert profile(join((2,4,1),(2,4,1)))[0]==9
    assert profile((32,1,35,7,15,26))[0]==19
    assert profile((32,17,3,7,13,30))[0]==37
    assert profile(join((2,4,1),join((6,12,8,1),(6,12,8,1))))[0]==75
    assert profile((2,4,33,23,39,75,31,127))[0]==35
    assert profile((8,17,11,50,69,23,47,127))[0]==77

    A,a=append((0,),0,TT(2),1)
    assert src(A,a)==(3,1,0)
    A,a=append(TT(2),0,TT(2),1);a-=1
    assert src(A,a)==(5,2,1)
    B=(2,4,33,23,39,75,31,127);b=0
    assert src(B,b)==(35,12,7)
    lab=(941,165); P=[]
    while lab!=(5,1):
        h,d,k=parent(*lab)
        P.append(k)
        lab=(h,d)
    assert P==[2,2,2]
    for k in P[::-1]:
        A,a=reset(A,a,k+1);B,b=reset(B,b,k+1)
        L=src(A,a);assert src(B,b)==tuple(7*x for x in L)
        if not fixtures: print('return through nonzero carry',L[0],7*L[0],L)
    assert L[0]==941 and all(941%d for d in range(2,isqrt(941)+1))

    T,a=(2,4,1),0
    for target in ((17,5,0),(97,29,0),(553,165,0)):
        T,b=append(T,a,(2,4,1),1);a=b-1
        assert src(T,a)==target
    assert profile(add((2,20,17,7,41,15,63,127),18))[0]==79
    if not fixtures:print('independent parents',79,553)

    # A second route starts with a rooted physical rotor pair.
    A=(16,385,131,327,14,159,55,345,1653,1279,255)
    B=(526,540,568,624,609,579,519,383,127,4480,1023,2047,3583)
    a,b=3,9
    if not fixtures:
        assert rooted(A,a)[:3][::2]==(5071,10842)
        assert rooted(B,b)[:3][::2]==(35497,37947)
    B=join(B,(0,))
    A,a=append(A,a,TT(2),1);B,b=append(B,b,TT(2),1)
    for _ in range(2):
        A=append(A,a,(0,),0)[0]
        B=append(B,b,(0,),0)[0]
    p=91807
    assert all(p%d for d in range(2,isqrt(p)+1))
    rot=lambda j:(j+1)%7 if j<7 else j
    assert all((B[i]>>j&1)==(B[rot(i)]>>rot(j)&1)
               for i in range(len(B)) for j in range(len(B)))
    if fixtures:
        for T,r in ((A,a),(B,b)):print(len(T),r,*T)
    else: print('rotational side output',p,7*p)
if __name__=='__main__':run('--fixtures' in sys.argv)
