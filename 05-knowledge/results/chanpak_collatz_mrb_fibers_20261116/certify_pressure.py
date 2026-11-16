"""Run from any working directory after eight_braid_root/lattice.py."""
import sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
for p in ['eight_braid_root','decimal_rank_root','chanpak_eight_braid_20261116','chanpak_decimal_rank_four_20261116']:
    sys.path.insert(0,str(root/p))
import numpy as np
from mixed_controller import load,calibrate,need
from check_return_growth import op,h,mul

def upper_integer(V,p,q,bits):
    D=10**12
    hist=np.bincount(V>>bits)
    rhs=D**q
    out=0
    for b in np.flatnonzero(hist):
        H=int(b)<<bits
        need(H>0)
        factor=H**p
        lo,hi=0,D
        while hi-lo>1:
            m=(hi+lo)//2
            if m**q*factor>=rhs: hi=m
            else:lo=m
        need((hi-1)**q*factor<rhs<=hi**q*factor)
        out+=int(hist[b])*hi
    return out,D

def main():
    D=load()
    A,_,_,_,U,_,_=calibrate(D)
    U=list({A[i]:i for i in U}.values())
    for name,J,p,q in [('all',range(168),37,43),('invertible',U,18,23)]:
        W=[tuple(D['words'][i])+(0,3,5,2,4) for i in J]
        M=np.stack([op(w) for w in W]);v=M[:,:,0].T
        f=h.reshape((1,-1))
        r=np.stack([mul(f,X)[0] for X in M])
        V=np.concatenate([mul(mul(r,X),v).ravel() for X in M])
        need(len(V)==len(W)**3 and np.min(V)==13460883)
        num,den=upper_integer(V,p,q,15)
        print(name,len(W),p,q,num,den)
        if not num<den: raise AssertionError("Pressure upper bound does not contract")
if __name__=='__main__':main()
