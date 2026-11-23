"""Classify all persistent parity faces of the nine-vertex switching fixture."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import lcm
from pathlib import Path
import subprocess
import sys

RESULTS=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(RESULTS/"chanpak_mod24_moonshine_ports_20261122"))
from check import CODE, values_at
from nine_ribbon import read


def faces(z,d):
    h,q,c,t,e,r,k=z
    dh,dq,dc,dt,de,dr,dk=d
    assert dh==0
    for side in ("left","right"):
        if side=="right":
            c,t=t,c
            e,r=r,e
            dc,dt=dt,dc
            de,dr=dr,de
        for p,b,step in ((p,b,s) for p in (0,1) for b in (0,1) for s in (-2,2)):
            offset=dk+b*de+p*dr-step
            if dq==dc==0:
                if dt==offset==0:
                    yield side,p,b,step,"plane"
                continue
            y=Fraction(-dt,dq) if dq else Fraction(-offset,dc)
            if y.denominator!=1 or y<0 or y%2!=b:
                continue
            y=int(y)
            if dq*y+dt or dc*y+offset:
                continue
            yield side,p,y,step,(h,q*y+t,k+c*y+e*b+r*p+h*read(y)+step//2)


def classify(binary):
    results=[]
    count=0
    with subprocess.Popen([str(binary.resolve())],stdout=subprocess.PIPE,text=True) as proc:
        for line in proc.stdout:
            v=list(map(int,line.split()))
            assert len(v)==18
            count+=1
            results.extend((tuple(v[:4]),*f) for f in faces(v[4:11],v[11:]))
        if proc.wait():
            raise RuntimeError("nine_switch failed")
    assert count==23180 and len(results)==4
    assert all(f[-1]==(15,388,370) for f in results)
    for f in results:
        print(f)
    print(count,"equal deletion-height markings among",256*9*36)


def selectors():
    hs=[values_at(CODE,c,False) for c in range(256)]
    heights=set(hs)
    assert (len(heights),min(heights),max(heights))==(125,59,3255)
    def smallest_options(bank,targets):
        for k in range(1,len(bank)+1):
            good=[p for p in combinations(bank,k)
                  if len({v%lcm(*p) for v in targets})==len(targets)]
            if good:
                return good
        return []
    for bank in ((3,5,11,17,41),(11,19,43,67,163),(1,2,5,8,20)):
        print("smallest exact height selectors from",bank,smallest_options(bank,heights))
    print("least unrestricted injective modulus",next(m for m in range(125,3256) if
                                                      len({h%m for h in heights})==125))
    for cut,h in enumerate(hs):
        sums=Counter()
        for i,value in enumerate(hs):
            d=(cut^i).bit_count()
            sums[min(d,9-d)]+=value
        assert 9*h==7*sums[1]-5*sums[2]+3*sums[3]-sums[4]
    return heights


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--native",type=Path)
    opts=parser.parse_args()
    selectors()
    if opts.native:
        classify(opts.native)
