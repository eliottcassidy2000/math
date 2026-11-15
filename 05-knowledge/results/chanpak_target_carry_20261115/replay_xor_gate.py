# Load the published predecessor packet (or its local working name).
from pathlib import Path
import sys
base=Path(__file__).resolve().parent.parent
for directory in ("chanpak_absolute_steering_20261115","chanpak_absolute_steering"):
    if (base/directory).is_dir():
        sys.path.insert(0,str(base/directory))
from random import Random
from replay_steering import p

def check(x):
    if not x:
        raise AssertionError('xor gate control failed')

def J(v,r):
    return (p(v)+p((1<<r)-1-v))//2

rng=Random(37117)
total=0
for r in range(1,11):
    for v in range(1<<r):
        lag=sum((1<<max(0,d-2))*sum(((v>>i)^(v>>(i+d)))&1 for i in range(r-d))
                for d in range(1,r))
        check(lag==J(v,r))
        if r>=2: check(lag%2==((v>>1)^(v>>(r-2)))&1)
        total+=1
gates=0
for t in range(1,30):
    for i in range(50):
        j=t+2+rng.randrange(50)
        r=j+t+3+rng.randrange(50)
        v=rng.getrandbits(r)&~(1<<j)
        for d in range(1,t+2):
            b=(v>>(j-d))&1
            v=v&~(1<<(j+d)) | ((b^(d<t+1))<<(j+d))
        check((J(v|(1<<j),r)-J(v,r))%(1<<(t+1))==1<<t)
        gates+=1
print('PASS: XOR and parity controls:',total,'deterministic guarded-bit controls:',gates)
