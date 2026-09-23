"""All candidates x = -p/3^j (j <= J, 1 <= |x| <= 1.6) lying in an exceptional class of Bad_M (thread dump),
classified exactly by the prover (collatz_procgen_20260922_dim_hostile_prover.py):
   |x| < 3/2 : prove(x)  -> 'hostile' (complete certificate) or 'descends' (explicit witness)
   |x| > 3/2 : counted only (route(x) in the prover module finds explicit descents for such points).
For a class c mod 2^M and a given j there is exactly one p in [0,2^M) with -p/3^j = c mod 2^M, so the
candidate list is exact; a random class hits the window [3^j, 1.6*3^j] with probability 0.6*3^j/2^M.
usage: python3 ..._dim_z13scan.py DUMPFILE M J [nodelimit]      (env JMIN=j0: only j0..J, output hostile_z13_extra_M*_J*.txt)
"""
import sys, time
sys.path.insert(0, __file__.rsplit('/',1)[0])
from fractions import Fraction as F
from collections import Counter, defaultdict
from collatz_procgen_20260922_dim_hostile_prover import prove, route

dump=sys.argv[1]; M=int(sys.argv[2]); J=int(sys.argv[3]); NL=int(sys.argv[4]) if len(sys.argv)>4 else 300000
import os
JMIN=int(os.environ.get('JMIN','0'))
MOD=1<<M
bad=[int(l.split()[0]) for l in open(dump)]
cands=defaultdict(set)
for j in range(JMIN,J+1):
    q=3**j
    for x in bad:
        p=(-x*q)%MOD
        if q<=p<=q*8//5 and p%2==1 and (j==0 or p%3):
            cands[j].add(F(-p,q))
hostile=[]; total=Counter()
print(f"# M={M} classes={len(bad)} J={J}; expected random hits for j: 0.6*3^j*classes/2^M")
for j in range(JMIN,J+1):
    c=Counter(); t0=time.time(); ex=[]
    for x in sorted(cands[j]):
        if -x<F(3,2):
            r=prove(x,nodelimit=NL)
            c[r[0]]+=1
            if r[0]=='hostile': hostile.append(x)
            elif r[0]=='undecided': ex.append(str(x))
        else:
            c['|x|>3/2 (not proved; see route())']+=1
    total+=c
    print(f"j={j:2d}: candidates {len(cands[j]):4d}  {dict(c)}  (random hits expected {0.6*3**j*len(bad)/MOD:.2e})  {time.time()-t0:.1f}s"
          + (f"  undecided: {ex[:3]}" if ex else ""), flush=True)
print("TOTAL",dict(total))
with open(dump.rsplit('/',1)[0]+(f'/hostile_z13_M{M}_J{J}.txt' if JMIN==0 else f'/hostile_z13_extra_M{M}_J{JMIN}-{J}.txt'),'w') as f:
    for x in sorted(hostile): f.write(f"{x}\n")
for L in (20,30,40,50):
    if L>M: continue
    proj={x%(1<<L) for x in bad}
    cov={(x.numerator*pow(x.denominator,-1,1<<L))%(1<<L) for x in hostile}
    print(f"level {L}: projected classes {len(proj)}; containing a PROVED hostile -p/3^j (j<={J}): {len(proj&cov)}")
