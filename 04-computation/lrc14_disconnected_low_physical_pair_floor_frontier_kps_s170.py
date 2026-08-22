#!/usr/bin/env python3
"""Exact finite scout and analytic constants for the disconnected-low pair floor.

This does not prove the universal floor.  It certifies the context atlas,
exhausts all high physical level pairs through ``--cap`` in exact arithmetic,
checks one low-channel zero hostile, and prints the proved large-ratio and
common-gcd reductions used in the accompanying reflection.
"""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd
from pathlib import Path
import argparse, os

ROOT=Path(os.environ.get('LRC_REPO',Path(__file__).resolve().parents[1])).resolve()
SELECTOR=ROOT/'04-computation/lrc14_j7_reflected_low_two_star_limit_selector_scout_20260812.py'
MASS_PATH=ROOT/'04-computation/lrc_general_reflected_pair_mass_thm3352.py'
EXPECTED_SELECTOR='cd8b08087f0f7e1e0c0c7d0be673629c7c2702c170c5c1e771e1d76df1d3cd1c'
EXPECTED_MASS='afd417297131401254769e1ef172d89c109ad2f9a843ea55e2badc3e7891435b'
EXPECTED_CONTEXT='9f7c7f24f81c409c09becfa921aa53387e5093710657bd3e5a0e935ecf4ea6c2'
EXPECTED_CAP80=(F(158,46397),3,5,336,174,12,3)
PAIR_TARGET=F(1,294)
DEBT_MAX=F(186636088362,11773143757375)
MASS=None;CONTEXTS=()

def req(x,d):
 if not x:raise RuntimeError(d)

def load(name,path):
 s=spec_from_file_location(name,path);req(s and s.loader,path)
 m=module_from_spec(s);s.loader.exec_module(m);return m

def init(contexts):
 global MASS,CONTEXTS
 CONTEXTS=contexts;MASS=load('pair_mass_'+str(os.getpid()),MASS_PATH)

def pair_min(pair):
 p,q=pair;best=None
 for L,cell,e,f in CONTEXTS:
  rec=(MASS.mass(L,cell,e,p,f,q),p,q,L,cell,e,f)
  if best is None or rec<best:best=rec
 return best

def contexts():
 sel=load('pair_selector',SELECTOR);out=set()
 for body,L in sel.MS.body_universe():
  cell,*_=sel.body_geometry(body,L)
  for e in body:
   for f in body:
    if e!=f:out.add((L,cell,e,f))
 return tuple(sorted(out))

def gamma(k):return F(1,2) if k%2==0 else F(k*k+1,2*k*k)

def midpoint_error(L,e,f,d,P,Q):
 C=Q*e-P*f
 return (gamma(P)*F(e*P,d*L*P-e)+gamma(Q)*F(f*Q,d*L*Q-f)
         +F(abs(C)*(abs(C)//L+1),2*d*d*L*P*Q))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--cap',type=int,default=80)
 ap.add_argument('--workers',type=int,default=10);a=ap.parse_args()
 req(sha256(SELECTOR.read_bytes().replace(b'\r\n',b'\n')).hexdigest()==EXPECTED_SELECTOR,'selector hash')
 req(sha256(MASS_PATH.read_bytes().replace(b'\r\n',b'\n')).hexdigest()==EXPECTED_MASS,'mass hash')
 ctx=contexts();req(len(ctx)==4044,len(ctx));sem=sha256()
 for x in ctx:sem.update((repr(x)+'\n').encode())
 req(sem.hexdigest()==EXPECTED_CONTEXT,sem.hexdigest())
 pairs=[]
 for p in range(1,a.cap+1):
  for q in range(p+1,a.cap+1):
   d=gcd(p,q)
   if p//d+q//d>=8:pairs.append((p,q))
 with ProcessPoolExecutor(max_workers=a.workers,initializer=init,initargs=(ctx,)) as pool:
  rows=tuple(pool.map(pair_min,pairs,chunksize=1))
 best=min(rows);below=sum(x[0]<PAIR_TARGET for x in rows)
 if a.cap==80:req(best==EXPECTED_CAP80,best)
 mass=load('pair_mass_control',MASS_PATH)
 hostile=mass.mass(168,90,1,10,3,2);req(hostile==0,hostile)

 # Large-ratio interval-discrepancy gate: if 5w>=36L(p+1), I>=1/294.
 # Inside its complement at d>=19, integrality gives Q/P<=29/4.  The two
 # shear summands are at most those of (P,Q)=(1,7), while the determinant
 # summand has the following ratio-only envelope.
 E19=midpoint_error(168,12,6,19,1,7)
 req(E19==F(673605911,139615031640),E19)
 req(E19<F(3,490),E19)
 shear19=F(1,265)+F(25,13027)
 det19=F(2945,45472*19*19)
 analytic19=shear19+det19
 req(shear19==F(19652,3452155),shear19)
 req(analytic19==F(2501969023,426078778720),analytic19)
 req(F(3,490)-analytic19==F(106676561,426078778720),analytic19)
 hunter_margin=5*PAIR_TARGET-DEBT_MAX
 req(hunter_margin==F(570672686921,494472037809750) and hunter_margin>0,hunter_margin)
 print('LRC14 DISCONNECTED-LOW PHYSICAL PAIR FLOOR FRONTIER')
 print('status=FINITE-EXACT scout + PROVED analytic reductions; universal pair floor OPEN')
 print('contexts',len(ctx),'context_semantic_sha256',sem.hexdigest())
 print('cap',a.cap,'high_pairs',len(pairs),'below_1_over_294',below,'minimum',best)
 print('low_channel_zero_control',(168,90,1,10,3,2),hostile)
 print('large_ratio_gate=5*(L*q-f)>=36*L*(p+1) implies mass>=1/294')
 print('bounded_ratio_common_gcd_gate=d>=19 implies mass>=1/294')
 print('d19_contextual_control_error',E19)
 print('d19_universal_analytic_error',analytic19,'error_margin_to_3_over_490',F(3,490)-analytic19)
 print('five_edge_Hunter_margin',hunter_margin)

if __name__=='__main__':main()
