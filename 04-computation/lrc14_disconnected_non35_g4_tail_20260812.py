#!/usr/bin/env python3
"""Exact scratch proof: all small-ruler non-3:5 high channels close for dilation g>=4.

Uses projective phase floor 1/49-12/(49PQ), THM-3350 midpoint error,
a one-variable monotone envelope for P>=N_g, and exact error enumeration below N_g.
Every residual error-bound failure is then checked by the literal mass engine.
"""
from fractions import Fraction as F
from importlib.util import spec_from_file_location,module_from_spec
from pathlib import Path
from math import gcd
from hashlib import sha256
import argparse
ROOT=Path(__file__).resolve().parents[1]
# This focused sequel obtains the feasible-context universe directly from the
# canonical THM-3350 tail, avoiding THM-3352's much larger head-bank build,
# and loads the canonical exact mass engine.
TAIL=ROOT/'04-computation/lrc14_connected_low_uniform_high_forest_tail_thm3350.py'
MASS=ROOT/'04-computation/lrc_general_reflected_pair_mass_thm3352.py'
EXPECTED_TAIL='78daaf73966d283c0c0bafa1c0975684e6167d2ef6375a3abeece4e00cdc87f9'
EXPECTED_MASS='afd417297131401254769e1ef172d89c109ad2f9a843ea55e2badc3e7891435b'
TARGET=F(1,294); CUT={4:12,5:7,6:6,7:5,8:5,9:5,10:5,11:5,12:5,13:4,14:4,15:4}
def req(x,d):
 if not x:raise RuntimeError(d)
def load(n,p):s=spec_from_file_location(n,p);m=module_from_spec(s);s.loader.exec_module(m);return m
def filehash(p):return sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def gamma(k):return F(k*k+1,2*k*k)
def triangle(rad,z):
 b=rad.numerator//rad.denominator+3
 return sum((max(F(),rad-abs(z+n)) for n in range(-b,b+1)),F())
def correction(a,b):
 A=F((a+b)%14,14);B=F((b-a)%14,14);ev={F(),F(1)}
 for rad in (A,B):
  for n in range(-3,4):
   for z in (-F(n),rad-n,-rad-n):
    if 0<=z<=1:ev.add(z)
 return min((triangle(A,z)-A*A)-(triangle(B,z)-B*B) for z in ev)
def phase(P,Q):return max(F(1,105),F(1,49)-F(12,49*P*Q))
def error(P,Q,g,L,e,f):
 C=abs(Q*e-P*f)
 return gamma(P)*F(e*P,g*L*P-e)+gamma(Q)*F(f*Q,g*L*Q-f)+F(C*(C//L+1),2*g*g*L*P*Q)
def envelope(P,g,L,e,f):
 # For Q in [P,8P], gamma_Q <= gamma_P; Q/(gLQ-f)<=P/(gLP-f).
 # Also |Qe-Pf|<=P max(e+f,8e+f), and floor(x)+1<=(x+L)/L.
 det=max((e+f)**2,F((8*e+f)**2,8))
 return (gamma(P)*F(e*P,g*L*P-e)+gamma(P)*F(f*P,g*L*P-f)
         +F(det,2*g*g*L*L)+F(e+f,2*g*g*L*P))
def env_margin(P,g,t):
 L,e,f=t;return F(1,49)-F(12,49*P*(P+1))-envelope(P,g,L,e,f)-TARGET
def main():
 req(filehash(TAIL)==EXPECTED_TAIL,('tail hash',filehash(TAIL)))
 req(filehash(MASS)==EXPECTED_MASS,('mass hash',filehash(MASS)))
 T=load('non35t',TAIL);M=load('non35m',MASS)
 contexts=set()
 for body,L in T.SEL.MS.body_universe():
  cell,*_=T.SEL.body_geometry(body,L)
  if L<4592:
   for e in body:
    for f in body:
     if e!=f:contexts.add((L,cell,e,f))
 contexts=tuple(sorted(contexts));lanes={}
 for L,j,e,f in contexts:lanes[e,f]=min(lanes.get((e,f),10**99),L)
 lanes=tuple(sorted((L,e,f) for(e,f),L in lanes.items()));req(len(contexts)==2530 and len(lanes)==176,(len(contexts),len(lanes)))
 # Complete linked residue audit for Phi>=1/49-12/(49PQ).
 corr=tuple(correction(a,b) for a in range(14) for b in range(14))
 req(min(corr)==-F(12,49),(min(corr),))
 # Monotone envelope.  For real x>0,
 # d/dx [(x^2+1)/(2x^2) * ex/(gLx-e)]
 #   = -e(2gLx+ex^2-e)/(2x^2(gLx-e)^2) < 0.
 # The other nonconstant losses are c/x and 12/[49x(x+1)], also decreasing.
 # Hence the cleared envelope margin is strictly increasing for every real
 # x>=1; one exact base check at N_g proves the full integer tail.
 envrows=[]
 for g,N in CUT.items():
  for t in lanes:
   v=env_margin(N,g,t);req(v>0,('envelope base',g,N,t,v))
  envrows.append((g,N,min((env_margin(N,g,t),t) for t in lanes)))
 # g>=16: the exact midpoint bound has no residual channel (checked below at g=16, P<4); envelope handles P>=4 and decreases in g.
 for t in lanes:req(env_margin(4,16,t)>0,('g16 envelope',t,env_margin(4,16,t)))
 for P in range(1,4):
  for Q in range(P+1,8*P):
   if gcd(P,Q)>1 or P+Q<8 or (P,Q)==(3,5):continue
   req(phase(P,Q)-max(error(P,Q,16,*t) for t in lanes)>TARGET,('g16 finite',P,Q))
 # Exact finite residual under actual midpoint error, P<N_g. Physical check every context.
 residual=[];checks=0;weak=None
 for g,N in CUT.items():
  for P in range(1,N):
   for Q in range(P+1,8*P):
    if gcd(P,Q)>1 or P+Q<8 or (P,Q)==(3,5):continue
    E,t=max((error(P,Q,g,*t),t) for t in lanes)
    if phase(P,Q)-E<=TARGET:
     residual.append((g,P,Q,t,phase(P,Q)-E))
     for L,j,e,f in contexts:
      v=M.mass(L,j,e,g*P,f,g*Q);checks+=1
      row=(v-TARGET,g,P,Q,L,j,e,f,v)
      if weak is None or row<weak:weak=row
      req(v>=TARGET,('literal failure',row))
 req(len(residual)==118,('residual count',len(residual)))
 sem=sha256(repr((tuple(envrows),tuple(residual),checks,weak)).encode()).hexdigest()
 print('LRC14 DISCONNECTED NON35 SMALL-RULER G>=4 CLOSURE')
 print('contexts',len(contexts),'endpoint_lanes',len(lanes),'phase_correction_min',min(corr))
 print('envelope_cutoffs',tuple(envrows))
 print('midpoint_residual_channels',len(residual),'literal_checks',checks)
 print('weakest_literal_margin',weak)
 print('semantic_sha256',sem)
 print('conclusion=for all small-ruler contexts, all primitive high P<Q<8P except3:5, and all g>=4, physical overlap>=1/294')
if __name__=='__main__':main()
