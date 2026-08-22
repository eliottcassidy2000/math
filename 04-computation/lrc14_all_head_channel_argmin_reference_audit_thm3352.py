#!/usr/bin/env python3
"""Recompute every fast channel argmin and audit it with reference mass."""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import spec_from_file_location,module_from_spec
from pathlib import Path
import os

ROOT=Path(os.environ.get('LRC_REPO',Path(__file__).resolve().parents[1])).resolve()
TAIL=ROOT/'04-computation/lrc14_connected_low_uniform_high_forest_tail_thm3350.py'
FAST_PATH=ROOT/'04-computation/lrc_general_reflected_pair_mass_thm3352.py'
REF_PATH=ROOT/'04-computation/lrc_general_reflected_pair_mass_reference_audit_thm3352.py'
EXPECTED_TAIL='78daaf73966d283c0c0bafa1c0975684e6167d2ef6375a3abeece4e00cdc87f9'
EXPECTED_FAST='afd417297131401254769e1ef172d89c109ad2f9a843ea55e2badc3e7891435b'
EXPECTED_REF='da941a4267147d5442be81ae81880742d2f6b901bfc1d20fb667822402a2950e'
EXPECTED_CHANNEL='40472e196369baf0db300d96c028c2b7f78836bf8a06373b2170d4769ee9eff4'

def req(x,d):
 if not x:raise RuntimeError(d)

def load(name,p):
 s=spec_from_file_location(name,p);m=module_from_spec(s);s.loader.exec_module(m);return m

for path,want in ((TAIL,EXPECTED_TAIL),(FAST_PATH,EXPECTED_FAST),(REF_PATH,EXPECTED_REF)):
 req(sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()==want,('hash',path))
T=load('arg_tail',TAIL)
CONTEXTS=()

def init(ctx):
 global CONTEXTS,FAST,REF
 CONTEXTS=ctx
 FAST=load('fast'+str(os.getpid()),FAST_PATH)
 REF=load('ref'+str(os.getpid()),REF_PATH)

def channel_min_audit(ch):
 scale,d,P,Q=ch;best=None
 for L,c,e,f in CONTEXTS:
  v=FAST.mass(L,c,e,scale*d*P,f,scale*d*Q);rec=(v,L,c,e,f)
  if best is None or rec<best:best=rec
 v,L,c,e,f=best;w=REF.mass(L,c,e,scale*d*P,f,scale*d*Q)
 if v!=w:return ch,best,('reference_mismatch',w,v-w)
 return ch,best,None

def main():
 contexts=set();channels=set();heads=0
 for body,L in T.SEL.MS.body_universe():
  c,*_=T.SEL.body_geometry(body,L)
  for e in body:
   for f in body:
    if e!=f:contexts.add((L,c,e,f))
 for row in T.rows():
  if row in T.DENSE:continue
  S,*_=T.threshold(row);edges=T.high_edges(row)
  for scale in range(1,S):
   heads+=1;channels.update((scale,d,P,Q) for i,j,d,P,Q in edges)
 contexts=tuple(sorted(contexts));channels=tuple(sorted(channels));sem=sha256()
 req((len(contexts),heads,len(channels))==(4044,261254,4148),(len(contexts),heads,len(channels)))
 print('inventory',len(contexts),heads,len(channels),flush=True)
 with ProcessPoolExecutor(max_workers=3,initializer=init,initargs=(contexts,)) as pool:
  for ix,(ch,best,error) in enumerate(pool.map(channel_min_audit,channels,chunksize=1),1):
   if error:raise RuntimeError((ch,best,error))
   sem.update((repr((ch,best))+'\n').encode())
   if ix%100==0:print('progress',ix,'of',len(channels),flush=True)
 print('ALL-HEAD CHANNEL ARGMIN REFERENCE AUDIT')
 print('contexts',len(contexts),'heads',heads,'channels',len(channels),'argmins_checked',len(channels))
 print('channel_semantic_sha256',sem.hexdigest())
 req(sem.hexdigest()==EXPECTED_CHANNEL,sem.hexdigest())
 print('failures 0')

if __name__=='__main__':main()
