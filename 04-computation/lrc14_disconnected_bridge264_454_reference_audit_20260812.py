#!/usr/bin/env python3
"""Independent exact audit of the raw 264--454 C bridge ledger."""
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec,spec_from_file_location
from math import gcd
from pathlib import Path
from random import Random
import ast

ENGINE=Path('/tmp/canonical_reference_engine_thm3352.py')
CONTEXTS=Path('/tmp/disconnected_head263_contexts.txt')
LEDGER=Path('/tmp/disconnected_bridge264_454_exact_scan.ledger')
CSOURCE=Path('/tmp/disconnected_bridge264_454_exact_scan.c')
EXPECTED_ENGINE='da941a4267147d5442be81ae81880742d2f6b901bfc1d20fb667822402a2950e'
EXPECTED_CONTEXT='efea6bd97522fe1c31a5a88ca9f3223f9e7a8c08e3be85c493e9f62fdfaf06e4'
EXPECTED_LEDGER='2865e79a9aed02b293329cc7ed219ca0f428a376175f9b4fbb73552b43127053'
EXPECTED_CSOURCE='0da8ff6e5a4caf7ee396265f0c50a1c6587f2adcf98343ca32f1c8712e40d04f'
EXPECTED_TASKS='e967030c1566b49a580692c712d0ab1def306db52e5f598009a6198074893f36'
TARGET=F(186636088362,58865718786875)

def req(x,d):
 if not x: raise RuntimeError(d)
def fh(p): return sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def load(n,p):
 s=spec_from_file_location(n,p);m=module_from_spec(s);s.loader.exec_module(m);return m

def main():
 for path,want in ((ENGINE,EXPECTED_ENGINE),(CONTEXTS,EXPECTED_CONTEXT),
                   (LEDGER,EXPECTED_LEDGER),(CSOURCE,EXPECTED_CSOURCE)):
  req(fh(path)==want,('hash',path,fh(path)))
 tree=ast.parse(Path(__file__).read_text(),filename=__file__)
 req(not any(isinstance(x,ast.Assert) for x in ast.walk(tree)),'assert node')
 M=load('bridge_reference_mass',ENGINE)
 ctx=tuple(tuple(map(int,x.split())) for x in CONTEXTS.read_text().splitlines())
 rows=[]
 for line in LEDGER.read_text().splitlines():
  p,q,n,d,L,j,e,f=map(int,line.split());rows.append((p,q,F(n,d),L,j,e,f))
 expected=tuple((p,q) for p in range(264,455) for q in range(p+1,8*p) if gcd(p,q)<=3)
 payload=''.join(f'{p} {q}\n' for p,q in expected).encode()
 req(sha256(payload).hexdigest()==EXPECTED_TASKS,'task digest')
 req(tuple(x[:2] for x in rows)==expected,('universe/order',len(rows)))
 req(all(x[2]>TARGET for x in rows),'threshold failure')
 # Exactness at every C-reported argmin.
 for p,q,v,L,j,e,f in rows:
  req(M.mass(L,j,e,p,f,q)==v,('argmin mismatch',p,q,v,L,j,e,f))
 # Independent global minima on 20 weakest plus deterministic 200 channels.
 ranked=sorted(rows,key=lambda x:x[2]); chosen={x[:2] for x in ranked[:20]}
 rng=Random(3354)
 while len(chosen)<220: chosen.add(rows[rng.randrange(len(rows))][:2])
 lookup={x[:2]:x for x in rows};controls=[]
 for p,q in sorted(chosen):
  best=min((M.mass(L,j,e,p,f,q),L,j,e,f) for L,j,e,f in ctx)
  claimed=lookup[p,q]
  req(best==(claimed[2],*claimed[3:]),('global min',p,q,best,claimed))
  controls.append(((p,q),best))
 weak=min(rows,key=lambda x:x[2])
 print('DISCONNECTED RAW BRIDGE 264--454 INDEPENDENT AUDIT')
 print('rows',len(rows),'contexts',len(ctx),'argmin_exact_checks',len(rows),'full_context_controls',len(controls),'control_mass_checks',len(controls)*len(ctx))
 print('failures',0,'weakest',weak,'margin',weak[2]-TARGET)
 print('task_semantic_sha256',sha256(payload).hexdigest())
 print('engine_sha256',fh(ENGINE),'context_sha256',fh(CONTEXTS),'ledger_sha256',fh(LEDGER),'c_source_sha256',fh(CSOURCE))
 print('control_semantic_sha256',sha256(repr(tuple(controls)).encode()).hexdigest(),'python_assert_nodes',0)

if __name__=='__main__': main()
