#!/usr/bin/env python3
"""Independent exact audit of the raw 455--678 C bridge ledger."""
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec,spec_from_file_location
from math import gcd
from pathlib import Path
from random import Random
import ast

ENGINE=Path('/tmp/canonical_reference_engine_thm3352.py')
CONTEXTS=Path('/tmp/disconnected_head263_contexts.txt')
LEDGER=Path('/tmp/disconnected_bridge455_678_exact_scan.ledger')
CSOURCE=Path('/tmp/disconnected_bridge455_678_exact_scan.c')
EXPECTED_ENGINE='da941a4267147d5442be81ae81880742d2f6b901bfc1d20fb667822402a2950e'
EXPECTED_CONTEXT='efea6bd97522fe1c31a5a88ca9f3223f9e7a8c08e3be85c493e9f62fdfaf06e4'
EXPECTED_LEDGER='f12b9cc843fd3d4dd8e300a2461df00721778f8665ab096d322a31dcc53b908a'
EXPECTED_CSOURCE='805b27801e1d4f12fa3ea1847d174b141ab8d5b07908e82ecc3fb965bc27ad7b'
EXPECTED_TASKS='3eea5e3450a6f188c9ac8574f57c8ef5f04e64595af148096df42cfe3bbd9561'
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
 expected=tuple((p,q) for p in range(455,679) for q in range(p+1,8*p) if gcd(p,q)<=3)
 payload=''.join(f'{p} {q}\n' for p,q in expected).encode()
 req(sha256(payload).hexdigest()==EXPECTED_TASKS,'task digest')
 req(tuple(x[:2] for x in rows)==expected,('universe/order',len(rows)))
 req(all(x[2]>TARGET for x in rows),'threshold failure')
 # Exactness/global minimality on the 20 weakest plus a deterministic spread
 # of 20 channels, each checked over every context by canonical reference.
 ranked=sorted(rows,key=lambda x:x[2]); chosen={x[:2] for x in ranked[:20]}
 rng=Random(3354)
 while len(chosen)<40: chosen.add(rows[rng.randrange(len(rows))][:2])
 lookup={x[:2]:x for x in rows};controls=[]
 for p,q in sorted(chosen):
  best=min((M.mass(L,j,e,p,f,q),L,j,e,f) for L,j,e,f in ctx)
  claimed=lookup[p,q]
  req(best==(claimed[2],*claimed[3:]),('global min',p,q,best,claimed))
  controls.append(((p,q),best))
 weak=min(rows,key=lambda x:x[2])
 print('DISCONNECTED RAW BRIDGE 455--678 INDEPENDENT AUDIT')
 print('rows',len(rows),'contexts',len(ctx),'full_context_controls',len(controls),'control_mass_checks',len(controls)*len(ctx))
 print('failures',0,'weakest',weak,'margin',weak[2]-TARGET)
 print('task_semantic_sha256',sha256(payload).hexdigest())
 print('engine_sha256',fh(ENGINE),'context_sha256',fh(CONTEXTS),'ledger_sha256',fh(LEDGER),'c_source_sha256',fh(CSOURCE))
 print('control_semantic_sha256',sha256(repr(tuple(controls)).encode()).hexdigest(),'python_assert_nodes',0)

if __name__=='__main__': main()
