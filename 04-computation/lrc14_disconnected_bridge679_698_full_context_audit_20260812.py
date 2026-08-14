#!/usr/bin/env python3
"""Independent exact audit of the raw 679--698 bridge."""
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec,spec_from_file_location
from math import gcd
from pathlib import Path
from random import Random
import ast

ROOT=Path(__file__).resolve().parents[1]
ENGINE=ROOT/'04-computation/lrc_general_reflected_pair_mass_reference_audit_thm3352.py'
CONTEXTS=ROOT/'04-computation/lrc14_disconnected_head263_contexts_20260812.txt'
LEDGER=Path('/tmp/disconnected_bridge679_698_exact_scan.ledger')
CSOURCE=ROOT/'04-computation/lrc14_disconnected_bridge679_698_exact_scan_20260812.c'
SCANOUT=ROOT/'05-knowledge/results/lrc14_disconnected_bridge679_698_exact_scan_20260812.out'
TARGET=F(186636088362,58865718786875)
EXPECTED={'engine':'da941a4267147d5442be81ae81880742d2f6b901bfc1d20fb667822402a2950e',
'context':'efea6bd97522fe1c31a5a88ca9f3223f9e7a8c08e3be85c493e9f62fdfaf06e4',
'ledger':'0915c0c3a6eaf9fa4577e4449900865992d8fa646184e4ba6600801e89a55421',
'source':'3453ff5e6fbb9df7b4a47c74bce73867bf9efd4cf9f607fed5181cf2334946f7',
'scanout':'23aa26a7bee6a910df5970fd0a28e757d2fc2701d04ab3e23a1bb4b1764f9573',
'tasks':'9b2087d4a05e1ed453b05594d3c4d6ec748ff9c4f5cab95ea615d8238e002fa3'}
def req(x,d):
 if not x:raise RuntimeError(d)
def fh(p):return sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def load(n,p):s=spec_from_file_location(n,p);m=module_from_spec(s);s.loader.exec_module(m);return m

def main():
 for p,k in ((ENGINE,'engine'),(CONTEXTS,'context'),(LEDGER,'ledger'),(CSOURCE,'source'),(SCANOUT,'scanout')):req(fh(p)==EXPECTED[k],(k,fh(p)))
 req(not any(isinstance(x,ast.Assert) for x in ast.walk(ast.parse(Path(__file__).read_text()))),'assert')
 ctx=tuple(tuple(map(int,x.split())) for x in CONTEXTS.read_text().splitlines())
 rows=[]
 for line in LEDGER.read_text().splitlines():
  p,q,n,d,L,j,e,f=map(int,line.split());rows.append((p,q,F(n,d),L,j,e,f))
 tasks=tuple((p,q) for p in range(679,699) for q in range(p+1,8*p) if gcd(p,q)<=3)
 req(tuple(x[:2] for x in rows)==tasks,'universe/order')
 req(sha256(''.join(f'{p} {q}\n' for p,q in tasks).encode()).hexdigest()==EXPECTED['tasks'],'task digest')
 req(all(x[2]>TARGET for x in rows),'threshold')
 # Canonical physical engine checks all contexts for weakest 10 plus 10 spread.
 M=load('bridge679_reference',ENGINE);ranked=sorted(rows,key=lambda x:x[2]);chosen={x[:2] for x in ranked[:10]};rng=Random(3354)
 while len(chosen)<20:chosen.add(rows[rng.randrange(len(rows))][:2])
 lookup={x[:2]:x for x in rows};controls=[]
 for p,q in sorted(chosen):
  best=min((M.mass(L,j,e,p,f,q),L,j,e,f) for L,j,e,f in ctx);claim=lookup[p,q]
  req(best==(claim[2],*claim[3:]),('global min',p,q,best,claim));controls.append(((p,q),best))
 weak=min(rows,key=lambda x:x[2])
 print('DISCONNECTED RAW BRIDGE 679--698 INDEPENDENT AUDIT')
 print('rows',len(rows),'contexts',len(ctx),'full_context_controls',len(controls),'control_mass_checks',len(controls)*len(ctx),'failures',0)
 print('weakest',weak,'margin',weak[2]-TARGET)
 print('task_semantic_sha256',EXPECTED['tasks'],'ledger_sha256',fh(LEDGER))
 print('engine_sha256',fh(ENGINE),'context_sha256',fh(CONTEXTS),'c_source_sha256',fh(CSOURCE),'scan_output_sha256',fh(SCANOUT))
 print('control_semantic_sha256',sha256(repr(tuple(controls)).encode()).hexdigest(),'python_assert_nodes',0)
if __name__=='__main__':main()
