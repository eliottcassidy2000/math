#!/usr/bin/env python3
"""Mutation controls: reject damaged finite universe, algebra, CM bound, and control."""
import json,subprocess,sys,tempfile
from pathlib import Path
source=Path(__file__).with_name('check.py').read_text()
cases=[
 ('missing_conductor', '72,84,120,168]', '72,84,120]', 'complete group screen'),
 ('damaged_rational_map', 'rel=J*(2-t)**2*(2+t)-64*(10+3*t)**3','rel=J*(2-t)**2*(2+t)-63*(10+3*t)**3','field substitution'),
 ('false_exponential_threshold', 'exp_lower>trace+2079','exp_lower>trace+10**30','exponential lower bound'),
 ('false_positive_control', 'J-8000,\'positive CM control\'', 'J-8001,\'positive CM control\'', 'positive CM control'),
]
records=[]
with tempfile.TemporaryDirectory() as td:
 for name,old,new,expected in cases:
  if source.count(old)!=1:raise RuntimeError('mutation target not unique: '+name)
  path=Path(td)/(name+'.py');path.write_text(source.replace(old,new))
  result=subprocess.run([sys.executable,*(['-O'] if not __debug__ else []),str(path)],capture_output=True,text=True,timeout=45)
  if result.returncode==0 or expected not in result.stderr:raise RuntimeError('failed to reject '+name+'\n'+result.stderr)
  records.append(dict(case=name,rejected=True,reason=expected))
print(json.dumps(records,indent=2))
