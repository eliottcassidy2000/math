#!/usr/bin/env python3
"""Exact dimension/trace controls for the completion theorem; not a Lie proof."""
import json
from pathlib import Path
original={0:133,1:56,2:1}
paired={-2:1,-1:56,0:133,1:56,2:1}
e8={-2:1,-1:56,0:134,1:56,2:1}
parabolic={0:134,1:56,2:1}
cotangent={j:original.get(j,0)+original.get(-j,0) for j in range(-2,3)}
assert sum(original.values())==190
assert sum(paired.values())==247
assert sum(e8.values())==248
assert sum(parabolic.values())==191
assert sum(cotangent.values())==380
trace=lambda profile:sum(k*v for k,v in profile.items())
assert trace(parabolic)==58
assert trace(e8)==trace(cotangent)==trace(paired)==0
# Hostile: zero grading trace does NOT imply a realizable Lie bracket or metric.
result={'original':original,'paired_but_impossible_under_theorem':paired,'minimum_realized_E8':e8,'nonunimodular_parabolic':parabolic,'cotangent_positive_control':cotangent,'parabolic_grading_trace':58,'logical_scope':'Numerical checks only; no structure constants or Jacobi verification. The proof excludes247 despite zero trace.'}
Path(__file__).with_name('completion_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
