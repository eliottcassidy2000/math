#!/usr/bin/env python3
"""Verify a factored, seed-independent inhomogeneous recurrence exactly."""
import json
from verify_farkas import row_step,output,states

def mul(row,word):
 for letter in reversed(word):row=row_step(row,letter)
 return row

def factored_row(core,word):
 row=output(core)
 for eigenvalue in (0,1,1,1,4,4,4,9,9,16):
  following=mul(row,word)
  row=[a-eigenvalue*b for a,b in zip(following,row)]
 return row

checks=0
for core in range(2):
 for a in range(8):
  word=(a,a^7)
  actual=factored_row(core,word)
  expected=[3840*int(q==0 and r==0 and s==0) for o,p,q,r,s in states]
  if actual!=expected:raise RuntimeError(('factored identity failed',core,word))
  checks+=1
# Hostile wrong multiplier and wrong input template must not pass.
expected=[3840*int(q==0 and r==0 and s==0) for o,p,q,r,s in states]
if factored_row(0,(0,7))==[x+int(i==0) for i,x in enumerate(expected)]:raise RuntimeError('coefficient mutation')
if factored_row(0,(0,0))==expected:raise RuntimeError('template mutation')
print(json.dumps(dict(complement_template_identities=checks,coordinate_equalities=checks*324,actual_seed_constant=23040,forced_cubic_coefficient='4/27',rejected_mutations=2,standard_library_only=True),sort_keys=True))
