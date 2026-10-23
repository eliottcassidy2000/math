import json
checks=0
def ck(v,s):
 global checks
 checks+=1
 if not v:raise RuntimeError(s)
for L in range(2,101):
 for D in range(L-1):
  for a in range(D+1):
   for b in range(D-a+1):
    for target in(-1,1):
     if (a-b-target)%L==0:ck(a-b==target,'alias below bound')
 ck((-(L-1)-1)%L==0,'shear cyclic symmetry')
 ck(-(L-1)!=1,'shear continuous defect')
for q in range(-10000,10001):
 for target in(-1,1):ck(((q-target)%8==0 and(q-target)%9==0)==((q-target)%72==0),'CRT equivariance')
for a in range(71):
 for b in range(71-a):
  for target in(-1,1):
   if(a-b-target)%72==0:ck(a-b==target,'degree70')
ck((-71-1)%8==0 and(-71-1)%9==0,'degree71 shear')
print(json.dumps(dict(status='PASS',checks=checks,scope='exact cyclic equivariance plus bounded total degree; no arbitrary Keller claim'),indent=2))
