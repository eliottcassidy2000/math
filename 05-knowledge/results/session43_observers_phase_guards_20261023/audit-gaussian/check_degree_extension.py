import hashlib,json,pathlib
count=0
def need(ok,msg):
 global count
 if not ok:raise RuntimeError(msg)
 count+=1
# Find the FIRST noninteger-weight congruence monomial, independently by degree layers.
for L in range(2,121):
 for target in [-1,1]:
  first=None
  for degree in range(L):
   charges=[2*a-degree for a in range(degree+1)]
   bad=[q for q in charges if q!=target and (q-target)%L==0]
   if bad:
    first=degree;break
  need(first==L-1,'first alias degree')
# Simultaneous guards: construct degree layers, rather than iterating arbitrary charges.
for degree in range(75):
 for a in range(degree+1):
  q=2*a-degree
  for target in [-1,1]:
   simultaneous=(q-target)%8==0 and (q-target)%9==0
   need(simultaneous==((q-target)%72==0),'combined action')
   if degree<=70 and simultaneous:need(q==target,'low-degree lift')
need((-71-1)%8==0 and (-71-1)%9==0 and -71!=1,'boundary shear')
parent=pathlib.Path(__file__).resolve().parent.parent
root=next((p for p in (parent/'session43-root', parent) if (p/'FINITE_GUARDS_AND_DEGREE.md').is_file()), None)
if root is None: raise RuntimeError('Root degree source missing in workspace and publication layouts')
pins={n:hashlib.sha256((root/n).read_bytes()).hexdigest() for n in ['FINITE_GUARDS_AND_DEGREE.md','check_degree_guards.py']}
print(json.dumps({'status':'PASS','checks':count,'source_pins':pins},sort_keys=True,indent=2))
