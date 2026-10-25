import pathlib,hashlib,json,itertools,subprocess,sys
P=pathlib.Path(__file__).resolve().parent;S=P.parent/'padic-lifts'
def ck(x):
 if not x:raise RuntimeError('audit failure')
def magic(c,a,b):return (a,3*c-a-b,b,c-a+b,c,c+a-b,2*c-b,a+b-c,2*c-a)
out={}
for m in (2,3,4,8,9,16,27,32,64,128):
 Q={x*x%m for x in range(m)};count=0
 for c,a,b in itertools.product(Q,repeat=3):
  v=tuple(x%m for x in magic(c,a,b))
  if len(set(v))==9 and set(v)<=Q:count+=1
 ck(bool(count)==(m in (27,128)));out[str(m)]={'square_residues':len(Q),'labelled_magic_arrays':count}
L=(8,1,6,3,5,7,4,9,2)
for p,a,b,k0 in ((2,8,-7,7),(3,3,-2,3)):
 v=tuple(a*x+b for x in L)
 for k in range(k0,k0+5):
  m=p**k;Q={x*x%m for x in range(m)};ck(set(v)<=Q and len({x%m for x in v})==9)
# Exhaust every lift of fixed root vector (1,...,1) mod3 to mod9.
count=0
for r in itertools.product((1,4,7),repeat=9):
 v=tuple(x*x%9 for x in r)
 if v==tuple(x%9 for x in magic(v[4],v[0],v[2])):count+=1
ck(count==27);out['fixed_mod3_root_children_mod9']=count
out['source_sha256']={f:hashlib.sha256((S/f).read_bytes()).hexdigest() for f in ('PROOF.md','verify.py','certificate.json')}
print(json.dumps(out,sort_keys=True,indent=2))
