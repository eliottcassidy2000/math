import importlib.util,itertools
spec=importlib.util.spec_from_file_location('bank',str(__import__('pathlib').Path(__file__).resolve().parent / 'bank.py'))
b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
def strong_reg(w,ab):
 # ab = 1 iff 0 -> 1; types bits denote base -> apex
 a,c=(0,1) if ab else (1,0)
 if w and w[0]==3:return False
 if w and w[-1]==0:return False
 if any(all(x & (1<<c) for x in w[:s]) and all(not (x & (1<<a)) for x in w[s:]) for s in range(len(w)+1)):
  return False
 return True
for m in range(9):
 n=0
 for w in itertools.product(range(4),repeat=m):
  for ab in (0,1):
   if b.strong(b.rows(w,ab))!=strong_reg(w,ab): raise AssertionError((w,ab))
   n+=1
 print('m',m,'checked',n,flush=True)
