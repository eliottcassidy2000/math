"""Exact finite inverse for the four generalized physical families, not all C2."""
from functools import lru_cache
from core import constants,glue,gadget,birth,require
from pair_families import scalar,parameter_limit,construct

@lru_cache(None)
def coefficients(L):return constants(L)

def source_bank(bound):
 """All strong start-one response quadruples up to bound, plus two transients."""
 if bound<1:return {}
 base=(1,1,1,1);recipes={base:()};queue=[base]
 for h,e,R,d in queue:
  for L in range(1,bound//R+1):
   g,Q,B=coefficients(L)
   out=(g*h-(L%2)*e+L*R,h-int(L==1)*e,Q*h+B*R,h)
   if out[0]<=bound and out not in recipes:
    recipes[out]=recipes[h,e,R,d]+(L,);queue.append(out)
 return recipes

def candidate_sources(target):
 """One recipe/tail per relevant quadruple, with all erased d values retained."""
 found={}
 for s,recipe in source_bank(target//7).items():
  if s==(1,0,2,1):continue # Canonically supplied by singleton plus one sink.
  D,E,U,d=s
  found.setdefault(s,(recipe,0))
  k=1
  while D+6*(U<<k)<=target:
   found.setdefault((D,0,U<<k,d),(recipe,k));k+=1
 return found

def support(bound):
 values={};inputs=candidate_sources(bound)
 for s,(recipe,k)in inputs.items():
  for kind in ('N','M','P','Q'):
   if s==(1,1,1,1)and kind!='P':continue
   for v in range(2 if kind=='N'else 1,parameter_limit(kind,s,bound)+1):
    h=scalar(kind,s,v)
    if h<=bound:
     values.setdefault(h,dict(family=kind,source=list(s),core_recipe=list(recipe),tail=k,parameter=v))
 return values,inputs

def inverse(target):
 """All numerical presentations; identical source quadruples use one witness."""
 out=[]
 if target<1:return out
 for s,(recipe,k)in candidate_sources(target).items():
  for kind in ('N','M','P','Q'):
   if s==(1,1,1,1)and kind!='P':continue
   for v in range(2 if kind=='N'else 1,parameter_limit(kind,s,target)+1):
    if scalar(kind,s,v)==target:
     out.append(dict(family=kind,source=list(s),core_recipe=list(recipe),tail=k,parameter=v))
 return out

def build_source(recipe,tail=0):
 T=[0];r=0
 for L in recipe:T,r=glue(T,{r},gadget(L))
 for _ in range(tail):T,_=birth(T,[r],0,1)
 return T,r

def realize(record):
 T,r=build_source(record['core_recipe'],record['tail'])
 return construct(record['family'],T,r,record['parameter'])
