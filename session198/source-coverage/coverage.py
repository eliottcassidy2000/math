"""Complete scalar spectra using physical recipes and proved pruning only."""
from source_oracle import source_bank
import c3

ONE=(1,1,1,1)
WEIGHT=(1,1,2,2)
BASIS=((1,0,1,0),(0,0,1,0),(2,1,2,1),(1,0,2,1))
LOWER={'A':(7,1),'B':(7,1),'C':(14,1),'J0':(6,1),'J1':(29,4)}
EXPECTED={
'A':((14,10,21,20),(10,7,16,15),(21,16,33,31),(20,15,31,29)),
'B':((14,10,26,23),(11,8,20,18),(20,15,35,32),(19,14,34,31)),
'C':((23,18,44,40),(18,14,34,31),(44,34,80,73),(40,31,73,67)),
'J0':((12,12,24,24),(12,12,21,21),(15,12,30,27),(15,12,30,27)),
'J1':((12,12,21,21),(12,12,18,18),(18,15,32,30),(18,15,31,29))}

def ycoords(s):
 D,E,U,d=s
 return D-E-d,U-D-d+E,E,d-E

def record(kind,v,w,bank):
 return dict(kind=kind,left=v,right=w,left_recipe=bank[v],right_recipe=bank[w])

def c3_support(bound):
 bank=source_bank(2*bound//15)
 sources=sorted((v for v in bank if v!=ONE),key=lambda s:(s[2],s))
 out={};tests={k:0 for k in LOWER}
 for kind,(num,den) in LOWER.items():
  # Handle the one admitted singleton side explicitly.
  if kind in ('C','J0') and ONE in bank:
   for s in sources:
    v,w=(ONE,s) if kind=='C' else (s,ONE)
    h=c3.value(kind,v,w);tests[kind]+=1
    if h<=bound:out.setdefault(h,record(kind,v,w,bank))
  for i,v in enumerate(sources):
   for j in range(i if kind in ('A','C') else 0,len(sources)):
    w=sources[j]
    if num*v[2]*w[2]>den*bound:break
    tests[kind]+=1;h=c3.value(kind,v,w)
    if h<=bound:out.setdefault(h,record(kind,v,w,bank))
 return out,bank,tests

def c1_support(bound):
 bank=source_bank(bound);out={}
 for v,recipe in bank.items():
  if v==ONE:continue
  D,E,U,d=v;h=D-E+U
  if h<=bound:out.setdefault(h,dict(kind='C1',source=v,recipe=recipe))
 return out,bank
