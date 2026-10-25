import json
from pathlib import Path
from fractions import Fraction as F

def req(x,m):
 if not x:raise RuntimeError(m)
def det(a):
 a=[row[:] for row in a];sgn=1;prev=1
 for k in range(len(a)-1):
  q=next((i for i in range(k,len(a)) if a[i][k]),None)
  if q is None:return 0
  if q!=k:a[k],a[q]=a[q],a[k];sgn=-sgn
  pivot=a[k][k]
  for i in range(k+1,len(a)):
   for j in range(k+1,len(a)):
    z=a[i][j]*pivot-a[i][k]*a[k][j];req(z%prev==0,'Bareiss exact division');a[i][j]=z//prev
   a[i][k]=0
  prev=pivot
 return sgn*a[-1][-1]
def main():
 data=json.loads(Path(__file__).with_name('direct_certificate.json').read_text());A=data['matrix'];cols=data['columns'];basis={};chosen=[];vectors=[]
 def add(v,label):
  w=[F(x) for x in v]
  for p,b in sorted(basis.items()):
   if w[p]:
    c=w[p];w=[x-c*y for x,y in zip(w,b)]
  p=next((i for i,x in enumerate(w) if x),None)
  if p is None:return False
  c=w[p];basis[p]=[x/c for x in w];chosen.append(label);vectors.append(v);return True
 for j,c in enumerate(cols):add([row[j] for row in A],['linear',c])
 req(len(basis)==32,'linear rank')
 for j,c in enumerate(cols):
  add([row[j]**2 for row in A],['square',c])
  if len(basis)==34:break
 if len(basis)<34:
  for j,c in enumerate(cols):
   for k,d in enumerate(cols[j+1:],j+1):
    add([row[j]*row[k] for row in A],['product',c,d])
    if len(basis)==34:break
   if len(basis)==34:break
 req(len(basis)==34,'quadratic full rank')
 minor=[[v[i] for v in vectors] for i in range(34)];D=det(minor);req(D!=0,'nonzero exact minor')
 w={int(k):v for k,v in data['signed_kernel'].items()};contrast=[]
 for label,v in zip(chosen,vectors):
  if label[0]!='linear':contrast.append({'feature':label,'signed_mixture_numerator':sum(w[r]*z for r,z in zip(data['representatives'],v))})
 out={'linear_rank':32,'quadratic_rank':34,'selected_features':chosen,'selected_integer_minor':minor,'determinant':str(D),'prime':1000000007,'determinant_mod_prime':D%1000000007,'mixture_contrasts':contrast}
 Path(__file__).with_name('quadratic_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='selected_integer_minor'},indent=2))
if __name__=='__main__':main()
