"""All-order ordinal-triangle zero classification: bounded exact controls."""
import itertools,json,hashlib
from quotient_kernel import kernel

def coefficients(k,r):
 if k==1:return 1,2,{}, {}, {}
 h=3**(k-1);L=2*3**(k-2); others=[i for i in range(k) if i!=r]
 root={};cross={}
 for i in others:
  d=abs(i-r);big=h*4**(d-1);small=big+(3**(k-2) if d==1 else 0)
  root[i]=(2*big,small) if r<i else (small,2*big)
 for i,j in itertools.combinations(others,2):
  cross[i,j]=h*4**(j-i-1)+(3**(k-3) if j==i+1 else 0)
 return h,2*h,{i:L for i in others},root,cross

def predicted(k,r):
 t=[0 if i<r else 3 for i in range(k)];t[r]=1;z={tuple(t)}
 if r:
  t=[0 if i<r else 3 for i in range(k)];t[r]=2;t[r-1]=1;z.add(tuple(t))
 if r+1<k:
  t=[0 if i<r else 3 for i in range(k)];t[r]=2;t[r+1]=2;z.add(tuple(t))
 return z

def value(k,r,t,coeff):
 h,pair,selfc,root,cross=coeff;a=t[r]&1;b=(t[r]>>1)&1;v=h+pair*a*(1-b)
 v+=sum(c for i,c in selfc.items() if t[i] in (1,2))
 for i,(u,w) in root.items():
  v+=(u*a+w*b)*(3-t[i]) if r<i else (u*(1-a)+w*(1-b))*t[i]
 for (i,j),c in cross.items():v+=c*t[i]*(3-t[j])
 return v

def literal_matrix(k,r,coeff):
 h,pair,selfc,root,cross=coeff;sizes=[3]*k;sizes[r]=2;offset=[];s=0
 for n in sizes:offset.append(s);s+=n
 C=[[0]*s for _ in range(s)];C[offset[r]][offset[r]+1]=C[offset[r]+1][offset[r]]=pair
 for i,c in selfc.items():
  for a,b in itertools.combinations(range(3),2):C[offset[i]+a][offset[i]+b]=C[offset[i]+b][offset[i]+a]=c
 for i,(u,w) in root.items():
  for a in range(3):
   for j,c in enumerate((u,w)):C[offset[i]+a][offset[r]+j]=C[offset[r]+j][offset[i]+a]=c
 for (i,j),c in cross.items():
  for a in range(3):
   for b in range(3):C[offset[i]+a][offset[j]+b]=C[offset[j]+b][offset[i]+a]=c
 return C

if __name__=='__main__':
 rows=[];checks=0
 for k in range(1,9):
  q=tuple(sum(1<<j for j in range(i+1,k)) for i in range(k))
  for r in range(k):
   cc=coefficients(k,r)
   if k<=6:
    h,C=kernel(q,r)
    if h!=cc[0] or C!=literal_matrix(k,r,cc):raise RuntimeError(('formula/kernel mismatch',k,r))
   zeros={t for t in itertools.product(range(4),repeat=k) if value(k,r,t,cc)==3**k};checks+=4**k
   if zeros!=predicted(k,r):raise RuntimeError(('classification',k,r,zeros,predicted(k,r)))
   weight=sum(3**sum(t[i] in (1,2) for i in range(k) if i!=r) for t in zeros)
   rows.append({'k':k,'root':r,'zero_profiles':sorted(zeros),'labelled_zero_count':weight})
  print(json.dumps({'k':k,'roots':k,'profile_checks':k*4**k,'labelled_zero_counts':[row['labelled_zero_count'] for row in rows if row['k']==k]},sort_keys=True),flush=True)
 print('PASS '+json.dumps({'rows':len(rows),'profile_checks':checks,'coefficient_audit_through_k':6,'zero_classification_through_k':8,'sha256':hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest()},sort_keys=True))
