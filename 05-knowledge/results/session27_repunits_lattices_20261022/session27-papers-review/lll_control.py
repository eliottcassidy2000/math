from pathlib import Path
import numpy as np,json
G=np.array([[int(x) for x in row.split(',')] for row in open(str(Path(__file__).with_name('downloaded_G96_reduced.txt'))).read().strip().strip('[]').split(';')],dtype='int64');n=len(G);V=np.eye(n,dtype='int64');L=np.linalg.cholesky(G.astype(float));B=np.diag(L)**2;mu=L/np.diag(L)[None,:];k=1;steps=0
while k<n and steps<200000:
 steps+=1
 for j in range(k-1,-1,-1):
  q=round(mu[k,j])
  if q:V[k]-=q*V[j];mu[k,:j+1]-=q*mu[j,:j+1]
 if B[k]>=(.75-mu[k,k-1]**2)*B[k-1]:k+=1
 else:
  nu=mu[k,k-1];bn=B[k]+nu*nu*B[k-1];mu[k,k-1]=nu*B[k-1]/bn;B[k]=B[k]*B[k-1]/bn;B[k-1]=bn
  mu[[k,k-1],:k-1]=mu[[k-1,k],:k-1]
  for i in range(k+1,n):
   t=mu[i,k];mu[i,k]=mu[i,k-1]-nu*t;mu[i,k-1]=t+mu[k,k-1]*mu[i,k]
  V[[k,k-1]]=V[[k-1,k]];k=max(1,k-1)
H=V@G@V.T; norms=H.diagonal();i=int(np.argmin(norms));x=V[i].tolist();exact=sum(int(x[a])*int(G[a,b])*int(x[b]) for a in range(n) for b in range(n));print('steps',steps,'min basis norm',exact,'max coordinate',max(map(abs,x)));print('vector',x)
json.dump({'norm':exact,'vector':x,'source_sha256':'7a6a3b8f761a9f02e1a687be8040aaf403c7ac75c33d92f010f2dd131f17bc24'},open(str(Path(__file__).with_name('lll_witness.json')),'w'))
