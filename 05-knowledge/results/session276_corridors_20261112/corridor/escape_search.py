from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'session275_exact_selection_20261111'/'switches'))
from carrier import profile,matrix,mul,from_matrix,fuse,hp
I=[[int(i==j) for j in range(4)] for i in range(4)]
T43=(52,13,24,49,34,6);T45=(52,13,24,33,42,6)
A=matrix(profile(T43));B=matrix(profile(T45))
def context(n,a):
 return (a<<1,)+tuple((((1<<n)-1)^((1<<(i+1))-1))<<1|(1-(a>>i&1)) for i in range(n))

def run():
 out=[]
 for n in range(1,9):
  for cut in range(1<<n):
   C=context(n,cut);M=matrix(profile(C));P=I
   seq=[[],[]];base=[[],[]]
   for k in range(6):
    for side in range(2):
     pa=mul(P,A) if side==0 else mul(A,P);pb=mul(P,B) if side==0 else mul(B,P)
     seq[side].append(from_matrix(pb)['H']-from_matrix(pa)['H']);base[side].append(from_matrix(pa)['H'])
    P=mul(P,M)
   for side in range(2):
    if seq[side]==[2]*6 and len(set(base[side]))>1:out.append(dict(n=n,cut=cut,rows=C,side=side,H=base[side],difference=seq[side],profile=profile(C)))
 print(json.dumps({'T43':profile(T43),'T45':profile(T45),'renewable':out},indent=2))
if __name__=='__main__':run()
