"""Exact checks for the two diagrams; deliberately uses no assert statements."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,xml.etree.ElementTree as ET
P=Path(__file__).resolve().parent
def need(ok,msg):
 if not ok:raise ValueError(msg)
def rank(a):
 a=[[Fraction(x) for x in row] for row in a]; r=0
 for c in range(len(a[0])):
  k=next((k for k in range(r,len(a)) if a[k][c]),None)
  if k is None:continue
  a[r],a[k]=a[k],a[r];p=a[r][c];a[r]=[x/p for x in a[r]]
  for k in range(r+1,len(a)):
   if a[k][c]:
    q=a[k][c];a[k]=[x-q*y for x,y in zip(a[k],a[r])]
  r+=1
  if r==len(a):break
 return r
def det(a):
 a=[list(r) for r in a];n=len(a);last=1;sign=1
 for k in range(n-1):
  p=next((i for i in range(k,n) if a[i][k]),None)
  if p is None:return 0
  if p!=k:a[k],a[p]=a[p],a[k];sign=-sign
  pivot=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    num=a[i][j]*pivot-a[i][k]*a[k][j]
    need(num%last==0,'Bareiss nonexact division');a[i][j]=num//last
   a[i][k]=0
  last=pivot
 return sign*a[-1][-1]
def main():
 # Coefficients of E,x,z, in magic placement.
 values={'A':(1,1,0),'B':(1,-1,-1),'C':(1,0,1),'D':(1,-1,1),'E':(1,0,0),'F':(1,1,-1),'G':(1,0,-1),'H':(1,1,1),'I':(1,-1,0)}
 def add(names):return tuple(sum(values[n][k] for n in names) for k in range(3))
 h={'ABC','GHI','AEI','BEH','CEG'};v={'ADG','CFI','AEI','DEF','CEG'}
 need(len(h|v)==8 and len(h&v)==2,'line overlap')
 need(all(add(line)==(3,0,0) for line in h|v),'symbolic magic sums')
 hs=set('ABCEGHI');vs=set('ACDEFGI');need(hs&vs==set('ACEGI'),'cell seam')
 grid=['B','I','D','G','E','C','F','A','H'];tree=[]
 for i in [-1,0,1]:
  for j in [-1,0,1]:
   name=grid[(i+1)*3+j+1];need(values[name]==(1,i,j),'grid-to-magic label mapping')
   if i==0 or j==0:tree.append((('A',i),('B',j)));need(name in hs&vs,'gray tree')
   elif i==j:need(name in hs-vs,'blue chords')
   else:need(name in vs-hs,'orange chords')
 seen={('A',0)}
 while True:
  new=seen|{b for a,b in tree if a in seen}|{a for a,b in tree if b in seen}
  if new==seen:break
  seen=new
 need(len(tree)==5 and len(seen)==6,'spanning tree')
 need(7-6+1==2 and 9-6+1==4,'cycle ranks')
 need(rank([[0,1,0,0,-2,0,0,1,0],[0,0,0,1,-2,1,0,0,0]])==2,'AP independent constraints')
 texts=[t.text or '' for t in ET.parse(P/'hourglass-carrier.svg').iter('{http://www.w3.org/2000/svg}text')]
 labels={'A':'E + x','B':'E − x − z','C':'E + z','D':'E − x + z','E':'E','F':'E + x − z','G':'E − z','H':'E + x + z','I':'E − x'}
 # E appears both as cell name and as its displayed value; only expected sequence is tested.
 expected=[]
 for n in list('ABCDEFGHI')+grid:expected.extend([n,labels[n]])
 actual=[]
 for start in [texts.index('A'),texts.index('B',texts.index('B₁ = E + z')+1)]:actual.extend(texts[start:start+18])
 need(actual==expected,'SVG displayed cell names and values')
 fills=[r.attrib['fill'] for r in ET.parse(P/'hourglass-carrier.svg').iter('{http://www.w3.org/2000/svg}rect') if r.attrib.get('stroke')=='#ffffff']
 need(fills==[('#e2e8ee' if n in hs&vs else '#d7eafe' if n in hs-vs else '#ffe5c8') for n in list('ABCDEFGHI')+grid],'SVG category colors')
 d=json.loads((P/'observer-data.json').read_text());a=d['matrix'];cols=d['columns'];features=d['selected_features']
 need(len(a)==34 and len(cols)==78 and all(len(r)==78 for r in a),'profile matrix shape')
 need(len({tuple(r) for r in a})==34,'individual profiles')
 need(rank(a)==32,'exact linear rank')
 need(features[-2:]==[['square',3],['square',15]],'specified squared coordinates')
 need(len(features)==34 and all(t=='linear' for t,q in features[:32]),'feature types')
 minor=[[row[cols.index(q)]**(1 if t=='linear' else 2) for t,q in features] for row in a]
 need(minor==d['minor'],'minor reconstructed from original profile matrix')
 determinant=det(minor)
 need(determinant==211991592542375584355844096000000,'exact determinant')
 need(determinant%1000000007==967372386,'determinant residue')
 ot=' '.join(t.text or '' for t in ET.parse(P/'observer-profile.svg').iter('{http://www.w3.org/2000/svg}text'))
 need(str(determinant) in ot,'displayed determinant');need('No individual-profile collision' in ot and 'No claim' in ot,'scope guards')
 for name,digest in json.loads((P/'SHA256SUMS.json').read_text()).items():need(hashlib.sha256((P/name).read_bytes()).hexdigest()==digest,'artifact hash: '+name)
 print('PASS: exact magic sums; 5-edge tree; 2+2 chords; displayed label map; 34 distinct profiles; rank 32; reconstructed rank-34 minor; artifact hashes.')
if __name__=='__main__':main()
