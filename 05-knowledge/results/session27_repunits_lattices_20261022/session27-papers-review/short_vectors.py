from pathlib import Path
import numpy as np,json
p=str(Path(__file__).with_name('downloaded_G96_reduced.txt'));G=np.array([[int(x) for x in row.split(',')] for row in open(p).read().strip().strip('[]').split(';')],dtype='int64')
v,U=np.linalg.eigh(G); best=10000;bx=None
# Exact-evaluated integer vectors from rounded multiples of short real eigenvectors.
for j in range(4):
 w=U[:,j]/max(abs(U[:,j]))
 for scale in np.arange(1,20000,.25):
  x=np.rint(scale*w).astype('int64');q=int(x@G@x)
  if q>0 and q<best:
   best=q;bx=x;print('best',best,'scale',scale,'eigen',j,flush=True)
   if q<=8:break
 if best<=8:break
print('RESULT',best,bx.tolist(),flush=True)
open(str(Path(__file__).with_name('short_vector_witness.json')),'w').write(json.dumps({'norm':best,'vector':bx.tolist()}))
