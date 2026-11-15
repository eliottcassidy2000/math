from score import read, mark_w
import json,sys,glob
directory=sys.argv[1] if len(sys.argv)>1 else "."
targets=[int(x) for x in sys.argv[2:]] or [32825,34111,35099,35395,35557,36833]
# For each fixed target use the exact necessary interval derived from 0<=W<=M.
W=[0]*59999
for m in range(1,len(W)):
    # Shift existing marks; adding the bit at zero contributes floor(m/2)
    # plus one extra if its adjacent bit is set.
    k=m>>1
    W[m]=W[k]+((k+(k&1)) if m&1 else 0)
answers={N:[] for N in targets}
for path in [f'{directory}/pairline10.tsv']+sorted(glob.glob(f'{directory}/pairline_1_*.tsv')):
 for n,x in read(path).items():
  for h,q,z,st in x:
   if q>=59999: continue
   for N in targets:
    lo=max(1,(N-h+q-1)//q)
    hi=(N-h)//(q-h)
    for m in range(lo,hi+1):
     if h+m*q-W[m]*h==N:
      answers[N].append([n,h,q,z,m,W[m]])
print(json.dumps(answers,indent=2))
