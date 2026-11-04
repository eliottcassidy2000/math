"""Independent subset-DP review of the modulo-sixteen boundary."""
from verify import q_partitions,hp,blow,need
from pathlib import Path
import hashlib,json
A=(48,17,59,19,32,10);B=(50,16,59,19,0,26)
for T in(A,B):
 need(all(not(T[i]>>i&1) for i in range(6)),'loop')
 need(all(((T[i]>>j)&1)+((T[j]>>i)&1)==1 for i in range(6) for j in range(i+1,6)),'not tournament')
qa,qb=q_partitions(A)[1:],q_partitions(B)[1:]
need(qa==(11,178,828,1752,1800,720) and qb==(3,90,612,1608,1800,720),'complete core profiles')
x,y=hp(blow((2,4,1),(A,(0,),(0,)))),hp(blow((2,4,1),(B,(0,),(0,))))
need((x,y)==(211,99) and x-y==3*(qa[0]-qb[0])+qa[1]-qb[1],'literal ambient counts')
need((qa[0]-qb[0])%16==8 and (x-y)%16==0,'mod16 sharpness')
r=Path(__file__).resolve().parent.parent/'mod8_transfer'
print(json.dumps({'status':'PASS','method':'Independent recursive endpoint DP and canonical unordered-cover partitions, not producer literal-word enumeration.',
                  'q_A':qa,'q_B':qb,'ambient_counts':[x,y],'core_gap':qa[0]-qb[0],'ambient_gap':x-y,
                  'source_pins':{f:hashlib.sha256((r/f).read_bytes()).hexdigest() for f in ['SHARPNESS.md','sharpness.py','sharpness.json','sharpness.optimized.json']}},indent=2,sort_keys=True))
