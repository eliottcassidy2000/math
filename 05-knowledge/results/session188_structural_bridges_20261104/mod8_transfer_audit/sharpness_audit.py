"""Independent endpoint-DP check of the modulo-sixteen boundary."""
from pathlib import Path
import hashlib,json,subprocess,sys
from audit import prof,root,H,need
HERE=Path(__file__).resolve().parent;SOURCE=HERE.parent/'mod8_transfer'
A=(48,17,59,19,32,10);B=(50,16,59,19,0,26)
p,q=prof(A),prof(B);a,b=H(root((2,4,1),0,A)),H(root((2,4,1),0,B))
need(p==[0,11,178,828,1752,1800,720] and q==[0,3,90,612,1608,1800,720],'independent full profiles')
need((a,b)==(211,99) and (p[1]-q[1])%16==8 and (a-b)%16==0 and a!=b,'actual modular boundary')
normal=subprocess.check_output([sys.executable,str(SOURCE/'sharpness.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(SOURCE/'sharpness.py')]);need(normal==optimized,'producer boundary modes agree')
x={'status':'PASS','A':A,'B':B,'qA':p,'qB':q,'ambient_counts':[a,b],'core_gap':p[1]-q[1],'ambient_gap':a-b,'source_sha256':{str(v.relative_to(HERE.parent)):hashlib.sha256(v.read_bytes()).hexdigest() for v in [SOURCE/'SHARPNESS.md',SOURCE/'sharpness.py',HERE/'audit.py',Path(__file__)]},'producer_normal_optimized_equal':True,'producer_receipt_sha256':hashlib.sha256(normal).hexdigest()}
print(json.dumps(x,sort_keys=True,indent=2))
