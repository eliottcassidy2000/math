"""Independent sparse terminal witness replay, no worker imports."""
import argparse,gzip,hashlib,json
from pathlib import Path
P=Path(__file__).resolve().parent.parent
ap=argparse.ArgumentParser()
ap.add_argument('--edges',type=Path,default=P/'session38-wiring/edges.tsv.gz')
ap.add_argument('--witness',type=Path,default=P/'session39-linking/n8_terminal_lift.json')
args=ap.parse_args()
edge=args.edges;witness=args.witness
z=json.loads(witness.read_text());v=z['lift'];current=[0]*len(v)
energy=0;edge_count=0
with gzip.open(edge,'rt') as f:
 if next(f).strip()!='source\ttarget\tblue\tblack':raise ValueError('schema')
 for line in f:
  a,b,blue,black=map(int,line.split());weight=blue+black
  current[a]+=weight*(v[a]-v[b]);current[b]+=weight*(v[b]-v[a])
  energy+=weight*(v[a]-v[b])**2;edge_count+=1
if any(a%256 for a in current):raise ValueError('nonintegral torsion representative')
if current!=[256*x for x in z['A_lift_divided_by_modulus']]:raise ValueError('current mismatch')
if energy!=sum(x*y for x,y in zip(v,current)):raise ValueError('energy mismatch')
if energy!=14960072960 or v[-1]!=0 or not any(x%2 for x in v):raise ValueError('terminal witness')
u=(energy//256)%256
if (169*u)%256!=1:raise ValueError('normalization')
print(json.dumps({'edge_count':edge_count,'vertices':len(v),'energy':energy,'unit':u,'normalized_unit':169*u%256,'scaled_normalized_unit':169*3*u%256,'edge_sha256':hashlib.sha256(edge.read_bytes()).hexdigest(),'witness_sha256':hashlib.sha256(witness.read_bytes()).hexdigest()},indent=2))
