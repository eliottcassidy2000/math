"""Read source set from pinned manuscript; independently verify recovered witness."""
import argparse, ast, hashlib, json, re
from pathlib import Path
from audit import adjacency, shifted_from_formula, require

parser=argparse.ArgumentParser()
parser.add_argument('--repo',type=Path,required=True)
args=parser.parse_args()
pins=json.loads(Path(__file__).with_name('SOURCE_PINS_AUDITED.json').read_text())
verified=[]
for bundle in [pins,pins['recovered_concurrent_source']]:
    for pin in bundle['files']:
        source=args.repo/pin['path']
        digest=hashlib.sha256(source.read_bytes()).hexdigest()
        require(digest==pin['sha256'],'source pin '+pin['path'])
        verified.append({'path':pin['path'],'sha256':digest})
path='05-knowledge/results/unprotected_ports_cancellation_20261025/counting/POSITIVE_JET.md'
text=(args.repo/path).read_text()
block=re.search(r'S=(\{.*?\})\.\s*\(J7\)',text,re.S)
require(block is not None,'source J7 connection set recovered literally')
S=ast.literal_eval(block.group(1))
M=adjacency()
xy=[divmod(i,5) for i in range(25)]
source=[[int(((b[0]-a[0])%5,(b[1]-a[1])%5) in S) for b in xy] for a in xy]
phi=[5*((-x)%5)+(-y)%5 for x,y in xy]
require(all(M[i][j]==source[phi[i]][phi[j]] for i in range(25) for j in range(25)),'core isomorphism all625')
def lift(A,i,j):
    return int(i%2<j%2) if i//2==j//2 else A[i//2][j//2]
phi50=[2*phi[i//2]+i%2 for i in range(50)]
require(all(lift(M,i,j)==lift(source,phi50[i],phi50[j]) for i in range(50) for j in range(50)),'inflation isomorphism all2500')
a=shifted_from_formula(5,2,1)
Q=[(-1)**(j+1)*v for j,v in enumerate(a)]
require(all(v>=0 for v in Q),'positive-source-jet convention')
require([v%5 for v in Q]==[0]*21+[4,3,3,4,1],'source mod5 polynomial')
print(json.dumps({'verdict':'PASS','source_commit':pins['recovered_concurrent_source']['commit'],
 'classification':'Recovered concurrent theorem and exact witness; not a new witness or new regular-spectrum theorem',
 'source_parsed_from':path,'source_connection_set':sorted(S),'source_hashes_verified':verified,
 'isomorphism':'(u,b) -> (-u,b)','core_entries':625,'inflated_entries':2500,
 'positive_source_shift':'Q(u) = -P(-1-u)','source_mod5_nonzero_coefficients':[4,3,3,4,1]},indent=2,sort_keys=True))
