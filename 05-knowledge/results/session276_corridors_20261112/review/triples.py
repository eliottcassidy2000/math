import contextlib,io,json
from pathlib import Path
with contextlib.redirect_stdout(io.StringIO()):
 from independent import graph,subset_count
 from sandwich import one,fuse,converse
 from carries import g

def F(A,B,c=641):return c+175*(g(A)+g(B))+1988*A*B+1036*A+777*B+100*(A&1)+75*(B&1)
Gs=[graph(7,set(range(7))-A) for A in ({0,2,4},{0,1,3})]
# Literal receipt uses the converse cyclic-core orientation (code0), same root cut.
for G in Gs:
 for i in range(1,8):
  for j in range(i+1,8):G[i][j],G[j][i]=G[j][i],G[i][j]
assert [subset_count(G) for G in Gs]==[641,643]
receipt=json.loads((Path(__file__).resolve().parents[1]/'digits/OVERLAP.json').read_text())
for row in receipt['literal']:
 G=Gs[(row['coreH']-641)//2];T=fuse(fuse(one(row['A']),G),converse(one(7)))
 rows=[sum(1<<j for j,v in enumerate(r) if v) for r in T]
 assert rows==row['T'];assert F(row['A'],7,row['coreH'])==row['literalH']
examples=[]
for t in [0,1,2,3,99,10**30,10**100]:
 A=170+568*t;B=7+25*t;N=F(A,B)
 assert F(A+1,B)==N+2
 examples.append({'t':t,'N':N})
print(json.dumps({'status':'PASS','independent_source_counts':[641,643],'independent_literal_graph_constructions_matched':len(receipt['literal']),'parameter_examples':examples,'scope':'independent adjacency and exact integer checks; nineteen-vertex DP result independently compared with scalar, not rerun by this script'},indent=2))
