import json,itertools,sys
from pathlib import Path
def step(s,c):
 D,E,U,d=s
 if c=='O':return D,0,2*U,d
 if c=='Z':return D-E+U,d,2*U-d,d
 if c=='J':return D-E+U,D-E,D+U,D
 raise ValueError(c)
def F(s):return s[0]-s[1]+s[2]

def need(v,msg):
 if not v:raise RuntimeError(msg)

def graph(w):
 n=len(w)+2
 # Begin with all chronological forward edges, then reverse only the recorded root edge.
 rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
 root=0
 for newborn,c in enumerate(w,2):
  if c in 'ZJ':
   rows[root]^=1<<newborn
   rows[newborn]|=1<<root
  if c=='J':root=newborn
 return rows

def paths(rows):
 n=len(rows);d=[{} for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=1
 for mask in range(1,1<<n):
  for v,count in d[mask].items():
   left=rows[v]&~mask
   while left:
    b=left&-left;left-=b;w=b.bit_length()-1
    target=d[mask|b];target[w]=target.get(w,0)+count
 return sum(d[-1].values())

def state(w):
 s=(1,0,2,1)
 for c in w:s=step(s,c)
 return s

def strong(rows):
 n=len(rows);full=(1<<n)-1
 for i in range(n):
  reached=1<<i;front=reached
  while front:
   b=front&-front;front-=b
   new=rows[b.bit_length()-1]&~reached
   reached|=new;front|=new
  if reached!=full:return False
 return True
checks=0
for length in range(6):
 for tup in itertools.product('OZJ',repeat=length):
  w=''.join(tup);s=state(w)
  need(paths(graph(w))==s[0],f'D literal {w}')
  need(paths(graph(w+'J'))==F(s),f'F literal {w}')
  need(strong(graph(w+'J')),f'strong {w}')
  checks+=2
# A frozen producer adjacency artifact is input, not authority for reconstruction.
if len(sys.argv)!=2: raise SystemExit('Usage: python3 verify_graphs_independent.py PATH_TO_graph_witnesses.json')
src=Path(sys.argv[1])
for item in json.loads(src.read_text()):
 w=item['source_word']+'J';rows=graph(w)
 need(rows==list(map(int,item['outgoing_rows_decimal'])),'full witness edge mismatch')
 need(len(rows)==item['n'] and strong(rows),'full witness order/connectivity')
 need(F(state(item['source_word']))==int(item['H']),'full witness certified H')
out={'status':'PASS','independent_DP_checks':checks,'source_word_lengths':[0,5],'max_graph_order':8,'full_witnesses':'all edges reconstructed independently; strong connectivity checked; recurrence H matched, not full DP at36/72'}
Path(__file__).with_name('graph_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
