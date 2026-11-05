import collections,hashlib,itertools,json,os,pathlib,subprocess,tempfile
from graphs import require,run_dp,strong
BASE=pathlib.Path(__file__).resolve().parent

def rank(rows,c):
 used={0}
 for row in rows:
  used|={mask|1<<j for mask in list(used)for j in range(c)if row>>j&1 and not mask>>j&1}
 return max(mask.bit_count()for mask in used)

def patterns(r,c):
 maps=[[sum(1<<p[j]for j in range(c)if x>>j&1)for x in range(1<<c)]for p in itertools.permutations(range(c))]
 representatives=set()
 for rows in itertools.combinations_with_replacement(range(1<<c),r):
  if rank(rows,c)==2:representatives.add(min(tuple(sorted(t[x]for x in rows))for t in maps))
 return sorted(representatives)

def canonical(rows):
 groups=collections.defaultdict(list)
 for v in range(8):groups[(bool(rows[8]>>v&1),rows[v].bit_count())].append(v)
 best=None
 for orders in itertools.product(*(itertools.permutations(group)for _,group in sorted(groups.items()))):
  order=sum((list(x)for x in orders),[])+[8]
  code=tuple(sum(1<<j for j,w in enumerate(order)if rows[v]>>w&1)for v in order)
  if best is None or code<best:best=code
 return best

def literal_runs(rows):
 result=[0]*5;full=255
 def walk(last,used,k):
  if used==full:result[k]+=1
  for v in range(8):
   if not used>>v&1 and rows[last]>>v&1:walk(v,used|1<<v,k)
  if last!=8 and k<4 and rows[last]>>8&1:walk(8,used,k+1)
 for v in range(8):walk(v,1<<v,0)
 walk(8,0,1)
 return result

cases=[];lines=['3']
for r in [2,3,4]:
 c=8-r;ps=patterns(r,c);require(len(ps)=={2:42,3:38,4:37}[r],'cross-pattern inventory')
 require(all(rank(p,c)==2 for p in ps),'pattern rank')
 lines.append(f'{r} {c} {len(ps)}');lines.extend(' '.join(map(str,p))for p in ps);cases.append((r,c,ps))
text='\n'.join(lines)+'\n'
with tempfile.TemporaryDirectory(prefix='clock-check-',dir=BASE)as tmp:
 exe=pathlib.Path(tmp)/'census'
 flags=['-O3'] if os.environ.get('CLOCK_CPP_MODE')!='debug'else['-O0','-D_GLIBCXX_ASSERTIONS']
 subprocess.run(['g++','-std=c++17',*flags,str(BASE/'census.cpp'),'-o',str(exe)],check=True)
 output=subprocess.run([str(exe)],input=text,text=True,capture_output=True,check=True).stdout
catalog={};receipts=[]
for line in output.splitlines():
 words=line.split()
 if words[0]=='CASE':receipts.append(list(map(int,words[1:])));continue
 require(words[0]=='UNIT'and len(words)==15,'census output shape')
 r,pattern,mask,A,beta=map(int,words[1:6]);rows=list(map(int,words[6:]));versions=[rows]
 if r<4:versions.append([511^row^(1<<i)for i,row in enumerate(rows)])
 for T in versions:
  code=canonical(T);a=run_dp(list(code),8)
  require(a[1:4]==[3*A+beta,3*A,A]and not any(a[4:]),'independent positive run check')
  catalog[code]={'rows':list(code),'root':8,'A':A,'beta':beta,'run':a[:4],'strong':strong(list(code))}
require(receipts==[[2,6,42,2752512,0],[3,5,38,311296,40],[4,4,37,151552,72]],'complete census inventories')
by_profile=collections.Counter((v['A'],v['beta'])for v in catalog.values())
actual={'rooted_classes':len(catalog),'by_profile':[{'A':a,'beta':b,'classes':n}for(a,b),n in sorted(by_profile.items())],'catalog':[catalog[k]for k in sorted(catalog)]}
require(actual==json.loads((BASE/'catalog.json').read_text()),'complete rooted catalog')
for T in actual['catalog']:
 a=literal_runs(T['rows']);require(a[:4]==T['run']and a[4]==0,'literal root-word enumeration')
 O=[v for v in range(8)if T['rows'][8]>>v&1];I=[v for v in range(8)if T['rows'][v]>>8&1]
 require(rank([sum(1<<j for j,w in enumerate(I)if T['rows'][v]>>w&1)for v in O],len(I))==2,'catalog matching')
print(json.dumps({'status':'PASS','cases':receipts,'candidate_graphs':sum(x[3]for x in receipts),'rooted_classes':len(catalog),'strong_rooted_classes':sum(v['strong']for v in catalog.values()),'profiles':actual['by_profile'],'pattern_input_sha256':hashlib.sha256(text.encode()).hexdigest(),'census_output_sha256':hashlib.sha256(output.encode()).hexdigest(),'catalog_sha256':hashlib.sha256((BASE/'catalog.json').read_bytes()).hexdigest(),'literal_run_controls':len(catalog)},sort_keys=True,indent=2))
