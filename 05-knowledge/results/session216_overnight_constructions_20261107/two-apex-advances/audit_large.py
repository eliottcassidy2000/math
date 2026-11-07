import struct,random,subprocess,json,sys
from pathlib import Path
p=Path(__file__).parent
file=p/(sys.argv[1] if len(sys.argv)>1 else 'census14.bin')
b=file.read_bytes();a=struct.unpack('<'+'I'*(len(b)//4),b)
found=[h for h in range(1,len(a),2) if a[h]!=0xffffffff]
r=random.Random(216);sample=sorted(set(found[:100]+found[-100:]+r.sample(found,min(1000,len(found)))))
lines=[]
for h in sample:
 x=a[h];n=x>>28;w=[x>>(2*i)&3 for i in range(n)]
 while w and w[0]==0:w.pop(0)
 while w and w[-1]==3:w.pop()
 lines.append(str(h)+' '+(''.join(map(str,w)) or '-'))
r=subprocess.run(['/tmp/session216-inverse-audit'],input='\n'.join(lines)+'\n',text=True,capture_output=True)
if r.returncode:raise RuntimeError(r.stderr)
print(json.dumps({'file':file.name,'found':len(found),'checked':len(sample),'seed':216,'audit_stdout':r.stdout.strip(),'largest_checked':max(sample)},indent=2))
