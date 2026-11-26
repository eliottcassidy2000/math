"""Optional search driver; exact static replay uses replay.py and check.cpp."""
from pathlib import Path
import subprocess,sys
here=Path(__file__).parent
root=Path(sys.argv[1])
i=int(sys.argv[2]);cap=int(sys.argv[3]);n=int(sys.argv[4])
f=here/f'extended{i}.b64'
if not f.exists():f.write_text((here/f'after{i}.b64').read_text())
for r in range(n):
 p=subprocess.run([str(here/'check'),str(root/f'b{i}.txt'),str(f),str(cap)],
                  capture_output=True,text=True,check=True)
 print(p.stdout,end='',flush=True)
 gain=next((line for line in p.stdout.splitlines() if line.startswith('GAIN ')),None)
 if not gain:break
 move=gain.removeprefix('GAIN ')
 rm,add=map(str.split,move.split('>'))
 s=f.read_text().strip()
 rows=[s[j:j+6] for j in range(0,len(s),6)]
 for z in rm:rows.remove(z)
 rows+=add
 assert len(add)==len(rm)+1
 with (here/f'new{i}.txt').open('a') as out:
  out.write(str(i)+' '+' '.join(rm)+' > '+' '.join(add)+'\n')
 f.write_text(''.join(rows)+'\n')
