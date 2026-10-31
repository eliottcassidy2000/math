#!/usr/bin/env python3
"""Frozen strong representatives; invokes the installed nauty gentourng."""
import subprocess,sys
maximum=int(sys.argv[1]) if len(sys.argv)>1 else 9
print('1 -')
for n in range(3,maximum+1):
    p=subprocess.Popen(['gentourng','-cq',str(n)],stdout=subprocess.PIPE,text=True)
    for line in p.stdout:print(n,line.strip())
    if p.wait()!=0:raise RuntimeError('gentourng failed')
