#!/usr/bin/env python3
"""Mutation checks are run as subprocesses in normal and optimized modes."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parent
CHECK=ROOT/'check_independent.py'


def need(c,msg):
    if not c: raise ValueError(msg)


def main():
    baseline=json.loads(subprocess.check_output([sys.executable,str(CHECK)],text=True))
    outcomes=[]
    with tempfile.TemporaryDirectory() as d:
        tmp=Path(d)
        (tmp/'sources').mkdir()
        (tmp/'sources/session143_tournaments.json').write_bytes((ROOT/'sources/session143_tournaments.json').read_bytes())
        for mode in ([],['-O']):
            clean=tmp/'clean.json'; clean.write_text(json.dumps(baseline))
            proc=subprocess.run([sys.executable,*mode,str(CHECK),'--certificate',str(clean)],capture_output=True,text=True)
            need(proc.returncode==0,'positive certificate control')
            outcomes.append({'mode':'optimized' if mode else 'normal','control':'unmodified certificate','exit':proc.returncode})
            for label,path,value in [
                ('wrong N3 count',['models',0,'N3'],999),
                ('wrong magic line',['magic_grid',0,0],109),
                ('wrong elliptic fiber partition',['trace_fiber_partition'],[3,2]+[1]*7),
                ('false twist fixed point',['twist_fixed_points'],1),
                ('false symmetric cube trace',['operator_traces',0,'symmetric_cube_trace'],999),
                ('missing tournament self-converse',['tournament_control','self_converse'],0),
            ]:
                j=copy.deepcopy(baseline); p=j
                for step in path[:-1]:p=p[step]
                p[path[-1]]=value
                bad=tmp/'bad.json';bad.write_text(json.dumps(j))
                proc=subprocess.run([sys.executable,*mode,str(CHECK),'--certificate',str(bad)],capture_output=True,text=True)
                need(proc.returncode!=0 and 'certificate mismatch' in proc.stderr,'hostile certificate passed: '+label)
                outcomes.append({'mode':'optimized' if mode else 'normal','control':label,'exit':proc.returncode,'reason':'certificate mismatch'})
            source=CHECK.read_text()
            for label,old,new,reason in [
                ('reducible cubic modulus','Field([1,2,0,1])','Field([0,2,0,1])','reducible extension modulus'),
                ('omit infinity point','n = 1+sum(sq[r] for r in rhs)','n = sum(sq[r] for r in rhs)','point count histogram'),
                ('replace power trace by symmetric trace','ns[2]==126-t**3+15*t','ns[2]==126-t**3+10*t','extension recurrence comparison'),
                ('replace nonsquare twist by square twist','(4*r[0]%5,3*r[1]%5)','(r[0]%5,4*r[1]%5)','free quadratic twist involution'),
                ('replace symmetric trace by power trace','sym=tr(symcube(m))','sym=tr(mm(mm(m,m),m))','operator traces'),
            ]:
                need(source.count(old)==1,'mutation target uniqueness: '+label)
                mutant=tmp/'mutant.py';mutant.write_text(source.replace(old,new))
                proc=subprocess.run([sys.executable,*mode,str(mutant)],capture_output=True,text=True)
                need(proc.returncode!=0 and reason in proc.stderr,'hostile code passed: '+label)
                outcomes.append({'mode':'optimized' if mode else 'normal','control':label,'exit':proc.returncode,'reason':reason})
    print(json.dumps({'status':'PASS','positive_controls':2,'negative_controls':22,'outcomes':outcomes},indent=2,sort_keys=True))

if __name__=='__main__':main()
