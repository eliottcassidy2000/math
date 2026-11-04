"""Compile all one-start classes through a certified parent-order bound."""
from pathlib import Path
import json,sys,subprocess,tempfile,time
HERE=Path(__file__).resolve().parent
def need(b,s):
    if not b:raise ValueError(s)
def main():
    exe=Path(tempfile.gettempdir())/'funnel-bank18979.exe'
    subprocess.run(['g++','-std=c++20','-O3',str(HERE/'funnel_bank.cpp'),'-o',str(exe)],check=True)
    start=time.monotonic();p=subprocess.Popen([str(exe),'18'],stdout=subprocess.PIPE,text=True)
    rows=[]
    for line in p.stdout:
        row=json.loads(line);rows.append(row);print(row,flush=True)
    need(p.wait()==0,'complete graded grammar run')
    report=dict(status='FINITE_EXACT_FULL_ONE_START_GRAMMAR_TO_PARENT18',rows=rows,seconds=time.monotonic()-start)
    (HERE/'one_start_bank.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
