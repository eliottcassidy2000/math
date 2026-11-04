from pathlib import Path
import json,tempfile,subprocess
HERE=Path(__file__).resolve().parent
def main():
    exe=Path(tempfile.gettempdir())/'feedback56941.exe'
    subprocess.run(['g++','-std=c++20','-O3',str(HERE/'feedback.cpp'),'-o',str(exe)],check=True)
    p=subprocess.Popen([str(exe),'3'],stdout=subprocess.PIPE,text=True);rows=[]
    for line in p.stdout:rows.append(json.loads(line));print(line.strip(),flush=True)
    if p.wait():raise ValueError('feedback search failed')
    (HERE/'feedback.json').write_text(json.dumps(dict(status='FINITE_EXACT_ALL_CUTS_WITH_UP_TO_THREE_SPLIT_PREFIX_VERTICES',rows=rows),indent=2)+'\n')
if __name__=='__main__':main()
