from pathlib import Path
import json,subprocess,tempfile
HERE=Path(__file__).resolve().parent
def main():
    rows=json.loads((HERE.parent/'return18979_palindrome_20261103/returns.json').read_text())
    exe=Path(tempfile.gettempdir())/'births56941.exe';subprocess.run(['g++','-std=c++20','-O3',str(HERE/'births.cpp'),'-o',str(exe)],check=True)
    text=''.join(' '.join(map(str,[len(x['next_rows']),x['next_H'],*x['next_rows']]))+'\n' for x in rows)
    r=subprocess.run([str(exe)],input=text,text=True,capture_output=True,check=True);out=[]
    for line in r.stdout.splitlines():
        H,b,*cuts=map(int,line.split());out.append(dict(H=H,offset=b,incoming_cuts=cuts));print(H,b,cuts,flush=True)
    (HERE/'births.json').write_text(json.dumps(dict(status='FINITE_EXACT_ALL_ONE_BIRTH_CUTS_OF_FOUR_ACTUAL_SOURCES',branches=out),indent=2)+'\n')
if __name__=='__main__':main()
