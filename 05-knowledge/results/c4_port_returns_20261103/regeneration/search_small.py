from pathlib import Path
import json,subprocess,tempfile
HERE=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='regen-small-') as temp:
    exe=Path(temp)/'search.exe'
    subprocess.run(['g++','-O3','-std=c++20',str(HERE/'small_bank.cpp'),'-o',str(exe)],check=True)
    for mode in ('local','full','two','complete7'):
        value=json.loads(subprocess.check_output([str(exe),mode],text=True))
        for hit in value['hits']:
            if mode=='two':hit['birth_cuts']=[hit['cut']%32,hit['cut']//32]
            elif mode=='local':hit['birth_cuts']=[hit['cut']]
        (HERE/(mode+'_small.json')).write_text(json.dumps(value,indent=2)+'\n')
        print(json.dumps({'mode':mode,'counts':value['counts'],'witnesses':len(value['hits'])}))
