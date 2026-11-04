"""Exact local menu; --all-prefixes terminates by positive endpoint count."""
from pathlib import Path
import argparse,json,subprocess,tempfile

HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('depth',type=int,nargs='?',default=1)
parser.add_argument('moving',type=int,nargs='?',default=0)
parser.add_argument('--seed',type=int,help='Index in realized_depth3.json')
parser.add_argument('--all-prefixes',action='store_true')
args=parser.parse_args()
if args.seed is None:
    old=HERE.parents[1]/'triangle_return_170807_20261103'
    record=json.loads((old/'zero_ready.json').read_text())[1]
    source_H=170807;inport=record['inport'];minimum=2
    source='triangle_return_170807_20261103/zero_ready.json record1'
    name=f'suffix_depth{args.depth}_moving{args.moving}'
else:
    record=json.loads((HERE/'realized_depth3.json').read_text())[args.seed]
    source_H=record['target_H'];inport=record['target_inport']
    minimum=len(record['prefix_names'])
    source=f'realized_depth3.json index{args.seed}'
    name=f'generation2_seed{args.seed}_depth{args.depth}'
maximum=60 if args.all_prefixes else 5
if args.all_prefixes:name+='_all_prefixes'
data=' '.join(map(str,[len(record['suffix_rows']),record['suffix_names'].index(inport),
    args.depth,args.moving,source_H,minimum,maximum,*record['suffix_rows']]))+'\n'
with tempfile.TemporaryDirectory(prefix='strong-suffix-') as temp:
    exe=Path(temp)/'search.exe'
    subprocess.run(['g++','-O3','-std=c++20',str(HERE/'search.cpp'),'-o',str(exe)],check=True)
    value=json.loads(subprocess.check_output([str(exe)],input=data,text=True))
value.update(suffix_names=record['suffix_names'],source_rows=record['suffix_rows'],
             source_inport=inport,input_source=source,moving_inport=bool(args.moving),
             all_prefixes=args.all_prefixes,
             prefix_termination='h + 2^s F is strictly increasing since F > 0; stop after first overshoot')
(HERE/(name+'.json')).write_text(json.dumps(value,indent=2)+'\n')
print(json.dumps(value))
