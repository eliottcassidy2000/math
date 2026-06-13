#!/usr/bin/env python3
"""
tt_suite.py -- Tournament Theory Suite: Unified Tool Launcher
kind-pasteur-2026-03-24-S20cq

One command to access ALL tournament-theory-powered tools.

USAGE:
  python tt_suite.py status            # show all tools and versions
  python tt_suite.py compress FILE     # compress a file
  python tt_suite.py fingerprint FILE  # data fingerprint
  python tt_suite.py hash --text "..." # tournament hash
  python tt_suite.py rank --demo       # pairwise ranking demo
  python tt_suite.py predict --demo    # time series prediction
  python tt_suite.py sort --bench      # adaptive sort benchmark
  python tt_suite.py codec --demo      # staircase codec
  python tt_suite.py all               # run ALL demos

LICENSE: MIT
"""

import sys
import os
import argparse

__version__ = "1.0.0"
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

TOOLS = [
    ('tpress_v5',         'Universal Compressor v5',  '5.0.0'),
    ('tpress_cli',        'Compression CLI (.tp)',     '5.0.0'),
    ('dataprint',         'Data Fingerprinter',        '1.0.0'),
    ('tournament_hash',   'LSH via Tournaments',       '1.0.0'),
    ('pairwise_db',       'Pairwise Comparison DB',    '1.0.0'),
    ('tournament_predict','Time Series Predictor',     '1.0.0'),
    ('tournament_sort',   'Adaptive Sort',             '1.0.0'),
    ('staircase_codec',   'Tournament Theory Codec',   '1.0.0'),
    ('formalrank_production', 'FormalRank Engine',      '2.0.0'),
]


def cmd_status(args):
    print(f"Tournament Theory Suite v{__version__}")
    print("=" * 70)
    print(f"\n  {'Module':>25} {'Description':>28} {'Version':>8} {'Status':>8}")
    print("  " + "-" * 68)
    ok = 0
    for name, desc, _ in TOOLS:
        try:
            mod = __import__(name)
            ver = getattr(mod, '__version__', '?')
            ok += 1
            print(f"  {name:>25} {desc:>28} {ver:>8} {'OK':>8}")
        except Exception as e:
            print(f"  {name:>25} {desc:>28} {'---':>8} {'MISSING':>8}")
    print(f"\n  {ok}/{len(TOOLS)} tools available")
    print(f"  Path: {SCRIPT_DIR}")


def cmd_all(args):
    print("=" * 80)
    print("  TOURNAMENT THEORY SUITE -- Running All Demos")
    print("=" * 80)
    for name, desc, _ in [
        ('tpress_v5', 'COMPRESSION', ''),
        ('dataprint', 'DATA FINGERPRINT', ''),
        ('tournament_hash', 'TOURNAMENT HASHING', ''),
        ('pairwise_db', 'PAIRWISE RANKING', ''),
        ('tournament_predict', 'SEQUENCE PREDICTION', ''),
        ('staircase_codec', 'STAIRCASE CODEC', ''),
    ]:
        print(f"\n{'='*80}")
        print(f"  {desc}")
        print(f"{'='*80}")
        try:
            mod = __import__(name)
            # Find demo function
            for fn_name in ['demo', 'benchmark', '_run_demo', 'main']:
                fn = getattr(mod, fn_name, None)
                if fn and callable(fn):
                    fn()
                    break
        except Exception as e:
            print(f"  Error: {e}")


def main():
    parser = argparse.ArgumentParser(
        description=f'Tournament Theory Suite v{__version__}',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Commands:
  status     Show all tools and versions
  all        Run all demos
  compress   Compress a file (delegates to tpress_cli)
  fingerprint Analyze data structure (delegates to dataprint)
  hash       Tournament LSH (delegates to tournament_hash)
  rank       Pairwise ranking (delegates to pairwise_db)
  predict    Time series analysis (delegates to tournament_predict)
  sort       Adaptive sort benchmark (delegates to tournament_sort)
  codec      Staircase codec demo (delegates to staircase_codec)
""")
    sub = parser.add_subparsers(dest='command')
    sub.add_parser('status')
    sub.add_parser('all')

    for cmd_name, module_name in [
        ('compress', 'tpress_cli'),
        ('fingerprint', 'dataprint'),
        ('hash', 'tournament_hash'),
        ('rank', 'pairwise_db'),
        ('predict', 'tournament_predict'),
        ('sort', 'tournament_sort'),
        ('codec', 'staircase_codec'),
    ]:
        p = sub.add_parser(cmd_name)
        p.add_argument('rest', nargs='*')

    args = parser.parse_args()

    if args.command == 'status':
        cmd_status(args)
    elif args.command == 'all':
        cmd_all(args)
    elif args.command:
        # Delegate to the specific tool's CLI
        module_map = {
            'compress': 'tpress_cli',
            'fingerprint': 'dataprint',
            'hash': 'tournament_hash',
            'rank': 'pairwise_db',
            'predict': 'tournament_predict',
            'sort': 'tournament_sort',
            'codec': 'staircase_codec',
        }
        mod_name = module_map.get(args.command)
        if mod_name:
            # Re-invoke the module's main with remaining args
            sys.argv = [mod_name + '.py', '--demo'] + (args.rest or [])
            mod = __import__(mod_name)
            if hasattr(mod, 'main'):
                mod.main()
            elif hasattr(mod, 'demo'):
                mod.demo()
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
