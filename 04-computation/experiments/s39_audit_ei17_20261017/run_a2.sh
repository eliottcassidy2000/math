#!/bin/sh
# Runs the independent recounts of a2_embedding_search.py; outputs next to this file.
cd "$(dirname "$0")"
S1=331,661,881,991,1321,1871,2311,2531,2861,2971,3191,3301,3631,3851,4621,4951,5171,5281,5501,6271,6491,7151,7481,7591,8581,9241,9461,9791,9901,10781,10891,11551
T54=109,439,659,769,1319,1429,1759,1979,2089,2309,2749,2969,3079,3299,3739,4289,4729,5059,5279,5939
M10=43,109,131,197,241,263,307,373,439,461,571,593,659,769,857,967,1033,1187,1231,1297,1319,1429,1451,1583,1627,1693,1759,1847,1913,1979,2089,2111,2221,2243,2287,2309,2441,2551
python -B a2_embedding_search.py split table 3,5,7-700 > a2_split_controls_table.out 2>&1
python -B a2_embedding_search.py twisted table 3,5,7-600,2441,2551 > a2_twisted_controls_table.out 2>&1
python -B a2_embedding_search.py twisted table $T54 > a2_twisted_54mod55_table.out 2>&1
python -B a2_embedding_search.py split table $S1 > a2_split_1mod55_table.out 2>&1
python -B a2_embedding_search.py split brute 3,5,7-400 > a2_split_brute.out 2>&1
python -B a2_embedding_search.py twisted brute 3,5,7-400 > a2_twisted_brute.out 2>&1
python -B a2_embedding_search.py split C 103,139,331,661,677,881,991,2861 > a2_split_methodC.out 2>&1
python -B a2_embedding_search.py twisted C 109,127,149,383,439,659 > a2_twisted_methodC.out 2>&1
python -B a2_embedding_search.py split table $M10 > a2_split_10mod11_table.out 2>&1
python -B a2_embedding_search.py twisted table $M10 > a2_twisted_10mod11_table.out 2>&1
echo done > a2_done.flag
