#!/bin/sh
# usage: sh dedupe_triples.sh n   -> data/triples_n.txt  (columns: h e B strong count), needs nauty's gentourng and ./allpairs
n=$1
cc -O2 -o allpairs allpairs.c
gentourng -q "$n" | ./allpairs "$n" | sort | uniq -c | awk '{print $2, $3, $4, $5, $1}' > "../data/triples_$n.txt"
