#!/bin/sh
# c1 driver: direct Teichmueller-rigidity tests (second-pass audit of THM-4750 5'(a)).  BIN = compiled c1_rigidity_direct.c
BIN=$1
$BIN 7 3 selftest; $BIN 13 2 selftest; $BIN 3 4 selftest
$BIN 3 9 lifts 2
for a in 2 3 4; do $BIN 5 6 lifts $a; done
for a in 2 3 5 6; do $BIN 7 5 lifts $a; done
for a in 2 3 5 10; do $BIN 11 4 lifts $a; done
$BIN 13 4 lifts 2
for a in 3 5 12; do $BIN 13 3 lifts $a; done
for p in 17 19 23; do $BIN $p 3 lifts 2; done
for p in 29 31 37 41 43; do $BIN $p 2 lifts 2; done
