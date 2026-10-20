#include <cstdint>
#include <iostream>
#include <fstream>
#include <array>
using U=__uint128_t;
uint64_t rng=20261020;
uint64_t next(){rng+=0x9e3779b97f4a7c15ULL; uint64_t z=rng;z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27))*0x94d049bb133111ebULL;return z^(z>>31);}
int pc(U x){return __builtin_popcountll(uint64_t(x))+__builtin_popcountll(uint64_t(x>>64));}
int main(int argc,char **argv){uint64_t limit=argc>1?std::stoull(argv[1]):1000000000ULL;std::ifstream f(argc>2?argv[2]:"cubic8.txt");int b,n;f>>b>>n;U bas[64],tab[8][256];uint64_t cof[64],ctab[8][256],lo,hi;for(int i=0;i<b;i++){f>>lo>>hi>>cof[i];bas[i]=U(lo)|(U(hi)<<64);}int canon[256];for(int i=0;i<256;i++)f>>canon[i];for(int k=0;k<8;k++){tab[k][0]=0;ctab[k][0]=0;for(int j=1;j<256;j++){int bit=__builtin_ctz(j);tab[k][j]=tab[k][j^(1<<bit)]^(8*k+bit<b?bas[8*k+bit]:0);ctab[k][j]=ctab[k][j^(1<<bit)]^(8*k+bit<b?cof[8*k+bit]:0);}}
uint64_t weight=0,landing=0,bestcoef=0;U bestcur=0;int best=65;for(uint64_t j=0;j<limit;j++){uint64_t c=next()&((1ULL<<b)-1); if(j%3>=1)c&=next(); if(j%3>=2)c&=next();U cur=0;for(int k=0;k<8;k++)cur^=tab[k][(c>>(8*k))&255];if(pc(cur)!=60)continue;weight++;uint32_t seen=0;bool good=true;for(int k=0;k<n;k++){int v=canon[(cur>>(8*k))&255];if(v<0){good=false;break;}uint32_t bit=1U<<v;if(seen&bit){good=false;break;}seen|=bit;}if(good){uint64_t coeff=0;for(int k=0;k<8;k++)coeff^=ctab[k][(c>>(8*k))&255];landing++; if(__builtin_popcountll(coeff)<best){best=__builtin_popcountll(coeff);bestcoef=coeff;bestcur=cur; std::cerr<<"new best "<<best<<" at "<<j<<"\n";}}}
std::cout<<"{\"trials\":"<<limit<<",\"weight60\":"<<weight<<",\"witnesses\":"<<landing<<",\"best_weight\":"<<best<<",\"coefficient_mask\":"<<bestcoef<<",\"images\":[";for(int k=0;k<n;k++){if(k)std::cout<<",";std::cout<<int((bestcur>>(8*k))&255);}std::cout<<"]}\n";}
