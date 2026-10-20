#include <cstdint>
#include <iostream>
#include <fstream>
#include <array>
#include <map>
using U=__uint128_t;
int main(int argc,char **argv){std::ifstream f(argc>1?argv[1]:"input8.txt"); int b,n;if(!(f>>b>>n)||b!=29||n!=16){std::cerr<<"Expected 29-basis 16-orbit input\n";return 1;}U basis[32];uint64_t lo,hi;for(int i=0;i<b;i++){f>>lo>>hi;basis[i]=U(lo)|(U(hi)<<64);}int canon[256];for(int i=0;i<256;i++)f>>canon[i];U cur=0;uint64_t hist[17]={};uint32_t first=0;bool have=false;uint64_t total=1ULL<<b;for(uint64_t j=0;j<total;j++){if(j)cur^=basis[__builtin_ctzll(j)];uint32_t seen=0;bool good=true;for(int k=0;k<n;k++){int v=canon[(cur>>(8*k))&255];if(v<0){good=false;break;}seen|=1U<<v;}if(good){int cnt=__builtin_popcount(seen);hist[cnt]++;if(cnt==n&&!have){have=true;first=j^(j>>1);}}}std::cout<<"{\"maps\":"<<total<<",\"first_coefficient_mask\":"<<(have?std::to_string(first):"null")<<",\"target_valued_type_coverage\":{";bool sep=false;for(int i=0;i<=n;i++)if(hist[i]){if(sep)std::cout<<",";sep=true;std::cout<<"\""<<i<<"\":"<<hist[i];}std::cout<<"}}\n";}
