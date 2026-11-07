#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdint>
using namespace std;
using State=array<uint32_t,54>;
int dest[4][54][3], degree[4][54];bool accept[54];
int maxlen,cap;vector<uint32_t>witness;uint64_t counts[15]={};
void init(){for(int pi=0;pi<2;pi++){int p[2]={pi,1-pi};for(int i=0;i<27;i++){int s[3]={i/9,(i/3)%3,i%3};bool ok=true;for(int g=0;g<3;g++)if(s[g]!=2&&!(s[g]==0&&(g!=1||p[0]<p[1])))ok=false;accept[27*pi+i]=ok;for(int m=0;m<4;m++)for(int g=0;g<3;g++){if(s[g]==0&&g>0&&!(m>>p[g-1]&1))continue;int t[3]={s[0],s[1],s[2]};t[g]=(g==2||!(m>>p[g]&1))?2:1;dest[m][27*pi+i][degree[m][27*pi+i]++]=27*pi+9*t[0]+3*t[1]+t[2];}}}}
void visit(const State &s,int len,uint32_t word){counts[len]++;uint64_t h=0;for(int i=0;i<54;i++)if(accept[i])h+=s[i];if(h>uint64_t(cap))return;if(witness[h]==UINT32_MAX)witness[h]=(uint32_t(len)<<28)|word;if(len==maxlen)return;for(int m=0;m<4;m++){State q={};for(int i=0;i<54;i++)if(s[i])for(int j=0;j<degree[m][i];j++)q[dest[m][i][j]]+=s[i];visit(q,len+1,word|(uint32_t(m)<<(2*len)));}}
int main(int argc,char**argv){maxlen=argc>1?atoi(argv[1]):12;cap=argc>2?atoi(argv[2]):100000;if(maxlen>14||cap<1)return 2;init();witness.assign(cap+1,UINT32_MAX);State s={};s[0]=s[27]=1;visit(s,0,0);int missing=0,found=0;cout<<"{\"max_chain_length\":"<<maxlen<<",\"cap\":"<<cap<<",\"counts\":[";for(int i=0;i<=maxlen;i++){if(i)cout<<",";cout<<counts[i];}cout<<"],\"missing_odds\":[";for(int h=1;h<=cap;h+=2)if(witness[h]==UINT32_MAX){if(missing++)cout<<",";cout<<h;}else found++;cout<<"],\"found\":"<<found<<"}\n";ofstream out(argc>3?argv[3]:"census.bin",ios::binary);out.write((char*)witness.data(),witness.size()*sizeof(uint32_t));return 0;}
