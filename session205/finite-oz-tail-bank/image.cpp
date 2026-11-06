#include <vector>
#include <iostream>
#include <cstdint>
int main(){const int M=100000000;std::vector<uint32_t> q(M+1);std::vector<bool>s(M+1);s[1]=true;for(int m=1;m<=M;m++){q[m]=q[m/2]+((m&1)?((m/2)+(m/2%2)):0);auto h=2LL*m-q[m]+1;if(h<=M)s[h]=true;}long long num=0,run=0,best=0;for(int h=1;h<=M;h+=2){if(s[h]){num++;run=0;}else best=std::max(best,++run);if(h==99999||h==999999||h==9999999||h==99999999)std::cout<<h<<" "<<num<<" "<<best<<"\n";}}
