#include <iostream>
#include <cstdint>
using namespace std;
uint64_t joint(uint64_t f){uint64_t g=0,j=0;int e=0;while(f){int b=f&1;f>>=1;if(b)j+=g+e;g=2*g+b;e=b;}return j;}
int main(){const uint64_t N=4611841;uint64_t hits=0;for(int t=1;t<=23;t++){uint64_t a=1ULL<<t;for(uint64_t f=1;f<=(N-1)/(a-1);f++){auto j=joint(f);if(1+a*f-j==N){cout<<"HIT "<<t<<" "<<f<<" "<<j<<"\n";hits++;}}}cout<<"hits "<<hits<<"\n";for(uint64_t f=1;f<=100;f++){auto d=3*(int64_t)joint(f)-(int64_t)joint(3*f)-2;if(d==2||d==4)cout<<"digit "<<f<<" "<<d<<" "<<joint(f)<<" "<<joint(3*f)<<"\n";}}
