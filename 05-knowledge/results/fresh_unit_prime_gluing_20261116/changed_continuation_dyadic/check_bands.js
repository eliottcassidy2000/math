function run(log){
  let M0,root;
  const check=(f,msg)=>{if(!f)throw Error(msg);};
  const bit=x=>31-Math.clz32(x), pop=x=>{let p=0;while(x){p++;x&=x-1}return p};
  function validate(a){
    for(let i=0;i<a.length;i++){
      check(!(a[i]>>i&1)&&a[i]<(1<<a.length),"bad row");
      for(let j=0;j<i;j++)check((a[i]>>j&1)+(a[j]>>i&1)===1,"not a tournament");
    }
    const F=(1<<a.length)-1;
    for(let c=0;c<2;c++){let seen=1,b=1;
      while(b){const i=bit(b&-b);b&=b-1;const x=(c?(F^a[i]^(1<<i)):a[i])&~seen;seen|=x;b|=x;}
      check(seen===F,"non-strong");
    }
  }
  function dp(a,root=-1){
    const n=a.length,F=(1<<n)-1,D=new Int32Array((F+1)*n),H=new Int32Array(F+1);H[0]=1;
    for(let i=0;i<n;i++)if(root<0||root===i)D[(1<<i)*n+i]=1;
    for(let S=1;S<=F;S++)for(let b=S;b;b&=b-1){
      const i=bit(b&-b),x=D[S*n+i]; if(!x)continue; H[S]+=x;
      for(let c=a[i]&~S;c;c&=c-1){const v=c&-c;D[(S|v)*n+bit(v)]+=x;}
    }
    return H;
  }
  function state(a,r){
    const H=dp(a),P=dp(a,r),F=H.length-1;
    let E=0,O=0;
    for(let S=1;S<=F;S++){const z=P[S]*H[F^S];if(pop(S)&1)O+=z;else E+=z;}
    return [H[F],E,O,H[F^(1<<r)]];
  }
  function skeletons(a,removed=0,first=-1){
    const F=((1<<a.length)-1)^((1<<M0)|removed);
    const d=new Int32Array((1<<a.length)*3);
    if(first<0)d[0]=1;
    function visit(S,k,v){
      d[3*S+k]++;
      for(let b=a[v];b;b&=b-1){
        const j=bit(b&-b);
        if(j===M0){if(k<2)visit(S,k+1,j);}
        else if((F>>j&1)&&!(S>>j&1))visit(S|(1<<j),k,j);
      }
    }
    for(let v=0;v<a.length;v++)if(first<0||v===first){
      if(v===M0)visit(0,1,v);else if(F>>v&1)visit(1<<v,0,v);
    }
    return [F,d];
  }
  function table(a){
    const [F,H]=skeletons(a),[K,G]=skeletons(a,1<<root),[,P]=skeletons(a,0,root);
    let E=[0,0,0],O=[0,0,0];
    for(let S=0;S<=F;S++)if((F&S)===S)for(let i=0;i<3;i++)if(P[3*S+i])for(let j=0;j<3;j++)if(H[3*(F^S)+j]&&i+j){
      check(i+j<=2&&i<=1&&j<=1,"unanticipated rooted pattern");
      const z=P[3*S+i]*H[3*(F^S)+j], idx=i+j===2?2:i;
      E[idx]+=z;O[idx]+=(pop(S)&1?-1:1)*z;
    }
    return [[H[3*F+1],H[3*F+2]],[G[3*K+1],G[3*K+2]],E,O];
  }
  function substitute(a,M){
    let w=[];
    for(let i=0;i<a.length;i++)for(let j=0;j<(i===M0?M.length:1);j++)w.push([i,j]);
    return [w.map(([i,k])=>w.reduce((r,[j,l],t)=>r+(((a[i]>>j&1)||(i===M0&&j===M0&&(M[k]>>l&1)))?1<<t:0),0)),w.findIndex(([i])=>i===root)];
  }
  const TT=n=>Array.from({length:n},(_,i)=>((1<<n)-1)^((1<<(i+1))-1));
  const eq=(a,b)=>JSON.stringify(a)===JSON.stringify(b);

  const bank=[
    // actual base Hamilton count, common unsigned coefficient (per optional R),
    // deletion coefficient (per u), multiplier, marked site, terminal row masks.
    [1479,1107,15,3,8,[510,1004,1016,496,482,448,896,256,512,185]],
    [1481,1107,15,1,10,[2014,4092,4056,4048,4064,4045,1920,1792,3584,1024,2048,705]],
    [1483,1107,15,1,10,[2042,4092,4089,4016,4064,1984,1928,1792,3584,1024,2048,737]],
    [1469,1107,15,1,9,[2014,1020,984,976,2016,973,896,1792,512,1024,366]],
    [715,545,5,5,7,[254,252,504,240,224,448,128,256,91]],
    [717,545,5,1,9,[1022,2044,2040,2032,1888,960,1920,784,512,1024,417]],
    [719,545,5,1,9,[2038,1020,1016,2033,2016,960,896,1792,512,1024,358]]
  ];
  for(const [base,K,g,m,mark,a] of bank){
    M0=mark;root=a.length-1;
    validate(a);
    check(a[mark]===(1<<root),"unique exit");
    const [H,G,T,D]=table(a);
    // All scalar identities after multiplication by a real ordinal prefix.
    check(eq(H.map(x=>m*x),[base,K])&&eq(G.map(x=>m*x),[g,0]),"count and deletion coefficients");
    check(m*T[2]===K&&m*(T[0]+T[1]-2*T[2])===g,"unsigned cover coefficients");
    check(D[1]===D[2]&&m*(D[0]-D[2])===-g,"signed cover coefficients");
    const sample=[];
    for(const B of [TT(1),TT(2),TT(3),[2,4,1]]){
      const J=dp(B),F=J.length-1,u=J[F],R=J.reduce((p,h,S)=>p+h*J[F^S],0);
      const got=state(...substitute(a,B)).map(x=>m*x);
      check(eq(got,[K*R+(base-2*K)*u,K*R/2,K*R/2+g*u,g*u]),"independent DP");
      sample.push(got);
    }
    log({base,mark,m,H,G,T,D,sample});
  }
}
run(x=>console.log(JSON.stringify(x)));
