#!/usr/bin/env python3
"""Generator-independent exact certificate for THM-3393.

The truth-bearing object is an embedded normalized 668 by 668 sign matrix,
not the obfuscated source that led to it.  Each row is stored as 84 little-
endian bytes: bit j is one exactly when column j is -1, and the top four
padding bits are zero.  The concatenated bytes are zlib-compressed and
base85-encoded below.

Every check uses ``require`` and therefore remains active under ``python -O``.
Only Python's standard library is used.
"""

from __future__ import annotations

import base64
import hashlib
import zlib


ORDER = 668
HALF_ORDER = 334
CORE_ORDER = 166
ROW_BYTES = 84
PAIR_COUNT = ORDER * (ORDER - 1) // 2

SEED_SHA256 = "6671d64e5118226a7b2b0e7f7be6fe52a371860984fd6152a925a0d73c5919a7"
UNNORMALIZED_TEXT_SHA256 = (
    "bdeb5059d77e2703211082627b60441b8c888c928a55cc6f295e011941a387b0"
)
NORMALIZED_TEXT_SHA256 = (
    "73f1de1539849e1dc7e6085cc69c563fd2965c44970263e8203384bd1a46aa63"
)
RAW_BITSET_SHA256 = (
    "47e5b8b061401c8cf18c2dee97f581a0c31b3ca7ad76e2930e43f1e4f18b50ca"
)
COMPRESSED_SHA256 = (
    "6847f4d26b3e9ad284e57b9284f4297856a38d74065ba2b9c972f6fed70738d8"
)
BASE85_SHA256 = (
    "c77ec3b02b7356942fe250758c1b6ac36e9bfa22202253091014dd393a003243"
)

# Four length-166 sign sequences.  This seed is used only for the independent
# bordered Goethals--Seidel reconstruction; the explicit matrix certificate
# below can be decoded and checked without it.
SEED = (
    "+++-+-+-+-+--+--+++-++----++--+----+-----+++---++-+---+++++++-+--+-+-+--++-+--+--+-"
    "---+-+-+-+-++-++---+--++++--++-++++-+++++--+---++-+---+++++++++--+-+-+--++-+--+--+-"
    "----+++-+++-++-+-++++-++-+--+--+------++--+-++-+-+---+--+-+-+-+--++---+++++-----"
    "+++++++---+---+--+-+--+-+--+-++-++-++++++-+-++-++-+-+---+--+-+-+-+--++---+++++-----"
    "+++----+-+-++---+-++----+---+++++--+---+--++-+++++++-+--+---+----++-+++----++---"
    "+++---++++-+-+--+++-+--++++-+++-----++-+++-+++-++++++++-+--+---+----++-+++----++---"
    "+++---+++---+-+--++----+--+-+++-+--++------+++--++++++-+--++--+-++---+--++---+---"
    "++--+-+----+++-+-++--++++-++-+---+-++--++++++---++++++++-+--++--+-+++--+--++---+---"
    "++--+-+-"
)

# zlib level-9 compression of 668 consecutive 84-byte row bitsets.
CERTIFICATE_BASE85 = (
    "c-"
    "pPqdtA?V<Nxu``P}ED17BG*B|<`1CZ}o%QH$kNa%@^*8cAg=r$!DdQi?gYgeX)bnL~?OT}~BiQKNEbt0+b$`T2BRyRK~>e_V"
    "h3`kQXuUhn%YKA+d;aeE2`|Np;Q-~6v%|NHm+o&LI*{v`cDC;Xnz|L5;XUl-_a>!MS?9=_cCpm^Fxw~xfW*HbdA>xKCh)oU~"
    "7PxX3#URtHleAndMlGam^Mok6FYAZ#qwU#zh9IjeE?zE-"
    "y?jFzY?EHdvv^8Z0UK<)`Xk7ek$?P+cy|vvYzqjnUqo+J;=c%bygFO_XLs#9*YP#rivTyG#S+^Qj>8_^Em>1#evLU=?cyo}@"
    "qS>%hU$ci>!_<XAi^d;}FddlV@utNmE1+vp)CT|GhYXmJXCk*v{h+!2%U#;QK?{$K3rg%zA0jw1@73L-74ar+(-"
    "wqXkO$mMw?E%+zw7YD-Iiqw3|}fXb&P&mm@Az-`1HYv=Qk%#Pg&V&Q<16l&(F7eZI);|XQO)Nq=mIVy>jS!eSXEv_&Wo-"
    "m)8lVltwI?(lAcdBkkm{yZTeu`SbP(7Tp|T&iby3S`m3k)qV;Y@7*aqymGXMZK+d4rv{((c7A?4BD5X}11?F=9jYkIe6e`wz"
    "MEO!^~ssza=cV{G-"
    "ZUtkx(mX?5&EfPxPm7<#Tsy!<b!Ts`@#*I_CL()qV;;88p_FnOpWus55uk^!(^Z(}5>Eoa8C)iPfs0k56t5cw7_Wd}QA8$-"
    "&=8&ORvi{z(5|Vnt&YG-zA&r|=-"
    "#DlV!(|6ty?FrRZGPk#!l>SZrNeu**MR{Shx@2KUs_HRz~PT4CiwenajtxB42d9A9TZq=01=c%>plovha$9HPqKA5qqlTU`q"
    "SLjdS>20HN`8E0n(_w4OaqoNjQ*aJlyL{+x1yMiv|6Xt{X~eG+Mtj(lURzyflyX{hZq&`ud%lY)6{e*$O#Z>Zc+GC@>_O7E4"
    "<^~UWLmD+@$D%bc)3t%Ro$e2Fj0kLypJx@pThgO7f0;={geD!yVo@VH})*Gw(PmR?%{<H<;3aApiKaD4e)FMp05W2&vZ>yCj"
    "fc@c&-4RZ4SV5sw}-90KEr19f9Y{UBI)i*wYt)E&-"
    "nAz|+YSc*YC&oC2Ujfv1jm_5z+o3de;2bT{yn6VHDDPk%{J1pp<UYl&wL@U#<EWdTrs;CT^v*0cgoQ(^i%0Ga_j?*h-"
    "oGlA!3#mrIwN<0@4PaW_atw~P-"
    "pvQox1MyrBJS}B2KO&&Oa~JTu`U&tf6rbA(Ko0@WtHAT2A@Gb+1`PtBlYnOc@O&BtJo6=2=L66~z%vtgE~^Eeexm4`0JIEv`"
    "T)=2&A@Y@QhNk|1_IBcz_Z&j;F+k=P641<z;g=lj5!NDJ+#`R0JI!<9s!<XW&+Qia;-"
    "N2tp}cyfv07E;JHn$Jq|#Z0?$<7xiuSj{vdt(ujIQxiRVq=S@;Tg){44a0HCvg=X~Nh19%229ajNR;(3yI+5*oUO@sh|P6nP"
    "8#4{UsPM0+b0B9fJ*+4vffv2swN(eyXfM+Q1tV#o(3CeUa038cFbAjjMTfnneqW^E)J3xB@&l$iotQdIutM(AkZNM`fc;;RL"
    "o^}R~T><Dg;5ij|R^|cE9omSY0Q3Ox><c`DE&@*vU1JCUod`VRiRT62nWC$z1E76?=Q!ed0eF_`&XodC;%Nmu`#A&8I$h8~0"
    "IC3<vB0y_`@r*=?kWM50MD_&Q{4f0{w|N+4L}9JGZ}bFqk!jYv33{$C7xcuQ<M)p+l1Os0O|rfdjn5lHSm0+=vEIviKhwh6x"
    ";xwD>b&|0F-"
    "#Dfaj`0;OV6GJPtsKr<8c+15a~Z^lSiX1w2K>(+_yY=(HCBDDhlPJPm<ofv(#c0O|%jorz}@@T`&BN&x5q;Q2oAbgc!RO=3?"
    "404)WcTY#rUFz{?q?&$zP=L1hS;AxlxJX<x6-2o`^>;XKZRshckS&$0=^#-1?z|++ccsi-CMgq{Iz;i6{v@i#rd!^Az0D2R6"
    "CIio{F~IYhNNWs0iDwD$EF1$o1BBgr0?<>ya~|>h7I<bV96bQ25Ae(ap3TJ5NAky8iGKll40uiho-"
    "f^ir;Vy=008|Fc%~80NZ{GY;M^PlIv;r6BA$7`bCWFMdjOgRJQonpd%$y)x^W!<{Rw!UB%V)zXOh8H0|43xJf{&)Z{T@DyXQ"
    "Rs8Uj4ifTyJm@C?-"
    "14hEn(!1EUHG~5n6v*exxv<i4G0G^6s;8`Kwa|eJn0MC=ea}@AwP}=4I&`{u6K|GUy=X1@>U;sJ<cs2mfIg5eka#?yU0PO}m"
    "LxJb<MBq6>J@W?uY5+WQiDy0V+#^ko1)wtESxG!^0nclynN|Q)2s|5sXQu_gGr%By8~_yq&oJP*sUz^r)c*G}=1)M^0?%Q<^"
    "O6d9`pD1q0ifRl&#u7JPYOJbiX$ci&}!hR1)e9?0MA>(wp0Mx1U%)y(|Zo^{8i968Gyb5o?_s6^f>Tbs;D{wKqG*s5O_`r2c"
    "APE=cWSC5x{d5@a$^_JiDrbjsZ~Od4YIF6HjRs0o4M}yTr2qc&-<n`woCgf#(z8`PC=DbEGn8IRIS)JX?Tg-az2FUsLrQfL;"
    "KeD}d)c2jDqDmfiqB?*UIo;JIiQ@U#|t-"
    "UgscfTua|bnpb8al$>P0O(NQsUx1ffT#ZM41n$io^s;(58$c4I|HD^b1m`A0iOE1GXUxjJTC&znpWVczdHk<8Nl-{@LW6-"
    "c<S%Y04VWXL_BrCQ-"
    "5~`K#u`W2jaONc<S%Y2q^H}1w5~Q0zBLA&fW!j7<gU<o(~OyXZzjRyFjM^&j8^0Q~^BmB~|kQ=po>l2|Sl215ZCu`b_{@20V"
    "R$r~dBj9i9V)o<{&^An-g2Ji9Feo{5U+DF8GJcuoPHF=v6Nhen$QK+A#W5#TvyCh+Vj(@q1R^}us7@U-"
    "j?JdM@bEC9L`c%}l+t=YhHgY>T-"
    "2>Lfr;&~Hz7QO<W`J!&u0q88?IUjh|&H$c%O2>Twlz5&bp0>dAgeJlnfKCRU6~r?eczVkkR{_vIz_WpP`T|cIag`8&#sSYz;"
    "8~RhJhv;;#Q=0H@XQ6Ck8c4_{oNS=JqSE!0MD>u;Hkem1E8hAGaY#5UIL!_yE6bfA9zj$o|SpPQ-"
    "5~`K#6Bx;2Crgc<S%Y0B9ocj3=JXz%xa6wGMz1&vC?406a@|(WL<NA@Do^Jp1hfo^|qW2LUMYECHUKt^?0!Vn+h{3V6-"
    "~p6Wv2`MWY=7XW<&Jney}bT06Gt!Z-"
    "rpu}@0@D$kr&j@YXVgULb@Vo*%g{HvMN#0lkK#AuA;JL~Mc<vSd{Uep_K#Av5;OYDU@Vq8G7Xv_v=S$%E{xIPAP!M4bKr4V}"
    "EAZU18+a}ev=LC^89_W7f#(oG<Cg$58+e{4o?*b#Oi}d;fKCIRKLgK<_lT$F+zS9|13W(mo<Ds}JY_+T0cbq%%m$wIqlu^ZD"
    "gh;)w}EHqTHv`_`Oas4^+0<8&xOD<_agATAaPs{K#hUtDd1Uo7kJKAMLY+fD&Xl0JcAYi&jAK)^#GK3_5+@SPXo_HZR1Y>^h"
    "4meop??Mo*we6nE*5ic)9~m1LB#Y-"
    "V*>oiKi9tlw||YQmO4V07^WKi05+PSts%&pz**{NjygY&k$wwZ~*EKJR^x`H{f|j(rp<4?Fl?xh^G#C&JZ~^15iWYxu19*1f"
    "I6SpjrUB33yHbo_9)tXM&)r5P*&bp4Py#aUSq=SDf1lK)VCaIN<rh9(Y=5f-"
    "C{35%3&GJkJ16BU#mG0J;HqS`p6~#8aFOSx?es;AsRr_e2BFNaak(dPd!&um^aK4*;Hdl61&=HZ=oJSK?_5Jm;upLe|r18Su"
    "0qp1pzR0BJg8J^dVkry=mHivymCqM4BO_?l7J<3wQ(@bnO-L)LR@Gw{qer?3ZjrYJli>sfdScvb_?-"
    "CqFDiJCR^BhAu3E{9Gl58Ei)U(=pw8Sg*#TTPwaqKF-"
    "%M|<4cnH?C>S{!2M*Jr+cey@!4uH|#gj!mz+*w|*3r>zm>(T|j)e_V@f+Pr*>+F#T9HSH~m{dir7sl#~>OUsl+SA6ZA+Wd!_"
    "4jeS7OSSP&vn=WdL{GJSJUi@-(VX>*1(?Ex=pWZz@74nuJ=$N>=2tt#TUn1hquDhpWvk)D(+~RXG=F2advoHME(T2-"
    ">iZjpb}s0f5OzFw<%Z+mI|!O7g?rTgxK6aV`;TpZO}pXJ%0BLK-Dk+I_Af1r>OK8wrz<flcX*EWnDO<x7T4zqieZ~;;){osS"
    "GJ#mpqWy*lJ>{-"
    "NmD{m>;%0v<4tcSZ7<qlD}I<TZ*IOx`m#+A3eJ1<vP_8{J#vM`EBDB*mz#?2cdc*KpMq8qOex&K_TMK3)&1kG_Uo<n>h$KYj"
    "lWyM>Zir_c7CcK&yIRpGb?4A;lwfDc6Z$RTK{0e?>B2hU+7OkF0rE&&b|HjN%W;zecWC2*2dekP4ux#F23;6zl-"
    "U>uNG7#y=>}V`a{&=Z*O$iTiB+5FlU38$#Y-nPeH6Qq7-"
    "gh`^~i3rt_1wN9teG?)7dxm|<Gtf49}{>&=OIgFe5}+LBN=H^1ygVAS}z?GI*V&e^q<ul1)Ol!~Y|>NNv^0?%Ipf#)ifh*G$"
    "*X8<VhTxtS5oeeOBizA@GbEpb<?vv$FYh+UgK!Io11Hf~FxSCQpZvsj@gMsHEVKB8u$CCjl@!SDC%M^B$!hJdxfD+Gt0?+x9"
    "SW4lFV*x1e{CG6*9Hh#l6fVIEfC5izG4SkS;72K3Z+8F+Jmb0p&#z^e!cE!^K!Imb2jDqc>_sVDDggza{!zg5TVWEVa1(t1D"
    "DbqK3p{TKyeNe`kO4q}r>PO}3{)gh3g_t$K!K;GH}JHVcu{M#y9j^+&vjYAvx_K+QaHzW017-"
    "u{!ToFCDa;y*$03U&r0CAPT@x>TpIxeo?9A$=SWQkrEsqg15n`U_9O7@t}USyu8n{K&t!Mt*+ExLDcqZb02FwZSOL!{T`;9^"
    "FUJE=;8~XeJoDull)^Rd1)#t)#0Pj*i@hj?t8oIL!1Ig-cs41MD1|Fn2|$783>oljk$6!Gx8V%{1)jERf#(X51*LHM|Hi!o6"
    "nG|l2s|ByvD6y%BcQ;uSO7ftDqJXqo7M_Ifv5jI;5lAmK`C7O>x{nx1)g?KfM+ie@l2oH1VDkOX*lpSR=QAYbh#RU63<=0bA"
    "u*{Qn>cl8Gi>#Jg)*zR~hk)@30Ag0#ETn;Ax?@pw_4-0R^6`p90UV21e8xePjbbf#-#nz_U;*q7-"
    "hAHvk2mcUytyT%CwoqY(rYcrGHI89L&byuuHF0#66x`IBx9rEo3z0F-"
    "#{1)lXfOyQmpP~!PH@O&<}pcJm&7=RMb{{x<{)UnhW-RT8DiRUKZ`C95iDcrg70F-"
    "$81J5>*1*LGy_5x7iX$L%4Dlvt7<_JK6XM`#6JTGyf)~NAI017<MZw8*TL`jsw{rjDozXAoGvpxl${e{GH%kdWg6nOSG2A&B"
    "DKT6>~eGEW>XF_k_=`P8j6z)nD00o}K_Q2Ch<V7jmzu&3(D^TLO6L=a4lPHCAS`I*o=M~_&K~O>|+#hcx{sk!T+%Xz>x+?rA"
    "h5O^J#J>Oqo*v%<PYX>3rEvENDDX_V0Xz+5UX;QWWdcy(Sz1Uu)k&1XZ3_UPz;oVQ;2CL9LMdGLYXFpZ+5yiztskXuqX{VR>"
    "|zQ$=g2cCh5O^J#J>Oqo?mN#=W%ffrEpycDDhkeJj;dEl)?=rpup3z26)yBf+>aj&wJf}1`0g4ehfTA6&aMm{pY>zKLaJ63g"
    "9_RQ$i_R_ErE&JSD)htE`$*xCNE~lz1i)PjN7%a5o7k@oWN~>xFid!le>W;JLXRcz&-iq7-"
    "iGdjJ%8j`jeab2K7K;f{O&K!K-"
    "aPvAK~E20!`G64mihQ`3Nvs^?e+#v!AJfjW+&&^^HwMPHFjQJB#;F(_rJTD8^PztwiJ^%%tenr5uMv+G;oHGFho&zrd&&3){"
    ";UWnr@Jw_8o(|e(O5qff0Vwg@4?N9u&6L9Z>9+8{K!In<BH$S<52h5Z*CYT6JSUz8o&{<<O5xlIDDdn(9eCDAjVOiNMnHk5i"
    "7oJKQi-TF>a`z$0#8)}@C=uVD1{sA0ziT1SHA+!vm#94;s_}4%=-;^&Q#`6YxL*q?Ee4-p7*8!&tsZuO5wZ-DDnIncuv&@Q)"
    "_g*2LL6Wp90Umayv@lJ{=1{iRV9nXS_O=Qn=z+017<QbAV@&G>=la1S<dvJb!Ejp7%w5l*09P2cW=nK^5=}R$>bG=j-"
    "hM00o|d9s|!DjTfbGsRR^w?hFH-"
    ")3r&I!cFu6puqFWH^9?Y=S3;pUp{Z~2T<VopaFO$$df3A8|M!|f#=g(z_VEGMXgbM5da0A%iMwIeS;)Q;cVjpDDXTh1)j~?5"
    "^9ZR^Z}s6(*k%dllxH$caDGp&(95k=UH(ErEp&!2B5$*+ZuRQ2}>x2`|AgS{tXm(`ep&o$AW4~;ocBX;MuPNc!mjrDTRCK4n"
    "Totr;mW=FhvHXaLow-"
    "6nJh*0iInoUX;Ss^ar58bJSDfDNCZ(s9+WV1)fPSiKp6&Qn(H00Vwdiu>yFmHL#%8=x?7``xj8)8F&tOUev}?3fGT-"
    "0#Exlz;m|Tg<7L&tpF5w?)({e9#<z(3ir29to;ip@Vqh}c$OOw&-"
    "B?%02FvWSPwkwwSLqZU9JY8#Ipu?hRW?Ih5Pyv041JHz%y4JORZ72O#l>lwzdGzN`pLVjfN3W;2F^hJR7yu)EX_e0ieLssSS"
    "9B>6$5p`};>K+kpa4^EbdVSJzCfQ62#Wo-xF;QrAo=+}ck7DDW&Go{hR-O5x-M0F-#{1)eW-"
    "8I;1QqX8)K{2X|`l9x~lXKV&QiRb?T&(~@{YK<%%04VX?1U%cM8I;0(xCnp}Pk-"
    "R~MpQy6+&ThEJnewzN~IsAaO#Ty6nL&Q1)k?68PpmXt^uIHv;FSu?>uLTyeNfx=QF>0pun^J?(FY8`wJ~7g&VH`pukgqclHj"
    "?1Vt>hMmstHP~fS*J9~$xyTpZBBfIVZ6nN_I&fei^C9<Ft?%fZizXeJ>&j3#&A@MwQArgQR&r0CAK_H?O?%fZizXb|BcX$F%"
    "SH&7i;okjF`dgsD)8ioUw9pXG?VSlI@JuNMo`y0^;p{yCDDW)(fq1IDD24kh5r6{Ey1BqJ(!hdJxT*mF6nKW%0na>bETwQc1"
    "Qd7<H3gn?<SvxLK{YWW9(eB70MFxM3rgYc6Hwy04tSOeiRb=80t!5P)d0_W!5T{8pqi*l1)keJ2A-"
    "h`7fRuvnkbn}JQcumn1*=9xdj1G;wb^1U1gZU#X;6Hm3Sr*Pw^T`;VA1#C!S5fbG<N+QaH+bW&qF4<-"
    "qfM1@YW)2(q3tODXIDo^v#q!cBs#Cj@v}_S{EdPcWr$sRR^w8X8Zau*Z&4ILdlTfM-"
    "<iLlpKHQ3`hqvYuq%nQu}?VNVFnZY|MUk?p+y;OZCp2cx}i^Fi1Y{NW9ICpQEf_8lZQUzB|-)_H7Y>W~(NJ=J-"
    "+dMn~9!Ak;0=pT$cu+ulg`cjir8)Z`RJhLoq=UUxq-vwj$B?sL!TcI(bI&YTV%IZ&Zax%N=AB@<3)0N#FsYx2^s+#x8=a(s0"
    "F3Rtuo$N8eD|^1yQ5#9KTgUWPoHOs-"
    "@rltt7~#%K53X*cCdtU+tE|^HOOHIbD{j0s&B(g_6m<DCyLI!el~4Bh1T`LQe=ve8eoF!_VKleZr6OYcnaN9@2*VciQ^x5}A"
    "x2(Hvs;00tqi<wGdS#4`|piapPtCctfnSulUw(NPVV{#vnu!GcF8#XDHMo9Xm-"
    "o~t(DG!+rJt1Yy0mF=gjguK8vYInv%Tb)ZPR72Xmovn!CtKe+o51dz#(S|0>hi7hl=^oBsEP**>4o8jqf)CdtOj?e?|y2Xi+"
    "l&1!Y5{uG)NaWuPi4R}tcCdm|dhD-"
    "8kcIzSWq$Vi?cn%lU((G0f@cfaQB!A%9P56vvx59zv0&0>}z%xb>LbF@vfaf4;lH!5q7>Pa2Zv6~AcT$ryA9z}dOlWqCcwV6"
    "<DGqoVDkEuji+DbuCTSk<jFPz0?3NGkd`eA{G4OO%^`O};8{oN&nxtO9(?S|Yvs)&>^DH$<I^elgluxr;8-"
    "ZsPHA#Dcr<;&?b|apjQ<Jn3cqR*;(d^bf;F(QLk`wUsQiRa#RtfO*r6%b-"
    "@a!$gqS>tlz_TAUNh^V;iKvWbw*~{xPShlwBc94yn%(*acy6L5sTFvBB?+O~E#f(fnxt^x>7vS_*)8IkL`~8O;Mv`vjAplp="
    "M8F-s)6Sg*<704dIdZKsYx0IJl)h;G`sa0c-mt$2RxGvCeiE`@!UyG(k<XwBI`l3Tg3AUHA(frbDr9SX1BTm&j-"
    "{bl>^T#0}0J;MFY>L)Fd?l&kC)CX16W@&zICBRRhn3Itk5g%?6&W)FinA&r>=H&29|<o)Od}DS_v0orGq$h^G@ZN!x&DgFKR"
    "Ew}_`XHAxSF=X0?u&2CA6XADMj!1I-"
    "|2hDB?fM)?UNqd3kYmEubZmj~IpHPz&0X!pQku<yI0z3y&lQaT&I;n}LU_bEeOifZQ@Z2lyL9<)LQ$tNs81TF%GNIWm;whsh"
    "=@{?~P)5@1)(5~-Oij{E;F&4GRk$6%b2T+d-vLh_RS%lo`X}&orY7kX@U)S}(d^b-"
    ";JJ^Qr02kMyDFb%w=#g|b!w761D@^%#M6*?22zu>19%?LKBL(!;%QG!(o^8+C4WY<Tf}oGHAw-"
    "$^RPOEX19pv6>5@VfM=P3J<V=C1fHK#lQa`}{vtD>*{vnOGn1O6v%qtyI+A9$a)GB0MsvXPjMSB8w<>|B4K+y_z_U`-"
    "gJ!or1D@NdNtz2hgQO-jyJZhNi?K-"
    "pp4lo1&2H@kp7*Iqx&b_WrML=5JcFr8`UH5Usq$%ds}^|XP?IDBo>L6w((KkE;CTn5IpBFjmPNB$r-"
    "5e>HA%OD=Vb9Dn%$}bo`b1LdJa5Ol|5*7>oM^BhMJ^>!1JcYgl4yv0?#YdBy|Cv^R<yQyEPPeKA<LP5bzu%&!^cf;<<#Hq&~"
    "oNr+6;SZv7K@=1`N=5qMq^X3^}{T;O?!nxs_V`9Ls<X16kcXCpO9H-"
    "YC8K{Cy5`2x={YLb2do<kIqXm%?Ncn+f`sSJ3UX?oD?7V+FoP0}OcsWqY5tvcX&6{9)gsguy`)^p$)KuuD2;JH>8NwZtafM+"
    "H(Nh08RQJzn;Tf>2;4>d^=;CWYEOS4;M!1E|INlM`P<X?AZf1f_s3_NdPGzUCe1okw$H5z#SiqRbKTp=)_*)8Jv8#PHAfTyD"
    "%l4iFi0?#j~Nh$%J=7M~h-RccIM^KY=pLi;2X?BZvcB3ZAop@?OXm+a*c*bBfM?AImG`saF@Ek);Qa$mM$I<N8KY-"
    "`A)Fk}`JXfjnX?E*9;Q28%Nx{JLf<Z0KZW#biYig1X0MEPHXEeJN2|VMdNh$)Ki{v3RyX68ri>OJO3_MS(?P+$)5_tMklavZ"
    "P?--cS?AA8mX-7@cB;eU7i=^2tci?GCP10fD8K!on*{uV>Q$tNs3GmF7_Mq9V>%emzHA()!vr-jDvs;gV=SXUjHUZBdgM6CZ"
    "3I(3ssY#j&JhNrQvs(`E+(J##G2rQ|enzuf(}AZOHAz1M&opTW&2ACTWNMOD0?%8jESlXSo+Z>Iod=!^49aMBi+I*ilN13wP"
    "s(a(b}ImQhES6<9C%JshtTZSK;U_XnxxafGfkRBvs+z&=L~9+76Z>)qB5G@(g9CfYLYB~=K|$in%$BE&jf0cz5<>nC0R7PB?"
    "O+u)FinB&uOAbG`qD1c=}V5)CqX@Q}&?Qt?z-S9W_Zdz;nCAgl4zy0Z&tEk~#xVcael<w}OGEhMJ^-"
    "z|%@7q1mk=z*9y|k|pqr6-a1yi+GBuNm>s)#|k7gyQKl1tEoxy1fD$vku<v{0iGAAN!kTGjTEjlyCndg*d&D$Pe~7&-"
    "4X!LMbsn>AfBo?n%yFv4%8&=Bc9THn%!ChJk6*{`Ve@o77<V3hrm-yO_B(BIxDkic8hq*sY!|ip8F(aG`m#-"
    "JjK)`xdG1!qFS2WS^zwS)Ff>Ip4Q3`n%x=%JOu=F4)BbVWYO#v@mxhs(i-"
    "47PE|&;Ti1Z6Gc`&1z|%@Pmu9yDfagAHl70Z5u_}9--8un0uTzsW5qOR@FrnG4X}~j(nxquqnXHx2>{fr^X-"
    "`ekKY^#0Ttc&39f9XgYLdKx=V5gu&2B{j&nwg<9S5Fe2Kh9*l@B~Wr6y?(@cczqOS4<EfoCQ)No#@UQgsN;ZXE}nKG-"
    "Ay&ok02n%$}Zo;K7ZeGNP-Rb@20^*QkDL`{+}@C=gH((Kk?;JJyKq-"
    "@}st$IeYTRVW~C2Eoyfv2xDgl4y{0#83`lBNUCG?hKgZaoB^C#Xpp0z9V}#L?^)@hqn%DFb*Ok-"
    "5_BRxa@TjGClc;5k{GOtV{+!1DuYl70Z5smd~%-Kqzk-"
    "%yj}1w3zRYH4;W6nI{tCdmqT&e!59+z{aTfSM!~@Ejz!r`fGNz;g*TNgaXbPH`N~Ze0VOIn*Te0iIWcIJ-"
    "qW?@*Hz4?G_Tl4*A91n_L6Cg~U8xkNCDX1Bb7XBcHYH-YC6MKaB9r324l)Fk~1Jk2!3bLeE?xtp?{ABm^7mS(s50MDxw%>@x"
    "p-810nNjw86>*)$S*Xo`D&x6D>leVyi0?&)`XTWn}0r2#pEvziyc~|@lcpk2yu!pjqI^g+4cnF5-"
    "HqO#_X_Bnv)kC!UM$vFzrSY|<)4rdNojq-9pPIE!76-CkS8SYQYJP`)c!@U*)m=)PUUe}@R57B-"
    "Y`wlwR9tV|@UZ1hT6T{CX>RdN7tObNRYWWtmAY9MO+UQS8iwksr`SH8Jy_VVCp>zjzEP9}zQ|wF+IY*?C~-"
    "=1QOn&Jx84iSpR$;ISw4oMIZ+}E)h+Iw5O(~VRiWdD7wp$JiX#0fw1xfJPwCOu%YVg_G0AbKX8Gz*VJ$UDl0q1&JH5SlSouf"
    "JxxMz(=o>{tq5c$B=Ki)_lG5?8onwoo*SOoq+E1Z7HA%`q7^<88!~L%HJ@pS}vA$7^66jCieB~GJqEh{XX*KM9py9{%Q&>dN"
    "oXFwb9mU@1wawa4WBr41&^L;%?WZs+=sT;`^Yjm<Eh?_$`2zhZEILh1lC<-"
    "?JBlV#_1UCS|6t7Y39Uu@DICkLioKAfe=sXu$IW|nQhy2#W)#h-"
    "qkw1MT;Mqsc$x)KG$(Tdp7Aq)XJ6o{&8BEhV+}luY=CDx@La!~qB&(};OQR+Jd1$m_mvdQiHN6N5%9bZJnyAbG$$<to~Fbzn"
    "0OAxXf6<VZhithbBL!eMRQUI;5j-"
    "RcuohNchV@DQxQ+gbHLLUcs5R<Xink=JPm&ao(aG+><C43N^9U5)c`z;foJYgisr<Lz%xGtc-"
    "|+TXQ)Y%7XnYep}?~lcm|EYXzmg497sHu0nhB26wT?vfoGxyc%CJmkEu!0od%vBvj3Im3yS7+-"
    "oUdb@q7$CZ{4OQNoNf_w|x&h!+__n4HV7k5`pIrbAabC;Q8AvisrP1z_T_3cy<MzUlda`rwIg}!9M{{E%5yAK1Fj12jH1QJm"
    "tW%s)V9BVIuIHPCVBF&(CqlOXLPTZHeba;F;ZnqB-dy;F&-?X9G`PBZ}tK<-"
    "oJ}GvIj~c=j7h(VXlN@bv!xc$Ndt?UB?ZX<C4%-45Vc4?Nvns7X?^0M8wsz%vwhT3w`QPS^rGJ$eDpT;LhIilRAD3-C-"
    "Ko|V9J?E4hWNm_tssTz0&5zhe_&4mNcIyvwhOgz0QniD$!&u3!b`3>+?Zl-8X-WhoQE(V?-"
    "0Z&O!YLawO!1MK5;Mo&+iY8Jtrz-"
    "@WZ65+pW8k^Ej+!LxBjEXlcq)PC1sw7+XaSxpiKhg3&aR+nPE9<W<^xX=@ob=IPTm4M&9i~05O^kTrzS}k4m@LS1J6~!(_;Z"
    "eb6N-BS@0NmIs?y?n-t9%Bm&QxFyOflc$OAXG^Z{Eo=wE_I`Eu#iJBx?An<G<o`JwKD~X~x&1vA-"
    "S`R$4fakQ2DVkGy1J8)(z_S8)rWH~&C$<KjPQ<eTcuq;CXinZ4c<v>hp}@1R5k+&lDByXGc;*7nj`Jv*)4BoAfRBJ@CGgx>M"
    "op6D5b(_G0X!Rl=cr#Onp2hoPoECJ^9As<xJl8Rh<Ms;1fI)*=hjqek|d{rXQxZRa|G~o>qya@q8xZ`BA$DQ=QxVygpYvdsH"
    "eblJn-~7K+&AI1$ZVc1D?Hr=V3dH=E8yJ4M*T<3_K^9Qj=uh06YWt0?!S=Gu4BlIdvlN%p#uoz;p6tisod6z_WsQRs+vI-"
    "%>PZ5C}XQ76Z>F;2FP|qB+%R;Q8DEc!mSdq84hBq~5@Dc{kuW9C-"
    "Sl#UZalz;i?l@Z3W@BPg1amIKc{mw@Lr;<=ZiIn^WJd950FJ_Med$5WFe4F{e9#B&MpTu)7s$^m$04h5b=fM?A}YLX<Kfv1l"
    "-@H7LSO{*!IQz(GvQ62D%CZ4k?niE6;&s$#s&jR4N;y6WfLO0;~t26Md0iKQy7|k65o=f)u&&9ygd^kpP-oSGx@pJ&5I!BD="
    "h-"
    "cSS;AsXt<tr(g6D0yq?M>jR1)giqQ#2<j1fJ`O=X&6I(VU_=WgzezNjyga&%5KPNfJ2#&;6Z&=YHV%WIaW5($2thLNxH406b"
    "ePDVkG90Z;2O!1ECBZ2g#;B$*rVjN1=9%YbKF6g5d2Yv5UQ9eB<so~{(lDLVsCf8sd^csf;6G$$gSb{_!GF2r*yMRU?Z;Av_"
    "KJii8>F$#+2RDr;Avj%vM2A*T)QZy%Z0G^{Gf#<ixGnS$`74fum0iHL2r{O$`<|J;w)6f!l1`^LWisqEoz%$Aac-"
    "j+B9P$!%2A=uEvkUOFI6%>yGzxh7T?U>8#M6kPIdviM9C#0St^=NK^C+5=2LjK;C%|(g@Jx=QXij$;czO`e?!eP)97S_F;@L"
    "9-"
    "cy<7uz3)>rr?Up0#!kR9ig;o)r%ME$8+HTFeBz1GoVE~n=6eFqYU25fqB%_<@bt3+o=w2htAe6Ag#++B(E)h20MFhH6wL_}f"
    "v5LI;JJc$o}y?@<OV!##sE)8;`ts$bJ9b=b9;B-"
    "xfggIz#%VnIq)p*06fPN&rpiyWRHNS|3={13wX}UrD#sm0zB<T0Z(J#ne`<_bBY$=xr2CaAf9a$%?Vq8r-"
    "w1{bOoLZ&r>uf4hNnoYT#)>Jk2SZlRE&<QaSK61fIA1Q#7Yd1fF$b;28xx8>Uh;rzr%U|GGQ-"
    "8&6l@`TSRk<`j>B=kG$`X#qT6{YKH8ARKtUCZ2}C^Yu=O<^&GFGr|RUDu8E1KZ@proq?y*x4=_EJU_r_P60gkmI6-"
    "^@SOD?MRO7b@Vqt;cnX2%v3?ZIDHOo-"
    "A@LLd&mX5yG$)7xo=eJs=PKa2;0Q%?!b0FVq#k%W1J9EksYwzA0#CCL;JFWYR;VbNlQ;lR?P=h79e6I3QZ%O|p7JW-"
    "83;U2`B5|{b_1T`$G|fSc=~RoXijzrc&`2pcvb+<GzCR-nsVTIAsl!%5YM?3%_$xM&)JT^GZc9CEum;m5Dq*C>;|5>#Pcvkb"
    "HdZWGtm=xRszqBCKSzyy@96(@oWU1qxw*jB(nydDa7*y@U-YiO_D((@GK>s%ZcZFissaAz_ac$@Eid=-"
    "ELBoB)0~hA@_ji9^je0k)k<mXW)6}3Gf^bJiYuWn$swN=Zr6aXD{N}Oihw93V7O%0G`IcbJ9-~&4~(uXTmPvxdC{lZl-"
    "8X5(qrq{|P+vf#+lyMRQ6A;AvG0Jgb3cpHC>76A@3NVBpyVJmXJLG$(Zfp2}gsGaPspIbbw*2zW+z1)jr+=L(AEq~5?YZ!_@"
    "R13c}{Q8cHr2A*>+Q`iGMO=nRwC+!S82iyam4~b_CHA$){;F;JAJeL5^%VQ{-"
    "lNJI`kC(u6h?c@0issahfM*KvGy|SZA5t_Y4+owTTY+aZ@NC&f(VUie_KpCa1;ld^MRNw;z|(|y)&Nh(7K-Ln{|}=O_CN"
)


Matrix = tuple[tuple[int, ...], ...]


def require(condition: bool, message: str) -> None:
    """Optimization-stable replacement for assert."""

    if not condition:
        raise RuntimeError(message)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def decode_certificate() -> tuple[tuple[int, ...], bytes, bytes]:
    encoded = CERTIFICATE_BASE85.encode("ascii")
    require(sha256(encoded) == BASE85_SHA256, "base85 payload hash changed")
    compressed = base64.b85decode(encoded)
    require(
        sha256(compressed) == COMPRESSED_SHA256,
        "compressed certificate hash changed",
    )

    decoder = zlib.decompressobj()
    raw = decoder.decompress(compressed) + decoder.flush()
    require(decoder.eof, "compressed certificate ended early")
    require(not decoder.unused_data, "compressed certificate has trailing data")
    require(not decoder.unconsumed_tail, "compressed certificate was not consumed")
    require(len(raw) == ORDER * ROW_BYTES, "wrong decoded certificate length")
    require(sha256(raw) == RAW_BITSET_SHA256, "raw bitset hash changed")

    rows = tuple(
        int.from_bytes(raw[i * ROW_BYTES : (i + 1) * ROW_BYTES], "little")
        for i in range(ORDER)
    )
    require(all(row >> ORDER == 0 for row in rows), "nonzero padding bit")
    return rows, raw, compressed


def row_text(row: int) -> str:
    return "".join("-" if row >> j & 1 else "+" for j in range(ORDER))


def verify_explicit_matrix(rows: tuple[int, ...]) -> tuple[int, int, int]:
    require(len(rows) == ORDER, "wrong row count")
    require(rows[0] == 0, "first row is not normalized positive")
    require(all(row & 1 == 0 for row in rows), "first column is not normalized positive")

    text = "".join(row_text(row) + "\n" for row in rows).encode("ascii")
    require(
        sha256(text) == NORMALIZED_TEXT_SHA256,
        "normalized text matrix hash changed",
    )

    checked = 0
    minimum = ORDER + 1
    maximum = -1
    for i, left in enumerate(rows):
        for right in rows[:i]:
            distance = (left ^ right).bit_count()
            minimum = min(minimum, distance)
            maximum = max(maximum, distance)
            require(
                distance == HALF_ORDER,
                f"row distance failed at pair {i},{checked}",
            )
            checked += 1
    require(checked == PAIR_COUNT, "wrong pair count")
    require((minimum, maximum) == (HALF_ORDER, HALF_ORDER), "distance range changed")
    return checked, minimum, maximum


def circulant(first_row: tuple[int, ...]) -> Matrix:
    size = len(first_row)
    return tuple(
        tuple(first_row[(j - i) % size] for j in range(size))
        for i in range(size)
    )


def transpose(matrix: Matrix) -> Matrix:
    return tuple(tuple(column) for column in zip(*matrix, strict=True))


def times_back_identity(matrix: Matrix) -> Matrix:
    """Right multiply by the anti-identity permutation matrix."""

    size = len(matrix)
    return tuple(tuple(row[size - 1 - j] for j in range(size)) for row in matrix)


def negate(matrix: Matrix) -> Matrix:
    return tuple(tuple(-entry for entry in row) for row in matrix)


def sign_text(entries: tuple[int, ...] | list[int]) -> str:
    return "".join("+" if entry == 1 else "-" for entry in entries)


def verify_sds(sequences: tuple[tuple[int, ...], ...]) -> int:
    require(len(sequences) == 4, "wrong sequence count")
    require(all(len(row) == CORE_ORDER for row in sequences), "wrong seed split")
    require(
        all(entry in (-1, 1) for row in sequences for entry in row),
        "non-sign seed entry",
    )
    supports = tuple(
        frozenset(i for i, entry in enumerate(row) if entry == -1)
        for row in sequences
    )
    require(
        tuple(map(len, supports)) == (82, 83, 83, 83),
        "SDS support sizes changed",
    )
    require(tuple(map(sum, sequences)) == (2, 0, 0, 0), "seed row sums changed")

    checked = 0
    for shift in range(1, CORE_ORDER):
        differences = sum(
            sum((x + shift) % CORE_ORDER in support for x in support)
            for support in supports
        )
        require(differences == 164, f"SDS difference count failed at {shift}")
        paf_sum = sum(
            sum(row[i] * row[(i + shift) % CORE_ORDER] for i in range(CORE_ORDER))
            for row in sequences
        )
        require(paf_sum == -4, f"PAF sum failed at {shift}")
        checked += 1
    return checked


def bordered_goethals_seidel(sequences: tuple[tuple[int, ...], ...]) -> tuple[str, ...]:
    """Direct block reconstruction, independent of the source renderer."""

    a, b, c, d = map(circulant, sequences)
    br = times_back_identity(b)
    cr = times_back_identity(c)
    dr = times_back_identity(d)
    btr = times_back_identity(transpose(b))
    ctr = times_back_identity(transpose(c))
    dtr = times_back_identity(transpose(d))

    blocks: tuple[tuple[Matrix, ...], ...] = (
        (a, br, cr, dr),
        (negate(br), a, dtr, negate(ctr)),
        (negate(cr), negate(dtr), a, btr),
        (negate(dr), ctr, negate(btr), a),
    )
    top_left = ("-++-", "+-+-", "++--", "----")
    top_right = ("---+", "--+-", "-+--", "+---")
    bottom_left = ("+++-", "--+-", "-+--", "+---")

    output: list[str] = []
    for i in range(4):
        output.append(
            top_left[i]
            + "".join(top_right[i][block] * CORE_ORDER for block in range(4))
        )
    for block_row in range(4):
        for i in range(CORE_ORDER):
            output.append(
                bottom_left[block_row]
                + "".join(sign_text(blocks[block_row][block][i]) for block in range(4))
            )
    require(len(output) == ORDER, "construction row count changed")
    require(all(len(row) == ORDER for row in output), "construction is not square")
    return tuple(output)


def normalize(rows: tuple[str, ...]) -> tuple[int, ...]:
    """Normalize by sign switches while preserving row and column order."""

    h00 = 1 if rows[0][0] == "+" else -1
    column_factors = tuple(1 if entry == "+" else -1 for entry in rows[0])
    output = []
    for row in rows:
        row_factor = 1 if row[0] == "+" else -1
        bits = 0
        for j, entry in enumerate(row):
            value = (1 if entry == "+" else -1) * row_factor * column_factors[j] * h00
            if value == -1:
                bits |= 1 << j
        output.append(bits)
    return tuple(output)


def main() -> None:
    rows, raw, compressed = decode_certificate()
    pairs, minimum, maximum = verify_explicit_matrix(rows)

    require(len(SEED) == 4 * CORE_ORDER, "wrong seed length")
    require(sha256(SEED.encode("ascii")) == SEED_SHA256, "seed hash changed")
    sequences = tuple(
        tuple(1 if entry == "+" else -1 for entry in SEED[i : i + CORE_ORDER])
        for i in range(0, len(SEED), CORE_ORDER)
    )
    sds_shifts = verify_sds(sequences)
    constructed = bordered_goethals_seidel(sequences)
    constructed_text = ("\n".join(constructed) + "\n").encode("ascii")
    require(
        sha256(constructed_text) == UNNORMALIZED_TEXT_SHA256,
        "unnormalized construction hash changed",
    )
    require(normalize(constructed) == rows, "construction/certificate mismatch")

    print("THM-3393 explicit Hadamard order-668 certificate")
    print(f"order={ORDER}")
    print(f"encoding=rowwise_little_endian_minus_bits row_bytes={ROW_BYTES}")
    print(f"raw_bitset_bytes={len(raw)} sha256={RAW_BITSET_SHA256}")
    print(f"zlib_bytes={len(compressed)} sha256={COMPRESSED_SHA256}")
    print(f"normalized_text_sha256={NORMALIZED_TEXT_SHA256}")
    print(f"unnormalized_text_sha256={UNNORMALIZED_TEXT_SHA256}")
    print("normalization=first_row_positive first_column_positive PASS")
    print(f"rows_checked={len(rows)}")
    print(f"off_diagonal_pairs_checked={pairs}")
    print(f"hamming_distance_min={minimum} max={maximum}")
    print("gram_identity=H*H^T=668*I PASS")
    print("sds_parameters=4-(166;82,83,83,83;164)")
    print(f"sds_nonzero_shifts_checked={sds_shifts} PAF_sum=-4 PASS")
    print("construction=4-row-bordered_Goethals-Seidel_core_order_166")
    print("construction_rebuild_matches_explicit_certificate=PASS")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
