;restore,'time_nosheath.sav'
;restore,'flux_radial_cal4l_DIIID153764.sav'

set_plot,'ps'
device,file='pf_tcoef2_rplot_DIIID153764.ps',/color,bits=8

psn=(g.psixy[*,38]-g.psi_axis)/(g.psi_bndry-g.psi_axis)
pf_avt=fltarr(260)
for x=0,259 do pf_avt[x]=mean(AV_PF_EXBI4L3[x,50:72])

ne0=g.neexp*1e14
e0=2*0.65/g.rmag
vthi=9.79e3*sqrt(ti0[*,38]*1000/2)
nu=4.8e-8/sqrt(2.)*ni0[*,38]*4.13e13*1.523e1*ti0[*,38]^(-1.5)
ri=1.02*sqrt(2.)*ti0[*,38]^0.5/g.bxy[*,38]/1e4
nue=2.91e-6*ne0[*,38]*15.23*(te0[*,38]*1000)^(-1.5)
re=2.38e-2*(te0[*,38]*1000)^0.5/g.bxy[*,38]/1e4

q0=-g.shiftangle/2/!pi
;dperp=(1+ti0/te0[*,38])*nu*ri^2.
dperp=(1+ti0/te0[*,38])*nue*re^2.
dneo=(1+1.6*q0^2.)*dperp

gamma_neo=-dneo*deriv(g.rxy[*,38],ni0[*,38])*4.13e19
gamma_cl=-dperp*deriv(g.rxy[*,38],ni0[*,38])*4.13e19
;pfn15=fft(pf_exbi4l,dim=3)
;pf_avtn15=fltarr(260)
;pf_avtn15t=fltarr(260,nnt)

;for x=0,259 do for t=0,nnt-1 do pf_avtn15t[x,t]=mean(pfn15[x,4:59,3,t])
;for x=0,259 do  pf_avtn15[x]=mean(pf_avtn15t[x,50:*])

plot,psn,-pf_avt/deriv(g.rxy[*,38],g.niexp[*,38]*1e20),xtitle='!7W!3!iN!n',ytitle='D!ir!n (m!e2!n/s)',chars=1.5,thick=3,xrange=[0.96,0.999]
;oplot,psn,-pf_avtn15/deriv(g.rxy[*,38],g.niexp[*,38]*1e20),thick=3,linestyle=2
oplot,psn,dneo,thick=3,col=2
oplot,psn,dperp,thick=3,col=4
;plot,tt4l*3.1e-4,AV_PF_EXBI4L[195,*],xtitle='t (ms)',ytitle='!7C!3!ii!n (m!e-2!ns!e-1!n)',chars=1.5,thick=3

yx=0.98
yy=0.035
dy=0.003
dx1=yx-0.001
dx2=yx-0.005
dy1=0.0003
xyouts,yx,yy,'simulation',chars=1.
oplot,[dx2,dx1],[yy+dy1,yy+dy1],thick=3;,col=6
yy=yy-dy
xyouts,yx,yy,'neoclassical',chars=1.
oplot,[dx2,dx1],[yy+dy1,yy+dy1],thick=3,col=2
yy=yy-dy
xyouts,yx,yy,'classical',chars=1.
oplot,[dx2,dx1],[yy+dy1,yy+dy1],thick=3,col=4
yy=yy-dy

device,/close
set_plot,'x'
