;restore,'time_nosheath.sav'
;restore,'flux_radial_cal4l_DIIID153764.sav'

set_plot,'ps'
device,file='hf_tcoef2_rplot_DIIID153764.ps',/color,bits=8

psn=(g.psixy[*,38]-g.psi_axis)/(g.psi_bndry-g.psi_axis)
hfi_avt=fltarr(260)
hfe_avt=fltarr(260)
for x=0,259 do hfi_avt[x]=mean(AV_hF_EXBI4L3[x,50:*])
for x=0,259 do hfe_avt[x]=mean(AV_hF_EXBe4L3[x,50:*])

nu=4.8e-8/sqrt(2.)*ni0[*,38]*4.13e13*1.523e1*ti0[*,38]^(-1.5)
ri=1.02*sqrt(2.)*ti0[*,38]^0.5/g.bxy[*,38]/1e4
nue=2.91e-6*ne0[*,38]*15.23*(te0[*,38]*1000)^(-1.5)
re=2.38e-2*(te0[*,38]*1000)^0.5/g.bxy[*,38]/1e4

e0=2*0.65/g.rmag
q0=-g.shiftangle/2/!pi
xii=q0^2*e0^(-1.5)*nu*ri^2.
xie=q0^2*e0^(-1.5)*nue*re^2.

;hfn15=fft(hf_exbi4l,dim=3)
;hf_avtn15=fltarr(260)
;hf_avtn15t=fltarr(260,nnt)
;hfn15e=fft(hf_exbe4l,dim=3)
;hf_avtn15e=fltarr(260)
;hf_avtn15te=fltarr(260,nnt)

;for x=0,259 do for t=0,nnt-1 do hf_avtn15t[x,t]=mean(hfn15[x,4:59,3,t])
;for x=0,259 do  hf_avtn15[x]=mean(hf_avtn15t[x,50:*])
;for x=0,259 do for t=0,nnt-1 do hf_avtn15te[x,t]=mean(hfn15e[x,4:59,3,t])
;for x=0,259 do  hf_avtn15e[x]=mean(hf_avtn15te[x,50:*])

plot,psn,-hfi_avt/deriv(g.rxy[*,38],g.Tiexp[*,38]),xtitle='!7W!3!iN!n',ytitle='!7v!3!ii!n (m!e2!n/s)',chars=1.5,thick=3,ystyle=8,xrange=[0.85,0.999],yrange=[0,0.4]
oplot,psn,xii,thick=3,linestyle=3
;oplot,psn,-hf_avtn15/deriv(g.rxy[*,38],g.Tiexp[*,38]),thick=3,linestyle=3
axis,yrange=[0,0.8],yaxis=1,/save,ytitle='!7v!3!ie!n (m!e2!n/s)',col=2,chars=1.5
oplot,psn,-hfe_avt/deriv(g.rxy[*,38],g.Teexp[*,38]),thick=3,col=2
oplot,psn,xie,thick=3,linestyle=3,col=2
;oplot,psn,-hf_avtn15e/deriv(g.rxy[*,38],g.Teexp[*,38]),thick=3,col=2,linestyle=2
;plot,tt4l*3.1e-4,AV_PF_EXBI4L[195,*],xtitle='t (ms)',ytitle='!7C!3!ii!n (m!e-2!ns!e-1!n)',chars=1.5,thick=3

yx=0.90
yy=0.7
dy=0.03
dx1=yx-0.005
dx2=yx-0.015
dy1=0.005
xyouts,yx,yy,'simulation',chars=1.
oplot,[dx2,dx1],[yy+dy1,yy+dy1],thick=3;,col=6
yy=yy-dy
xyouts,yx,yy,'neoclassical',chars=1.
oplot,[dx2,dx1],[yy+dy1,yy+dy1],thick=3,linestyle=2
yy=yy-dy

device,/close
set_plot,'x'
