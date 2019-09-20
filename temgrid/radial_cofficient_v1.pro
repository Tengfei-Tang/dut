
;g=file_import('cmod_08MA_1100223012_1150_36x32y_0.9psi1.05_v1.bout.nc')
;g=file_import('hl2a_24959_260_64_80_110_d8.nc')
g=file_import('EAST_daoyuan_90_110_260_64.nc')
;g=file_import('hl2a_24959_85_105.nc')
;g=file_import('east062585.03800_x260y64_psi085to105_x108.nc')
;filename="cmod_08MA_1100223012_1150_36x32y_0.9psi1.05_v1.bout.nc"
;filename="east062585.03800_x260y64_psi085to105_x108.nc"
;filename="hl2a_24959_85_105.nc"
;filename="hl2a_24959_260_64_80_110_d8.nc"
filename="EAST_daoyuan_90_110_260_64.nc"
mu = 4.e-7*3.1415926
Ni_x=1.
Te_x= 10.
Ti_x= 10.
density_unit = 1.e20                          ; 1/m^3
ee = 1.6022e-19                               ; elementary charge 
d_0=10 ;diff=d_0*diff0,as curvature drift is neglected in this calculation

g11=g.Rxy*g.Rxy*g.Bpxy*g.Bpxy/g.rmag/g.rmag/g.bmag/g.bmag
J=(g.hthe/g.rmag)/(g.Bpxy/g.bmag)
lbar=g.rmag
unit_psi=g.rmag*g.rmag*g.bmag
jx=g.IXSEPS1
max=max(g.rxy[jx,*],index)
jy=index
rxy=g.Rxy[*,jy]/lbar
bpxy=g.Bpxy[*,jy]/g.bmag
x=(g.psixy[*,jy]-g.psi_axis)/(g.psi_bndry-g.psi_axis)
xreal=g.psixy[*,jy]/unit_psi

NX=g.NX
NY=g.NY
;jx=g.IXSEPS1

ni=g.Niexp/Ni_x
ti=g.Tiexp/Ti_x
te=g.Teexp/Te_x

dnidx=deriv(xreal,ni[*,jy])
dtidx=deriv(xreal,ti[*,jy])
dtedx=deriv(xreal,te[*,jy])

R_major=1.75
r_minor=0.5
P_input=2.e6
Elgt=1.67
q_input=P_input/4/3.14/3.14/R_major/r_minor/sqrt(Elgt)/2
chi_coef=q_input/(g.Niexp[0,jy]*1e20)/(abs(dtedx[0])/lbar/lbar/g.bmag*g.rxy[0,jy]*g.bpxy[0,jy]*Te_x*1e-19*1.6)
print, 'chi_coef for constant flux', chi_coef


print,'Boundary values', '    Core','     separatrix','       Edge'
print,'ni : ', ni[0,jy],'  ---',ni[jx,jy],'  ---',ni[NX-1,jy]
print,'ti : ', ti[0,jy],'  ---',ti[jx,jy],'  ---',ti[NX-1,jy]
print,'te : ', te[0,jy],'  ---',te[jx,jy],'  ---',te[NX-1,jy]

print,'Gradient', '    Core','     separatrix','       Edge'
print,'ni : ', dnidx[0],'  ---',dnidx[jx],'  ---',dnidx[NX-1]
print,'ti_no_source : ', dtidx[0],'  ---',dtidx[jx],'  ---',dtidx[NX-1]
print,'te_no_source : ', dtedx[0],'  ---',dtedx[jx],'  ---',dtedx[NX-1]

diff0 = dnidx[0]*J[0,jy]*g11[0,jy]/(dnidx*J[*,jy]*g11[*,jy])
chi0 = dtidx[0]*J[0,jy]*g11[0,jy]/(dtidx*J[*,jy]*g11[*,jy])
che0 = dtedx[0]*J[0,jy]*g11[0,jy]/(dtedx*J[*,jy]*g11[*,jy])
for i=jx+4,NX-1 do diff0[i] = diff0[jx+4]
;for i=jx+4,NX-1 do che0[i] = che0[jx+4]
;for i=jx+4,NX-1 do chi0[i] = chi0[jx+4]

V_conv = diff0/ni[*,jy]*dnidx
window,1
plot,V_conv,chars=2,title='velorcity V_conv'

term_conv_e = te[*,jy]*deriv(xreal,V_conv)+3./2.*V_conv*dtedx
term_conv_i = ti[*,jy]*deriv(xreal,V_conv)+3./2.*V_conv*dtidx

q_inner=q_input-abs(term_conv_e[0]*ni[0,jy]/lbar/lbar/g.bmag*g.rxy[0,jy]*g.bpxy[0,jy]*Te_x*1e-19*1.6*1e20*(xreal[1]-xreal[0]))
print, 'q (energy flux)  at the inner boundary', q_inner

Si=fltarr(NX)
Se=fltarr(NX)
for i=0,NX-1 do Si[i] = 0.0
for i=0,NX-1 do Se[i] = 0.0

for i=1,NX-1 do Si[i] = Si[i-1]-(term_conv_i[i]+term_conv_i[i-1])*(xreal[i]-xreal[i-1])/2.
for i=1,NX-1 do Se[i] = Se[i-1]-(term_conv_e[i]+term_conv_e[i-1])*(xreal[i]-xreal[i-1])/2.

window,2
plot,Si,chars=2,title='Source'
oplot,Se,color=2

chi = dtidx[0]*J[0,jy]*g11[0,jy]/(dtidx*J[*,jy]*g11[*,jy]) + Si*J[0,jy]*g11[0,jy]/(dtidx*J[*,jy]*g11[*,jy])
che = dtedx[0]*J[0,jy]*g11[*,jy]/(dtedx*J[*,jy]*g11[*,jy]) + Se*J[0,jy]*g11[0,jy]/(dtedx*J[*,jy]*g11[*,jy])
for i=jx+4,NX-1 do che[i] = che[jx+4]
for i=jx+4,NX-1 do chi[i] = chi[jx+4]

diff=fltarr(NX,NY)
chi_i=fltarr(NX,NY)
chi_e=fltarr(NX,NY)

for i= 0,NY-1 do diff[*,i]=diff0
;for i= 0,NY-1 do chi_i[*,i]=chi0
;for i= 0,NY-1 do chi_e[*,i]=che0
for i= 0,NY-1 do chi_i[*,i]=chi
for i= 0,NY-1 do chi_e[*,i]=che

print,'coefficient', '  separatrix'
print,'diff0: ', diff0[jx]
;print,'chi0: ', chi[jx]
print,'chi: ', chi0[jx]
;print,'che0: ', che0[jx]
print,'che: ', che[jx]

;window,1
;surface,diff,az=-20,chars=3,title='$D_{\perp}$'
;window,2
;surface,chi_i,az=-20,chars=3,title='$\chi_{i\perp}$'
;window,3
;surface,chi_e,az=-20,chars=3,title='$\chi_{e\perp}$'

window,4
plot,diff0,chars=2,title='perpendicular transport coefficient D'
window,5
plot,chi,chars=2,title='chi_i'
;oplot,chi0,color=2
window,6
plot,che,chars=2,title='chi_e'
;oplot,che0,color=2


handle = file_open(filename, /write)
s = file_write(handle, 'diff',diff)
s = file_write(handle, 'chi_i',chi_i)
s = file_write(handle, 'chi_e',chi_e)
file_close, handle
