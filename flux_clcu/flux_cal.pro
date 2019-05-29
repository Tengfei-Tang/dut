;restore,'vexb_nosheath.sav'
;restore,'vbtild_nosheath.sav'
;restore,'ni_nosheath.sav'
;restore,'ti_nosheath.sav'
;restore,'te_nosheath.sav'
;restore,'vepar_nosheath.sav'
;restore,'vpar_nosheath.sav'

;g=file_import("data/DIIID153764.01359_x260y64_psi085to105_fac_ohthe_imp3.nc")

;ni0=collect(path='data',var='n0')
;ti0=collect(path='data',var='ti0')
;te0=collect(path='data',var='te0')

;kappai=collect(path='data',var='kappa_par_i')
;kappae=collect(path='data',var='kappa_par_e')

nnt=n_elements(te4l3[0,0,0,*])

tet=te4l3
tit=ti4l3
for z=0,63 do for t=0,nnt-1 do tet[*,*,z,t]=te0+te4l3[*,*,z,t]
for z=0,63 do for t=0,nnt-1 do tit[*,*,z,t]=ti0+ti4l3[*,*,z,t]

vbtpi4l3=vp4l3*vbtild4l3
vbtpe4l3=vep4l3*vbtild4l3
dtidr=ti0
dtedr=te0
for y=0,63 do for t=0,nnt-1 do dtidr[*,y]=deriv(g.rxy[*,y],ti0[*,y])
for y=0,63 do for t=0,nnt-1 do dtedr[*,y]=deriv(g.rxy[*,y],te0[*,y])

qbri=kappai
qbre=kappae
for z=0,63 do for t=0,nnt-1 do qbri[*,*,z,t]=kappai[*,*,z,t]*dtidr*vbtpi4l3[*,*,z,t]
for z=0,63 do for t=0,nnt-1 do qbre[*,*,z,t]=kappae[*,*,z,t]*dtedr*vbtpe4l3[*,*,z,t]

vthi=9.79e3*sqrt(ti0*1000/2)
vthe=4.19e5*sqrt(te0*1000)

vttri=vthi*vbtpi4l3
vttre=vthe*vbtpe4l3
vtoti=kappai
vtote=kappae
for z=0,63 do for t=0,nnt-1 do vtoti[*,*,z,t]=(vthi/7.412e6+vp4l3[*,*,z,t])*vbtpi4l3[*,*,z,t]
for z=0,63 do for t=0,nnt-1 do vtote[*,*,z,t]=(vthe/7.412e6+vep4l3[*,*,z,t])*vbtpe4l3[*,*,z,t]
vptri=vp4l3*vbtpi4l3
vptre=vep4l3*vbtpe4l3

;moment_xyzt,vttri,rms=rmsvttri,dc=dcvttri
;moment_xyzt,vttre,rms=rmsvttre,dc=dcvttre
;moment_xyzt,vtoti,rms=rmsvtoti,dc=dcvtoti
;moment_xyzt,vtote,rms=rmsvtote,dc=dcvtote
;moment_xyzt,vptri,rms=rmsvptri,dc=dcvptri
;moment_xyzt,vptre,rms=rmsvptre,dc=dcvptre
;moment_xyzt,qbri,rms=rmsqbri,dc=dcqbri
;moment_xyzt,qbre,rms=rmsqbre,dc=dcqbre

pf_exbi4l3=ni4l3*vexbx4l3*4.13e19*7.412e6
pf_btri4l3=ni4l3*vtoti*4.13e19*7.412e6
hf_exbi4l3=fltarr(260,64,64,nnt)
hf_exbe4l3=hf_exbi4l3
hf_toti4l3=hf_exbi4l3
hf_tote4l3=hf_exbi4l3

for t=0, nnt-1 do begin $
   for z=0, 63 do begin $ 
   hf_exbi4l3[*,*,z,t]=(ti0+ti4l3[*,*,z,t])*vexbx4l3[*,*,z,t]*7.412e6*1000*2.3 &$
   hf_exbe4l3[*,*,z,t]=(te0+te4l3[*,*,z,t])*vexbx4l3[*,*,z,t]*7.412e6*1000*2.3 &$
   hf_toti4l3[*,*,z,t]=(ti0+ti4l3[*,*,z,t])*vtoti[*,*,z,t]*7.412e6*1000*2.3 &$
   hf_tote4l3[*,*,z,t]=(te0+te4l3[*,*,z,t])*vtote[*,*,z,t]*7.412e6*1000*2.3 &$
   end &$
   end

hf_qbri4l3=qbri*7.412e6*2.3
hf_qbre4l3=qbre*7.412e6*2.3

moment_xyzt,pf_exbi4l3,dc=dcpf_exbi4l3
moment_xyzt,pf_btri4l3,dc=dcpf_btri4l3
moment_xyzt,hf_exbi4l3,dc=dchf_exbi4l3
moment_xyzt,hf_exbe4l3,dc=dchf_exbe4l3
moment_xyzt,hf_toti4l3,dc=dchf_toti4l3
moment_xyzt,hf_tote4l3,dc=dchf_tote4l3
moment_xyzt,hf_qbri4l3,dc=dchf_qbri4l3
moment_xyzt,hf_qbre4l3,dc=dchf_qbre4l3

av_pf_exbi4l3=fltarr(260,nnt)
av_pf_btri4l3=av_pf_exbi4l3
av_hf_exbi4l3=av_pf_exbi4l3
av_hf_exbe4l3=av_pf_exbi4l3
av_hf_toti4l3=av_pf_exbi4l3
av_hf_tote4l3=av_pf_exbi4l3
av_hf_qbri4l3=av_pf_exbi4l3
av_hf_qbre4l3=av_pf_exbi4l3

for x=0,259 do begin $
   for t=0, nnt-1 do begin $
   av_pf_exbi4l3[x,t] = mean(dcpf_exbi4l3[x,4:59,t]) &$
   av_pf_btri4l3[x,t] = mean(dcpf_btri4l3[x,4:59,t]) &$
   av_hf_exbi4l3[x,t] = mean(dchf_exbi4l3[x,4:59,t]) &$
   av_hf_exbe4l3[x,t] = mean(dchf_exbe4l3[x,4:59,t]) &$
   av_hf_toti4l3[x,t] = mean(dchf_toti4l3[x,4:59,t]) &$
   av_hf_tote4l3[x,t] = mean(dchf_tote4l3[x,4:59,t]) &$
   av_hf_qbri4l3[x,t] = mean(dchf_qbri4l3[x,4:59,t]) &$
   av_hf_qbre4l3[x,t] = mean(dchf_qbre4l3[x,4:59,t]) &$
   end &$
   end 

save,av_pf_exbi4l3,av_pf_btri4l3,av_hf_exbi4l3,av_hf_exbe4l3,av_hf_toti4l3,av_hf_tote4l3,av_hf_qbri4l3,av_hf_qbre4l3,f="flux_radial_cal4l3_DIIID153764.sav"
