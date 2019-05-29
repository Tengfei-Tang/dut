
ni4ln=pcollect(path='.',var='ni' )
moment_xyzt,ni4ln,rms=rmsni4ln,dc=dcni4ln
save,ni4ln,rmsni4ln,dcni4ln,filename='ni_nosheath.sav'
ti4ln=pcollect(path='.',var='ti' )
moment_xyzt,ti4ln,rms=rmsti4ln,dc=dcti4ln
save,ti4ln,rmsti4ln,dcti4ln,filename='ti_nosheath.sav'
te4ln=pcollect(path='.',var='te' )
moment_xyzt,te4ln,rms=rmste4ln,dc=dcte4ln
save,te4ln,rmste4ln,dcte4ln,filename='te_nosheath.sav'
ps4ln=pcollect(path='.',var='psi' )
moment_xyzt,ps4ln,rms=rmsps4ln,dc=dcps4ln
save,ps4ln,rmsps4ln,dcps4ln,filename='psi_nosheath.sav'
vp4ln=pcollect(path='.',var='vipar' )
moment_xyzt,vp4ln,rms=rmsvp4ln,dc=dcvp4ln
save,vp4ln,rmsvp4ln,dcvp4ln,filename='vpar_nosheath.sav'
vep4ln=pcollect(path='.',var='vepar' )
moment_xyzt,vep4ln,rms=rmsvep4ln,dc=dcvep4ln
save,vep4ln,rmsvep4ln,dcvep4ln,filename='vepar_nosheath.sav'
ph4ln=pcollect(path='.',var='phi' )
moment_xyzt,ph4ln,rms=rmsph4ln,dc=dcph4ln
save,ph4ln,rmsph4ln,dcph4ln,filename='phi_nosheath.sav'
p4ln=pcollect(path='.',var='p' )
moment_xyzt,p4ln,rms=rmsp4ln,dc=dcp4ln
save,p4ln,rmsp4ln,dcp4ln,filename='p_nosheath.sav'
jp4ln=pcollect(path='.',var='jpar' )
moment_xyzt,jp4ln,rms=rmsjp4ln,dc=dcjp4ln
save,jp4ln,rmsjp4ln,dcjp4ln,filename='jpar_nosheath.sav'
u4ln=pcollect(path='.',var='u' )
moment_xyzt,u4ln,rms=rmsu4ln,dc=dcu4ln
save,u4ln,rmsu4ln,dcu4ln,filename='u_nosheath.sav'

kappae4ln=pcollect(path='.',var='kappa_par_e' )
kappai4ln=pcollect(path='.',var='kappa_par_i' )
save,kappae4ln,kappai4ln,filename='kappa_nosheath.sav'

heatfi4ln=pcollect(path='.',var='heatflux_par_i' )
moment_xyzt,heatfi4ln,rms=rmshfi4ln,dc=dchfi4ln
heatfe4ln=pcollect(path='.',var='heatflux_par_e' )
moment_xyzt,heatfe4ln,rms=rmshfe4ln,dc=dchfe4ln
save,heatfi4ln,heatfe4ln,rmshfi4ln,dchfi4ln,rmshfe4ln,dchfe4ln,filename='heatflux_nosheath.sav'

heatfi4ln_fl=pcollect(path='.',var='heatflux_par_flutter_i' )
moment_xyzt,heatfi4ln_fl,rms=rmshfi4ln_fl,dc=dchfi4ln_fl
heatfe4ln_fl=pcollect(path='.',var='heatflux_par_flutter_e' )
moment_xyzt,heatfe4ln_fl,rms=rmshfe4ln_fl,dc=dchfefr4ln_fl
save,heatfi4ln_fl,heatfe4ln_fl,rmshfi4ln_fl,dchfi4ln_fl,rmshfe4ln_fl,dchfe4ln_fl,filename='heatflux_flutter_nosheath.sav'

vexbx4ln=pcollect(path='.',var='vexb_x' )
moment_xyzt,vexbx4ln,rms=rmsvexb4ln,dc=dcvexb4ln
save,vexbx4ln,rmsvexb4ln,dcvexb4ln,f="vexb_nosheath.sav"
vbtild4ln=pcollect(path='.',var='vbtild_x' )
moment_xyzt,vbtild4ln,rms=rmsvbtild4ln,dc=dcvbtild4ln
save,vbtild4ln,rmsvbtild4ln,dcvbtild4ln,f="vbtild_nosheath.sav"


;qse=pcollect(path='.',var='q_se')
;qsi=pcollect(path='.',var='q_si')
;save,qse,qsi,filename='q_sie_nosheath.sav'

tt4ln=pcollect(path='.',var='t_array')
save,tt4ln,filename='time_nosheath.sav'