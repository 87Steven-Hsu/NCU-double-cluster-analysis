!p.background=255
!p.color=0
!p.charsize=1.2
!p.thick=1.3
!p.charthick=1.0
device,decompose = 0
!p.font=-1
red = [255,0,0]
green = [0,255,0]
blue = [0,0,255]
TVLCT,red,green,blue,1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;讀取資料
;test = read_csv('D:\學校資料\天文所陳老師\h and x\GAIA DR3 data\h_and_x_3deg_DR3.csv')
;
;;Transfer CSV transfer to sav.
;SAVE, test, FILE='D:\學校資料\天文所陳老師\h and x\GAIA DR3 data\h_and_x_3deg.dat'
;restore,'D:\學校資料\天文所陳老師\h and x\GAIA DR3 data\h_and_x_3deg.dat'
;
;;Transfer structure to array
;RA = [float(test.FIELD01)]    ;J2000 RA (deg)
;DEC = [float(test.FIELD02)]   ;J2000 DEC (deg)
;Plx = [float(test.FIELD03)]  ;PLX, Absolute steeGmagar paraeGmagax (paraeGmagax)
;ePlx = [float(test.FIELD04)] ;PLX error
;PM = [float(test.FIELD05)] ;PLX error
;pmRA = [float(test.FIELD06)] ;pmRA (mas/yr), Proper motion in right ascension
;epmRA = [float(test.FIELD07)];pmRA error (mas/yr)
;pmDEC = [float(test.FIELD08)]  ;pmDEC (mas/yr), Proper motion in declination
;epmDEC = [float(test.FIELD09)] ;pmDEC error (mas/yr)
;Gmag = [float(test.FIELD10)] ;Gmag (mag), Gband mean magnitude = brightness distribution
;eGmag = [float(test.FIELD11)];Gmag error (mag)
;BPmag = [float(test.FIELD12)] ;BPmag (mag), Integrated BP mean magnitude
;eBPmag = [float(test.FIELD13)];BPmag error
;RPmag = [float(test.FIELD14)]    ;RPmag (mag), Integrated RP mean magnitude
;eRPmag = [float(test.FIELD15)]   ;RPmag error

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ra0 = MEAN(RA)
de0 = MEAN(DEC)

;;;Coordiante of equatorial system ASCC16
window,1,xsize=900,ysize=900
plot,RA,DEC,psym=3,xtitle='RA (deg)',ytitle='DEC (deg)',title='Coordiante of Stars Observed by GAIA',yrange = [55,60]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; NGC869 = h
radcen869 = 14.4/60. ;~0.24 deg (radius)
ra869 = 34.75
dec869 = 57.13

;; NGC884 = x
radcen884 = 10.5/60. ;~0.175 deg (radius)
ra884 = 35.60
dec884 = 57.13

;; circle out a region of cluster core
window,2,xsize=900,ysize=900
PLOT, RA, DEC, psym=3, xtit='RA', ytit='DEC',title='NGC869 and NGC884',yrange = [55,60]
TVCIRCLE, 0.24, ra869, dec869, color=1,THICK = 2
TVCIRCLE, 0.175, ra884, dec884, color=3,THICK = 2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; calculate distance of all stars to (ra869,869)
rad869 = SQRT( (RA-ra869)*(RA-ra869) + (DEC-dec869)*(DEC-dec869) )
in869 = WHERE(rad869 LT radcen869)

rad884 = SQRT( (RA-ra884)*(RA-ra884) + (DEC-dec884)*(DEC-dec884) )
in884 = WHERE(rad884 LT radcen884)


;; filed star selection
;fieldR869 = !pi*radcen869*radcen869 ;;
field869rad = SQRT( (RA-33.5)*(RA-33.5) + (DEC-58)*(DEC-58) )
field869 = WHERE(field869rad LT radcen869)

field884rad = SQRT( (RA-33.5)*(RA-33.5) + (DEC-58)*(DEC-58) )
field884 = WHERE(field884rad LT radcen884)

;; circle out the inner and outer regionof Double Cluster
window,3,xsize=900,ysize=900
PLOT, RA, DEC, psym=3, xtit='RA', ytit='DEC',title='NGC869(h Persei), NGC884('+ cgGreek('chi') + ' Persei) and field star, GAIA DR3, N = 560569',yrange = [55,60]

TVCIRCLE, radcen869, ra869, dec869, color=1,THICK = 2
TVCIRCLE, radcen869, 33.5, 58, color=1,THICK = 2

TVCIRCLE, radcen884, ra884, dec884, color=3,THICK = 2
TVCIRCLE, radcen884, 33.5, 58, color=3,THICK = 2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; parallax of NGC869
window,4,xsize=900,ysize=900
plothist,Plx[in869],bin = 0.2,color = 1,xrange=[-3,5],title = 'Parallax of NGC869',xtitle='Plx' ;; central part of cluster (red)
plothist,Plx[field869],bin = 0.2,color = 3,/over ;; outer part of cluster (blue)
xyouts,3.,1500,'Red: cluster members',color = 1
xyouts,3.,1400,'Blue: filed stars',color = 3


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; parallax of NGC884
window,5,xsize=900,ysize=900
plothist,Plx[in884],bin = 0.2,color = 1,xrange=[-3,5],title = 'Parallax of NGC884',xtitle='Plx' ;; central part of cluster (red)
plothist,Plx[field884],bin = 0.2,color = 3,/over ;; outer part of cluster (blue)
xyouts,3.,1500,'Red: cluster members',color = 1
xyouts,3.,1400,'Blue: filed stars',color = 3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

px869 = where(Plx[in869] and Plx GT 0.37 and PLx LT 0.41)
py869 = where(Plx[in869] and Plx GT 0.37 and PLx LT 0.41)

;; proper motion of NGC869
window,6,xsize=900,ysize=900
plot,pmRA,pmDEC,psym = 3,xtitle = 'pmRA (mas/yr)',ytitle = 'pmDEC (mas/yr)',title = 'Proper Motion of NGC869',xrange=[-10,10],yrange=[-10,10],/nodata
oplot,pmRA,pmDEC,psym = 3,color = cgcolor('gray')
oplot,pmRA[in869],pmDEC[in869],psym = 3,color = cgcolor('black') ;; all star in cluster region
oplot,pmRA[px869],pmDEC[py869],psym = 1,color = cgcolor('red') ;; plx fit certain criterial 
xyouts,3.,8,'Black: all star in cluster region'
xyouts,3.,7,'Red: 0.37 < plx < 0.41',color = 1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

px884 = where(Plx[in884] and Plx GT 0.37 and PLx LT 0.41)
py884 = where(Plx[in884] and Plx GT 0.37 and PLx LT 0.41)

window,7,xsize=900,ysize=900
plot,pmRA,pmDEC,psym = 3,xtitle = 'pmRA (mas/yr)',ytitle = 'pmDEC (mas/yr)',title = 'Proper Motion of NGC884',xrange=[-10,10],yrange=[-10,10],/nodata
oplot,pmRA,pmDEC,psym = 3,color = cgcolor('gray')
oplot,pmRA[in884],pmDEC[in884],psym = 3,color = cgcolor('black')
oplot,pmRA[px884],pmDEC[py884],psym = 1,color = cgcolor('red')
xyouts,3.,8,'Black: all star in cluster region'
xyouts,3.,7,'Red: 0.37 < plx < 0.41',color = 1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

BPRP = BPmag - RPmag

;; color-color magnitude plot of ASCC21
window,8,xsize=900,ysize=900
plot,BPRP,Gmag,xtitle = 'G-RP [mag]',ytitle = 'G [mag]',title = 'color-color magnitude of NGC869',xrange=[-0.5,2.0],yrange=[max(BPRP),5.],ystyle=1,/nodata
oplot,BPRP,Gmag,psym = 3,color = cgcolor('gray')
oplot,BPRP[px869],Gmag[py869],psym = 1,color = cgcolor('red')
xyouts,1.3,7,'Black: all stars'
xyouts,1.3,8,'Red: 0.37 < plx < 0.41',color = 1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; color-color magnitude plot of ASCC21
window,9,xsize=900,ysize=900
plot,BPRP,Gmag,xtitle = 'G-RP [mag]',ytitle = 'G [mag]',title = 'color-color magnitude of NGC884',xrange=[-0.5,2.0],yrange=[max(BPRP),5.],ystyle=1,/nodata
oplot,BPRP,Gmag,psym = 3,color = cgcolor('gray')
oplot,BPRP[px884],Gmag[py884],psym = 1,color = cgcolor('red')
xyouts,1.3,7,'Black: all stars'
xyouts,1.3,8,'Red: 0.37 < plx < 0.41',color = 1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

rad869 = SQRT( (RA-ra869)*(RA-ra869) + (DEC-dec869)*(DEC-dec869) )
in869 = WHERE(rad869 LT radcen869)

rad884 = SQRT( (RA-ra884)*(RA-ra884) + (DEC-dec884)*(DEC-dec884) )
in884 = WHERE(rad884 LT radcen884)

;; circle out the core region and field star region
window,10,xsize=900,ysize=900
PLOT, RA, DEC, psym=3, xtit='RA', ytit='DEC',title='NGC869 and NGC884',yrange = [55,60]
TVCIRCLE,0.01,ra0,de0, color=2,THICK = 2

TVCIRCLE, radcen869, ra869, dec869, color=1,THICK = 2
TVCIRCLE, radcen884, ra884, dec884, color=3,THICK = 2
;TVCIRCLE, radcen, 33.5, 56.0, color=2,THICK = 2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;plot cluster star members with 0.37<Plx<0.41
out869 = WHERE(rad869 GT radcen869)

out884 = WHERE(rad884 GT radcen884)

NPlx869 = Plx
nanplx = where(NPlx869 LT 0.37 or NPLx869 GT 0.41)
NPlx869[nanplx] = !VALUES.F_NAN
NPlx869[out869] = !VALUES.F_NAN
index869 = where(Finite(Nplx869))

NPlx884 = Plx
nanplx = where(NPlx884 LT 0.37 or NPLx884 GT 0.41)
NPlx884[nanplx] = !VALUES.F_NAN
NPlx884[out884] = !VALUES.F_NAN
index884 = where(Finite(Nplx884))

window,11,xsize=900,ysize=900
PLOT, RA, DEC, psym=3, xtit='RA (J2000)', ytit='DEC (J2000)',title='NGC869(h Persei), NGC884('+ cgGreek('chi') + ' Persei)', $
  xrange=[34.5,35.75],yrange = [56.5,58],ystyle=1,/ISOTROPIC
OPLOT, RA[index869], DEC[index869],psym = 2,color = 1
OPLOT, RA[index884], DEC[index884],psym = 2,color = 3

TVCIRCLE, radcen869, ra869, dec869, color=1,THICK = 2
TVCIRCLE, radcen884, ra884, dec884, color=3,THICK = 2

xyouts,34.5,57.85,'Red: NGC869(h Persei)',color = 1
xyouts,34.5,57.80,'Blue: NGC884('+ cgGreek('chi') + ' Persei)',color = 3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; tangential velocity calculation
vtRA = 4.75*pmRA/Plx
vtDEC = 4.75*pmDEC/Plx

Q1 = where(vtRA GT 0 and vtDEC GT 0)
Q2 = where(vtRA LT 0 and vtDEC GT 0)
Q3 = where(vtRA LT 0 and vtDEC LT 0)
Q4 = where(vtRA GT 0 and vtDEC LT 0)

NPlx = Plx
nanplx = where(NPlx LT 0.37 or NPLx GT 0.41)
NPlx[nanplx] = !VALUES.F_NAN
index = where(Finite(NPlx)) ;; index of all srars with 0.37<PLx<0.41

window,12,xsize=900,ysize=900
PLOT, RA[index], DEC[index], psym=3, xtit='RA (J2000)', ytit='DEC (J2000)',title='NGC869(h Persei), NGC884('+ cgGreek('chi') + ' Persei), 0.37<Plx<0.41, N=22701', $
  xrange=[34.5,35.75],yrange = [56.5,58],ystyle=1,/ISOTROPIC
OPLOT, RA[index869], DEC[index869],psym = 2,color = 1
OPLOT, RA[index884], DEC[index884],psym = 2,color = 3

PRA = RA
PDEC = DEC
PRA[nanplx] = !VALUES.F_NAN
PDEC[nanplx] = !VALUES.F_NAN

A = FIndGen(16) * (!PI*2/16.)
UserSym, cos(A), sin(A), /fill

OPLOT, PRA[Q1], PDEC[Q1],psym = 8,color = cgcolor('Brown')
OPLOT, PRA[Q2], PDEC[Q2],psym = 8 ,color = cgcolor('green')
OPLOT, PRA[Q3], PDEC[Q3],psym = 8,color = cgcolor('pink')
OPLOT, PRA[Q4], PDEC[Q4],psym = 8,color = cgcolor('purple')

TVCIRCLE, radcen869, ra869, dec869, color=1,THICK = 2
TVCIRCLE, radcen884, ra884, dec884, color=3,THICK = 2

xyouts,34.5,57.85,'Red: NGC869(h Persei)',color = 1
xyouts,34.5,57.80,'Blue: NGC884('+ cgGreek('chi') + ' Persei)',color = 3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; plot histofram of cluster star members with 0.37<Plx<0.41
window,13,xsize=900,ysize=900
!P.Multi = [0, 2, 2]
PLOT, RA[index], DEC[index], psym=3, xtit='RA (J2000)', ytit='DEC (J2000)', $
  xrange=[34.4,35.8],yrange = [56.5,58],position = [0.07,0.07,0.70,0.70],xstyle=1,ystyle=1
OPLOT, RA[index869], DEC[index869],psym = 2,color = 1
OPLOT, RA[index884], DEC[index884],psym = 2,color = 3

TVCIRCLE, radcen869, ra869, dec884, color=1,THICK = 2
TVCIRCLE, radcen884, ra884, dec884, color=3,THICK = 2

plothist,RA[index],bin = 0.1,position = [0.07,0.70,0.70,0.95], $
  title='NGC869(h Persei), NGC884('+ cgGreek('chi') + ' Persei), 0.37<Plx<0.41, N=22701',xtitle=' ',yrange=[400,550],xrange=[34.4,35.8],XTICKFORMAT="(A1)"
plothist,DEC[index],bin = 0.1,position = [0.70,0.07,0.95,0.70],ytitle = ' ',/rotate,xrange=[600,1300],yrange = [56.5,58],YTICKFORMAT="(A1)" ;,XTICKFORMAT="(A1)"






































end