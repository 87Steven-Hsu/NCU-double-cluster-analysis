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
test = read_csv('D:\學校資料\天文所陳老師\h and x\data\H_and_X_3deg.csv')

;Transfer CSV transfer to sav.
SAVE, test, FILE='D:\學校資料\天文所陳老師\h and x\data\H_and_X_3deg.dat'
restore,'D:\學校資料\天文所陳老師\h and x\data\H_and_X_3deg.dat'

;Transfer structure to array
RA = [float(test.FIELD01)]    ;RA (deg)
DEC = [float(test.FIELD02)]   ;DEC (deg)
RA_ICRS = [float(test.FIELD03)]   ;
eRA_ICRS = [float(test.FIELD04)]  ;
DEC_ICRS = [float(test.FIELD05)]   ;
eDEC_ICRS = [float(test.FIELD06)]  ;
Plx = [float(test.FIELD07)]  ;PLX, Absolute steeGmagar paraeGmagax (paraeGmagax)
ePlx = [float(test.FIELD08)] ;PLX error
pmRA = [float(test.FIELD09)] ;pmRA (mas/yr), Proper motion in right ascension
epmRA = [float(test.FIELD10)];pmRA error (mas/yr)
pmDEC = [float(test.FIELD11)]  ;pmDEC (mas/yr), Proper motion in declination
epmDEC = [float(test.FIELD12)] ;pmDEC error (mas/yr)
Gmag = [float(test.FIELD13)] ;Gmag (mag), Gband mean magnitude = brightness distribution
eGmag = [float(test.FIELD14)];Gmag error (mag)
BPmag = [float(test.FIELD15)] ;BPmag (mag), Integrated BP mean magnitude
eBPmag = [float(test.FIELD16)];BPmag error
RPmag = [float(test.FIELD17)]    ;RPmag (mag), Integrated RP mean magnitude
eRPmag = [float(test.FIELD18)]   ;RPmag error

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ra0 = MEAN(RA)
de0 = MEAN(DEC)

;;;Coordiante of equatorial system ASCC16
window,1,xsize=900,ysize=900
plot,RA,DEC,psym=3,xtitle='RA (deg)',ytitle='DEC (deg)',title='Coordiante of Stars Observed by GAIA',yrange = [55,60]
TVCIRCLE,0.01,ra0,de0, color=2,THICK = 3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
radcen = 30./60.
ra869 = 34.75
dec869 = 57.13

ra884 = 35.60
dec884 = 57.13

;; circle out a region which is in 30' of ASCC16 and ASCC21
window,2,xsize=900,ysize=900
PLOT, RA, DEC, psym=3, xtit='RA', ytit='DEC',title='NGC869 and NGC884',yrange = [55,60]
TVCIRCLE,0.01,ra0,de0, color=2,THICK = 2
TVCIRCLE, radcen, ra869, dec869, color=1,THICK = 2
TVCIRCLE, radcen, ra884, dec884, color=3,THICK = 2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

rad869 = SQRT( (RA-ra869)*(RA-ra869) + (DEC-dec869)*(DEC-dec869) )
in869 = WHERE(rad869 LT radcen)

rad884 = SQRT( (RA-ra884)*(RA-ra884) + (DEC-dec884)*(DEC-dec884) )
in884 = WHERE(rad884 LT radcen)

radout1 = 0.74 ;;inner circle = 44.90 arcmin
radout2 = 0.9 ;; outer circle = 54 arcmin
out869 = where(rad869 gt radout1 and rad869 lt radout2)
out884 = where(rad884 gt radout1 and rad884 lt radout2)

;; circle out the inner and outer region which is in 30' of ASCC16 and ASCC21
window,3,xsize=900,ysize=900
PLOT, RA, DEC, psym=3, xtit='RA', ytit='DEC',title='ASCC16 and ASCC21',yrange = [55,60]
TVCIRCLE,0.01,ra0,de0, color=2,THICK = 2

TVCIRCLE, radcen, ra869, dec869, color=1,THICK = 2
TVCIRCLE, radout1, ra869, dec869, color=3,THICK = 2
TVCIRCLE, radout2, ra869, dec869, color=3,THICK = 2  ;; green

TVCIRCLE, radcen, ra884, dec884, color=3,THICK = 2
TVCIRCLE, radout1, ra884, dec884, color=1,THICK = 2
TVCIRCLE, radout2, ra884, dec884, color=1,THICK = 2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; parallax of NGC869
window,4,xsize=900,ysize=900
plothist,Plx[in869],bin = 0.2,color = 1,xrange=[-3,5],title = 'Parallax of NGC869',xtitle='Plx' ;; central part of cluster (red)
plothist,Plx[out869],bin = 0.2,color = 3,/over ;; outer part of cluster (blue)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; parallax of NGC884
window,5,xsize=900,ysize=900
plothist,Plx[in884],bin = 0.2,color = 1,xrange=[-3,5],title = 'Parallax of NGC884',xtitle='Plx' ;; central part of cluster (red)
plothist,Plx[out884],bin = 0.2,color = 3,/over ;; outer part of cluster (blue)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

px869 = where(Plx[in869] and Plx GT 0.38 and PLx LT 0.42)
py869 = where(Plx[in869] and Plx GT 0.38 and PLx LT 0.42)

;; proper motion of NGC869
window,6,xsize=900,ysize=900
plot,pmRA,pmDEC,psym = 3,xtitle = 'pmRA (mas/yr)',ytitle = 'pmDEC (mas/yr)',title = 'Proper Motion of NGC869',xrange=[-10,10],yrange=[-10,10],/nodata
oplot,pmRA,pmDEC,psym = 3,color = cgcolor('gray')
oplot,pmRA[in869],pmDEC[in869],psym = 3,color = cgcolor('black')
oplot,pmRA[px869],pmDEC[py869],psym = 1,color = cgcolor('red')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

px884 = where(Plx[in884] and Plx GT 0.38 and PLx LT 0.42)
py884 = where(Plx[in884] and Plx GT 0.38 and PLx LT 0.42)

window,7,xsize=900,ysize=900
plot,pmRA,pmDEC,psym = 3,xtitle = 'pmRA (mas/yr)',ytitle = 'pmDEC (mas/yr)',title = 'Proper Motion of NGC884',xrange=[-10,10],yrange=[-10,10],/nodata
oplot,pmRA,pmDEC,psym = 3,color = cgcolor('gray')
oplot,pmRA[in884],pmDEC[in884],psym = 3,color = cgcolor('black')
oplot,pmRA[px884],pmDEC[py884],psym = 1,color = cgcolor('red')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

BPRP = BPmag - RPmag

;; color-color magnitude plot of ASCC21
window,8,xsize=900,ysize=900
plot,GRP,Gmag,xtitle = 'G-RP [mag]',ytitle = 'G [mag]',title = 'color-color magnitude of NGC869',xrange=[-0.5,2.0],yrange=[max(GRP),5.],ystyle=1,/nodata
oplot,GRP,Gmag,psym = 3,color = cgcolor('gray')
oplot,GRP[px869],Gmag[py869],psym = 1,color = cgcolor('red')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; color-color magnitude plot of ASCC21
window,9,xsize=900,ysize=900
plot,BPRP,Gmag,xtitle = 'G-RP [mag]',ytitle = 'G [mag]',title = 'color-color magnitude of NGC884',xrange=[-0.5,2.0],yrange=[max(GRP),5.],ystyle=1,/nodata
oplot,BPRP,Gmag,psym = 3,color = cgcolor('gray')
oplot,BPRP[px884],Gmag[py884],psym = 1,color = cgcolor('red')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

rad869 = SQRT( (RA-ra869)*(RA-ra869) + (DEC-dec869)*(DEC-dec869) )
in869 = WHERE(rad869 LT radcen)

rad884 = SQRT( (RA-ra884)*(RA-ra884) + (DEC-dec884)*(DEC-dec884) )
in884 = WHERE(rad884 LT radcen)


;; circle out the inner and outer region which is in 30' of ASCC16 and ASCC21
window,10,xsize=900,ysize=900
PLOT, RA, DEC, psym=3, xtit='RA', ytit='DEC',title='NGC869 and NGC884',yrange = [55,60]
TVCIRCLE,0.01,ra0,de0, color=2,THICK = 2

TVCIRCLE, radcen, ra869, dec869, color=1,THICK = 2
TVCIRCLE, radcen, ra884, dec884, color=3,THICK = 2
TVCIRCLE, radcen, 33.5, 56.0, color=2,THICK = 2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

radref = SQRT( (RA-33.5)*(RA-33.5) + (DEC-56.0)*(DEC-56.0) )
inref = WHERE(radref LT radcen)

;; parallax of NGC869
window,11,xsize=900,ysize=900
plothist,Plx[in869],bin = 0.2,color = 1,xrange=[-3,5],title = 'Parallax of NGC869',xtitle='Plx' ;; central part of cluster (red)
plothist,Plx[inref],bin = 0.2,color = 3,/over ;; reference field (blue)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; parallax of NGC884
window,12,xsize=900,ysize=900
plothist,Plx[in884],bin = 0.2,color = 1,xrange=[-3,5],title = 'Parallax of NGC884',xtitle='Plx' ;; central part of cluster (red)
plothist,Plx[inref],bin = 0.2,color = 3,/over ;; reference field (blue)











































end