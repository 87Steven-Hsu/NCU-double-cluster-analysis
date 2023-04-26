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
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;讀取資料
;file = read_csv('D:\學校資料\天文所陳老師\Double cluster\data\double cluster xmatch - mod.csv')
;
;;Transfer CSV transfer to sav.
;SAVE, file, FILE='D:\學校資料\天文所陳老師\Double cluster\data\double cluster xmatch - mod.dat'
;restore,'D:\學校資料\天文所陳老師\Double cluster\data\double cluster xmatch - mod.dat'
;
;;Transfer structure to array
;;;2 MASS
;angDist = [float(file.FIELD01)]    
;sigDist = [float(file.FIELD02)]   
;RAJ2000 = [float(file.FIELD03)]   ;
;DECJ2000 = [float(file.FIELD04)]  ;
;errHalfMaj = [float(file.FIELD05)]   ;
;errHalfMin = [float(file.FIELD06)]  ;
;errPosAng = [float(file.FIELD07)]  ;
;Jmag = [float(file.FIELD08)] ;
;Hmag = [float(file.FIELD09)] ;
;Kmag = [float(file.FIELD10)];
;eJmag = [float(file.FIELD11)]  ;
;eHmag = [float(file.FIELD12)] ;
;eKmag = [float(file.FIELD13)] ;
;;;GAIA
;RA_epoch2000 = [float(file.FIELD14)];
;DEC_epoch2000 = [float(file.FIELD15)] ;
;errHalfMaj_1 = [float(file.FIELD16)] ;
;errHalfMin_1 = [float(file.FIELD17)]    ;
;errPosAng_1 = [float(file.FIELD18)]   ;
;RA_c = [float(file.FIELD19)]
;eRA_c = [float(file.FIELD20)]
;DEC_c = [float(file.FIELD21)]   ;
;eDEC_c = [float(file.FIELD22)]  ;
;Plx = [float(file.FIELD23)]   ;
;ePlx = [float(file.FIELD24)]  ;
;pmRA = [float(file.FIELD25)]  ;
;epmRA = [float(file.FIELD26)]
;pmDEC = [float(file.FIELD27)]
;epmDEC = [float(file.FIELD28)]   ;
;Gmag = [float(file.FIELD29)]  ;
;BPmag = [float(file.FIELD30)]   ;
;RPmag = [float(file.FIELD31)]  ;

ra0 = MEAN(RA_epoch2000)
de0 = MEAN(DEC_epoch2000)

;; coordinate and radius of ASCC16 and ASCC21
cen16 = 39.2/60.
ra16 = 81
de16 = 1.78
inc16 = 44.9/60 ;;inner circle = 44.90 arcmin
outc16 = 59.6/60 ;; outer circle = 54 arcmin

cen21 = 48.0/60
ra21 = 82
de21 = 3.65
inc21 = 60/60 ;;inner circle = 44.90 arcmin
outc21 = 76.84/60 ;; outer circle = 54 arcmin

;;2 MASS
RA5deg = WHERE(RAJ2000 LT 84 and RAJ2000 GT 79)
MRA = RAJ2000[RA5deg]
DEC5deg = WHERE(DECJ2000 LT 5 and DECJ2000 GT 0)
MDEC = DECJ2000[DEC5deg]

;; Magnitude filter
jmaglt16 = where(Jmag LT 16)
hpmaglt15 = where(Hmag lt 15)
kpmaglt15 = where(Kmag lt 15)

;;;Coordiante of equatorial system ASCC16 and ASCC21
;window,1,xsize=900,ysize=900
;plot,RAJ2000,DECJ2000,psym=3,xtitle='RA (deg)',ytitle='DEC (deg)',title='Coordiante of Stars Observed by 2 MASS'
;TVCIRCLE,0.01,ra0,de0, color=2,THICK = 3
;TVCIRCLE,radcen,ra16,de16, color=1,THICK = 3
;TVCIRCLE,radcen,ra21,de21, color=3,THICK = 3

;;;Coordiante of equatorial system ASCC16 and ASCC21
window,2,xsize=900,ysize=900
!P.Multi = [0, 2, 2]
plot,MRA[jmaglt16],MDEC[jmaglt16],psym=3,position = [0.07,0.07,0.70,0.70],xtitle='RA (deg)',ytitle='DEC (deg)',/isotropic
TVCIRCLE,cen16,ra16,de16, color=1,THICK = 3
TVCIRCLE,cen21,ra21,de21, color=3,THICK = 3
plothist,MRA[jmaglt16],bin = 0.1,yrange = [2700,4200],position = [0.07,0.70,0.70,0.95], title = '2MASS Jmag<16, N=291178'
plothist,MDEC[jmaglt16],bin = 0.1,xrange=[2800,3800],position = [0.70,0.07,0.95,0.70],/ROTATE

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; GAIA
RA5deg = WHERE(RA_epoch2000 LT 84 and RA_epoch2000 GT 79)
gaiaRA = RA_epoch2000[RA5deg]
DEC5deg = WHERE(DEC_epoch2000 LT 5 and DEC_epoch2000 GT 0)
gaiaDEC = DEC_epoch2000[DEC5deg]

;; Magnitude filter
gmaglt16 = where(Gmag LT 16)
bpmaglt16 = where(BPmag lt 16)
rpmaglt16 = where(rpmag lt 16)

;; parallax filter
plx16 = Plx[gmaglt16]

;;;boundary of cluster, inner circle and outer circle 
rad16 = SQRT( (gaiaRA - ra16)*(gaiaRA - ra16) + (gaiaDEC - de16)*(gaiaDEC - de16) )
in16 = WHERE(Plx[gmaglt16] and rad16 LT cen16)
out16 = where(Plx[gmaglt16] and rad16 gt inc16 and rad16 lt outc16)

rad21 = SQRT( (gaiaRA - ra21)*(gaiaRA - ra21) + (gaiaDEC - de21)*(gaiaDEC - de21) )
in21 = WHERE(Plx[gmaglt16] and rad21 LT cen21)
out21 = where(Plx[gmaglt16] and rad21 gt inc21 and rad21 lt outc21)

;; proper motion filter
pmRA16 = pmRA[gmaglt16]
pmDEC16 = pmDEC[gmaglt16]
px16 = where(Plx[gmaglt16] and Plx[in16] and Plx GT 2.5 and PLx LT 3.0)
py16 = where(Plx[gmaglt16] and Plx[in16] and Plx GT 2.5 and Plx Lt 3.0)
px21 = where(Plx[gmaglt16] and Plx[in21] and Plx GT 2.5 and PLx LT 3.0)
py21 = where(Plx[gmaglt16] and Plx[in21] and Plx GT 2.5 and Plx Lt 3.0)

;;;
window,3,xsize=900,ysize=900
!P.Multi = [0, 2, 2]
plot,gaiaRA[gmaglt16],gaiaDEC[gmaglt16],psym=3,position = [0.07,0.07,0.70,0.70],xtitle='RA (deg)',ytitle='DEC (deg)',/isotropic
TVCIRCLE,cen16,ra16,de16, color=1,THICK = 3
TVCIRCLE,cen21,ra21,de21, color=3,THICK = 3
plothist,gaiaRA[gmaglt16],bin = 0.1,yrange = [900,1600],position = [0.07,0.70,0.70,0.95], title = 'GAIA Gmag<16, N=104536'
plothist,gaiaDEC[gmaglt16],bin = 0.1,xrange=[900,1500],position = [0.70,0.07,0.95,0.70],/ROTATE

;;;ASCC16 Plx
window,4,xsize=900,ysize=900
!P.Multi = [0, 1, 1]
plothist,Plx[in16],bin = 0.2,color = 1,title = 'Parallax of ASCC16, Gmag<16',xrange=[-1,5],xtitle='Plx' ;; central part of cluster (red)
plothist,Plx[out16],bin = 0.2,color = 3,/over ;; outer part of cluster (blue)

;;;ASCC21 Plx
window,5,xsize=900,ysize=900
!P.Multi = [0, 1, 1]
plothist,Plx[in21],bin = 0.2,color = 1,title = 'Parallax of ASCC21, Gmag<16',xrange=[-1,5],yrange=[0,1100],xtitle='Plx' ;; central part of cluster (red)
plothist,Plx[out21],bin = 0.2,color = 3,/over ;; outer part of cluster (blue)

;;;ASCC16 proper motion
window,6,xsize=900,ysize=900
plot,pmRA[gmaglt16],pmDEC[gmaglt16],psym = 3,xtitle = 'pmRA (mas/yr)',ytitle = 'pmDEC (mas/yr)',title = 'Proper Motion of ASCC16, Gmag<16',xrange=[-10,10],yrange=[-10,10],/nodata
oplot,pmRA[gmaglt16],pmDEC[gmaglt16],psym = 3,color = cgcolor('gray')
oplot,pmRA[in16],pmDEC[in16],psym = 3,color = cgcolor('black')
oplot,pmRA[px16],pmDEC[py16],psym = 1,color = cgcolor('red')

;;;ASCC21 proper motion
window,7,xsize=900,ysize=900
!P.Multi = [0, 1, 1]
plot,pmRA[gmaglt16],pmDEC[gmaglt16],psym = 3,xtitle = 'pmRA (mas/yr)',ytitle = 'pmDEC (mas/yr)',title = 'Proper Motion of ASCC21, Gmag<16',xrange=[-10,10],yrange=[-10,10],/nodata
oplot,pmRA[gmaglt16],pmDEC[gmaglt16],psym = 3,color = cgcolor('gray')
oplot,pmRA[in21],pmDEC[in21],psym = 3,color = cgcolor('black')
oplot,pmRA[px21],pmDEC[py21],psym = 1,color = cgcolor('red')

;;magnitude filter
G16 = Gmag[gmaglt16]
GRP = Gmag - RPmag

in16 = WHERE(rad16 LT cen16) ;Plx[gmaglt16] and 
px16 = where(Gmag[gmaglt16] and Plx[in16] and Plx GT 2.5 and PLx LT 3.0)
py16 = where(Gmag[gmaglt16] and Plx[in16] and Plx GT 2.5 and Plx Lt 3.0)

in21 = WHERE(rad21 LT cen21)
px21 = where(Gmag[gmaglt16] and Plx[in21] and Plx GT 2.5 and PLx LT 3.0)
py21 = where(Gmag[gmaglt16] and Plx[in21] and Plx GT 2.5 and Plx Lt 3.0)


;; color-color magnitude plot of ASCC16
window,8,xsize=900,ysize=900
plot,GRP,Gmag,xtitle = 'G-RP [mag]',ytitle = 'G [mag]',title = 'color-color magnitude of ASCC16, Gmag<16',xrange=[-0.5,2.0],yrange=[max(Gmag),5.],ystyle=1,/nodata
oplot,GRP,Gmag,psym = 3,color = cgcolor('gray')
oplot,GRP[px16],Gmag[py16],psym = 1,color = cgcolor('red')

;;; color-color magnitude plot of ASCC21
window,9,xsize=900,ysize=900
plot,GRP,Gmag,xtitle = 'G-RP [mag]',ytitle = 'G [mag]',title = 'color-color magnitude of ASCC21, Gmag<16',xrange=[-0.5,2.0],yrange=[max(GRP),5.],ystyle=1,/nodata
oplot,GRP,Gmag,psym = 3,color = cgcolor('gray')
oplot,GRP[px21],Gmag[py21],psym = 1,color = cgcolor('red')

;device,decompose = 0
;JKmag = Jmag - Kmag
;;;; Colorcolor magnitude of 2 MASS
;window,5,xsize=900,ysize=900
;plot,JKmag,Kmag,psym=3,xtitle = 'J - K [mag]',ytitle = 'K [mag]' 
;;oplot,JKmag[plx16],Kmag[plx16],psym=3,color=cgcolor('red')
;oplot,JKmag[plx21],Kmag[plx21],psym=3,color=cgcolor('Cyan')
































end