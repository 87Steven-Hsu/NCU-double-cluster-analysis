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

;讀取資料 ASCC 16
test = read_csv('D:\學校資料\天文所陳老師\Double cluster\data\ASCC16.csv') 

;Transfer CSV transfer to sav.
SAVE, test, FILE='D:\學校資料\天文所陳老師\Double cluster\data\ASCC16.dat'
restore,'D:\學校資料\天文所陳老師\Double cluster\data\ASCC16.dat'

;Transfer structure to array
RA16 = [float(test.FIELD01)]    ;RA (deg)
DEC16 = [float(test.FIELD02)]   ;DEC (deg)
RA_ICRS16 = [float(test.FIELD03)]   ;
eRA_ICRS16 = [float(test.FIELD04)]  ;
DEC_ICRS16 = [float(test.FIELD05)]   ;
eDEC_ICRS16 = [float(test.FIELD06)]  ;
Plx16 = [float(test.FIELD07)]  ;PLX, Absolute steeGmagar paraeGmagax (paraeGmagax)    
ePlx16 = [float(test.FIELD08)] ;PLX error
pmRA16 = [float(test.FIELD09)] ;pmRA (mas/yr), Proper motion in right ascension
epmRA16 = [float(test.FIELD10)];pmRA error (mas/yr)
pmDEC16 = [float(test.FIELD11)]  ;pmDEC (mas/yr), Proper motion in declination
epmDEC16 = [float(test.FIELD12)] ;pmDEC error (mas/yr)
Gmag16 = [float(test.FIELD13)] ;Gmag (mag), Gband mean magnitude = brightness distribution
eGmag16 = [float(test.FIELD14)];Gmag error (mag)
BPmag16 = [float(test.FIELD15)] ;BPmag (mag), Integrated BP mean magnitude
eBPmag16 = [float(test.FIELD16)];BPmag error
RPmag16 = [float(test.FIELD17)]    ;RPmag (mag), Integrated RP mean magnitude
eRPmag16 = [float(test.FIELD18)]   ;RPmag error

;讀取資料 ASCC 21
test = read_csv('D:\學校資料\天文所陳老師\Double cluster\data\ASCC21.csv')

;Transfer CSV transfer to sav.
SAVE, test, FILE='D:\學校資料\天文所陳老師\Double cluster\data\ASCC21.dat'
restore,'D:\學校資料\天文所陳老師\Double cluster\data\ASCC21.dat'

;Transfer structure to array
RA21 = [float(test.FIELD01)]    ;RA (deg)
DEC21 = [float(test.FIELD02)]   ;DEC (deg)
RA_ICRS21 = [float(test.FIELD03)]   ;
eRA_ICRS21 = [float(test.FIELD04)]  ;
DEC_ICRS21 = [float(test.FIELD05)]   ;
eDEC_ICRS21 = [float(test.FIELD06)]  ;
Plx21 = [float(test.FIELD07)]  ;PLX, Absolute steeGmagar paraeGmagax (paraeGmagax)
ePlx21 = [float(test.FIELD08)] ;PLX error
pmRA21 = [float(test.FIELD09)] ;pmRA (mas/yr), Proper motion in right ascension
epmRA21 = [float(test.FIELD10)];pmRA error (mas/yr)
pmDEC21 = [float(test.FIELD11)]  ;pmDEC (mas/yr), Proper motion in declination
epmDEC21 = [float(test.FIELD12)] ;pmDEC error (mas/yr)
Gmag21 = [float(test.FIELD13)] ;Gmag (mag), Gband mean magnitude = brightness distribution
eGmag21 = [float(test.FIELD14)];Gmag error (mag)
BPmag21 = [float(test.FIELD15)] ;BPmag (mag), Integrated BP mean magnitude
eBPmag21 = [float(test.FIELD16)];BPmag error
RPmag21 = [float(test.FIELD17)]    ;RPmag (mag), Integrated RP mean magnitude
eRPmag21 = [float(test.FIELD18)]   ;RPmag error

;;;Coordiante of equatorial system ASCC16
window,1,xsize=900,ysize=900
plot,RA16,DEC16,yrange=[0.5,3.0],psym=3,xtitle='RA (deg)',ytitle='DEC (deg)',title='Coordiante of Stars Observed by GAIA',/ISOTROPIC   ;,yrange=[1.7,1.9]


;;;Coordiante of equatorial system ASCC21
window,2,xsize=900,ysize=900
plot,RA21,DEC21,yrange=[2.5,5.0],psym=3,xtitle='RA (deg)',ytitle='DEC (deg)',title='Coordiante of Stars Observed by GAIA',/ISOTROPIC ;,yrange=[0.5,3.0]

;coordinate of galactic coordinate
;window,2 
;b = findgen(n_elements(ra))   ;gra
;L1 = findgen(n_elements(dec)) ;
;L2 = findgen(n_elements(dec)) ;gdec
;for i = 0,n_elements(ra)-1 do begin
;  b[i] = 180/(!pi*asin(cos(dec[i])*cos(192.85948)*cos(ra[i]-192.85948)+sin(dec[i]*cos(27.12825))))   ;degree
;  L1[i] = -(asin((cos(dec[i])*sin(ra[i]-192.85948))/cos(b[i]))-122.93192) ;radius
;  L2[i] = 180/(!pi*L1[i])                                                 ;degree
;endfor
;plot,b,l2,yrange=[0.45,0.48],psym=3

;;;Coordiante of equatorial system
;window,2,xsize=800,ysize=800
;plot,RA,DEC,psym=3,xrange=[81.,81.2],yrange=[1.0,3.0],xtitle='RA (deg)',ytitle='DEC (deg)',title='Coordiante of Stars Observed by GAIA' ;,/ISOTROPIC  ;

;;Proper Motion
;window,3
;plot,pmRA,pmDEC,psym=3,ystyle=1,xstyle=1,xtitle='pmRA (deg)',ytitle='pmDEC (deg)',title='Coordinate of Proper Motion' ;,xrange=[-60,100],yrange=[-150,100]

;;Gmag
;window,5   
;x1 = findgen(22) & x2 = findgen(22) & ans = findgen(22)
;for i = 0,ceil(max(gmag))-1 do begin ;
;  x1[i] = n_elements(where(gmag lt i))
;  x2[i] = n_elements(where(gmag lt i+1))
;  ans[i] = x2[i] - x1[i]
;endfor
;red = [255,0,0] & green = [0,255,0] & blue = [0,0,255] & white = [0,0,0] & black = [255,255,255]
;tick = ['6', '7', '8', '9', '10', '11', '12' ,'13', '14', '15', '16', '17', '18', '19', '20', '21']
;box = barplot(ans,fill_color='black',xstyle=1,xsubticklen=0.5,xrange=[6,22],xtitle='Gmag',ytitle='count') ;,xtickname=["6", "7", '8', '9', '10', '11', '12' ,'13', '14', '15', '16', '17', '18', '19', '20', '21']) ;,colors=1 ;,psym=7,xrange=[5.5,21.5],xstyle=1 
;w = window(DIMENSIONS=[800,700])
;hisGmag = HISTOGRAM(Gmag, BINSIZE=binsize, LOCATIONS=binvals)
;Gmag = PLOT(binvals,hisGmag,/CURRENT,/STAIRSTEP)

;;Color-color Diagram
;window,6
;xx = BPmag-RPmag
;yy = Gmag-RPmag
;plot,xx,yy,xrange=[-2,8],yrange=[-0.5,2.5],psym=3,xtitle='BP-RP (mag)',ytitle='G-RP (mag)',title='Color-color Diagram'

;;;Color-color Diagram
;window,7
;xx = BPmag-RPmag
;yy = Gmag-BPmag
;plot,xx,yy,xrange=[-2,8],yrange=[-8,2],psym=3,xtitle='BP-RP (mag)',ytitle='C-RP (mag)',title='Color-color Diagram'


;;;;;; cluster identification of ASCC16
;; plot RA and DEC
;w = window(DIMENSIONS=[1000,900])
;coordinate = PLOT(RA16,Dec16,'.',/current,position=[0.07,0.07,0.70,0.70],xtitle='RA',ytitle='DEC')
;
;;;RA histogram
;hisRA = HISTOGRAM(RA16, BINSIZE=0.01, LOCATIONS=binvals)
;;RA = barplot(binvals,hisRA,/CURRENT,/HISTOGRAM,position=[0.07,0.70,0.70,0.95],xrange=[80,82.5],xtickname=[' '])
;RA = plot(binvals,hisRA,/CURRENT,/STAIRSTEP,position=[0.07,0.70,0.70,0.95],xrange=[80,82.5],xtickname=[' '])
;
;;; DEC histogram
;hisDec = HISTOGRAM(Dec16,BINSIZE=0.01, LOCATIONS=binvals)
;Dec = PLOT(hisDec,binvals,/CURRENT,/STAIRSTEP,yrange=[0.5,3],position=[0.70,0.07,0.95,0.70],ytickname=[' '])


;;;;;; cluster identification of ASCC21
;;; plot RA and DEC
;w = window(DIMENSIONS=[1000,900])
;coordinate = PLOT(RA21,Dec21,'.',/current,position=[0.07,0.07,0.70,0.70],xtitle='RA',ytitle='DEC')
;
;;;RA histogram
;hisRA = HISTOGRAM(RA21, BINSIZE=0.1, LOCATIONS=binvals)
;;RA = barplot(binvals,hisRA,/CURRENT,/HISTOGRAM,position=[0.07,0.70,0.70,0.95],xrange=[80,82.5],xtickname=[' '])
;RA = plot(binvals,hisRA,/CURRENT,/STAIRSTEP,position=[0.07,0.70,0.70,0.95],xrange=[80.5,83.5],xtickname=[' '])
;
;;; DEC histogram
;hisDec = HISTOGRAM(Dec21,BINSIZE=0.1, LOCATIONS=binvals)
;Dec = PLOT(hisDec,binvals,/CURRENT,/STAIRSTEP,yrange=[2.5,5],position=[0.70,0.07,0.95,0.70],ytickname=[' '])


;; circle out a region which is in 10' of ASCC16
ra0=MEAN(RA16)
de0=MEAN(DEC16)

rad = SQRT( (RA16-ra0)*(RA16-ra0) + (DEC16-de0)*(DEC16-de0) )

radcen = 30./60. ; 30 arcmin
in16 = WHERE(rad LT radcen)

window,3,xsize=900,ysize=900
PLOT, RA16, DEC16, psym=3, /ynozero, xtit='RA', ytit='DEC',title='ASCC16 All Star in 2 deg'
TVCIRCLE, radcen, ra0, de0, color=1,THICK = 2

radout1 = 0.74 ;;inner circle = 44.90 arcmin
radout2 = 0.9 ;; outer circle = 54 arcmin
out16 = where(rad gt radout1 and rad lt radout2)

TVCIRCLE, radout1, ra0, de0, color=3,THICK = 2 ;; blue 
TVCIRCLE, radout2, ra0, de0, color=2,THICK = 2  ;; green


;; parallax of ASCC16
window,4,xsize=900,ysize=900
plothist,Plx16[in16],bin = 0.2,color = 1,title = 'Parallax of ASCC16',xrange=[-1,5],xtitle='Plx' ;; central part of cluster (red)
plothist,Plx16[out16],bin = 0.2,color = 3,/over ;; outer part of cluster (blue)


;; proper motion of ASCC16
px16 = where(Plx16[in16] and Plx16 GT 2.5 and PLx16 LT 3.0)
py16 = where(Plx16[in16] and Plx16 GT 2.5 and Plx16 Lt 3.0)

window,5,xsize=900,ysize=900
plot,pmRA16,pmDEC16,psym = 3,xtitle = 'pmRA (mas/yr)',ytitle = 'pmDEC (mas/yr)',title = 'Proper Motion of ASCC16',xrange=[-10,10],yrange=[-10,10],/nodata
oplot,pmRA16,pmDEC16,psym = 3,color = cgcolor('gray')
oplot,pmRA16[in16],pmDEC16[in16],psym = 3,color = cgcolor('black')
oplot,pmRA16[px16],pmDEC16[py16],psym = 1,color = cgcolor('red')

;; color-color magnitude plot of ASCC16
GRP16 = Gmag16 - RPmag16
window,6,xsize=900,ysize=900
plot,GRP16,Gmag16,xtitle = 'G-RP [mag]',ytitle = 'G [mag]',title = 'color-color magnitude of ASCC16',xrange=[-0.5,2.0],yrange=[max(GRP16),5.],ystyle=1,/nodata
oplot,GRP16,Gmag16,psym = 3,color = cgcolor('gray')
oplot,GRP16[px16],Gmag16[py16],psym = 1,color = cgcolor('red')

;; circle out a region which is in 10' of ASCC21
ra0=MEAN(RA21)
de0=MEAN(DEC21)

rad = SQRT( (RA21-ra0)*(RA21-ra0) + (DEC21-de0)*(DEC21-de0) )

radcen = 30./60. ; = 30 arcmin
in21 = WHERE(rad LT radcen)

window,7,xsize=900,ysize=900
PLOT, RA21, DEC21, psym=3,/ynozero,xtit='RA',ytit='DEC',title='ASCC21 All Star in 2 deg'
TVCIRCLE, radcen, ra0, de0, color=1,THICK = 2

radout1 = 0.74 ;; inner circle = 44.90 arcmin
radout2 = 0.9 ;; outer circle = 54 arcmin
out21 = where(rad gt radout1 and rad lt radout2)

TVCIRCLE, radout1, ra0, de0, color=3,THICK = 2 ;; blue 
TVCIRCLE, radout2, ra0, de0, color=2,THICK = 2 ;; green


;; parallax of ASCC21
window,8,xsize=900,ysize=900
plothist,Plx21[in21],bin = 0.2,color = 1,title = 'Parallax of ASCC21',xrange=[-1,5],xtitle='Plx';; central part of cluster (red)
plothist,Plx21[out21],bin = 0.2,color = 3,/over ;; outer part of cluster (blue)


;; proper motion scatter plot of ASCC21
px21 = where(Plx21[in21] and Plx21 GT 2.5 and PLx21 LT 3.0)
py21 = where(Plx21[in21] and Plx21 GT 2.5 and Plx21 Lt 3.0)

window,9,xsize=900,ysize=900
plot,pmRA21,pmDEC21,psym = 3,xtitle = 'pmRA (mas/yr)',ytitle = 'pmDEC (mas/yr)',title = 'Proper Motion of ASCC21',xrange=[-10,10],yrange=[-10,10],/nodata
oplot,pmRA21,pmDEC21,psym = 3,color = cgcolor('gray')
oplot,pmRA21[in21],pmDEC21[in21],psym = 3,color = cgcolor('black')
oplot,pmRA21[px21],pmDEC21[py21],psym = 1,color = cgcolor('red')


;; color-color magnitude plot of ASCC21
GRP21 = Gmag21 - RPmag21
window,10,xsize=900,ysize=900
plot,GRP21,Gmag21,xtitle = 'G-RP [mag]',ytitle = 'G [mag]',title = 'color-color magnitude of ASCC21',xrange=[-0.5,2.0],yrange=[max(GRP21),5.],ystyle=1,/nodata
oplot,GRP21,Gmag21,psym = 3,color = cgcolor('gray')
oplot,GRP21[px21],Gmag21[py21],psym = 1,color = cgcolor('red')

;;; color-color magnitude plot of ASCC16 and  ASCC21
;window,11,xsize=900,ysize=900
;plot,GRP16,Gmag16,xtitle = 'G-RP [mag]',ytitle = 'G [mag]',title = 'color-color magnitude of ASCC16',xrange=[-0.5,2.0],yrange=[max(GRP16),5.],ystyle=1,/nodata
;oplot,GRP16,Gmag16,psym = 3,color = cgcolor('gray')
;oplot,GRP16[px16],Gmag16[py16],psym = 1,color = cgcolor('red')
;oplot,GRP21,Gmag21,psym = 3,color = cgcolor('blue')
;oplot,GRP21[px21],Gmag21[py21],psym = 1,color = cgcolor('green')


































end