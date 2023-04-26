!p.background=255
!p.color=0
!p.charsize=2.0
!p.thick=1.5
!p.charthick=1.0
device,decompose = 0
!p.font=-1
red = [255,0,0]
green = [0,255,0]
blue = [0,0,255]
TVLCT,red,green,blue,1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CD, 'D:\學校資料\天文所陳老師\h and x'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;read Xmatch data: GAIA eDR3 & B-J distance
;match = READ_ASCII('D:\學校資料\天文所陳老師\h and x\h and chi cross match Apr-13\Xmat_result.txt', DATA_START=0)
;
;;;;Transfer structure to array
;RA = [float(match.FIELD01[0,*])] ;
;DEC = [float(match.FIELD01[1,*])]      ;
;ID = [float(match.FIELD01[2,*])]     ;
;Plx = [float(match.FIELD01[3,*])]      ;[mas]
;ePlx = [float(match.FIELD01[4,*])]   ;[mas]
;PM = [float(match.FIELD01[5,*])]    ;[mas/year]
;pmRA = [float(match.FIELD01[6,*])]   ;[mas/yr]
;pmDEC = [float(match.FIELD01[7,*])]   ;[mas/yr]
;Gmag = [float(match.FIELD01[8,*])]    ;
;eGmag = [float(match.FIELD01[9,*])]  ;
;BPmag = [float(match.FIELD01[10,*])] ;
;eBPmag = [float(match.FIELD01[11,*])] ;
;RPmag = [float(match.FIELD01[12,*])] ;
;eRPmag = [float(match.FIELD01[13,*])] ;
;rgeo = [float(match.FIELD01[14,*])] ; [unit: pc]
;rgep16 = [float(match.FIELD01[15,*])] ;
;rgeo84 = [float(match.FIELD01[16,*])] ;
;rpgeo = [float(match.FIELD01[17,*])] ;
;rpgeo16 = [float(match.FIELD01[18,*])] ;
;rpgeo84 = [float(match.FIELD01[19,*])] ;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;IDL window 4: cluster member selection with 2291<d<2399
;IDL window 8: cluster member selection with 2291<d<2399, and pm in circle (win. 6 )
;IDL window 17: cluster member selection with 2291<d<2399, and pm in circle (win. 6 ) + $
;filed star has the same criteria with cluster member selection
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ra0 = MEAN(RA)
de0 = MEAN(DEC)

;; REF: https://www.aanda.org/articles/aa/pdf/2014/08/aa22720-13.pdf
;; NGC869 = h
radcen869 = 14.4/60. ;~0.24 deg (radius)
ra869 = 34.75
dec869 = 57.13

;; NGC884 = x
radcen884 = 10.5/60. ;~0.175 deg (radius)
ra884 = 35.50 ;35.60
dec884 = 57.13

rad869 = SQRT( (RA-ra869)*(RA-ra869) + (DEC-dec869)*(DEC-dec869) )
in869 = WHERE(rad869 LT radcen869)
out869 = WHERE(rad869 GT radcen869)

rad884 = SQRT( (RA-ra884)*(RA-ra884) + (DEC-dec884)*(DEC-dec884) )
in884 = WHERE(rad884 LT radcen884)
out884 = WHERE(rad884 GT radcen884)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; set a criteria of cluster region 
!P.Multi = [0, 1, 1]

ind = where(RA gt 34.0 and RA lt 36 and DEC gt 55.5 and DEC lt 57.5) ;

window,1,xsize=900,ysize=900
plot,RA[ind],DEC[ind],psym=3,xtitle='RA (deg)', ytitle='DEC (deg)', ystyle=1 ;,xstyle=1 ;,title='Xmatch result, N=1210298' ;,xrange=[33.5,36.5],yrange=[55.5,58.5],title='Coordiante of Stars Observed by GAIA eDR3'
TVCIRCLE, 0.24, ra869, dec869, color=1,THICK = 2
TVCIRCLE, 0.175, ra884, dec884, color=3,THICK = 2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; circle out a region of cluster core
window,2,xsize=900,ysize=900
PLOT,RA,DEC, psym=3, xtit='RA', ytit='DEC',title='NGC869 and NGC884',xrange=[34,36],yrange=[55.5,57.5],ystyle = 1 ;
TVCIRCLE, 0.24, ra869, dec869, color=1,THICK = 2
TVCIRCLE, 0.175, ra884, dec884, color=3,THICK = 2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; tangential velocity calculation (ref: https://apsu.edu/physics/astronomy/stellar-labs/proper-motion-of-a-star.pdf)
vtRA = 4.75*pmRA*0.001*rgeo ;pm unit: mas, rgeo unit: pc
vtDEC = 4.75*pmDEC*0.001*rgeo

;ref of distance of h and x:https://iopscience.iop.org/article/10.1086/341865/pdf (PAGE: 8)
;d = 7600 ly or 2344+55-53 pc = 2291-2399 pc

Nrgeo869 = rgeo
Nrgeo869[out869] = !VALUES.F_NAN
dind869 = where(Nrgeo869 gt 2291 and Nrgeo869 lt 2399) ; and rgeo[in869]
vtRA869 = vtRA[dind869]
vtDEC869 = vtDEC[dind869]

Nrgeo884 = rgeo
Nrgeo884[out884] = !VALUES.F_NAN
dind884 = where(Nrgeo884 gt 2291 and Nrgeo884 lt 2399)
vtRA884 = vtRA[dind884]
vtDEC884 = vtDEC[dind884]

window,3,xsize=900,ysize=900
plot,vtRA,vtDEC,psym = 3,xtitle = 'Vt_RA (km/yr)',ytitle = 'Vt_DEC (km/yr)',title = 'Tengential velocity', $
  xstyle=1,ystyle=1,xrange=[-300,300],yrange=[-300,300]
oplot,vtRA869,vtDEC869,psym = 1,color = cgcolor('red')
oplot,vtRA884,vtDEC884,psym = 1,color = cgcolor('blue')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; cluster position check
window,4,xsize=900,ysize=900
PLOT,RA,DEC, psym=3, xtit='RA', ytit='DEC',title='NGC869 and NGC884',xrange=[34,36],yrange=[55.5,57.5],ystyle = 1,xstyle = 1,/isotropic
TVCIRCLE, 0.24, ra869, dec869, color=1,THICK = 2
TVCIRCLE, 0.175, ra884, dec884, color=3,THICK = 2
oplot,RA[dind869],DEC[dind869],psym = 1,color = cgcolor('red')
oplot,RA[dind884],DEC[dind884],psym = 1,color = cgcolor('blue')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
BPRP = BPmag - RPmag

;; CMD
window,5,xsize=900,ysize=900
plot,BPRP,Gmag,xtitle = 'G-RP [mag]',ytitle = 'G [mag]',title = 'CMD',/nodata,xrange=[-0.5,3.0],yrange=[max(Gmag),5.],ystyle=1
oplot,BPRP[0,*],Gmag[0,*],psym = 3,color = cgcolor('gray')
oplot,BPRP[dind869],Gmag[dind869],psym = 1,color = cgcolor('red')
oplot,BPRP[dind884],Gmag[dind884],psym = 1,color = cgcolor('blue')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Proper motion 
window,6,xsize=900,ysize=900
plot,pmRA,pmDEC,psym = 3,xtitle = 'pmRA (mas/yr)',ytitle = 'pmDEC (mas/yr)',title = 'Proper Motion',xrange=[-5,5],yrange=[-5,5],/nodata
oplot,pmRA,pmDEC,psym = 3,color = cgcolor('gray')
oplot,pmRA[dind869],pmDEC[dind869],psym = 1,color = cgcolor('red')
oplot,pmRA[dind884],pmDEC[dind884],psym = 1,color = cgcolor('blue') 
TVCIRCLE, 0.5, -0.7, -1.1, color=2,THICK = 2 ;(-0.6, -1.2)

;calculate proper motion or every star which distance to (-0.6,-1.2)
;; this calculation include red and blue cluster member candidates and FIELD stars(gray dots)
pmind = SQRT( (pmRA+0.7)*(pmRA+0.7) + (pmDEC+1.1)*(pmDEC+1.1) ) 
pmin = WHERE(pmind LT 0.5)

;; Proper motion circle
window,7,xsize=900,ysize=900
plot,pmRA,pmDEC,psym = 3,xtitle = 'pmRA (mas/yr)',ytitle = 'pmDEC (mas/yr)',title = 'Proper Motion',/isotropic,xrange=[-5,5],yrange=[-5,5],/nodata
oplot,pmRA,pmDEC,psym = 3,color = cgcolor('gray')
oplot,pmRA[pmin],pmDEC[pmin],psym = 1,color = cgcolor('green')
TVCIRCLE, 0.5, -0.7, -1.1, color=1,THICK = 2
TVCIRCLE, 0.001, -0.7, -1.1, color=1,THICK = 2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; coordinate + proper motion + distance
Nrgeo869 = rgeo
Nrgeo869[out869] = !VALUES.F_NAN
dpin869 = where(Nrgeo869 gt 2291 and Nrgeo869 lt 2399 and pmind LT 0.5) ;pmind LT 0.4

Nrgeo884 = rgeo
Nrgeo884[out884] = !VALUES.F_NAN
dpin884 = where(Nrgeo884 gt 2291 and Nrgeo884 lt 2399 and pmind LT 0.5) ;pmind LT 0.4

window,8,xsize=900,ysize=900
PLOT,RA,DEC, psym=3, xtit='RA', ytit='DEC',title='NGC869 and NGC884',xrange=[34,36],yrange=[55.5,57.5],ystyle = 1,xstyle = 1,/isotropic
TVCIRCLE, 0.24, ra869, dec869, color=1,THICK = 2
TVCIRCLE, 0.175, ra884, dec884, color=3,THICK = 2
TVCIRCLE, 14.4/60., 35, 56.25, color=1,THICK = 2 ;field star region
TVCIRCLE, 10.5/60., 35, 56.25, color=3,THICK = 2
oplot,RA[dpin869],DEC[dpin869],psym = 1,color = cgcolor('red')
oplot,RA[dpin884],DEC[dpin884],psym = 1,color = cgcolor('blue')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; proper motion histogram
; field star number 
rad869 = SQRT( (RA-35.)*(RA-35.) + (DEC-56.25)*(DEC-56.25) )
field869 = WHERE(rad869 LT (14.4/60.) )
fout869 = WHERE(rad869 GT (14.4/60.) )

rad884 = SQRT( (RA-35.)*(RA-35.) + (DEC-56.25)*(DEC-56.25) )
field884 = WHERE(rad884 LT (10.5/60.) )
fout884 = WHERE(rad884 GT (10.5/60.) )

;; distance between 2291-2399 pc
window,9,xsize=900,ysize=900
plothist,PM[field869],bin = 0.2,xrange=[0,5],title='NGC869 proper motion vector' ;; central part of cluster (red),yrange=[0,200]
plothist,PM[dind869],bin = 0.2,color = 1,/over

window,10,xsize=900,ysize=900
plothist,PM[field884],bin = 0.2,xrange=[0,5],title='NGC884 proper motion vector' 
plothist,PM[dind884],bin = 0.2,color = 3,/over

;; proper motion vecotr of PM + distance 
window,11,xsize=900,ysize=900
plothist,PM[field869],bin = 0.1,xrange=[0,5],title='NGC869 proper motion vector',xticklen='-0.01' ;,yrange=[0,1000]
plothist,PM[dpin869],bin = 0.1,color = 1,/over
xyouts,3.5,190,'Red: NGC869(h Persei)',color = 1

window,12,xsize=900,ysize=900
plothist,PM[field884],bin = 0.1,xrange=[0,5],title='NGC884 proper motion vector',xticklen='-0.01' ;,yrange=[0,1000]
plothist,PM[dpin884],bin = 0.1,color = 3,/over
xyouts,3.5,110,'Blue: NGC884('+ cgGreek('chi')+ ' Persei)',color = 3

;; proper motion in RA and DEC of PM + distance
window,13,xsize=900,ysize=900
plothist,pmRA[field869],bin = 0.1,xrange=[-5,5],yrange=[0,180],title='NGC869 proper motion in RA. direction',xticklen='-0.01' ;
plothist,pmRA[dpin869],bin = 0.1,color = 1,/over
xyouts,2.5,150,'Red: NGC869(h Persei)',color = 1

window,14,xsize=900,ysize=900
plothist,pmDEC[field869],bin = 0.1,xrange=[-5,5],yrange=[0,180],title='NGC869 proper motion in DEC. direction',xticklen='-0.01' ;
plothist,pmDEC[dpin869],bin = 0.1,color = 1,/over
xyouts,2.5,150,'Red: NGC869(h Persei)',color = 1

window,15,xsize=900,ysize=900
plothist,pmRA[field884],bin = 0.1,xrange=[-5,5],yrange=[0,100],title='NGC884 proper motion in RA. direction',xticklen='-0.01' ;
plothist,pmRA[dpin884],bin = 0.1,color = 3,/over
xyouts,2.5,90,'Blue: NGC869(h Persei)',color = 3

window,16,xsize=900,ysize=900
plothist,pmDEC[field884],bin = 0.1,xrange=[-5,5],yrange=[0,105],title='NGC884 proper motion in DEC. direction',xticklen='-0.01' ;
plothist,pmDEC[dpin884],bin = 0.1,color = 3,/over
xyouts,2.5,90,'Blue: NGC869(h Persei)',color = 3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; cluster region histogram (window,1)
!P.Multi = [0, 2, 2]

NRA = RA[ind]
NDEC = DEC[ind]

window,17,xsize=900,ysize=900
PLOT,NRA,NDEC, psym=3, xtit='RA', ytit='DEC',xrange=[34,36],yrange=[55.5,57.5],ystyle = 1,xstyle = 1, $
  /isotropic,position = [0.07,0.07,0.70,0.70]
TVCIRCLE, 0.24, ra869, dec869, color=1,THICK = 2
TVCIRCLE, 0.175, ra884, dec884, color=3,THICK = 2
TVCIRCLE, 14.4/60., 35, 56.25, color=1,THICK = 2
TVCIRCLE, 10.5/60., 35, 56.25, color=3,THICK = 2
oplot,RA[dpin869],DEC[dpin869],psym = 1,color = cgcolor('red')
oplot,RA[dpin884],DEC[dpin884],psym = 1,color = cgcolor('blue')

plothist,NRA,bin = 0.1,position = [0.07,0.70,0.70,0.95],BOXPLOT=0, $
  title='NGC869(h Persei), NGC884('+ cgGreek('chi') + ' Persei)',xtitle=' ',yrange=[5800,6600],xrange=[34.,36.],XTICKFORMAT="(A1)"
plothist,NDEC,bin = 0.1,position = [0.70,0.07,0.95,0.70],BOXPLOT=0,ytitle = ' ',/rotate,xrange=[5500,8500],yrange = [55.5,57.5], $
  YTICKFORMAT="(A1)" ;,XTICKFORMAT="(A1)"
oPLOT,NRA,NDEC, psym=3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; find """filed star""" whcih satisfied the same condition with double cluster
Nrgeo869 = rgeo
Nrgeo869[fout869] = !VALUES.F_NAN
fin869 = where(Nrgeo869 gt 2291 and Nrgeo869 lt 2399 and pmind LT 0.5) ;pmind LT 0.4

Nrgeo884 = rgeo
Nrgeo884[fout884] = !VALUES.F_NAN
fin884 = where(Nrgeo884 gt 2291 and Nrgeo884 lt 2399 and pmind LT 0.5) ;pmind LT 0.4

!P.Multi = [0, 1, 1]

window,18,xsize=900,ysize=900
PLOT,RA,DEC, psym=3, xtit='RA', ytit='DEC',title='NGC869 and NGC884',xrange=[34,36],yrange=[55.5,57.5],ystyle = 1,xstyle = 1,/isotropic
TVCIRCLE, 0.24, ra869, dec869, color=1,THICK = 2
TVCIRCLE, 0.175, ra884, dec884, color=3,THICK = 2
TVCIRCLE, 14.4/60., 35, 56.25, color=1,THICK = 2 ;field star region
TVCIRCLE, 10.5/60., 35, 56.25, color=3,THICK = 2
oplot,RA[dpin869],DEC[dpin869],psym = 1,color = cgcolor('red') ; plot cluster member candidates
oplot,RA[dpin884],DEC[dpin884],psym = 1,color = cgcolor('blue')
oplot,RA[fin869],DEC[fin869],psym = 1,color = cgcolor('red') ; plot field star with condition satisfied with cluster
oplot,RA[fin884],DEC[fin884],psym = 1,color = cgcolor('green')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; cluster shape fitting
deg = findgen(360)
tilt = -30 ;[deg]
a = 0.24
b = 0.15
xcen = ra869
ycen = dec869
xx869 = xcen + (a*cos(!dtor*deg)*cos(!dtor*tilt)-b*sin(!dtor*deg)*sin(!dtor*tilt))
yy869 = ycen + (a*cos(!dtor*deg)*sin(!dtor*tilt)+b*sin(!dtor*deg)*cos(!dtor*tilt))

tilt = 20 ;[deg]
a = 0.17
b = 0.09
xcen = ra884
ycen = dec884
xx884 = xcen + (a*cos(!dtor*deg)*cos(!dtor*tilt)-b*sin(!dtor*deg)*sin(!dtor*tilt))
yy884 = ycen + (a*cos(!dtor*deg)*sin(!dtor*tilt)+b*sin(!dtor*deg)*cos(!dtor*tilt))

window,19,xsize=900,ysize=900
PLOT,RA,DEC, psym=3, xtit='RA', ytit='DEC',title='NGC869 and NGC884',xrange=[34,36],yrange=[55.5,57.5],ystyle = 1,xstyle = 1,/isotropic
TVCIRCLE, 0.24, ra869, dec869, color=1,THICK = 2
TVCIRCLE, 0.175, ra884, dec884, color=3,THICK = 2
TVCIRCLE, 14.4/60., 35, 56.25, color=1,THICK = 2 ;field star region
TVCIRCLE, 10.5/60., 35, 56.25, color=3,THICK = 2
oplot,RA[dpin869],DEC[dpin869],psym = 1,color = cgcolor('red')
oplot,RA[dpin884],DEC[dpin884],psym = 1,color = cgcolor('blue')
oplot,RA[fin869],DEC[fin869],psym = 1,color = cgcolor('red') ; plot field star with condition satisfied with cluster
oplot,RA[fin884],DEC[fin884],psym = 1,color = cgcolor('green')

oplot,xx869,yy869,color=2,THICK = 2
oplot,xx884,yy884,color=2,THICK = 2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; CMD
;BPRP = BPmag - RPmag

;; CMD of NGC869 
window,20,xsize=900,ysize=900
plot,BPRP[ind],Gmag[ind],xtitle = 'BP-RP [mag]',ytitle = 'G [mag]',title = 'CMD of NGC869',/nodata,xrange=[-0.5,3.0],yrange=[max(Gmag),5.],ystyle=1
oplot,BPRP[0,ind],Gmag[0,ind],psym = 3,color = cgcolor('gray')
oplot,BPRP[dpin869],Gmag[dpin869],psym = 1,color = cgcolor('red')

;;; CMD of NGC884
window,21,xsize=900,ysize=900
plot,BPRP[ind],Gmag[ind],xtitle = 'BP-RP [mag]',ytitle = 'G [mag]',title = 'CMD of NGC884',/nodata,xrange=[-0.5,3.0],yrange=[max(Gmag),5.],ystyle=1
oplot,BPRP[0,ind],Gmag[0,ind],psym = 3,color = cgcolor('gray')
oplot,BPRP[dpin884],Gmag[dpin884],psym = 1,color = cgcolor('blue')

;;;CMD of field stars
window,22,xsize=900,ysize=900
plot,BPRP[ind],Gmag[ind],xtitle = 'BP-RP [mag]',ytitle = 'G [mag]',title = 'CMD of Reference field',/nodata,xrange=[-0.5,3.0],yrange=[max(Gmag),5.],ystyle=1
oplot,BPRP[0,ind],Gmag[0,ind],psym = 3,color = cgcolor('gray')
oplot,BPRP[fin869],Gmag[fin869],psym = 6,color = cgcolor('red')

;; cymbol size calculation
A = FIndGen(5) * (!PI*2/5.)
UserSym, cos(A), sin(A), /fill
oplot,BPRP[fin884],Gmag[fin884],psym = 8,color = cgcolor('blue')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; absolute magnitude calculation
array = fltarr(n_elements(Gmag)) + 5*ALOG10(rgeo) - 5 + 1.55 ;reddening A_{G} = 1.55
MG = fltarr(n_elements(Gmag))
MG = Gmag[0,*]-array

array = fltarr(n_elements(BPmag)) + 5*ALOG10(rgeo) - 5 + 1.91 ;reddening A_{BP} = 1.91
MBP = fltarr(n_elements(BPmag))
MBP = BPmag[0,*]-array

array = fltarr(n_elements(RPmag)) + 5*ALOG10(rgeo) - 5 + 1.17  ;reddening A_{RP} = 1.17
MRP = fltarr(n_elements(RPmag))
MRP = RPmag[0,*]-array

MBPRP = MBP - MRP

;; luminosity function calculation
L = fltarr(n_elements(Gmag))
Ls = fltarr(n_elements(Gmag))+4.83 ; 4.83 = solar absolute magnitude
L = 10^(0.4*(Ls-MG)) ;luminosity ratio

;; H-R diagram
window,23,xsize=1000,ysize=900
plot,MBPRP[dpin869],MG[dpin869],xtitle = 'M!iBP!n-M!iRP!n [mag]',ytitle = 'M!iG!n',title = 'H-R Diagram',/nodata, $
  xrange=[-1.5,3.0],yrange=[8,-6],ystyle=8,xstyle=1,position = [0.1,0.1,0.90,0.90];,/isotropic ;,ystyle=16 ;,YTICKFORM="(A1)" ;
oplot,MBPRP[0,dpin869],MG[0,dpin869],psym = 1,color = cgcolor('red')
Axis, YAxis=1,YLog=1, YRange=[min(L[ind]), max(L[ind])], /Save,YMinor=10,ytitle='Luminosity (L!i*!n/L!i'+sunsymbol()+'!n)' ;;, ystyle=1 ;

window,24,xsize=1000,ysize=900
plot,MBPRP[dpin884],MG[dpin884],xtitle = 'M!iBP!n-M!iRP!n [mag]',ytitle = 'M!iG!n',title = 'H-R Diagram',/nodata, $
  xrange=[-1.5,3.0],yrange=[8,-6],ystyle=8,xstyle=1,position = [0.1,0.1,0.90,0.90];,/isotropic ;,ystyle=16 ;,YTICKFORM="(A1)" ;
oplot,MBPRP[0,dpin884],MG[0,dpin884],psym = 1,color = cgcolor('blue')
Axis, YAxis=1,YLog=1, YRange=[min(L[ind]), max(L[ind])], /Save,YMinor=10,ytitle='Luminosity (L!i*!n/L!i'+sunsymbol()+'!n)' ;;, ystyle=1 ;

window,25,xsize=1000,ysize=900
plot,MBPRP[fin869],MG[fin869],xtitle = 'M!iBP!n-M!iRP!n [mag]',ytitle = 'M!iG!n',title = 'H-R Diagram',/nodata, $
  xrange=[-1.5,3.0],yrange=[8,-6],ystyle=8,xstyle=1,position = [0.1,0.1,0.90,0.90];,/isotropic ;,ystyle=16 ;,YTICKFORM="(A1)" ;
oplot,MBPRP[0,fin869],MG[0,fin869],psym = 6,color = cgcolor('red')
Axis, YAxis=1,YLog=1, YRange=[min(L[ind]), max(L[ind])], /Save,YMinor=10,ytitle='Luminosity (L!i*!n/L!i'+sunsymbol()+'!n)' ;;, ystyle=1 ;

window,26,xsize=1000,ysize=900
plot,MBPRP[fin884],MG[fin884],xtitle = 'M!iBP!n-M!iRP!n [mag]',ytitle = 'M!iG!n',title = 'H-R Diagram',/nodata, $
  xrange=[-1.5,3.0],yrange=[8,-6],ystyle=8,xstyle=1,position = [0.1,0.1,0.90,0.90];,/isotropic ;,ystyle=16 ;,YTICKFORM="(A1)" ;
oplot,MBPRP[0,fin884],MG[0,fin884],psym = 8,color = cgcolor('blue')
Axis, YAxis=1,YLog=1, YRange=[min(L[ind]), max(L[ind])], /Save,YMinor=10,ytitle='Luminosity (L!i*!n/L!i'+sunsymbol()+'!n)' ;;, ystyle=1 ;

;; absolute magnitude histogram of all stars + NGC869 + NGC884
window,27,xsize=900,ysize=900
plothist,MG[dpin869],xhist,yhist,bin = 0.5,title='NGC869 Histogram',xtitle='Absolute Magnitude ', ytitle='count',yrange=[0,40],xrange=[-6,7], $
  color = cgcolor('red'), xticklen = '-0.02'
plothist,MG[fin869],xhist,yhistf,bin = 0.5,xtitle=' ', /OVERPLOT,color = cgcolor('green')

;; absolute magnitude histogram of NGC869 + NGC884 + field stars
window,28,xsize=900,ysize=900
plothist,MG[dpin884], bin=0.5, title='NGC884 Histogram',xtitle='Absolute Magnitude ', ytitle='count',yrange=[0,20],xrange=[-8,8],color = cgcolor('blue') ;,BOXPLOT=0
plothist,MG[fin884],bin = 0.5,xtitle=' ', /OVERPLOT,color = cgcolor('green') ;,BOXPLOT=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Mass function calculation
Mass = (L/1.4)^(1/3.5)

;;; Mass function hitogram of NGC869
window,29,xsize=900,ysize=900
plothist,Mass[dpin869],bin=0.5,yrange=[0,70],xrange=[0,12],xtitle='Mass (M!i*!n/M!i'+sunsymbol()+'!n)', ytitle='count',title='NGC869 Mass Function Histogram'

;;; Mass function histogram of NGC884
window,30,xsize=900,ysize=900
plothist,Mass[dpin884],bin=0.5,yrange=[0,40],xrange=[0,12],xtitle='Mass (M!i*!n/M!i'+sunsymbol()+'!n)', ytitle='count',title='NGC884 Mass Function Histogram'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
xx = findgen(28)/2-5.25
res = [1,0,2,2,0,2,5,1,4,5,12,19,23,30,39,28,21,17,19,7,5,6,2,2,3]
window,2,xsize=900,ysize=900
plot,xx,res,psym=10,xstyle=1,xrange=[-6.0,7],xticklen = '-0.02'



end