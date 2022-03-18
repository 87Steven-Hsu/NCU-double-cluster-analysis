!p.background=255
!p.color=0
!p.charsize=1.8
!p.thick=1.5
!p.charthick=1.0
device,decompose = 0
!p.font=-1
red = [255,0,0]
green = [0,255,0]
blue = [0,0,255]
TVLCT,red,green,blue,1
loadct,3
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CD, 'D:\學校資料\天文所陳老師\h and x'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;read Xmatch data: GAIA eDR3 & B-J distance
;match = READ_ASCII('D:\學校資料\天文所陳老師\h and x\h and chi cross match Jun-2\Xmat_result_Jun-4.txt', DATA_START=0)
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; REF: https://www.aanda.org/articles/aa/pdf/2014/08/aa22720-13.pdf
;;NGC869 = h
radcen869 = 14.4/60. ;~0.24 deg (radius)
ra869 = 34.75
dec869 = 57.13

;; NGC884 = x
radcen884 = 10.5/60. ;~0.175 deg (radius)
ra884 = 35.50 ;35.60
dec884 = 57.13
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
window,0,xsize=900,ysize=900
PLOT,RA,DEC,psym=3, xtit='RA', ytit='DEC',title='NGC869 and NGC884', xrange=[max(RA),min(RA)],xstyle = 1,ystyle = 1,/isotropic

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
MRA = mean(RA[0,1:*])

NRA = (RA-MRA)/15.
RANew = NRA*cos(DEC*!dtor)*15

;; set a criteria of a rectangle range of star candidate
ind = where(RANew gt (34.0-MRA) and RANew lt (36.25-MRA) and DEC gt 56.0 and DEC lt 59) ;(RA gt 33.5 and RA lt 36.5 and DEC gt 56.0 and DEC lt 59)
NRA = RANew[ind]
NDEC = DEC[ind]

Nind  = where(NRA gt (34.-MRA) and NRA lt (36.25-MRA) and NDEC gt 56.5 and NDEC lt 57.75)

window,1,xsize=900,ysize=900
PLOT,RANew,DEC,psym=3, xtit='RA difference (deg)', ytit='DEC (deg)',title='NGC869 and NGC884', xrange=[max(RANew),min(RANew)],xstyle = 1,ystyle = 1,/isotropic
oplot,NRA[Nind],NDEC[Nind],psym = 1,color = cgcolor('red') ; plot cluster member candidates

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; coordinate criteria + distance (2291-2399)
dis1 = 2300 ;2291
dis2 = 2400 ;2399
Nind  = where(RANew gt (34.-MRA) and RANew lt (36.25-MRA) and DEC gt 56.5 and DEC lt 57.75 and rgeo gt dis1 and rgeo lt dis2)

;; Coordinate + distance
window,2,xsize=900,ysize=900
PLOT,RANew,DEC, psym=3, xtit='RA difference (deg)', ytit='DEC (deg)',title='NGC869 and NGC884',xrange=[max(RANew),min(RANew)], $
  yrange=[55.,59.0],ystyle = 1,xstyle = 1,/isotropic
oplot,RANew[Nind],DEC[Nind],psym = 1,color = cgcolor('red') ; plot cluster member candidates
TVCIRCLE, 0.24, ra869, dec869, color=1,THICK = 2
TVCIRCLE, 0.175, ra884, dec884, color=3,THICK = 2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Proper Motion
window,3,xsize=900,ysize=900
plot,pmRA,pmDEC,psym = 3,xtitle = 'pmRA (mas/yr)',ytitle = 'pmDEC (mas/yr)',title = 'Proper Motion',xrange=[-5,5],yrange=[-5,5],/nodata
oplot,pmRA,pmDEC,psym = 3,color = cgcolor('gray')
oplot,pmRA[Nind],pmDEC[Nind],psym = 1,color = cgcolor('red')
;RofPM = 0.5
;TVCIRCLE, RofPM, -0.6, -1.15, color=2,THICK = 2 ;(-0.6, -1.2)

WRITE_JPEG,"PM.png",TVRD(/TRUE),QUALITY=100, /true

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Gaussian fit of proper motion
window,4,xsize=900,ysize=900
plothist,pmRA[Nind],xhist,yhist,bin=0.1,xticklen=-0.01,/noplot
plot,xhist,yhist,psym=10,xtitle='pmRA of candidates',ytitle='Count',xticklen=-0.01,xrange=[-2,1]
;gau =  GAUSSFIT(xhist, yhist, coeffRA, NTERMS=6)
;oplot,xhist,gau,color=cgcolor('red')
;oplot,[coeffRA[1],coeffRA[1]],!Y.CRANGE,linestyle=2,color=cgcolor('blue'),thick=2
;oplot,[coeffRA[1]-coeffRA[2],coeffRA[1]-coeffRA[2]],!Y.CRANGE,linestyle=2,color=cgcolor('green'),thick=2
;oplot,[coeffRA[1]+coeffRA[2],coeffRA[1]+coeffRA[2]],!Y.CRANGE,linestyle=2,color=cgcolor('green'),thick=2
;oplot,[coeffRA[1]-2*coeffRA[2],coeffRA[1]-2*coeffRA[2]],!Y.CRANGE,linestyle=2,color=cgcolor('green'),thick=2
;oplot,[coeffRA[1]+2*coeffRA[2],coeffRA[1]+2*coeffRA[2]],!Y.CRANGE,linestyle=2,color=cgcolor('green'),thick=2

;;158 = [ -1( 起算的地方之值(右邊數值) ) - (-16.8( bin的最小值 ) ) ]*10(bin大小) + 1 (種樹問題) -1 (IDL sart from 0)
;; 158 =  [-1-(-16.8)]*10+1-1
;gau1 =  GAUSSFIT(xhist[158:167], yhist[158:167], NTERMS=6) ;xhist[158:167], yhist[158:167]
;oplot,xhist[158:167],gau1,color=cgcolor('blue')

;; fit background
plothist,pmRA[Nind],xhist,yhist,bin=0.1,xticklen=-0.01,/noplot
yhist = float(yhist)
yhist[158:167] = !VALUES.F_NAN  ;F_NAN
gau2 =  MPFITPEAK(xhist[149:178], yhist[149:178], coeff2, NTERMS=6,/NAN) ;; fit background (dot line)
oplot,xhist[149:178],gau2,linestyle= 1 ;; fit background (dot line)

;; fit background + peak
plothist,pmRA[Nind],xhist,yhist,bin=0.1,xticklen=-0.01,/noplot
cc = fltarr(30)
cc[9:18] = yhist[158:167]-gau2[9:18]
gau3 =  MPFITPEAK(xhist[149:178], cc, coeff3, NTERMS=6,/NAN) ;; fit peak (dashed line)
;oplot,xhist[149:178],gau3,linestyle= 2 ;; fit peak (dashed line)
;oplot,xhist[149:178],gau2+gau3,color=cgcolor('red') ;; combine peak and background

xx = findgen(90)/30.7-1.95

yhis = INTERPOL(yhist[158:167],30) ;; interpol peak to 90 elements
fie = INTERPOL(gau2[9:18],30) ;; interpol background peak to 90 elements
cc = fltarr(90)
cc[27:56] = yhis - fie ;;cc[30:59] = yhis - fie subtrack peak and background peak
gau4 =  MPFITPEAK(xx, cc, coeff4, NTERMS=6,/NAN) ;; fit subtrack peak and background peak
oplot,xx,gau4,linestyle= 2 ;,color=cgcolor('blue') 

gau2int = INTERPOL(gau2,90) ;; interpol background to 90 elements
;oplot,xx,gau2int,color=cgcolor('green')

oplot,xx,gau2int+gau4,color=cgcolor('red') ;; combine inter. background and fit subtrack peak and background peak


WRITE_JPEG,"pmRA fit_ver2.png",TVRD(/TRUE),QUALITY=100, /true

window,5,xsize=900,ysize=900
plothist,pmDEC[Nind],xhist,yhist,bin=0.1,xticklen=-0.01,/noplot
plot,xhist,yhist,psym=10,xtitle='pmDEC of candidates',ytitle='Count',xticklen=-0.01,xrange=[-2.6,0.4],xstyle=1
;gau =  GAUSSFIT(xhist, yhist, coeffDec, NTERMS=6)
;oplot,xhist,gau,color=cgcolor('red')
;oplot,[coeffDec[1],coeffDec[1]],!Y.CRANGE,linestyle=2,color=cgcolor('blue'),thick=2
;oplot,[coeffDec[1]-coeffDec[2],coeffDec[1]-coeffDec[2]],!Y.CRANGE,linestyle=2,color=cgcolor('green'),thick=2
;oplot,[coeffDec[1]+coeffDec[2],coeffDec[1]+coeffDec[2]],!Y.CRANGE,linestyle=2,color=cgcolor('green'),thick=2
;oplot,[coeffDec[1]-2*coeffDec[2],coeffDec[1]-2*coeffDec[2]],!Y.CRANGE,linestyle=2,color=cgcolor('green'),thick=2
;oplot,[coeffDec[1]+2*coeffDec[2],coeffDec[1]+2*coeffDec[2]],!Y.CRANGE,linestyle=2,color=cgcolor('green'),thick=2

;;158 =  [-1-(-32.3)]*10+1-1
;gau =MPFITPEAK(xhist[310:314], yhist[310:314], NTERMS=3,/NAN)  
;oplot,xhist[310:314],gau,color=cgcolor('blue')

yhist = float(yhist)
yhist[310:314] = !VALUES.F_NAN  ;F_NAN
gau2 =  MPFITPEAK(xhist[298:327], yhist[298:327],coeff2, NTERMS=6,/NAN)
oplot,xhist[298:327],gau2,linestyle= 1

plothist,pmDEC[Nind],xhist,yhist,bin=0.1,xticklen=-0.01,/noplot
cc = fltarr(30)
cc[12:16] = yhist[310:314]-gau2[12:16]
;gau3 = GAUSSFIT(xhist[298:327], cc, coeff3, NTERMS=6)
gau3 =  MPFITPEAK(xhist[298:327], cc,coeff3, NTERMS=6,/NAN,/GAUSSIAN)
;oplot,xhist[298:327],gau3,linestyle= 2
;oplot,xhist[298:327],gau2+gau3,color=cgcolor('red')

xx = findgen(90)/31.22-2.5
yhis = INTERPOL(yhist[310:314],15) ;; interpol peak to 90 elements
fie = INTERPOL(gau2[12:16],15) ;; interpol background peak to 90 elements
cc = fltarr(90)
cc[35:49] = yhis - fie ;;cc[30:59] = yhis - fie subtrack peak and background peak
gau4 =  MPFITPEAK(xx, cc, coeff4, NTERMS=6,/NAN) ;; fit subtrack peak and background peak
oplot,xx,gau4,linestyle= 2 ;,color=cgcolor('blue') 

gau2int = INTERPOL(gau2,90) ;; interpol background to 90 elements
;oplot,xx,gau2int,color=cgcolor('green')

oplot,xx,gau2int+gau4,color=cgcolor('red')

WRITE_JPEG,"pmDEC fit_ver2.png",TVRD(/TRUE),QUALITY=100, /true

;PMV = sqrt(pmRA[Nind]*pmRA[Nind] + pmDEC[Nind]*pmDEC[Nind])
;theta = atan(pmDEC[Nind]/pmRA[Nind]) ;/!dtor => transfer to deg

;sele = where(pmRA[Nind] gt (coeffRA[1]-coeffRA[2]) and pmRA[Nind] lt (coeffRA[1]+coeffRA[2]) and $
;  pmDec[Nind] gt (coeffDec[1]-coeffDec[2]) and pmDec[Nind] lt (coeffDec[1]+coeffDec[2]))

;sele = where(pmRA[Nind] gt (coeffRA[1]-2*coeffRA[2]) and pmRA[Nind] lt (coeffRA[1]+2*coeffRA[2]) and $
;  pmDec[Nind] gt (coeffDec[1]-2*coeffDec[2]) and pmDec[Nind] lt (coeffDec[1]+2*coeffDec[2]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; plot PM density map
;window,6,xsize=900,ysize=900
;plot,pmRA,pmDEC,psym = 3,xtitle = 'pmRA (mas/yr)',ytitle = 'pmDEC (mas/yr)',title = 'Proper Motion',xrange=[-5,5],yrange=[-5,5],/nodata
;oplot,pmRA,pmDEC,psym = 3,color = cgcolor('gray')
;oplot,pmRA[Nind],pmDEC[Nind],psym = 1,color = cgcolor('red')
;oplot,pmRA[Nind[sele]],pmDEC[Nind[sele]],psym = 1,color = cgcolor('green')

pmn = intarr(24,24) ;儲存每行有的亂數數量

for i = 0,23 do begin
  for k = 0,23 do begin
    pmn[i,k] = n_elements( where( pmRA[Nind] gt (-6+i*0.5) and pmRA[Nind] lt (-6+(i+1)*0.5) and $
       pmDEC[Nind] gt (-6+k*0.5)and pmDEC[Nind] lt (-6+(k+1)*0.5) )) 
  endfor
endfor

xx = indgen(24)/2.-6 & yy = indgen(24)/2.-6

nn = reverse(pmn)

window,7,xsize=900,ysize=900
contour,pmn,xx,yy,xtitle='pmRA',ytitle='pmDEC',xrange=[-6,6],yrange=[-6,6],xticks=6,yticks=6, $
  xstyle=1,ystyle=1,/nodata,xtickFORMAT='(f5.2)',ytickFORMAT='(f5.2)' ;,NLEVELS = 3,,position=[0.13,0.11,0.9,0.56]
oplot,pmRA[Nind],pmDEC[NInd],psym = 1
contour,pmn,xx,yy,NLEVELS=10,color = cgcolor('red'),/overplot
plot,[-6,6],[-6,6],/nodata,/noerase,xticks=9,yticks=10,XTickLen=1.0,YTickLen=1.0, $
  XGridStyle=10,YGridStyle=10,xtickformat="(A1)",ytickformat="(A1)" ;position=[0.13,0.11,0.9,0.56]

WRITE_JPEG,"PM contour.png",TVRD(/TRUE),QUALITY=100, /true
  
window,7,xsize=900,ysize=900
contour,pmn,xx,yy,xtitle='pmRA',ytitle='pmDEC',xrange=[-6,6],yrange=[-6,6],xticks=6,yticks=6, $
  xstyle=1,ystyle=1,/nodata,position=[0.13,0.11,0.9,0.9],xtickFORMAT='(f5.2)',ytickFORMAT='(f5.2)';,NLEVELS = 3 ,xrange=[34,36.25],yrange=[56.5,57.75]
oplot,pmRA[Nind],pmDEC[NInd],psym = 1
tvscale,bytscl(pmn),position=[0.13,0.11,0.9,0.9],MAXVALUE=( mean(bytscl(pmn))+2*stdev(bytscl(pmn)) ),MINVALUE= (min(bytscl(pmn)) )
colorbar,ticknames=['0.00','0.25','0.50','0.75','1.00'],divisions=4,/vertical,position=[0.93,0.1,0.94,0.9],title="Normolization Power",/right
contour,pmn,xx,yy,NLEVELS=7,color = cgcolor('red'),/overplot
plot,[-6,6],[-6,6],/nodata,/noerase,xticks=9,yticks=10,XTickLen=1.0,YTickLen=1.0, $
  XGridStyle=10,YGridStyle=10,xtickformat="(A1)",ytickformat="(A1)" ;position=[0.13,0.11,0.9,0.56]

WRITE_JPEG,"PM flase color map.png",TVRD(/TRUE),QUALITY=100, /true

;;; plot PM vector histogram and fit Gauaaian
;window,6,xsize=900,ysize=900
;plothist,PMV*cos(theta),xhist,yhist,bin=0.1,xticklen=-0.01,/noplot
;plot,xhist,yhist,psym=10,xtitle='Proper motion vector of candidates',ytitle='Count',xticklen=-0.01,xrange=[0,5]
;gau =  GAUSSFIT(xhist, yhist, coeff, NTERMS=6)
;oplot,xhist,gau,color=cgcolor('red')
;oplot,[coeff[1],coeff[1]],!Y.CRANGE,linestyle=2,color=cgcolor('green'),thick=2
;oplot,[coeff[1]-coeff[2],coeff[1]-coeff[2]],!Y.CRANGE,linestyle=2,color=cgcolor('green'),thick=2
;oplot,[coeff[1]+coeff[2],coeff[1]+coeff[2]],!Y.CRANGE,linestyle=2,color=cgcolor('green'),thick=2
;
;sele = where( PMV*cos(theta) gt (coeff[1]-coeff[2]) and PMV*cos(theta) lt (coeff[1]+coeff[2]) )
;
;window,7,xsize=900,ysize=900
;plot,PMV,theta,psym=3,xstyle=4,ystyle=4,xrange=[0,5],yrange=[-5,5],/polar ;
;AXIS, 0, 0, XAX=0
;AXIS, 0, 0, YAX=0
;oplot,PMV[sele],theta[sele],psym=3,color=cgcolor('red'),/polar
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

window,8,xsize=900,ysize=900
PLOT,RANew,DEC, psym=3, xtit='RA difference (deg)', ytit='DEC (deg)',title='NGC869 and NGC884',xrange=[1.0,-1.0], $
  yrange=[56,58],ystyle = 1,xstyle = 1,/isotropic,xticklen=-0.01
TVCIRCLE, 1.5/60, -0.13, 57.13, color=cgcolor('red'),THICK = 2 ;; NGC869
TVCIRCLE, 1.5/60., 0.28, 57.10, color=cgcolor('blue'),THICK = 2  ;;NGC884 
;oplot,RANew[cen869],DEC[cen869],psym = 1,color = cgcolor('red')

WRITE_JPEG,"PM center.png",TVRD(/TRUE),QUALITY=100, /true

rad869 = SQRT( (RANew+0.13)*(RANew+0.13) + (DEC-57.13)*(DEC-57.13) )
cen869 = WHERE(rad869 LT (1.5/60.))

rad884 = SQRT( (RANew-0.28)*(RANew-0.28) + (DEC-57.10)*(DEC-57.10) )
cen884 = WHERE(rad884 LT (1.5/60.))

window,9,xsize=900,ysize=900
plothist,rgeo[cen869],color = cgcolor('red'), bin = 100, xtit = 'Distance (pc)', ytit = 'count', title = 'NGC869',xticklen = -0.01
;plothist,rgeo[cen884],color = cgcolor('blue'), bin = 100, xtit = 'DIstance (pc)', ytit = 'count', title = 'NGC884',xticklen = -0.01 ;/over,

window,10,xsize=900,ysize=900
plot,pmRA,pmDEC,psym = 3,xtitle = 'pmRA (mas/yr)',ytitle = 'pmDEC (mas/yr)',title = 'Proper Motion', $
  xstyle=1,ystyle=1,/nodata ;,xrange=[-5,5],yrange=[-5,5]
oplot,pmRA,pmDEC,psym = 3,color = cgcolor('gray')
oplot,pmRA[cen869],pmDEC[cen869],psym = 1,color = cgcolor('red')
oplot,pmRA[cen884],pmDEC[cen884],psym = 1,color = cgcolor('blue')


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Distance + Proper motion
pmind = SQRT( (pmRA+0.6)*(pmRA+0.6) + (pmDEC+1.15)*(pmDEC+1.15) )
pmin = WHERE(pmind LT 0.5)

;dis1 = 2300 ;2291
;dis2 = 2400 ;2399
;dpin = where(RANew gt (34.-MRA) and RANew lt (36.25-MRA) and DEC gt 56.5 and DEC lt 57.75 and rgeo gt dis1 and rgeo lt dis2 and pmind LT RofPM)

dpin = where(RANew gt (34.-MRA) and RANew lt (36.25-MRA) and DEC gt 56.5 and DEC lt 57.75 and rgeo gt dis1 and rgeo lt dis2 and $
  pmRA gt (coeffRA[1]-coeffRA[2]) and pmRA lt (coeffRA[1]+coeffRA[2]) and $
  pmDec gt (coeffDec[1]-coeffDec[2]) and pmDec lt (coeffDec[1]+coeffDec[2]))

;;; filed star number calculation (35.0, 56.25)
fco869 = SQRT( (RANew-(35-MRA))*(RANew-(35-MRA)) + (DEC-56.25)*(DEC-56.25) )
coin = WHERE(fco869 LT (14.4/60))
fRA869 = RANew[coin]
fDEC869 = DEC[coin]
;fdin869 = where(rgeo gt dis1 and rgeo lt dis2 and pmind LT RofPM) ;pmind LT 0.4

fco884 = SQRT( (RANew-(35-MRA))*(RANew-(35-MRA)) + (DEC-56.25)*(DEC-56.25) )
coin = WHERE(fco884 LT (10.5/60))
fRA884 = RANew[coin]
fDEC884 = DEC[coin]
;fdin884 = where(rgeo gt dis1 and rgeo lt dis2 and pmind LT RofPM) ;pmind LT 0.4

window,9,xsize=900,ysize=900
PLOT,RANew,DEC, psym=3, xtit='RA difference (deg)', ytit='DEC (deg)',title='NGC869 and NGC884',xrange=[1.5,-1.5],yrange=[56.,59.0],ystyle = 1,xstyle = 1,/isotropic ;,xrange=[33.5,36.5],yrange=[56.,59.0]
;oplot,RA[Nind],DEC[Nind],psym = 1,color = cgcolor('red') ; plot cluster member candidates
oplot,RANew[dpin],DEC[dpin],psym = 1,color = cgcolor('red')
;TVCIRCLE, 0.24, ra869, dec869, color=1,THICK = 2
;TVCIRCLE, 0.175, ra884, dec884, color=3,THICK = 2
;TVCIRCLE, 14.4/60., 35, 56.25, color=1,THICK = 2 ;field star region
;TVCIRCLE, 10.5/60., 35, 56.25, color=3,THICK = 2
;oplot,fRA869[fdin869],fDEC869[fdin869],psym = 6,color = cgcolor('red')
;A = FIndGen(5) * (!PI*2/5.)
;UserSym, cos(A), sin(A), /fill
;oplot,fRA884[fdin884],fDEC884[fdin884],psym = 8,color = cgcolor('blue')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
n = intarr(18,10) ;儲存每行有的亂數數量

for i = 0,17 do begin
  for k = 0,9 do begin
    n[i,k] = n_elements( where( RANew[dpin] gt (34-MRA+i*(12.5/100)) and RANew[dpin] lt (34-MRA+(i+1)*(12.5/100.)) and DEC[dpin] gt (56.5+k*(12.5/100.))and DEC[dpin] lt (56.5+(k+1)*(12.5/100.)) ) )
  endfor
endfor

;; plot histogram of star count
window,10,xsize=900,ysize=900
plothist,n, bin=1,xtitle='pixel density',ytitle='Count',xticklen=-0.01,HALFBIN=0
;window,5,xsize=900,ysize=900
;plot,[34.0,36.25],[55.5,57.75],position=[0.13,0.11,0.9,0.56],/nodata,xtickformat="(A1)",ytickformat="(A1)" ;,/ISOTROPIC
;maxx = max(bytscl(n)) ;- 1.1*stdev(bytscl(n))
;minn = min(bytscl(n)) ;+ 1.1*stdev(bytscl(n))
;tvscale,bytscl(n),position=[0.13,0.11,0.9,0.56],MAXVALUE=maxx,MINVALUE=minn
;plot,[34.0,36.25],[55.5,57.75],/nodata,/noerase,position=[0.13,0.11,0.9,0.56],xticks=9,yticks=5,XTickLen=1.0,YTickLen=1.0,xtitle='RA (deg)',ytitle='DEC (deg)', $
;  xtickname = ['34.0','34.25','34.5','34.75','35.0','35.25','35.5','35.75','36.0','36.25'],ytickname = ['56.5','56.75','57.0','57.25','57.5','57.75']
;plot,[34,36.25],[55.5,57.75],/nodata,/noerase,position=[0.13,0.11,0.9,0.56],xticks=9,yticks=5,XTickLen=1.0,YTickLen=1.0, $
;  XGridStyle=10,YGridStyle=10,COLOR=255,xtickformat="(A1)",ytickformat="(A1)"
;colorbar,ticknames=['0.00','0.25','0.50','0.75','1.00'],divisions=4,/vertical,position=[0.93,0.1,0.94,0.56],title="Normolization Power",/right

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; plot contour
xx = indgen(18)/8.+(34.-MRA) & yy = indgen(10)/8.+56.5

nn = reverse(n)

window,11,xsize=900,ysize=900
plot,RANew[dpin],DEC[dpin],psym=1,xrange=[1.25,-1],xstyle=1,ystyle=1,xtit='RA difference (deg)', ytit='DEC (deg)' $
  ,/isotropic,xticklen=-0.01 ;,position=[0.13,0.11,0.9,0.56] ,xticks=4,yticks=10

window,12,xsize=900,ysize=900
contour,nn,xx,yy,xtitle='RA difference (deg)',ytitle='DEC (deg)',xrange=[1.25,-1],yrange=[56.5,57.75],xticks=9,yticks=5, $
  xstyle=1,ystyle=1,position=[0.13,0.11,0.9,0.56],/nodata,xtickFORMAT='(f5.2)' ;,NLEVELS = 3,xrange=[34,36.25],yrange=[56.5,57.75]
oplot,RANew[dpin],DEC[dpin],psym = 1
contour,nn,xx,yy,NLEVELS=7,color = cgcolor('red'),/overplot
plot,[1.25,-1],[55.5,57.75],/nodata,/noerase,position=[0.13,0.11,0.9,0.56],xticks=9,yticks=10,XTickLen=1.0,YTickLen=1.0, $
  XGridStyle=10,YGridStyle=10,xtickformat="(A1)",ytickformat="(A1)"

window,13,xsize=900,ysize=900
contour,nn,xx,yy,xtitle='RA difference (deg)',ytitle='DEC (deg)',xrange=[1.25,-1],yrange=[56.5,57.75],xticks=9,yticks=5, $
  xstyle=1,ystyle=1,position=[0.13,0.11,0.9,0.56],/nodata,xtickFORMAT='(f5.2)' ;,NLEVELS = 3 ,xrange=[34,36.25],yrange=[56.5,57.75]
oplot,RANew[dpin],DEC[dpin],psym = 1
contour,nn,xx,yy,NLEVELS=10,color = cgcolor('red'),/overplot
;plot,[1.25,-1],[55.5,57.75],/nodata,/noerase,position=[0.13,0.11,0.9,0.56],xticks=9,yticks=5,XTickLen=1.0,YTickLen=1.0, $
;  XGridStyle=10,YGridStyle=10,xtickformat="(A1)",ytickformat="(A1)"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; plot contour + false color
window,14,xsize=900,ysize=900
contour,nn,xx,yy,xtitle='RA (deg)',ytitle='DEC (deg)',xrange=[1.25,-1],yrange=[56.5,57.75],xticks=9,yticks=5, $
  xstyle=1,ystyle=1,position=[0.13,0.11,0.9,0.56],/nodata,xtickFORMAT='(f5.2)' ;,NLEVELS = 3 ,xrange=[34,36.25],yrange=[56.5,57.75]
oplot,RANew[dpin],DEC[dpin],psym = 1
tvscale,bytscl(n),position=[0.13,0.11,0.9,0.56],MAXVALUE=( mean(bytscl(nn))+2*stdev(bytscl(nn)) ),MINVALUE= (min(bytscl(nn)) )
colorbar,ticknames=['0.00','0.25','0.50','0.75','1.00'],divisions=4,/vertical,position=[0.93,0.1,0.94,0.56],title="Normolization Power",/right
;contour,nn,xx,yy,NLEVELS=7,color = cgcolor('white'),/overplot;
;plot,[1.25,-1],[56.5,57.75],/nodata,/noerase,position=[0.13,0.11,0.9,0.56],xticks=9,yticks=5,XTickLen=1.0,YTickLen=1.0, $
;    XGridStyle=10,YGridStyle=10,COLOR=255,xtickformat="(A1)",ytickformat="(A1)"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; plot 3D

;p = plot3D(RA[dpin],DEC[dpin],rgeo[dpin], '.',/SYM_FILLED,AXIS_STYLE=2, $
;  xtitle='RA (deg)',ytitle='DEC (deg)',ztitle='distance (pc)')
;ax = p.AXES
;ax[2].HIDE = 1
;ax[6].HIDE = 1
;ax[7].HIDE = 1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; plot all star spatial distribution
;dpinN = where(rgeo gt dis1 and rgeo lt dis2 and pmind LT RofPM) ;pmind LT 0.4
;
;window,15,xsize=900,ysize=900
;PLOT,RANew,DEC,psym=3, xtit='RA difference (deg)',ytit='DEC (deg)',title='NGC869 and NGC884', $
;   xrange=[max(RANew),min(RANew)],ystyle = 1,xstyle = 1,/isotropic ;,xrange=[-1.5,1.5],yrange=[56.,59.0]
;oplot,RANew[dpinN],DEC[dpinN],psym = 1,color = cgcolor('red')
;
;window,16,xsize=900,ysize=900
;plot,RANew[dpinN],DEC[dpinN],psym = 1,xrange=[max(RANew),min(RANew)],ystyle = 1,xstyle = 1,$
;   xtit='RA difference (deg)',ytit='DEC (deg)',/isotropic

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;dis = 1./Plx*1000
;
;sp = sort(dis)
;lo = min(where(dis(sp) gt 0))
;hi = max(where(dis(sp) gt 0))
;sr = sort(rgeo)
;
;xx = findgen(1e5)
;window,17,xsize=900,ysize=900
;plot,xx,dis[sp[lo:hi]],psym = 3, ystyle = 1, xstyle = 1, ytitle='pc',yrange=[100,2e3],xrange=[100,3e5],/ylog,/xlog,/isotropic ;,yrange=[100,2e3],xrange=[100,3e5]
;oplot,xx,rgeo[sr],psym = 3,color=cgcolor('red')



































end