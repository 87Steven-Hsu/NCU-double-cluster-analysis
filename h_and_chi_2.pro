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
;; set a criteria of a rectangle range of star candidate
ind = where(RA gt 33.5 and RA lt 36.5 and DEC gt 56.0 and DEC lt 59) ;RA gt 34.0 and RA lt 36 and DEC gt 55.5 and DEC lt 57.5
NRA = RA[ind]
NDEC = DEC[ind]

Nind  = where(NRA gt 34. and NRA lt 36.25 and NDEC gt 56.5 and NDEC lt 57.75)

window,1,xsize=900,ysize=900
PLOT,NRA,NDEC, psym=3, xtit='RA', ytit='DEC',title='NGC869 and NGC884',ystyle = 1,xstyle = 1,/isotropic
oplot,NRA[Nind],NDEC[Nind],psym = 1,color = cgcolor('red') ; plot cluster member candidates

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; coordinate criteria + distance (2291-2399)
dis1 = 2300 ;2291
dis2 = 2400 ;2399
Nind  = where(RA gt 34. and RA lt 36.25 and DEC gt 56.5 and DEC lt 57.75 and rgeo gt dis1 and rgeo lt dis2)

;; Coordinate + distance
window,2,xsize=900,ysize=900
PLOT,RA,DEC, psym=3, xtit='RA', ytit='DEC',title='NGC869 and NGC884',xrange=[33.5,36.5],yrange=[56.,59.0],ystyle = 1,xstyle = 1,/isotropic
oplot,RA[Nind],DEC[Nind],psym = 1,color = cgcolor('red') ; plot cluster member candidates
;TVCIRCLE, 0.24, ra869, dec869, color=1,THICK = 2
;TVCIRCLE, 0.175, ra884, dec884, color=3,THICK = 2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Proper Motion
window,3,xsize=900,ysize=900
plot,pmRA,pmDEC,psym = 3,xtitle = 'pmRA (mas/yr)',ytitle = 'pmDEC (mas/yr)',title = 'Proper Motion',xrange=[-5,5],yrange=[-5,5],/nodata
oplot,pmRA,pmDEC,psym = 3,color = cgcolor('gray')
oplot,pmRA[Nind],pmDEC[Nind],psym = 1,color = cgcolor('red')
RofPM = 0.5
TVCIRCLE, RofPM, -0.7,-1.2, color=3,THICK = 2 ;(-0.6, -1.2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Distance + Proper motion
pmind = SQRT( (pmRA+0.6)*(pmRA+0.6) + (pmDEC+1.15)*(pmDEC+1.15) )
pmin = WHERE(pmind LT 0.5)

;dis1 = 2300 ;2291
;dis2 = 2400 ;2399
dpin = where(RA gt 34. and RA lt 36.25 and DEC gt 56.5 and DEC lt 57.75 and rgeo gt dis1 and rgeo lt dis2 and pmind LT RofPM) ;pmind LT 0.4

;;; filed star number calculation (35.0, 56.25)
fco869 = SQRT( (RA-35)*(RA-35) + (DEC-56.25)*(DEC-56.25) )
coin = WHERE(fco869 LT (14.4/60))
fRA869 = RA[coin]
fDEC869 = DEC[coin]
fdin869 = where(rgeo gt dis1 and rgeo lt dis2 and pmind LT RofPM) ;pmind LT 0.4

fco884 = SQRT( (RA-35)*(RA-35) + (DEC-56.25)*(DEC-56.25) )
coin = WHERE(fco884 LT (10.5/60))
fRA884 = RA[coin]
fDEC884 = DEC[coin]
fdin884 = where(rgeo gt dis1 and rgeo lt dis2 and pmind LT RofPM) ;pmind LT 0.4

window,4,xsize=900,ysize=900
PLOT,RA,DEC, psym=3, xtit='RA', ytit='DEC',title='NGC869 and NGC884',xrange=[33.5,36.5],yrange=[56.,59.0],ystyle = 1,xstyle = 1,/isotropic ;,xrange=[33.5,36.5],yrange=[56.,59.0]
;oplot,RA[Nind],DEC[Nind],psym = 1,color = cgcolor('red') ; plot cluster member candidates
oplot,RA[dpin],DEC[dpin],psym = 1,color = cgcolor('red')
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
    n[i,k] = n_elements( where( RA[dpin] gt (34+i*(12.5/100)) and RA[dpin] lt (34+(i+1)*(12.5/100.)) and DEC[dpin] gt (56.5+k*(12.5/100.))and DEC[dpin] lt (56.5+(k+1)*(12.5/100.)) ) )
  endfor
endfor

xx = indgen(18)/8.+34. & yy = indgen(10)/8.+56.5

device,decompose = 0
!p.background=255
loadct,3
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
xx = indgen(18)/8.+34. & yy = indgen(10)/8.+56.5

window,5,xsize=900,ysize=900
plot,RA[dpin],DEC[dpin],psym = 1,xticks=18,yticks=10,xstyle=1,ystyle=1,position=[0.13,0.11,0.9,0.56]

window,6,xsize=900,ysize=900
contour,n,xx,yy,xtitle='RA (deg)',ytitle='DEC (deg)',xrange=[34,36.25],yrange=[56.5,57.75],xticks=18,yticks=10, $
  xstyle=1,ystyle=1,position=[0.13,0.11,0.9,0.56],/nodata ;,NLEVELS = 3
oplot,RA[dpin],DEC[dpin],psym = 1
contour,n,xx,yy,NLEVELS=5,color = cgcolor('red'),/overplot
;plot,[34,36.25],[55.5,57.75],/nodata,/noerase,position=[0.13,0.11,0.9,0.56],xticks=18,yticks=10,XTickLen=1.0,YTickLen=1.0, $
;  XGridStyle=10,YGridStyle=10,xtickformat="(A1)",ytickformat="(A1)"

window,7,xsize=900,ysize=900
contour,n,xx,yy,xtitle='RA (deg)',ytitle='DEC (deg)',xrange=[34,36.25],yrange=[56.5,57.75],xticks=9,yticks=5, $
  xstyle=1,ystyle=1,position=[0.13,0.11,0.9,0.56],/nodata ;,NLEVELS = 3
oplot,RA[dpin],DEC[dpin],psym = 1
contour,n,xx,yy,NLEVELS=10,color = cgcolor('red'),/overplot
;plot,[34,36.25],[55.5,57.75],/nodata,/noerase,position=[0.13,0.11,0.9,0.56],xticks=9,yticks=5,XTickLen=1.0,YTickLen=1.0, $
;  XGridStyle=10,YGridStyle=10,xtickformat="(A1)",ytickformat="(A1)"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; plot contour + false color
window,10,xsize=900,ysize=900
contour,n,xx,yy,NLEVELS=7,xrange = [34,36.25],yrange=[56.5,57.75],xtitle='RA (deg)',ytitle='DEC (deg)', $
  xstyle=1,ystyle=1,position=[0.13,0.11,0.9,0.56],/nodata ;,NLEVELS = 3
oplot,RA[dpin],DEC[dpin],psym = 1
;; option 2: star candidates + contour + false color
;; ( tvscale, contour) & plot must run separatelty 
tvscale,bytscl(n),position=[0.13,0.11,0.9,0.56],MAXVALUE=( mean(bytscl(n))+2*stdev(bytscl(n)) ),MINVALUE= (min(bytscl(n)) ) ;-3*stdev(bytscl(n))
colorbar,ticknames=['0.00','0.25','0.50','0.75','1.00'],divisions=4,/vertical,position=[0.93,0.1,0.94,0.56],title="Normolization Power",/right
;contour,n,xx,yy,NLEVELS=7,color = cgcolor('white'),/overplot;
plot,[34,36.25],[55.5,57.75],/nodata,/noerase,position=[0.13,0.11,0.9,0.56],xticks=9,yticks=5,XTickLen=1.0,YTickLen=1.0, $
    XGridStyle=10,YGridStyle=10,COLOR=255,xtickformat="(A1)",ytickformat="(A1)"
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; plot 3D

;p = plot3D(RA[dpin],DEC[dpin],rgeo[dpin], '.',/SYM_FILLED,AXIS_STYLE=2, $
;  xtitle='RA (deg)',ytitle='DEC (deg)',ztitle='distance (pc)')
;ax = p.AXES
;ax[2].HIDE = 1
;ax[6].HIDE = 1
;ax[7].HIDE = 1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; plot histogram of star count
window,11,xsize=900,ysize=900
plothist,n, bin=1,xtitle='pixel density',ytitle='Count',xticklen=-0.01,HALFBIN=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Gaussian fit of proper motion
window,12,xsize=900,ysize=900
plothist,pmRA[Nind],xhist,yhist,bin=0.1,xticklen=-0.01,/noplot
plot,xhist,yhist,psym=10,xtitle='pmRA of candidates',ytitle='Count',xticklen=-0.01,xrange=[-2,0.6]
gau =  GAUSSFIT(xhist, yhist, coeff, NTERMS=6)
oplot,xhist,gau,color=cgcolor('red')

window,13,xsize=900,ysize=900
plothist,pmDEC[Nind],xhist,yhist,bin=0.1,xticklen=-0.01,/noplot
plot,xhist,yhist,psym=10,xtitle='pmDEC of candidates',ytitle='Count',xticklen=-0.01,xrange=[-2,0.6]
gau =  GAUSSFIT(xhist, yhist, coeff, NTERMS=6)
oplot,xhist,gau,color=cgcolor('red')

PMV = sqrt(pmRA[Nind]*pmRA[Nind] + pmDEC[Nind]*pmDEC[Nind])
theta = atan(pmDEC[Nind]/pmRA[Nind]) ;/!dtor => transfer to deg

;; plot PM vertor on polar coordinate

;pol = POLARPLOT(PMV,theta,SYMBOL=3, LINESTYLE=6)

;; plot PM vector histogram and fit Gauaaian
window,14,xsize=900,ysize=900
plothist,PMV*cos(theta),xhist,yhist,bin=0.1,xticklen=-0.01,/noplot
plot,xhist,yhist,psym=10,xtitle='Proper motion vector of candidates',ytitle='Count',xticklen=-0.01,xrange=[0,5]
gau =  GAUSSFIT(xhist, yhist, coeff, NTERMS=6)
oplot,xhist,gau,color=cgcolor('red')
oplot,[coeff[1],coeff[1]],!Y.CRANGE,linestyle=2,color=cgcolor('green'),thick=2
oplot,[coeff[1]-coeff[2],coeff[1]-coeff[2]],!Y.CRANGE,linestyle=2,color=cgcolor('green'),thick=2
oplot,[coeff[1]+coeff[2],coeff[1]+coeff[2]],!Y.CRANGE,linestyle=2,color=cgcolor('green'),thick=2

sele = where( PMV*cos(theta) gt (coeff[1]-coeff[2]) and PMV*cos(theta) lt (coeff[1]+coeff[2]) )

window,15,xsize=900,ysize=900
plot,PMV,theta,psym=3,xstyle=4,ystyle=4,xrange=[0,5],yrange=[-5,5],/polar ;
AXIS, 0, 0, XAX=0
AXIS, 0, 0, YAX=0
oplot,PMV[sele],theta[sele],psym=3,color=cgcolor('red')


































end