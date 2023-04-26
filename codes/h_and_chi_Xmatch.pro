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
;test = read_csv('D:\學校資料\天文所陳老師\h and x\h and chi cross match\h and chi Xmatch.csv')
;
;;Transfer CSV transfer to sav.
;SAVE, test, FILE='D:\學校資料\天文所陳老師\h and x\h and chi cross match\h and chi Xmatch.dat'
;restore,'D:\學校資料\天文所陳老師\h and x\h and chi cross match\h and chi Xmatch.dat'
;;
;;;;Transfer structure to array
;angdist = [float(test.FIELD01)] ;
;RA = [float(test.FIELD02)]      ;J2000
;DEC = [float(test.FIELD03)]     ;J2000
;Plx = [float(test.FIELD04)]     ;
;pmotion = [float(test.FIELD05)]      ;
;pmRA = [float(test.FIELD06)]    ;
;pmDec = [float(test.FIELD07)]   ;
;Gmag = [float(test.FIELD08)]    ;
;BPmag = [float(test.FIELD09)]   ;
;RPmag = [float(test.FIELD10)]   ;
;rego = [float(test.FIELD11)]    ;Median of the geometric distance
;rego16 = [float(test.FIELD12)]  ;16th percentile of the geometric distance
;rego84 = [float(test.FIELD13)]  ;84th percentile of the geometric distance
;rpgeo = [float(test.FIELD14)]   ;Median of the photogeometric distance
;rpgeo16 = [float(test.FIELD15)] ;16th percentile of the photogeometric distance
;rpgeo84 = [float(test.FIELD16)] ;84th percentile of the photogeometric distance

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

RAindex = where(RA GT 34.25 and RA LT 36.25)
DECindex = where(DEC GT 56 and DEC LT 58)
NRA = RA(RAindex)
NDEC = DEC(DECindex)
;; NGC869 = h
radcen869 = 14.4/60. ;~0.24 deg (radius)
ra869 = 34.75
dec869 = 57.13

;; NGC884 = x
radcen884 = 10.5/60. ;~0.175 deg (radius)
ra884 = 35.60
dec884 = 57.13

;; circle out a region of cluster core
window,1,xsize=900,ysize=900
PLOT, NRA, NDEC, psym=3, xtit='RA', ytit='DEC',title='NGC869 and NGC884, GAIA eDR3, N=283173' ;,/isotropic,yrange=[56,58],xrange=[34,36.5]
TVCIRCLE, 0.24, ra869, dec869, color=1,THICK = 2
TVCIRCLE, 0.175, ra884, dec884, color=3,THICK = 2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Proper motion 
window,2,xsize=900,ysize=900
PLOT, pmRA, pmDEC, psym=3, xtit='RA', ytit='DEC',title='NGC869 and NGC884',yrange = [55,60]

window,3,xsize=900,ysize=900
PLOT, RA, DEC, psym=3, xtit='RA', ytit='DEC',title='NGC869 and NGC884',/isotropic,yrange=[56,58],xrange=[34,36.5]
TVCIRCLE, 0.24, ra869, dec869, color=1,THICK = 2
TVCIRCLE, 0.175, ra884, dec884, color=3,THICK = 2



end
