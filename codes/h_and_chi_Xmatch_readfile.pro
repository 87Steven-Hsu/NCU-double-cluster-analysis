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
test = read_csv('D:\學校資料\天文所陳老師\h and x\h and chi cross match\h and chi eDR3 Apr-1.csv')

;Transfer CSV transfer to sav.
SAVE, test, FILE='D:\學校資料\天文所陳老師\h and x\h and chi cross match\h and chi eDR3 Apr-1.dat'
restore,'D:\學校資料\天文所陳老師\h and x\h and chi cross match\h and chi eDR3 Apr-1.dat'
;
;;;Transfer structure to array
angdist = [float(test.FIELD01)] ;
RA = [float(test.FIELD02)]      ;J2000
DEC = [float(test.FIELD03)]     ;J2000
Plx = [float(test.FIELD12)]     ;
pmotion = [float(test.FIELD15)]      ;
pmRA = [float(test.FIELD16)]    ;
pmDec = [float(test.FIELD18)]   ;
Gmag = [float(test.FIELD33)]    ;
BPmag = [float(test.FIELD36)]   ;
RPmag = [float(test.FIELD38)]   ;
rego = [float(test.FIELD63)]    ;Median of the geometric distance
rego16 = [float(test.FIELD64)]  ;16th percentile of the geometric distance
rego84 = [float(test.FIELD65)]  ;84th percentile of the geometric distance
rpgeo = [float(test.FIELD66)]   ;Median of the photogeometric distance
rpgeo16 = [float(test.FIELD67)] ;16th percentile of the photogeometric distance
rpgeo84 = [float(test.FIELD68)] ;84th percentile of the photogeometric distance


data = TRANSPOSE([ [angdist], [RA], [DEC], [PLx], [pmotion], [pmRA], [pmDec], [Gmag], [BPmag], [RPmag], [rego], [rego16], [rego84], [rpgeo], [rpgeo16], [rpgeo84]])
WRITE_CSV, 'h and chi Xmatch.csv', data





end
