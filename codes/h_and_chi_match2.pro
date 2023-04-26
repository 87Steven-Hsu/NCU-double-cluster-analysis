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

;;;讀取資料 GAIA eDR3
;test = read_csv('D:\學校資料\天文所陳老師\h and x\h and chi cross match Apr-1\h and chi eDR3 Apr-1.csv')
;
;;Transfer CSV transfer to sav.
;SAVE, test, FILE='D:\學校資料\天文所陳老師\h and x\h and chi cross match Apr-1\h and chi eDR3 Apr-1.dat'
;restore,'D:\學校資料\天文所陳老師\h and x\h and chi cross match Apr-1\h and chi eDR3 Apr-1.dat'
;;
;;;;Transfer structure to array
;RA_J2000= [float(test.FIELD01)] ;
;DEC_J2000 = [float(test.FIELD02)]      ;
;RA_ICRS = [float(test.FIELD03)]     ;
;eRA_ICRS = [float(test.FIELD04)]     ;
;DEC_ICRS = [float(test.FIELD05)]      ;
;DEC_ICRS = [float(test.FIELD06)]    ;
;Source = [float(test.FIELD07)]   ;
;Plx = [float(test.FIELD08)]    ;
;ePlx = [float(test.FIELD09)]   ;
;pm = [float(test.FIELD10)]   ;
;pmRA = [float(test.FIELD11)]    ;
;epmRA = [float(test.FIELD12)]  ;
;pmDEC = [float(test.FIELD13)]  ;
;epmDEC = [float(test.FIELD14)]   ;
;Gmag = [float(test.FIELD15)] ;
;eGmag = [float(test.FIELD16)] ;
;BPmag = [float(test.FIELD17)] ; 
;eBPmag = [float(test.FIELD18)] ; 
;RPmag = [float(test.FIELD19)] ;
;eRPmag = [float(test.FIELD20)] ;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;讀取資料 Bailer- Jones
;test = read_csv('D:\學校資料\天文所陳老師\h and x\h and chi cross match Apr-1\h and chi B-J Apr-1.csv')
;
;;Transfer CSV transfer to sav.
;SAVE, test, FILE='D:\學校資料\天文所陳老師\h and x\h and chi cross match Apr-1\h and chi B-J Apr-1.dat'
;restore,'D:\學校資料\天文所陳老師\h and x\h and chi cross match Apr-1\h and chi B-J Apr-1.dat'
;
;RAJ2000BJ = [float(test.FIELD01)]      ;
;DECJ2000BJ = [float(test.FIELD02)]     ;
;;SourceBJ = [float(test.FIELD03)] ;
;RA_ICRSBJ = [float(test.FIELD04)]      ;
;DEC_ICRSBJ = [float(test.FIELD05)]     ;
;rgeo = [float(test.FIELD06)]    ;Median of the geometric distance
;rgeo16 = [float(test.FIELD07)]  ;16th percentile of the geometric distance
;rgeo84 = [float(test.FIELD08)]  ;84th percentile of the geometric distance
;rpgeo = [float(test.FIELD09)]   ;Median of the photogeometric distance
;rpgeo16 = [float(test.FIELD10)] ;16th percentile of the photogeometric distance
;rpgeo84 = [float(test.FIELD11)] ;84th percentile of the photogeometric distance

matchres = match_2d(RA_J2000, DEC_J2000, RAJ2000BJ, DECJ2000BJ, 5./3600.,MATCH_DISTANCE=md)

eqneg = where(matchres eq -1)

rgeo(eqneg) = !VALUES.F_NAN
rgeo16(eqneg) = !VALUES.F_NAN
rgeo84(eqneg) = !VALUES.F_NAN
rpgeo(eqneg) = !VALUES.F_NAN
rpgeo16(eqneg) = !VALUES.F_NAN
rpgeo84(eqneg) = !VALUES.F_NAN

;matchres(eqneg) = !VALUES.F_NAN

index = intarr(560569)

for i = 0,560568 do begin
  for j = i+1,560568 do begin
    if (matchres[i] eq matchres[j]) then begin
        index[j] = j
    endif
  endfor
endfor

nezero = where(index ne 0)


;data = TRANSPOSE([ [RA_J2000], [DEC_J2000], [Plx], [ePLx], [pm], [pmRA], [epmRA], [pmDec], [epmDEC], [Gmag], [BPmag], $
;  [RPmag], [rgeo], [rgeo16], [rgeo84], [rpgeo], [rpgeo16], [rpgeo84]])
;WRITE_CSV, 'h and chi Xmatch Apr-1.csv', data








end