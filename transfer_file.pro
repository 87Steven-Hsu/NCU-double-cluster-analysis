openr,1,'D:\學校資料\天文所陳老師\GAIA\asu.txt'

a = 1. ;RA (deg)
b = 1. ;RAerror (deg in mili arcsec)
cc = 1. ;DEC (deg)
d = 1. ;DEC error (deg in mili arcsec)
e = 1. ;PLX, Absolute steeGmagar paraeGmagax (paraeGmagax)
f = 1. ;PLX error
g = 1. ;pmRA (mas/yr), Proper motion in right ascension
h = 1. ;pmRA error (mas/yr)
i = 1. ;pmDEC (mas/yr), Proper motion in declination
j = 1. ;pmDEC error (mas/yr)
k = 1. ;Gmag (mag), Gband mean magnitude
l = 1. ;Gmag error (mag)
m = 1. ;BPmag (mag), Integrated BP mean magnitude
n = 1. ;BPmag error
o = 1. ;RPmag (mag), Integrated RP mean magnitude
p = 1. ;RPmag error
q = 1. ;RV (km/sec), Spectroscopic radial velocity in the solar barycentric reference frame
r = 1. ;RV error

emp = ' '

RA = findgen(152000) & eRA = findgen(152000) & DEC = findgen(152000) & eDEC = findgen(152000) & PLX = findgen(152000) & ePLX = findgen(152000) & pmRA = findgen(152000)
epmRA = findgen(152000) & pmDEC = findgen(152000) & epmDEC = findgen(152000) & Gmag = findgen(152000) & eGmag = findgen(152000) & BPmag = findgen(152000) & eBPmag = findgen(152000)
RPmag = findgen(152000) & eRPmag = findgen(152000) & RV = findgen(152000) & eRV = findgen(152000)

for y = 1,74 do readf,1,emp
                                                    ;,f,g,h,i,j,k,l,m,n,o,p,q,r
for z = 0,99999 do begin                                   ; a     b          c       d       e       f       g 
  readf,1,a,b,cc,d,e,format = '(f15.11,1x,f7.4,1x,f15.11,4x,f7.4,1x,f7.4) 
  ;,2x,f7.4,4x,f9.3,2x,f6.3,4x,f9.3,2x,f6.3,1x,f7.4,1x,f6.4,1x,f7.4,1x,f6.4,1x,f7.4,1x,f6.4,1x,f7.2,1x,f5.2)' ;,3x,f6.2,2x,f4.2

  if a eq !values.f_nan then a = !values.f_nan &　if b eq !values.f_nan then b = !values.f_nan　&　if cc eq !values.f_nan then cc = !values.f_nan　& if d eq !values.f_nan then d = !values.f_nan ;
  if e eq !values.f_nan then e = !values.f_nan & if f eq !values.f_nan then f = !values.f_nan & if g eq !values.f_nan then g = !values.f_nan & if h eq !values.f_nan then h = !values.f_nan
  if i eq !values.f_nan then i = !values.f_nan &　if j eq !values.f_nan then j = !values.f_nan　&　if k eq !values.f_nan then k = !values.f_nan　&　if l eq !values.f_nan then l = !values.f_nan
  if m eq !values.f_nan then m = !values.f_nan & if n eq !values.f_nan then n = !values.f_nan & if o eq !values.f_nan then o = !values.f_nan & if p eq !values.f_nan then p = !values.f_nan
  ;if q eq !values.f_nan then q = !values.f_nan &　if r eq !values.f_nan then r = !values.f_nan

  RA[z] = a & eRA[z] = b & DEC[z] = cc & eDEC[z] = d &　PLX[z] = e & ePLX[z] = f & pmRA[z] = g & epmRA[z] = h & pmDEC[z] = i ;
  epmDEC[z] = j & Gmag[z] = k & eGmag[z] = l & BPmag[z] = m &　eBPmag[z] = n & RPmag[z] = o & eRPmag[z] = p ; RV[z] = q & eRV[z] = r

  if  EOF(1) eq 1 then break
endfor

RA = RA[0:z] & eRA = eRA[0:z] & DEC = DEC[0:z] & eDEC = eDEC[0:z] &　PLZ = PLZ[0:z] & ePLZ = ePLZ[0:z] & pmRA = pmRA[0:z] & epmRA = epmRA[0:z]  
pmDEC = pmDEC[0:z] & epmDEC = epmDEC[0:z] & Gmag = Gmag[0:z] & eGmag = eGmag[0:z] & BPmag = BPmag[0:z] &　eBPmag = eBPmag[0:z] & RPmag = RPmag[0:z]
eRPmag = eRPmag[0:z] ; RV = RV[0:z] & eRV = eRV[0:z]

SAVE,RA,eRA,DEC,eDEC,PLX,ePLX,pmRA,epmRA,pmDEC,epmDEC,Gmag,eGmag,BPmag,eBPmag,RPmag,eRPmag,FILENAME = 'D:\學校資料\天文所陳老師\GAIA\asu.sav' ;,RV,eRV

close,1

end