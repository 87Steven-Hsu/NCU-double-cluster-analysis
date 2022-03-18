;!p.background=255
;!p.color=0
;!p.charsize=1.2
;!p.thick=1.3
;!p.charthick=1.0
;device,decompose = 0
;!p.font=-1
;red = [255,0,0]
;green = [0,255,0]
;blue = [0,0,255]
;TVLCT,red,green,blue,1
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;讀取資料
;test = read_csv('D:\學校資料\天文所陳老師\Double cluster\data\ASCC_double_cluster.csv')
;
;;Transfer CSV transfer to sav.
;SAVE, test, FILE='D:\學校資料\天文所陳老師\Double cluster\data\ASCC_double_cluster.dat'
;restore,'D:\學校資料\天文所陳老師\Double cluster\data\ASCC_double_cluster.dat'
;
;;Transfer structure to array
;RA = [float(test.FIELD01)]    ;RA (deg)
;DEC = [float(test.FIELD02)]   ;DEC (deg)
;RA_ICRS = [float(test.FIELD03)]   ;
;eRA_ICRS = [float(test.FIELD04)]  ;
;DEC_ICRS = [float(test.FIELD05)]   ;
;eDEC_ICRS = [float(test.FIELD06)]  ;
;Plx = [float(test.FIELD07)]  ;PLX, Absolute steeGmagar paraeGmagax (paraeGmagax)
;ePlx = [float(test.FIELD08)] ;PLX error
;pmRA = [float(test.FIELD09)] ;pmRA (mas/yr), Proper motion in right ascension
;epmRA = [float(test.FIELD10)];pmRA error (mas/yr)
;pmDEC = [float(test.FIELD11)]  ;pmDEC (mas/yr), Proper motion in declination
;epmDEC = [float(test.FIELD12)] ;pmDEC error (mas/yr)
;Gmag = [float(test.FIELD13)] ;Gmag (mag), Gband mean magnitude = brightness distribution
;eGmag = [float(test.FIELD14)];Gmag error (mag)
;BPmag = [float(test.FIELD15)] ;BPmag (mag), Integrated BP mean magnitude
;eBPmag = [float(test.FIELD16)];BPmag error
;RPmag = [float(test.FIELD17)]    ;RPmag (mag), Integrated RP mean magnitude
;eRPmag = [float(test.FIELD18)]   ;RPmag error
;
;
;
n = intarr(150,150) ;儲存每行有的亂數數量

for i = 0,149 do begin  
  for k = 0,149 do begin 
    n[i,k] = n_elements( where(RA gt (min(RA)+i*1/30.) and RA lt (min(RA)+(i+1)*1/30.)  and DEC gt (min(DEC)+k*1/30.)and DEC lt (min(DEC)+(k+1)*1/30.) )) 
  endfor
endfor

device,decompose = 0
!p.background=255
loadct,3
window,11,xsize=900,ysize=900
plot,[0,1],[0,1],/noerase,position=[0.05,0.1,0.89,0.95],xticks=10,yticks=10,xsty=1,ystyle=1,XGridStyle=10,YGridStyle=10,/ISOTROPIC 
;tvimage,bytscl(nfn),position=[0.05,0.1,0.92,0.95]
plot,!x.crange,[0,1],/nodata,xtickformat="(A1)",ytickformat="(A1)",xticks=10,yticks=10,xsty=1,ystyle=1,XGridStyle=10,YGridStyle=10,/ISOTROPIC 
;nfn是矩陣資料點正確數量, 但要用bytscl(fn), 顛倒的矩陣畫tvimage才可以跟nfn的數量對應 (原因未知 )
tvimage,bytscl(n)
colorbar,maxrange=1.,minrange=0.,divisions=4,ticknames=['0.00','0.25','0.50','0.75','1.00'],/vertical,position=[0.93,0.1,0.94,0.95],title="Normolization Power",/right


end