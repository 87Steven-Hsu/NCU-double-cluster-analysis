set_plot,'win'
device,decompose = 0
!p.background=255
!p.color='0'
!p.charsize=1.2
!p.thick=1.3
!p.charthick=1.0
!p.font=-1

r = (RANDOMU(Seed, 100,100)*10000)  ;uniform distribution
;r = (RANDOMN(Seed, 10000)*10000)  ;gaussian distribution
xx = indgen(10000)

;;畫亂數
window,0
plot,xx,r,psym=3,xticks=10,yticks=10,XTickLen=1.0,YTickLen=1.0,XGridStyle=10,YGridStyle=10,position=[0.05,0.1,0.96,0.98] $
 ,xtickname=[0,10,20,30,40,50,60,70,80,90,100],ytickname=[0,10,20,30,40,50,60,70,80,90,100]

n = intarr(10,10) ;儲存每行有的亂數數量
;;以亂數下標判斷每行有的亂數合
for i = 0,9 do begin  ;判斷x
   for k = 0,9 do begin ;判斷y
      n[i,k] = n_elements(where(r[i*1000:(i+1)*1000-1] lt (k+1)*1000))
   endfor
endfor

fn = intarr(10,10) ;儲存每格有的亂數數量, 儲存格子上下顛倒
;;判斷每格有的亂數數量
for i = 0,9 do begin  ;判斷y
   for k = 1,9 do begin ;下一個-前一個即為此格有的數量
    fn[i,k] = n[i,k] - n[i,k-1]
   endfor
endfor

fn[*,0] = n[*,0] ;將n矩陣的第0列儲存至fn矩陣的第0列
nfn = intarr(10,10) ;正確儲存每格有的亂數數量, 沒有上下顛倒
for i = 0,9 do begin ;將n矩陣上下反轉並儲存
  nfn[*,i] = fn[*,-(i+1)]
endfor

print,fn," ",nfn," " ;,bytscl(nfn) ;n," ",fn,
print,'max =',max(nfn) & print,'min =',min(nfn) & print,'STDEV =',stddev(nfn)

device,decompose = 0
!p.background=255
loadct,3

;;假色圖
window,2
plot,[0,1],[0,1],/noerase,position=[0.05,0.1,0.92,0.95],xticks=10,yticks=10,xsty=1,ystyle=1,XGridStyle=10,YGridStyle=10
;tvimage,bytscl(nfn),position=[0.05,0.1,0.92,0.95]
plot,!x.crange,[0,1],/nodata,position=[0.05,0.1,0.92,0.95],xtickformat="(A1)",ytickformat="(A1)",xticks=10,yticks=10,xsty=1,ystyle=1,XGridStyle=10,YGridStyle=10
;nfn是矩陣資料點正確數量, 但要用bytscl(fn), 顛倒的矩陣畫tvimage才可以跟nfn的數量對應 (原因未知 )
tvimage,bytscl(fn),position=[0.05,0.1,0.92,0.95]
plot,[0.5,10.],[0.5,10.],/nodata,/noerase,position=[0.05,0.1,0.92,0.95],xticks=10,yticks=10,XTickLen=1.0,YTickLen=1.0,XGridStyle=10,YGridStyle=10,background=255
colorbar,maxrange=1.,minrange=0.,divisions=4,ticknames=['0.00','0.25','0.50','0.75','1.00'],/vertical,position=[0.93,0.1,0.94,0.95],title="Normolization Power",/right

;;勿用
;--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
;;;contour
;normolize = fltarr(n_elements(nfn))
;normolize = nfn/float(max(nfn))
;
;x = indgen(10) & y = indgen(10)
;;!p.multi = [0,2,1]
;window,1
;;tv,normolize ;,xsty=1,ysty=1,xticklen=0.04,xminor=10,/noerase
;contour,normolize,x,y,xrange=[0,9],yrange=[0,9],xticks=9,yticks=9,/fill,position=[0.05,0.1,0.92,0.95],XTickLen=1.0,YTickLen=1.0,XGridStyle=10,YGridStyle=10
;colorbar,maxrange=1.,minrange=0,divisions=4,/vertical,ticknames=['0.00','0.25','0.50','0.75','1.00'],position=[0.93,0.1,0.95,0.95],title="Normolization Power",/right
;--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

;;;加入星團
;star = indgen(20)*5+5000 & yy = intarr(8)+5000  ;長方形星團
;window,3
;plot,xx,r,psym=3,xticks=10,yticks=10,XTickLen=1.0,YTickLen=1.0,XGridStyle=10,YGridStyle=10,position=[0.05,0.1,0.96,0.98],$ ;畫亂數
;  xtickname=[0,10,20,30,40,50,60,70,80,90,100],ytickname=[0,10,20,30,40,50,60,70,80,90,100]
;for i = 1,8 do begin ;將星團畫入亂數圖中
;  yyy = yy+(i+1)*10
;  oplot,star,yyy,psym=3
;  ;print,yyy
;  
;endfor
;
;fn[5,6] = fn[5,6]+64
;
;nfn = intarr(10,10)
;for i = 0,9 do begin
;  nfn[*,i] = fn[*,-(i+1)]
;endfor
;
;print,fn," ",nfn
;print,'max =',max(nfn) & print,'min =',min(nfn) & print,'STDEV =',stddev(nfn)
;
;device,decompose = 0
;!p.background=255
;loadct,3
;
;;;假色圖
;window,4
;plot,[0,1],[0,1],/noerase,position=[0.05,0.1,0.92,0.95],xticks=10,yticks=10,xsty=1,ystyle=1,XGridStyle=10,YGridStyle=10
;;tvimage,bytscl(nfn),position=[0.05,0.1,0.92,0.95]
;plot,!x.crange,[0,1],/nodata,position=[0.05,0.1,0.92,0.95],xtickformat="(A1)",ytickformat="(A1)",xticks=10,yticks=10,xsty=1,ystyle=1,XGridStyle=10,YGridStyle=10
;tvimage,bytscl(fn),position=[0.05,0.1,0.92,0.95]
;plot,[0.5,10.],[0.5,10.],/nodata,/noerase,position=[0.05,0.1,0.92,0.95],xticks=10,yticks=10,XTickLen=1.0,YTickLen=1.0,XGridStyle=10,YGridStyle=10,background=255
;colorbar,maxrange=1.,minrange=0.,divisions=4,ticknames=['0.00','0.25','0.50','0.75','1.00'],/vertical,position=[0.93,0.1,0.94,0.95],title="Normolization Power",/right

;加入球狀星團
star = (RANDOMU(Seed, 20,20)*200) ;(Seed, 30,30)表示星團星星(資料點)密度   ;;*1代表星團半徑大小
nums = findgen(n_elements(star))

window,5,xsize = 800,ysize = 800
plot,star,nums,psym=3,/polar

window,6,xsize = 800,ysize = 800
plot,xx,r,psym=3,xticks=10,yticks=10,XTickLen=1.0,YTickLen=1.0,XGridStyle=10,YGridStyle=10,position=[0.05,0.1,0.96,0.98],$
  xtickname=[0,10,20,30,40,50,60,70,80,90.0,100],ytickname=[0,10,20,30,40,50,60,70,80,90,100]
oplot,star,nums,/polar,psym=3

number = where(star )


end