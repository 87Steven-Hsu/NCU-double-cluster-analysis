set_plot,'win'
device,decompose = 0
!p.background=255
!p.color='0'
!p.charsize=1.2
!p.thick=1.3
!p.charthick=1.0
!p.font=-1

r = (RANDOMU(Seed, 30,30)*10) ;100,100表示星團星星(資料點)密度 *1代表星團半徑大小
xx = findgen(n_elements(r))

window,0,xsize = 800,ysize = 800
PLOT,r,xx,/polar,psym=3,xstyle=1,ystyle=1

window,1,xsize = 800,ysize = 800
PLOT,[-100,-50,0,1,2,3,50,100],[-100,-50,0,1,2,3,50,100],xstyle=1,ystyle=1
oplot,r,xx,/polar,psym=3
















end

