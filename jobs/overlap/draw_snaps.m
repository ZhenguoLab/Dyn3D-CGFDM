clc;
clear;
% close all
addpath("../../matlab/");
parfile = 'params.json';
par = get_params(parfile);
NY = par.NY;
NZ = par.NZ;
DT = par.DT;
TMAX = par.TMAX;
TSKP = par.EXPORT_TIME_SKIP;
% TSKP = 1;
nfault = par.num_fault;
j1 = par.Fault_grid(1:4:end-3);
j2 = par.Fault_grid(2:4:end-2);
k1 = par.Fault_grid(3:4:end-1);
k2 = par.Fault_grid(4:4:end);
NT = floor(TMAX/(DT*TSKP));
out = par.OUT;
% out = 'outC0.48'
its = 100:5:NT;
its = floor(its);
nt = length(its);

[x,y,z] = gather_coord(parfile);
x = x*1e-3;
y = y*1e-3;
z = z*1e-3;
figure;
filename1='Tn_40_E.gif';

tn0 = gather_snap(parfile,out,'tn',1);
ts1 = gather_snap(parfile,out,'ts1',1);
ts2 = gather_snap(parfile,out,'ts2',1);
ts0 = sqrt(ts1.^2+ts2.^2);
% v = Vs1;
% ts2 = gather_snap(parfile,out,'ts2',1);
% ts0 = sqrt(ts1.^2+ts2.^2);
for i = 1:nt 
it = its(i);
disp(it);

% out = './compareE';
% out = 'outC0.45'
% Vs1 = gather_snap(parfile,out,'tn',it);

Vs1 = gather_snap(parfile,out,'Vs1',it);
% v = Vs1/1e6;
Vs2 = gather_snap(parfile,out,'Vs2',it);
v = sqrt(Vs1.^2+Vs2.^2);
% v = (v - ts0)/1e6;
% v = (Vs1 - tn0)/1e6;
v = log10(v + 1e-90);
% v = v ./ tn;
% v = Vs1;
% if flag_image
% pcolor(y(:,:,1),z(:,:,1),Vs1(:,:,1));
for n = 1 : nfault
    x1 = squeeze(x(:,:, n));
    y1 = squeeze(y(:,:, n));
    z1 = squeeze(z(:,:, n));
%     surf(x1(j1(n)-30:j2(n)+30,k1(n)-30:k2(n)), y1(j1(n)-30:j2(n)+30,k1(n)-30:k2(n)), z1(j1(n)-30:j2(n)+30,k1(n)-30:k2(n)), squeeze(v(j1(n)-30:j2(n)+30,k1(n)-30:k2(n), n)));
    surf(x1(j1(n):j2(n),k1(n):k2(n)), y1(j1(n):j2(n),k1(n):k2(n)), z1(j1(n):j2(n),k1(n):k2(n)), squeeze(v(j1(n):j2(n),k1(n):k2(n), n)));
%     surf(x1, y1, z1, squeeze(v(:,:, n)));
    hold on;
end
hold off;
view([-60, 30]); 
axis equal; 
shading interp;
% ylim([40, 43]);
% zlim([-1,1]);
% ylim([-16, -14]);

caxis([-16 1.5]);
% caxis([-20, 20]);
colormap( 'jet' );
colorbar;
title(['Rate t=', num2str(it*DT*TSKP),'s'],'FontSize',12)
set(gca,'FontSize',12)
%axis([-1 1 -1 0]*15)
pause(0.01)
drawnow;
%     im=frame2im(getframe(gcf));
%     [imind,map]=rgb2ind(im,256);
%    if i==1
%         % ���Ѹ���ģʽд��ָ����gif�ļ�
%        imwrite(imind,map,filename1,'gif','LoopCount',Inf,'DelayTime',0.2);
%    else
%         % ����׷��ģʽ��ÿһ֡д��gif�ļ�
%        imwrite(imind,map,filename1,'gif','WriteMode','append','DelayTime',0.2);
%    end
end
