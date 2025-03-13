clc
clear
addpath('../../matlab');
par = get_params('params.json');
nx = par.NX;
ny = par.NY;
nz = par.NZ;
dh = par.DH;
OUT = par.OUT;
srci = par.src_i;
num_fault = par.num_fault;
f0 = par.RS_f0
fw = par.RS_fw
Vini = par.RS_Vini;
V0 = par.RS_V0;

nuc_seg = 2;
j0 = 10000;
k0 = -7500; %7500
R0 = 2000;
R = 2000;
tn0 = -40*1e6;
Te = 0.8;
tau0 = (Te*f0 + (1-Te)*fw) * tn0
Te = (abs(tau0)-fw*abs(tn0)) / ((f0-fw)*abs(tn0))

fnm = par.Fault_geometry;
x = ncread(fnm, 'x');
y = ncread(fnm, 'y');
z = ncread(fnm, 'z');
vec_n = ncread(fnm, 'vec_n');
vec_m = ncread(fnm, 'vec_s1');
vec_l = ncread(fnm, 'vec_s2');


Tx = zeros(ny, nz, num_fault);
Ty = zeros(ny, nz, num_fault);
Tz = zeros(ny, nz, num_fault);

Tn = zeros(ny, nz, num_fault);
Tm = zeros(ny, nz, num_fault);
Tl = zeros(ny, nz, num_fault);

dTx = zeros(ny, nz, num_fault);
dTy = zeros(ny, nz, num_fault);
dTz = zeros(ny, nz, num_fault);



F = zeros(ny, nz);
sigman_plot = zeros(nz, 1);
tau_plot = zeros(nz, 1);
for i = 1 : num_fault
    ix = srci(i);
%     x1 = squeeze(x(:,:, ix));
    y1 = squeeze(y(:,:, ix));
    z1 = squeeze(z(:,:, ix));
    for k = 1:nz
        if(~mod(k,50)) disp([num2str(k), '/', num2str(nz)]); end
%         depth = -z1(1, k);
        
        for j = 1:ny
        depth = abs(z1(1,k));
        if depth < 4.0e3 + dh/2 && depth > 4.0e3 - dh/2
            ks = k;
        end
        if depth < 11e3 + dh/2 && depth > 11e3 - dh/2
            kd1 = k;
        end
        if depth < 14e3 + dh/2 && depth > 14e3 - dh/2
            kd2 = k;
        end
        if depth < 5.0e3
            tn = tn0 * max(depth, dh/3.0)/5.e3;
            tau = tau0 * max(depth, dh/3.0)/5.e3;
            sigman_plot(k) = tn;
            tau_plot(k) = tau;
        else
            tn = tn0;
            tau = tau0;
            sigman_plot(k) = tn;
            tau_plot(k) = tau;
        end
            en = squeeze(vec_n(i,:,j,k));
            em = squeeze(vec_m(i,:,j,k));
            el = squeeze(vec_l(i,:,j,k));
%             Tn(j,k,i) = -120*1e6;
            Tn(j,k,i) = tn;
%             if depth < 5.0e3
%                Tn(j,k,i) = Tn(j,k,i)*max(depth, dh/3.0)/5.e3;
%             end
%             if (k == nz) 
%                 Tn(j,k,i) = -7378.0*dh/3.0; 
%             end
            Tl(j,k,i) = 0.0;
%             Tm(j,k,i) = -75*1e6;
            Tm(j,k,i) = tau;

            Tx(j,k,i)=Tn(j,k,i)*en(1)+Tm(j,k,i)*em(1)+Tl(j,k,i)*el(1);
            Ty(j,k,i)=Tn(j,k,i)*en(2)+Tm(j,k,i)*em(2)+Tl(j,k,i)*el(2);
            Tz(j,k,i)=Tn(j,k,i)*en(3)+Tm(j,k,i)*em(3)+Tl(j,k,i)*el(3);

            Txyz = [Tx(j, k, i), Ty(j, k, i), Tz(j, k, i)];

            r = sqrt((y1(j,k)-j0)^2+(z1(j,k)-k0)^2);
   
%     R = 1400;
%     R = 2000;
    if(r<R0 && i == nuc_seg)
      Fr = exp(r^2/(r^2-R0^2));
%       F(j,k) = Fr;
    else
      Fr = 0;
    end
%     tau_nuke = 100e6*Fr;
    tau_nuke = (0.87*abs(tn0)-abs(tau))*Fr; %0.86
    tn = dot(Txyz, en);
    ts_vec = Txyz-tn*en;
    ts = norm(ts_vec);
    ts = max(ts, 1);
    ts_vec_nuke = tau_nuke/ts * ts_vec;

    Txyz = 0*en+ts_vec_nuke;
    dTx(j,k,i) = Txyz(1);
    dTy(j,k,i) = Txyz(2);
    dTz(j,k,i) = Txyz(3);

            Tn(j, k, i)=Tx(j, k, i)*en(1)+Ty(j, k, i)*en(2)+Tz(j, k, i)*en(3);
            Tm(j, k, i)=Tx(j, k, i)*em(1)+Ty(j, k, i)*em(2)+Tz(j, k, i)*em(3);
            Tl(j, k, i)=Tx(j, k, i)*el(1)+Ty(j, k, i)*el(2)+Tz(j, k, i)*el(3);

        end
    end
end
if R<R0
Tnuc0 = sum(sum(dTy(:,:,nuc_seg)/1e6))

ix = srci(nuc_seg);
y1 = squeeze(y(:,:, ix));
z1 = squeeze(z(:,:, ix));
r = sqrt((y1-j0).^2+(z1-k0).^2);

F(r<R) = exp(r(r<R).^2./(r(r<R).^2-R^2));

Tnuc1 = abs(Tnuc0) / sum(sum(F))
tau_nuke = Tnuc1*F*1e6; %0.86
dTx(:,:,:) = 0;
dTy(:,:,:) = 0;
dTz(:,:,:) = 0;
for k = 1:nz
    if(~mod(k,50)) disp([num2str(k), '/', num2str(nz)]); end
    for j = 1:ny
        en = squeeze(vec_n(nuc_seg,:,j,k));
        Txyz = [Tx(j, k, nuc_seg), Ty(j, k, nuc_seg), Tz(j, k, nuc_seg)];
        tn = dot(Txyz, en);
        ts_vec = Txyz-tn*en;
        ts = norm(ts_vec);
        ts = max(ts, 1);
        ts_vec_nuke = tau_nuke(j,k)/ts * ts_vec;

        Txyz = 0*en+ts_vec_nuke;
        dTx(j,k,nuc_seg) = Txyz(1);
        dTy(j,k,nuc_seg) = Txyz(2);
        dTz(j,k,nuc_seg) = Txyz(3);
    end
end
Tnuc2 = sum(sum(dTy(:,:,nuc_seg)))
Tnuc2 - Tnuc0*1e6
end

j1 = par.Fault_grid(1:4:end-3);
j2 = par.Fault_grid(2:4:end-2);
k1 = par.Fault_grid(3:4:end-1);
k2 = par.Fault_grid(4:4:end);

%% rate state friction
b = zeros(ny, nz, num_fault);
a = zeros(ny, nz, num_fault);
Vw = zeros(ny, nz, num_fault);
L = zeros(ny, nz, num_fault);
State = zeros(ny, nz, num_fault);
B = ones(ny, nz, num_fault);
% 
for i = 1 : num_fault
    ix = srci(i);
%     x1 = squeeze(x(:,:, ix));
    y1 = squeeze(y(:,:, ix));
    z1 = squeeze(z(:,:, ix));
    
%     if i == 1
%         W = 25000; w = 3000;
%         
%     else
%         W = 25000; w = 3000;
%     end
%     By = Bfunc(y1-25000, W, 0);
%     Bz = Bfunc(z1+7.5e3, 7500, 0);
%     
%     B = (1-By.*Bz);
%  
%     if i == nuc_seg
%         B(j2(i)+1:end, :) = 1;
%     else
%         B(1:j1(i)-1, :) = 1;
%     end
    a(:,:,i) = 0.01;
    Vw(:,:,i) = 0.1;%+0.9*B;%(:,:,i);
    for j = 1 : ny
        a(j,ks:k2(i),i) = linspace(0.01, 0.02, length(ks:k2(i)));
        Vw(j,ks:k2(i),i) = linspace(0.1, 1, length(ks:k2(i)));
        a(j,kd2:kd1,i) = linspace(0.024, 0.01, length(kd2:kd1));
        Vw(j,kd2:kd1,i) = linspace(1, 0.1, length(kd2:kd1));
    end
%     a(:,:,i) = 0.01+0.01*B;%(:,:,i);
%     a(:,1:128,i) = 0.01+0.01*B2(:,1:128);%(:,:,i);
    a(:,1:kd2,i) = 0.024;
    a(:,k2(i):end,i) = 0.02;
    Vw(:,k2(i):end,i) = 1;
    Vw(:,1:kd2,i) = 1;
    b(:,:,i) = 0.014;
%     L(:,:,i) = 0.02015;  % AL
    L(:,:,i) = 0.2;        % SL
%     b = ones(size(a))*0.012;
%     L = ones(size(a))*0.02;

    Tau = sqrt(Tm(:,:,i).^2+Tl(:,:,i).^2);
    State(:,:,i) = a(:,:,i).*log(2.0*V0/Vini*sinh(Tau./abs(Tn(:,:,i))./a(:,:,i)));
end
figure;
subplot(1,3,1);
plot(a(j1(2)+1,:,2), z1(1,:)/1000, 'r','LineWidth',1.5);
hold on;
plot(b(j1(2)+1,:,2), z1(1,:)/1000, 'k','LineWidth',1.5);
hold on;
ylim([-15, 0]);
xlim([0, 0.03]);
ylabel('Depth (km)');
legend('a','b');
set(gca, 'FontSize', 12);
% ax1_pos = ax1.Position;
% set(ax1,'XColor','k','YColor','k');
% ax2=axes('Position',ax1_pos,...
%            'XAxisLocation','top',...
%            'XColor','g');
% hold on;
subplot(1,3,2);
plot(Vw(j1(2)+1,:,2), z1(1,:)/1000, 'b', 'LineWidth', 1.5);
xlim([0,1]);
ylim([-15, 0]);
xlabel('m/s');
% ylabel('Depth (km)');
legend('Vw');
set(gca, 'FontSize', 12);

% figure;
subplot(1,3,3);
plot(-sigman_plot/1e6, z1(1,:)/1000, 'r','LineWidth',1.5);
hold on;
plot(-tau_plot/1e6, z1(1,:)/1000, 'k','LineWidth',1.5);
hold on;
ylim([-15, 0]);
xlim([0, 50]);
% ylabel('Depth (km)');
xlabel('MPa');
legend('\sigma','\tau');
set(gca, 'FontSize', 12);

% 
% Tau = sqrt(Tm.^2+Tl.^2);
% State = a.*log(2.0*V0/Vini*sinh(Tau./abs(Tn)./a));

fnm_out = par.Fault_init_stress
ncid = netcdf.create(fnm_out, 'NETCDF4');

dimid(1) = netcdf.defDim(ncid,'ny',ny);
dimid(2) = netcdf.defDim(ncid,'nz',nz);
dimid(3) = netcdf.defDim(ncid,'num_fault',num_fault);

dimidxyz(1) = dimid(1);
dimidxyz(2) = dimid(2);
dimidxyz(3) = netcdf.defDim(ncid, 'nx', nx);
varid(1) = netcdf.defVar(ncid,'x','NC_DOUBLE',dimidxyz);
varid(2) = netcdf.defVar(ncid,'y','NC_DOUBLE',dimidxyz);
varid(3) = netcdf.defVar(ncid,'z','NC_DOUBLE',dimidxyz);
varid(4) = netcdf.defVar(ncid,'Tx','NC_DOUBLE',dimid);
varid(5) = netcdf.defVar(ncid,'Ty','NC_DOUBLE',dimid);
varid(6) = netcdf.defVar(ncid,'Tz','NC_DOUBLE',dimid);
varid(7) = netcdf.defVar(ncid,'dTx','NC_DOUBLE',dimid);
varid(8) = netcdf.defVar(ncid,'dTy','NC_DOUBLE',dimid);
varid(9) = netcdf.defVar(ncid,'dTz','NC_DOUBLE',dimid);
varid(10) = netcdf.defVar(ncid,'a','NC_DOUBLE',dimid);
varid(11) = netcdf.defVar(ncid,'b','NC_DOUBLE',dimid);
varid(12) = netcdf.defVar(ncid,'L','NC_DOUBLE',dimid);
varid(13) = netcdf.defVar(ncid,'Vw','NC_DOUBLE',dimid);
varid(14) = netcdf.defVar(ncid,'State','NC_DOUBLE',dimid);
% varid(15) = netcdf.defVar(ncid,'faultgrid', 'NC_INT',dimid);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid(1),x);
netcdf.putVar(ncid,varid(2),y);
netcdf.putVar(ncid,varid(3),z);
netcdf.putVar(ncid,varid(4),Tx);
netcdf.putVar(ncid,varid(5),Ty);
netcdf.putVar(ncid,varid(6),Tz);
netcdf.putVar(ncid,varid(7),dTx);
netcdf.putVar(ncid,varid(8),dTy);
netcdf.putVar(ncid,varid(9),dTz);
netcdf.putVar(ncid,varid(10),a);
netcdf.putVar(ncid,varid(11),b);
netcdf.putVar(ncid,varid(12),L);
netcdf.putVar(ncid,varid(13),Vw);
netcdf.putVar(ncid,varid(14),State);
% netcdf.putVar(ncid,varid(15),faultgrid);
netcdf.close(ncid);

figure;
for i = 1 : num_fault
    ix = srci(i);
    x1 = squeeze(x(:,:, ix));
    y1 = squeeze(y(:,:, ix));
    z1 = squeeze(z(:,:, ix));
    Stress = squeeze(Ty(j1(i):j2(i),k1(i):k2(i), i) + dTy(j1(i):j2(i),k1(i):k2(i), i));
    surf(x1(j1(i):j2(i),k1(i):k2(i)), y1(j1(i):j2(i),k1(i):k2(i)), z1(j1(i):j2(i),k1(i):k2(i)), Stress);
    hold on;
    view([-60,30]);
    shading interp;
    colorbar;
    axis image;
end
