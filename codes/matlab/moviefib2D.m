clear;

% Path to the directory where the data has been stored
DIRNAME = '../../../../fiberruns/data/';

IT = 0:200:50000;
dt = 5e-4;

figure(1), clf

fname = [DIRNAME sprintf('fiber%d.nc',IT(1))];
I = ncinfo(fname);
Ns = I.Dimensions(2).Length;
Data = ncread(fname,'Pos');
hpl = plot(Data(:,1),Data(:,2),'.-');
xlim([-2 2])
ylim([-2 2])
axis equal;
%%
V = zeros(length(IT)-1,1);
L = zeros(length(IT),1);
L(1) = sqrt(sum((Data(end,:)-Data(1,:)).^2));  
cnt = 1;
for it=IT(2:end)
    fname = [DIRNAME sprintf('fiber%d.nc',it)];
    Data = ncread(fname,'Pos');
    set(hpl,'XData',Data(:,1),'YData',Data(:,2))
    L(cnt+1) = sqrt(sum((Data(end,:)-Data(1,:)).^2));
    Data = ncread(fname,'Vel');
    V(cnt) = mean(Data(:,1));
    cnt = cnt+1;
    pause(0.1)
end