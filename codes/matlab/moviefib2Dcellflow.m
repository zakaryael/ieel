%clear;
%Store = [];
% Path to the directory where the data has been stored
DIRNAME = '../../../runs/test1/';

IT = 0:2000:2000000;

figure(1), clf

fname = [DIRNAME sprintf('fiber%d.nc',IT(1))];
I = ncinfo(fname);
Ns = I.Dimensions(2).Length;
Data = ncread(fname,'Pos');
hpl = plot(Data(:,1),Data(:,2),'-','Linewidth',2);
for ik = -10:10
    line([ik ik]-0.5,[-10 10])
    line([-10 10],[ik ik]-0.5)
end
axis equal;

xlim(mean(Data(:,1))+[-1 1])
ylim(mean(Data(:,2))+[-1 1])

%%
V = zeros(length(IT)-1,1);
L = zeros(length(IT),1);
L(1) = sqrt(sum((Data(end,:)-Data(1,:)).^2));  
cnt = 1;
for it=IT(2:end)
    fname = [DIRNAME sprintf('fiber%d.nc',it)];
    Data = ncread(fname,'Pos');
    set(hpl,'XData',Data(:,1),'YData',Data(:,2))
    xlim(mean(Data(:,1))+[-1 1])
    ylim(mean(Data(:,2))+[-1 1])

    L(cnt+1) = sqrt(sum((Data(end,:)-Data(1,:)).^2));
    Data = ncread(fname,'Vel');
    V(cnt) = mean(Data(:,1));
    cnt = cnt+1;
    pause(0.1)
end
%%
figure(2)
dt = 1e-3;
plot(IT(2:end)*dt,V)
Store = [Store -mean(V(50:249))];
