clear;

% Path to the directory where the data has been stored
DIRNAME = '../../../../fiberruns/data_swim/';

IT = 0:200:50000;

figure(1), clf

fname = [DIRNAME sprintf('fiber%d.ff',IT(1))];
fid = fopen(fname,'r');
Data = fread(fid,[3 inf],'double');
Ns = size(Data,2)/3;
fclose(fid);
hpl = plot3(Data(1,1:Ns),Data(2,1:Ns),Data(3,1:Ns),'.-');
axis equal;
xlim([-3 1])
ylim([-0.5 0.5])
zlim([-0.5 0.5])

%%
for it=IT(2:end)
    fname = [DIRNAME sprintf('fiber%d.ff',it)];
    fid = fopen(fname,'r');
    Data = fread(fid,[3 inf],'double');
    fclose(fid);
    set(hpl,'XData',Data(1,1:Ns),'YData',Data(2,1:Ns),'ZData',Data(3,1:Ns))
    ylim([-0.5 0.5])
    zlim([-0.5 0.5])
    pause(0.1)
end