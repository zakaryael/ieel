clear;

% Path to the directory where the data has been stored
DIRNAME = '../../../../fiberruns/test/';

IT = 0:200:50000;

figure(1), clf

fname = [DIRNAME sprintf('fiber%d.nc',IT(1))];
I = ncinfo(fname);
Ns = I.Dimensions(2).Length;
Data = ncread(fname,'Pos');
hpl = plot3(Data(1:Ns,1),Data(1:Ns,2),Data(1:Ns,3),'.-');
axis equal;
%xlim([-3 1])
ylim([-0.25 0.25])
zlim([-0.25 0.25])
%view([-90 0])
%%
for it=IT(2:end)
    fname = [DIRNAME sprintf('fiber%d.nc',it)];
    Data = ncread(fname,'Pos');
    set(hpl,'XData',Data(1:Ns,1),'YData',Data(1:Ns,2),'ZData',Data(1:Ns,3))
    ylim([-0.25 0.25])
    zlim([-0.25 0.25])
    pause(0.1)
end