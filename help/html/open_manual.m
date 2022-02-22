function open_manual()
%find the location of the asmita app and open the manual
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},'Asmita'));
fpath = [appinfo(idx(1)).location,'/Asmita/doc/ASMITA manual.pdf'];
open(fpath)
