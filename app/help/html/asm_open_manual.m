function asm_open_manual()
%find the location of the asmita app and open the manual
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},'Asmita'));
fpath = [appinfo(idx(1)).location,[filesep,'Asmita',filesep,'app',...
                            filesep,'doc',filesep,'ASMITA manual.pdf']];
open(fpath)
