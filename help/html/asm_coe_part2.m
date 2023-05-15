function asm_coe_part2()
%find the location of the asmita app and open the manual
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},'Asmita'));
fpath = [appinfo(idx(1)).location,[filesep,'Asmita',filesep,'doc',...
                                            filesep,'COE paper-p2.pdf']];
open(fpath)