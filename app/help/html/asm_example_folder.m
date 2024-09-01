function asm_example_folder()
%find the location of the asmita demo folder and open it
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},'Asmita'));
fpath = [appinfo(idx(1)).location,[filesep,'Asmita',filesep,'app',filesep,'example']];
try
    winopen(fpath)
catch
    msg = sprintf('The examples can be found here:\n%s',fpath);
    msgbox(msg)
end