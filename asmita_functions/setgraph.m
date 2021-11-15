function setgraph(mobj,src,~)
%
%-------function help------------------------------------------------------
% NAME
%   setgraph.m
% PURPOSE
%   update the Network or Flow tabs with directed graphs of network
% USAGE
%   menu/tab calbback function
% INPUT
%   mobj - handle to mui model UI
%   src - handle to graphical object making the call (eg a tab)
% OUTPUT
%   graph displayed  on tab with pan and zoom button
% SEE ALSO
%   setmatrix.m
% NOTE
%   if the pan/zoom button is not reset the other tabs will not work (do
%   not know why!)
%
% Author: Ian Townend
% CoastalSEA (c)Oct 2021
%--------------------------------------------------------------------------
%
    %update the Network or Flow tabs with directed graphs of network
    estobj = getClassObj(mobj,'Inputs','Element');
    advobj = getClassObj(mobj,'Inputs','Advection');
    if isempty(estobj) && isempty(advobj)
        cla(src.Children);
        return
    end

    switch src.Tag
        case 'Network'
            axtag = 'axNetwork';
            [g,nlabel] = Estuary.initialiseDispersionGraph(mobj);
        case 'Rivers'
            axtag = 'axRivers';
            [g,nlabel] = Advection.initialiseRiverGraph(mobj);
        case 'Drift'
            axtag = 'axDrift';
            [g,nlabel] = Advection.initialiseDriftGraph(mobj);  
        case 'TP Network'
            axtag = 'axTidalPump';
            [g,nlabel] = Advection.initialiseQtpGraph(mobj); 
        otherwise
            return;
    end

    %plot resulting graph on tab
    elabel = g.Edges.Weight;            
    h_ax = findobj('Tag',axtag);
    if isempty(h_ax)
        h_ax = axes('Parent',src,'Tag',axtag,'Color','none', ...
            'Color',[0.9,0.9,0.9], ...
            'Position', [0 0 1 1], ...
            'XColor','none', 'YColor','none','ZColor','none', ...
            'NextPlot','replacechildren');    
    end
    %LWidths = 5*elabel/max(elabel);
    plot(g,'Parent',h_ax,'EdgeLabel',elabel,'NodeLabel',nlabel);
        %'LineWidth',LWidths,

    hpan = findobj(src,'Style','pushbutton');
    if isempty(hpan)
        uicontrol('Parent',src, ...
            'Style', 'pushbutton', 'String', 'Pan',...
            'TooltipString','When active, right click for zoom menu',...
            'Units','normalized', ...
            'Position', [0.02 0.94 0.06 0.05], ...
            'Callback', @(src,evtdat)panButton(h_ax,src,evtdat));    
        uicontrol('Parent',src, ...
            'Style', 'pushbutton', 'String', 'Off',...
            'Units','normalized', ...
            'Position', [0.08 0.94 0.06 0.05], ...
            'Callback', 'zoom(gcbf,''off'');pan(gcbf,''off'')'); 
    end
end 
%%
function panButton(h_ax,~,~)            
    %Create push button to enable pan and zoom
    hCM = uicontextmenu;
    uimenu('Parent',hCM,'Label','Switch to zoom',...
        'Callback','zoom(gcbf,''on'')');
    hPan = pan;
    setAllowAxesPan(hPan,h_ax,true);
    hPan.UIContextMenu = hCM;
    pan('on')
end
        
