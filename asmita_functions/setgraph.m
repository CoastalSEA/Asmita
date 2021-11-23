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
    eleobj = getClassObj(mobj,'Inputs','Element');
    advobj = getClassObj(mobj,'Inputs','Advection');
    if isempty(eleobj) && isempty(advobj)
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
            %convert from drift rate to equivalent flow rate
            nele = length(eleobj);
            [q,exchIn,exchOut,nodetxt] = graph2matrix(g,nele);
            qIn = exchIn(:,2);
            qOut = exchOut(:,1);
            Element.setEqConcentration(mobj);
            [q,qIn,qOut] = getDriftFlow(advobj,mobj,q,qIn,qOut);
            equivFlow = matrix2graph(q,qIn,qOut,nodetxt);
        case 'TP Network'
            axtag = 'axTidalPump';
            [g,nlabel] = Advection.initialiseQtpGraph(mobj); 
            
        otherwise
            return;
    end

    %get element type to set the node colour
    elabel = g.Edges.Weight;  %magnitude of exchanges
    if isempty(g.Nodes)
        ntype = 0;
    else
        eletype = g.Nodes.Type;   %type of element
        ntype = zeros(size(eletype));
        for i=1:length(eletype)
            nt = find(strcmp(mobj.GeoType,eletype{i}));
            if ~isempty(nt)
                ntype(i) = nt;
            end
        end
    end
    
    %plot resulting graph on tab
    h_ax = findobj('Tag',axtag);
    if isempty(h_ax)
        h_ax = axes('Parent',src,'Tag',axtag,'Color','none', ...
            'Color',[0.9,0.9,0.9], ...
            'Position', [0 0 1 1], ...
            'XColor','none', 'YColor','none','ZColor','none', ...
            'NextPlot','replacechildren');    
    end
    %LWidths = 5*elabel/max(elabel);
    hg = plot(g,'Parent',h_ax,'EdgeLabel',elabel,'NodeLabel',nlabel);
        %'LineWidth',LWidths,
    hg.NodeCData = ntype;
    
    hpan = findobj(src,'Style','pushbutton','Tag','panbut');
    if isempty(hpan)
        uicontrol('Parent',src, ...
            'Style', 'pushbutton', 'String', 'Pan',...
            'TooltipString','When active, right click for zoom menu',...
            'Units','normalized', ...
            'Position', [0.02 0.94 0.06 0.05], ...
            'Callback', @(src,evtdat)panButton(h_ax,src,evtdat),...
            'Tag','panbutxt');    
        uicontrol('Parent',src, ...
            'Style', 'pushbutton', 'String', 'Off',...
            'Units','normalized', ...
            'Position', [0.08 0.94 0.06 0.05], ...
            'Callback', 'zoom(gcbf,''off'');pan(gcbf,''off'')', ...
            'Tag','panbut'); 
    end
    %
    if strcmp(src.Tag,'Drift') 
        hdrift = findobj(src,'Style','pushbutton','Tag','drift');
        h_ax.UserData = equivFlow.Edges.Weight;
        if isempty(hdrift)
            uicontrol('Parent',src, ...
                'Style', 'pushbutton', 'String', 'Qs',...
                'Units','normalized', ...
                'Position', [0.92 0.94 0.06 0.05], ...
                'Callback', @(src,evtdat)driftSwitch(h_ax,src,evtdat), ...
                'Tag','drift'); 
        else
            hdrift.String = 'Qs';
        end
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
%%
function driftSwitch(h_ax,src,~)
    %switch between drift rate and equivalent flow rate
    hg = findobj(h_ax,'Type','GraphPlot');
    
    if strcmp(src.String,'Qs')
        src.String = 'Fs';        
    else
        src.String = 'Qs';
    end
    temp = hg.EdgeLabel;
    hg.EdgeLabel = h_ax.UserData;
    h_ax.UserData = temp;
end
