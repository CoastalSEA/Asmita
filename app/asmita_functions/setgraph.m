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
    elseif isempty(advobj)
        advobj = Advection.initialiseTransients(mobj);
    elseif isempty(eleobj)
        getdialog('No elements defined'); return;
    end

    %check that graphs have been initialised
    ok = initialiseGraphs(mobj,advobj); 
    if ok<1, return; end

    switch src.Tag
        case 'Network'
            axtag = 'axNetwork';
            [g,nlabel] = Estuary.initialiseDispersionGraph(mobj);
        case 'Rivers'
            axtag = 'axRivers';
            [g,nlabel] = Advection.initialiseRiverGraph(mobj);
        case 'Drift'
            axtag = 'axDrift';
            Estuary.initialiseDispersionGraph(mobj);
            [g,nlabel] = Advection.initialiseDriftGraph(mobj);
            if height(g.Edges)>0
                %convert from drift rate to equivalent flow rate
                nele = length(eleobj);
                [q,exchIn,exchOut,nodetxt] = graph2matrix(g,nele);
                qIn = exchIn(:,2);
                qOut = exchOut(:,1);
                Element.setEqConcentration(mobj);
                [q,qIn,qOut] = getDriftFlow(advobj,mobj,q,qIn,qOut);
                equivFlow = matrix2graph(q,qIn,qOut,nodetxt);
            end
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
    
    %add control buttons
    hpan = findobj(src,'Style','pushbutton','-regexp','Tag','panbut');
    if isempty(hpan)
        uicontrol('Parent',src, ...
            'Style', 'pushbutton', 'String', 'Pan',...
            'TooltipString','When active, right click for zoom menu',...
            'Units','normalized', ...
            'Position', [0.02 0.94 0.06 0.05], ...
            'Callback', @(src,evtdat)panButton(h_ax,src,evtdat),...
            'Tag','panbut_on');    
        uicontrol('Parent',src, ...
            'Style', 'pushbutton', 'String', 'Off',...
            'TooltipString','Pan/zoom status. Turn "off" when finished',...
            'Units','normalized', ...
            'Position', [0.08 0.94 0.06 0.05], ...
            'Callback', @(src,evtdat)onoffButton(src,evtdat), ...
            'Tag','panbut_off'); 
    end
    %
    hb = findobj(src,'Style','pushbutton','Tag','FigButton');
    if isempty(hb)
        %button to create plot as stand-alone figure
        uicontrol('Parent',src,'Style','pushbutton',...
            'String','>Figure','Tag','FigButton',...
            'TooltipString','Create plot as stand alone figure',...
            'Units','normalized','Position',[0.88 0.94 0.10 0.05],...
            'Callback',@(src,evtdat)graphPlot(hg,src,evtdat));  
    else
        hb.Callback = @(src,evtdat)graphPlot(hg,src,evtdat);
    end
    %
    if strcmp(src.Tag,'Drift') && height(g.Edges)>0
        hdrift = findobj(src,'Style','pushbutton','Tag','drift');
        h_ax.UserData = equivFlow.Edges.Weight;
        if isempty(hdrift)
            uicontrol('Parent',src, ...
                'Style', 'pushbutton', 'String', 'Qs',...
                'TooltipString','Switch between drift rate, Qs, and equivalent flow rate, Fs',...
                'Units','normalized', ...
                'Position', [0.02 0.87 0.06 0.05], ...
                'Callback', @(src,evtdat)driftSwitch(h_ax,src,evtdat), ...
                'Tag','drift'); 
        else
            hdrift.String = 'Qs';
        end
    end
end 
%%
function ok = initialiseGraphs(mobj,advobj)
    %check that the dispersion and advection graphs have been initialised
    estobj = getClassObj(mobj,'Inputs','Estuary');
    rncobj = getClassObj(mobj,'Inputs','RunConditions');
    if isempty(estobj) || isempty(rncobj)
        getdialog('Estuary or RunConditions not defined')
        ok = 0; return; 
    end

    if isempty(estobj.DispersionGraph)
        Estuary.initialiseDispersionGraph(mobj);
    end

    if isempty(advobj.RiverGraph) && rncobj.IncRiver
        Advection.initialiseRiverGraph(mobj);
    end

    if isempty(advobj.DriftGraph)  && rncobj.IncDrift
        Advection.initialiseDriftGraph(mobj);
    end
    ok = 1;
end
%%
function panButton(h_ax,src,~)            %#ok<INUSD> 
    %Create push button to enable pan and zoom
    hCM = uicontextmenu;
    uimenu('Parent',hCM,'Label','Switch to zoom',...
        'Callback',@(src,evtdat)zoomButton(src,evtdat));
    hPan = pan;
    setAllowAxesPan(hPan,h_ax,true);
    hPan.UIContextMenu = hCM;    
    hon = findobj(h_ax.Parent.Children,'Tag','panbut_on');
    hon.String = 'Pan';
    pan('on')
    hoff = findobj(h_ax.Parent.Children,'Tag','panbut_off');
    hoff.String = 'On';
end
%% 
function zoomButton(src,~)
    hzm = findobj(src.Parent.Parent.Children,'Tag','panbut_on');
    hzm.String = 'Zoom';
    zoom('on');    
end
%%
function onoffButton(src,~)
    zoom(gcbf,'off');
    pan(gcbf,'off');
    src.String = 'Off';
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
%%
function graphPlot(hg,src,~)
    %generate plot for figure button callback
    hfig = figure('Tag','PlotFig','Color','w');    
    ax = axes('Parent',hfig,'Tag','PlotFig');
    copyobj(hg,ax)
    title(sprintf('Connectivity of %s',src.Parent.Tag))
    axis off
end
