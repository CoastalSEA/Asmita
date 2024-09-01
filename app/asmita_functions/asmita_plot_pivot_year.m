function asmita_plot_pivot_year(src,~)
%
%-------function help------------------------------------------------------
% NAME
%   asmita_plot_pivot_year.m
% PURPOSE
%   modify an asmita timeseries plot so that the volumes are set to the 
%   initial element volume at the defined pivot year.
% USAGE
%   asmita_plot_pivot_year - from command line
%   asmita_plot_pivot_year(src,evt) as a button function callback
% INPUTS
%   src - handle to calling button, parent is figure called from
%   evt - not used but required for call
% OUTPUT
%   prompts user to set initial value and pivot year and creates a new
%   figure with the adjusted data.
% USAGE
%   asmita_plot_pivot_year()
% NOTES
%   can also be called from command line with no arguments. Select figure
%   to be updated and then run function.
% SEE ALSO
%   button set in ASM_Plots.callPlotType and defined in ASM_Plots.addPivotButton
%
% Author: Ian Townend
% CoastalSEA (c) May 2023
%--------------------------------------------------------------------------
%
    if nargin<1
        ax1 = gca;
    else
        hf = src.Parent;  %figure called from
        ax1 = findobj(hf,'Type','axes');
    end

    %create a new figure for the updated lines
    hf = figure('Name','asmPlot','Tag','PlotFig');
    % Get handles for all children from ax1
    ax2 = copyobj(ax1,hf);
    clear ax1
    legend(ax2,'show')

    plines = findobj(ax2,'Type','line');
    nline = length(plines);
    xlen = length(plines(1).XData);
    x = NaT(xlen,nline); y = zeros(xlen,nline); 

    isdur = false;
    if isduration(plines(1).XData)
        %convert durations to datetime from year 0
        isdur = true;
    end

    for i=1:nline
        runtime = plines(i).XData;
        if isdur
            %convert durations to datetime from year 0
            runtime = datetime(0,1,1,0,0,0)+runtime;
            runtime.Format = 'dd-MMM-yyyy HH:mm:ss';  
        end
        x(:,i) = runtime;
        y(:,i) = plines(i).YData;
    end

    %defaults values from existing plot
    y0 = mean(y(1,:));     %average of y-values at x=0
    x0 = x(1,1);           %initial year as default pivot year
    xend = x(end,1);       %

    %check values to use with user
    promptxt = {'Initial value','Pivot year','Scale factor','X-start','X-end'};
    defaults = {num2str(y0),char(x0),'1',char(x0),char(xend)};
    answers = inputdlg(promptxt,'asmplot',1,defaults);
    if isempty(answers), return; end          %user cancelled
    y0 = str2double(answers{1});              %intial value
    x0 = datetime(answers{2});                %pivot year
    sf = str2double(answers{3});              %scale factor
    xst = datetime(answers{4});               %start year
    xnd = datetime(answers{5});               %end year
    %update lines based on defined initial value and pivot year
    for i=1:nline
        ypivot = interp1(x(:,i),y(:,i),x0,'Linear');
        yoffset = y0-ypivot;
        if sf==1
            plines(i).YData = y(:,i)+yoffset;
        else
            plines(i).YData = (y(:,i)-ypivot)*sf+y0;
        end
    end
    
    if ~isduration(ax2.XLim)
        ax2.XLim = [xst,xnd];
    end

    refresh(ax2.Parent)
end