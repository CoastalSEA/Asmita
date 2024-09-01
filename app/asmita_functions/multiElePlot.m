function multiElePlot(obj,mobj)
%
%-------function help------------------------------------------------------
% NAME
%   multiElePlot.m
% PURPOSE
%   generate a plot of all asmita elements for selected variable and
%   scaling
% USAGE
%   multiElePlot(obj,mobj,dtype)
% INPUT
%   obj - instance of ASM_plots
%   mobj - instance of mui model
% OUTPUT
%   plot of all elements for selected variable
% SEE ALSO
%   used in ASM_Plots 
%
% Author: Ian Townend
% CoastalSEA (c) Nov 2023
%--------------------------------------------------------------------------
%   
    dtype = 'array';
    muicat = mobj.Cases;
    ok = 1;
    dim = 1; %dimension to apply scaling function if matrix
    usescale = obj.UIset.scaleList{obj.UIsel(1).scale};

    time = getProperty(muicat,obj.UIsel(2),dtype);
    eleobj = getClassObj(mobj,'Inputs','Element');
    elename = getEleProp(eleobj,'EleName');
    nvar = length(elename);
    %initialise struct used in muiCatalogue.getProperty
    props(nvar) = setPropsStruct(muicat);
    %get the data and labels for each variable
    for i=1:nvar                
        obj.UIsel(1).dims(2).value = elename(i);
        props(i) = getProperty(muicat,obj.UIsel(1),dtype);
        if isempty(props(i).data), ok = ok-1; end
        if obj.UIsel(1).scale>1 %apply selected scaling to variable               
            data2use = scalevariable(props(i).data,usescale,dim);
            props(i).data = data2use;
            props(i).label = sprintf('%s-%s',usescale,props(i).label);
        end
    end
    if ok<=0, return; end
    
    %get an existing figure of create a new one
    getFigure(obj); 
    idx = obj.Plot.FigNum==obj.Plot.CurrentFig.Number;
    hfig = findobj('Number',obj.Plot.FigNum(idx));
    figax = axes('Parent',hfig,'Tag','PlotFigAxes'); 
    %generate plot            
    plot(figax,time.data,props(1).data,'DisplayName',elename{1},'ButtonDownFcn',@godisplay)
    hold on
    for j=2:nvar
        plot(figax,time.data,props(j).data,'DisplayName',elename{j},'ButtonDownFcn',@godisplay)
    end
    hold off
    xlabel(time.label)
    ylabel(props(1).label)
    legend
    titletxt = sprintf('%s for %s',props(1).desc,props(1).case);
    title(titletxt)
    %assign muiPlots instance to handle
    mobj.mUI.Plots = obj;
    if isvalid(obj.Plot.CurrentFig)
        obj.Plot.CurrentFig.Visible = 'on';
    end
end