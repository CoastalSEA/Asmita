    %inputs needed by get_sed_flux based on inputs defined in LiveScript
    inp.Volume = Vm0;          %element volume (m^3)  
    inp.SurfaceArea = S;       %element surface area (m^2)
    inp.Prism = P;             %tidal prism of channel (m^3)
    inp.EqScaleCoeff = alpha;  %equilibrium scale coeffient, alpha
    inp.EqShapeCoeff = 1;      %equilibrium shape coeffient, beta 
    inp.VerticalExchange = w;  %vertical exchange (m/s)
    inp.HorizontalExchange = d;%horizontal exchange (m/s)
    inp.RiverDischarge = 0;    %river discharge (m^3/s) +ve downstream
    inp.TransportCoeff = n;    %transport coefficient n (3-5)
    inp.EqConc = cE;           %equilibrium concentration (-)
    inp.RiverConc = 0;         %river load imported by advection (-)
    inp.BedConc = 1;           %concentration of bed (-)
    inp.y2s = y2s;             %conversion from years to seconds
    
    sedvol = zeros(1,100); watervol = sedvol; conc = sedvol; %inisialise
    for i=1:100
        [sedvol(i),watervol(i),conc(i)] = get_sed_flux(inp,dslr*accel(i)); %call function for sediment flux      
    end

    %report sedvol as +ve volume for sediment import
    sedvol = -sign(inp.TransportCoeff)*sedvol*dt; 
    watervol = watervol*dt;
    t = years(1:100);
    figure
    yyaxis left
    plot(t,sedvol)
    ylabel('Sediment volume (m^3)')
    yyaxis right
    plot(t,watervol)
    ylabel('Water volume (m^3)')
    xlabel ('Time')
    legend({'Sediment flux','Change in Water volume'},'Location','northwest')
    title('Single element model')
    
    