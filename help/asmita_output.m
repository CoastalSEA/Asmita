%% Output
% At the end of a model run the sediment and water balances are displayed 
% in a message window. The water balance considers only elements that are 
% within the estuary, i.e channel, tidal flats, saltmarsh and storage 
% element. 

%% Output Variables
% Output variables are held in two tables: one contains variables that apply to all
% elements and the other has variables that only apply to reaches
%%
% _*Element output variables*_
%%
% * *Moving Volume*, $m^3$, Volume of the element including any changes due to water level variations (e.g. due to slr and tidal range).
% * *Fixed Volume*	$m^3$,	Morphological volume of the element (excludes water level change)
% * *Equilibrium Volume*, $m^3$, Target volume defined by the specified equilibrium relationship
% * *Moving Depth*, $m$, Depth of the element including any changes due to water level variations (e.g. due to slr and tidal range).
% * *Fixed Depth*	$m$,	Depth of the element relative to a fixed surface (excludes water level change)
% * *Equilibrium Depth*, $m$, Target depth defined by the specified equilibrium relationship
% * *Biological Production*, $m^3$,	Volume of organic matter produced in the element
% * *Concentration*, ppm,	Sediment concentration in the element
% * *Water Level Change*, $m$, Water level change in the element
%%
% _*Reach output variables*_
%%
% * *Reach Prism*, $m^3$,	Volume of the tidal prism within each reach
% * *Tidal Prism*, $m^3$,	Cumulative tidal prism for all reach elements upstream (*)
% * *Upstream CSA*, $m^2$,	Cross-sectional area of the upstream interface of a reach
% * *River Flow*, $m/s$,	River flow velocity at the upstream interface of a reach
% * *Mean Water Level*, $mAD$,	Mean sea level in each reach
% * *High Water evel*, $mAD$,	High water level in each reach
% * *Low Water Level*, $mAD$,	Low water level in each reach
% * *Tidal Range*, $m$, Tidal range in each reach
%%
% [ mAD – metres above datum. (*) - includes reach prism of the element ]
% 
%% Interpretation when using water and sediment volumes (±n)
% <<output_note.png>>
%%
% *Elements using water volumes*
%
% The changes to the fixed volume, $V_f$, are due to sediment exchange, 
% ignoring any changes in the water level. It can be thought of as the change 
% in the water volume below a fixed initial surface (such as high water), 
% or the morphological change.  In contrast, the moving volume,  $V_m$, 
% includes both the morphological changes and the water volume changes, 
% and is equivalent to the total change in the water volume. 
% Another way of thinking of this is that positive changes of the water 
% volume, $+\Delta V_m$  , increases the accommodation space, whereas the import 
% of sediment into the element $(-\Delta V_f)$ reduces it. Hence for a linear sea 
% level rise the typical behaviour is for the fixed volume to reduce 
% linearly and the moving volume to initially increase to create an excess 
% water volume (model spin up in response to the imposition of the sea level 
% perturbation) and then, if there is sufficient sediment supply, to be 
% approximately constant. The adjustment of the moving volume, relative to 
% the initial volume, implies the development of a constant over-depth 
% (because the area of each element is constant). Any excess accommodation 
% space due sea level rise is infilled albeit with a lag that reflects the 
% morphological response time (Townend et al., 2007).
%
% *Elements using sediment volumes*
%
% Change in the fixed volume, $V_f$,  reflect import or export of sediment. 
% The moving volume,  $V_m$,   changes relative to the initial volume to 
% reflect changes in sediment volume and any changes in the water volume. 
% In this case it is only meaningful to think of this as changes in the 
% accommodation space, relative to some initial volume. Under rising sea 
% levels, changes in the moving volume,  $\Delta V_m$ , are negative (because of 
% the change in sign of n) and, if the element infills, there is a positive 
% change in the fixed volume,  $+\Delta V_f$. Thus, as with water volume elements, 
% an increase in water level increases  the accommodation space and import 
% of sediment reduces it. Hence, in an infilling element the fixed volume, 
% $V_f$, increases over time, whereas the moving volume initially reduces 
% and then becomes approximately constant. This again reflects the creation 
% of a constant excess water volume, or over-depth, once the initial excess 
% accommodation space, caused by the imposition of sea level rise, is 
% infilled.

%% See Also
% The <matlab:open_manual manual> provides further details of setup and 
% configuration of the model.
%%
% Townend IH, Wang ZB, Rees JG, 2007, Millennial to annual volume changes in the Humber Estuary, Proc.R.Soc.A, 463, 837-854.
% 
