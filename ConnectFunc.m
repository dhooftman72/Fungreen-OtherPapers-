% Matlab code belonging to Connectivity calculations
% 
% Copyright Lactuca; D.A.P. Hooftman 2018 for the FUNGREEN project
% 
% Note code may cross lines here whereas in real code does not, no additional line breaks (…) are added here. 
for site = 1:1:12
    cd ‘name folder’ % change per country to appropriate folder
    central_sites = x20; % Sweden: 120; Belgium 220; Germany 320
    alpha = 1;
 
    % Load all distance rasters in one go, assuming a strong computer
    % "cut_top_off_ascii" is restricted code
    file = ['Movement raster',int2str(site),'.asc'];
    [data_block, ~,~,~,~]   = cut_top_off_ascii(file);
    distance_raster = data_block;
    
    % Vegetation = LCM raster
    % Note this file is leading in ncols and nrows, check BEFORE whether
    % grids are identical, else redo extent clipping in ArcMAP
    file = ['Vegetation raster',int2str(site),'.asc'];
    [data_block, ncols,nrows,NODATA_value,cellsize]  = cut_top_off_ascii(file);
    lcm_raster = data_block;
    
    % Define GI, note the LCM raster is leading on NoDATA values, those 
    % rasters are circular in extent (so outside circle is NODATA) whereas distance rasters are not.
    GI_raster = zeros(nrows,ncols);
    GI_raster(lcm_raster == NODATA_value) = NODATA_value;
    GI_raster(lcm_raster == 6) = 1; % semi-natural
    GI_raster(lcm_raster == 7) = 1; % Islets
    GI_raster(lcm_raster == 12) = 1; % logged forest
    GI_raster(lcm_raster == 14) = 1; % wetlands
    GI_raster(lcm_raster == 16) = 1; % Historic semi-open forest
    GI_raster(lcm_raster == 17) = 1; % Historic grazed former arable
    GI_raster(lcm_raster > 100) = 1; % roads, hedges, complex borders, powerlines, all target sites
    cd ..
       
    % run the distance raster against the lcm        
        % Initiate values
    connectivity(type,site) = 0;
    main_size(site,1) = 0;
    semi_natural_size(site,1) = 0;
    linear_feature_size(site,1) = 0;
    inside_road_veg(site,1) = 0;
    other_veg(site,1) = 0;
    total_GI(site,1) = 0;
    Islets_veg(site,1) = 0;
    outside_road_veg(site,1) = 0;
    
    grid_per_ha = ((cellsize^2)./10000);
    
   % Run over all cells individually
    for x = 1:ncols
    for y = 1:nrows
    connectivity_added = 0;
    if GI_raster(y,x) == 1 
    if (lcm_raster(y,x) ~= (central_sites+site)) && (distance_raster(y,x) ~=  NODATA_value) 
    % Hanski IFM model per gridcell: area connected to weighted by distance"
    	  connectivity_added = exp(-alpha*(distance_raster(y,x)/1000)).*grid_per_ha; 
     % note in km (1000m) and in hectare (10000 square meters)
            total_GI(site) = total_GI(site)+grid_per_ha; 
    end
    end
    connectivity(type,site) = connectivity(type,site) + connectivity_added;
                
    calculate cover of all types for statistics
    if lcm_raster(y,x) == (central_sites+site)
       main_size(site) = main_size(site) + grid_per_ha;
    end                
    if lcm_raster(y,x) == 6 || (lcm_raster(y,x) > central_sites && lcm_raster(y,x) ~= (central_sites+site))
       semi_natural_size(site) = semi_natural_size(site) + grid_per_ha;
    end
    if lcm_raster(y,x) == 109 || lcm_raster(y,x) == 110 ||  lcm_raster(y,x) == 111 || lcm_raster(y,x) == 9
        linear_feature_size(site) =  linear_feature_size(site) + grid_per_ha;
    end
    if lcm_raster(y,x) == 180 || lcm_raster(y,x) == 80
    % split roads fully within or (partly outside urban
    % areas), testing whether a road is surrounded by only
    % urban (or other road and NODATA)
    if y ~= nrows && x~= ncols && y ~= 1 && x ~= 1
       surround_raster = unique([lcm_raster(y+1,x+1), lcm_raster(y+1,x),lcm_raster(y+1,x-1),...
lcm_raster(y,x+1),lcm_raster(y,x-1), lcm_raster(y-1,x+1),...
lcm_raster(y-1,x),lcm_raster(y-1,x-1)]);
test= find((surround_raster ~=180 | surround_raster == 80) & surround_raster ~=5 &…
 surround_raster~=NODATA_value);
  if isempty(test) == 0
     outside_road_veg(site) =  outside_road_veg(site) + grid_per_ha;
  else
      inside_road_veg(site) =inside_road_veg(site) + grid_per_ha;
  end
    end
    end
    if lcm_raster(y,x) == 7
        Islets_veg(site) =  Islets_veg(site) +  grid_per_ha;
    end
    if lcm_raster(y,x) == 12 || lcm_raster(y,x) == 14 || lcm_raster(y,x) == 16 || lcm_raster(y,x) == 17
        other_veg(site) =  other_veg(site) +  grid_per_ha;
    end
    end % end of nrows
    end % end of ncols
    connectivity_ratio(type,site) = connectivity(type,site)./ connectivity(1,site);
end % end of site/landscape
connectivity =  connectivity';
connectivity_ratio = connectivity_ratio';
save('Filename','connectivity','connectivity_ratio','total_GI','main_size',...
'semi_natural_size','linear_feature_size','outside_road_veg','Islets_veg'other_veg');
 clear all
 load('File name')
