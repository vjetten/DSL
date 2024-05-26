# Ethiopia erosion model
# water balance on daily basis 
# erosion calc from LISEM assuming 1 hour of runoff
# V. Jetten May 2024, ITC
# version 1.12

#! --matrixtable --radians --lddout 

binding

### input maps
  dem = dem.map;
  Grad = grad.map;
  LDD = ldd.map;
  landuse = landuse.map;
  mask = mask0.map;                 
  chanmask = chanmask.map;
  
  Ksat1 = ksat1.map;                # ksat map in mm/h
  theta_s = thetas1.map;            # porosity (fraction)
  thetainit = thetai1.map;          # initial moisture content
  theta_fc = fieldcap1.map;         # field capacity (fraction)
  theta_wp = wilting1.map;          # wilting point (fraction)
  Smax = smax.map;                  # max canopy storage, interception mm
  SD1 = soildep1.map;      # soildepth root zone 
  #SD2 = soildep2.map; 
  Kfactor = aggrstab.map;   # splash erosion g/J (EUROSEM)
  Coh = coh.map;      # coh soil kP
  CohVeg = cohadd.map;  # coh roots kPa
  D50 = d50.map;
  PH = ch.map;  #plant height m

# not used
  #ID = id.map;
  #outlet = outlet.map;
  #stations = rainstations.map;
  #Points = outlet.map;  #


### calibration factors
  coh_cal = 2.0;      # lower is more flow erosion 
  infil_cal = 0.75;    # lower is more runoff, calibration factor Ksat, influences infiltration  
                    # calibrate until the fraction runoff is acceptable

  rain_fact = 2.0; # CHECK if the rainfall images are in mm/day, if so set this factor to 1.0 

  dt = scalar(600);

### Output maps
  SoilMoisture = moist;                 # daily soil moisture maps (mm)
  lambda  = lambda.map;
  intccum = interc.map;                 # cumulatiuve interception (mm)
  Pcum = raintotal.map;

# output tables
  p_tss = rainavg.tss;
  pcum_tss = raincumavg.tss;
  ETpcum_tss = ETpcumavg.tss;      # average cumulative ETp mm
  ETpavg_tss = etpavg.tss;         # average daily ETp (mm)
  intccum_tss = interception.tss;        # average cumulative interception (mm)
  eta_tss = ETaavg.tss;            # average daily ETa (mm)
  etacum_tss = ETacumavg.tss;      # average cumulative ETa (mm)
  ETfact_tss = ETfactor.tss;       # average daily ratio ETa/ETp
  theta_tss = theta.tss;           # average daily theta  (-)
  moisture_tss = moisture.tss;     # average daily soil moisture (mm) 

areamap
  mask;
  
  
timer
  1 122 1;  # runs 365 days from the startday onwards

initial
  hoursrain = 1.0;    # nr hours it rains on average, influences infl

  count = 1*mask; # used to count rainy days
  nLU = ordinal(landuse); # output
  report nlu.map=nLU;

  #############################################  
  ### read soil and vegetation related data ###
  ############################################# 
  
  Kc = 1.0;     # FAO crop factor Kc
  day = scalar(1);  

  Ksat = Ksat1 * mask;           # saturated conductivity mm/day
  theta = thetainit;             # initialize soil moisture theta at half the wilting point (arbitrary)           
  SoilMoisture = theta*SD1;      # initial soil moisture in mm 

  ks = ln(min(1000.0,max(0.5,Ksat)));
  lambda = max(0.1, min(0.7, 0.0849*ks+0.159));	#Brooks Corey lambda needed for percolation

  # Givers transport capacity factors
  Cg = ((D50+5)/0.32)**-0.6;
  Dg = ((D50+5)/300)**0.25;
  
  Grad += 0.005;
  Grad = if(chanmask eq 1,0.01,Grad);

  Ksat = Ksat1 * mask;  

  #https://www.researchgate.net/publication/266227873_A_simple_formula_for_predicting_settling_velocity_of_sediment_particles
  GRAV = 9.81;
  Vs = 2*(2650.0-1000.0)*GRAV*((D50/2000000.0)**2)/(9*0.001);  # Stokes settling velocity if D50 < 100
  report vs.map=Vs;

  ########################################
  ### Initialize vegetation cover maps ###
  ########################################  

  ndvi0 = "$1/ndvi_1.map";
  ndvi1 = "$1/ndvi_2.map";
  ndvi2 = "$1/ndvi_3.map";
  ndvi3 = "$1/ndvi_4.map"; 
  ndvi4 = "$1/ndvi_5.map";
  ndvi5 = "$1/ndvi_6.map";
  ndvi6 = "$1/ndvi_7.map";
  ndvi7 = "$1/ndvi_8.map";
  ndvi8 = "$1/ndvi_8.map";

  daycov0 = 1;
  daycov1 = 16;
  daycov2 = 32;
  daycov3 = 48;
  daycov4 = 54;
  daycov5 = 70;
  daycov6 = 86;
  daycov7 = 102;
  daycov8 = 122;

# corrections: if a value in a ndvi map is less than 0.2 (suspiciopus/cloud)
# AND the values of the previous AND next images are higher, so it presents a dip, than take the average 
# of the previous and the next inmage, else take the recorded ndvi value
  thr = 0.2;
 
  ndvia = ndvi0;
  ndvib = ndvi1;
  ndvic = ndvi2;
  ndvi1 = if ((ndvib lt thr) and (ndvia gt thr) and (ndvic gt thr) and (ndvia gt ndvib) and (ndvic gt ndvib),(ndvia+ndvic)/2.0,ndvib);
  
  ndvia = ndvi1;
  ndvib = ndvi2;
  ndvic = ndvi3;
  ndvi2 = if ((ndvib lt thr) and (ndvia gt thr) and (ndvic gt thr) and (ndvia gt ndvib) and (ndvic gt ndvib),(ndvia+ndvic)/2.0,ndvib);
  
  ndvia = ndvi2;
  ndvib = ndvi3;
  ndvic = ndvi4;
  ndvi3 = if ((ndvib lt thr) and (ndvia gt thr) and (ndvic gt thr) and (ndvia gt ndvib) and (ndvic gt ndvib),(ndvia+ndvic)/2.0,ndvib);
  
  ndvia = ndvi3;
  ndvib = ndvi4;
  ndvic = ndvi5;
  ndvi4 = if ((ndvib lt thr) and (ndvia gt thr) and (ndvic gt thr) and (ndvia gt ndvib) and (ndvic gt ndvib),(ndvia+ndvic)/2.0,ndvib);
  
  ndvia = ndvi4;
  ndvib = ndvi5;
  ndvic = ndvi6;
  ndvi5 = if ((ndvib lt thr) and (ndvia gt thr) and (ndvic gt thr) and (ndvia gt ndvib) and (ndvic gt ndvib),(ndvia+ndvic)/2.0,ndvib);
  
  ndvia = ndvi5;
  ndvib = ndvi6;
  ndvic = ndvi7;
  ndvi6 = if ((ndvib lt thr) and (ndvia gt thr) and (ndvic gt thr) and (ndvia gt ndvib) and (ndvic gt ndvib),(ndvia+ndvic)/2.0,ndvib);

  ndvia = ndvi6;
  ndvib = ndvi7;
  ndvic = ndvi8;
  ndvi7 = if ((ndvib lt thr) and (ndvia gt thr) and (ndvic gt thr) and (ndvia gt ndvib) and (ndvic gt ndvib),(ndvia+ndvic)/2.0,ndvib);

  #########################
  ### initialize totals ###
  #########################

  Perc = 0*mask;
  Pcum = 0*mask;
  Covcum = 0*mask;
  kc = 0*mask;
  ETp = 0*mask;
  ETpcum = 0*mask;
  ETacum = 0*mask;
  FTfactor = 0*mask;
  nrCells = maptotal(mask);
  interception = 0*mask;
  intccum = 0*mask; 
  infilcum = 0*mask;
  perccum = 0*mask;
  rocumout = 0*mask;

  nrCells = maptotal(mask);
  nrCellsNC = maptotal(if(chanmask eq 1,0,mask));
  rocum = 0*mask;

  dx = celllength();

  Sed = 0*mask;
  erosion = 0*mask;
  deposition= 0*mask;
  Dscum =0*mask;
  Dfcum = 0*mask;
  Depcum = 0*mask;
  ETfcum = 0*mask;
  Conc = 0*mask;

 
dynamic

  ########################
  ### meteo data input ###
  ########################


  rain = timeinput("$1/CHIR")*mask/rain_fact;
  #CASE SENSITIVE!!!! is the filename is CHIR0000.001 then this should be CHIR, not chir
  
 # hoursrain = if (rain gt 10, 2.0, 1.0);

  report pdayavg.tss = maptotal(rain)/nrCells;
  report p_tss = maptotal(Pcum)/nrCells;
  # write a graph of the average daily rainfall
  
  report Pcum = Pcum + rain;
  #calculate cumulative P for outut
  report Pcum_tss = maptotal(Pcum)/nrCells;   
  # write a graph of the average cumulative rainfall

 # get the rainfall values 
 # ETo = timeinput(ERA)*mask;

 # ETo = timeinputscalar(ETo_tss, nominal(mask));    # read potential evapotranspiration from a file and give the whole area that value
  #ETo = timeinput(ETA)*mask/rain_fact;
 
   #Tahir, Mekete & Mengistu, Ashenafi & Mersso, Berhan. (2018). 
   #Evaluation of livestock feed balance under mixed crop–livestock production system in the central highlands of Ethiopia. Agriculture and Food Security. 7. 10.1186/s40066-018-0170-8. 
   ETo = 4.2*mask;
   ETo = if(day gt 30, 3.0, ETo);
   ETo = if(day gt 61, 3.1, ETo);
   ETo = if(day gt 92, 3.5, ETo);
   ETo = if(day gt 122, 3.9, ETo);
    
  Kc = 1.0; # 
  ETp = Kc * ETo;
  ETpcum = ETpcum + ETp;


  ######################################
  ### vegetation ndvi interpolation ###
  ######################################     
    
  #linear interpolation between two images
  Ndvi = ndvi0*mask;
  Ndvi = if (day ge daycov0, ndvi0 + (day - daycov0)/(daycov1-daycov0 +0.01)*(ndvi1-ndvi0), Ndvi); 
  Ndvi = if (day ge daycov1, ndvi1 + (day - daycov1)/(daycov2-daycov1 +0.01)*(ndvi2-ndvi1), Ndvi); 
  Ndvi = if (day ge daycov2, ndvi2 + (day - daycov2)/(daycov3-daycov2 +0.01)*(ndvi3-ndvi2), Ndvi); 
  Ndvi = if (day ge daycov3, ndvi3 + (day - daycov3)/(daycov4-daycov3 +0.01)*(ndvi4-ndvi3), Ndvi); 
  Ndvi = if (day ge daycov4, ndvi4 + (day - daycov4)/(daycov5-daycov4 +0.01)*(ndvi5-ndvi4), Ndvi); 
  Ndvi = if (day ge daycov5, ndvi5 + (day - daycov5)/(daycov6-daycov5 +0.01)*(ndvi6-ndvi5), Ndvi); 
  Ndvi = if (day ge daycov6, ndvi6 + (day - daycov6)/(daycov7-daycov6 +0.01)*(ndvi7-ndvi6), Ndvi); 
  Ndvi = if (day ge daycov7, ndvi7 + (day - daycov7)/(daycov8-daycov7+0.01)*(ndvi8-ndvi7), Ndvi); 
  Ndvi = if (day ge daycov8, ndvi8, Ndvi); 
 # Ndvi = if (day ge daycov8, ndvi8 + (day - daycov8)/(daycov9-daycov8+0.01)*(ndvi9-ndvi8), Ndvi); 
 # Ndvi = if (day ge daycov9, ndvi9 + (day - daycov9)/(daycov10-daycov9+0.01)*(ndvi10-ndvi9), Ndvi); 
 # Ndvi = if (day ge daycov10, ndvi10 + (day - daycov10)/(daycov11-daycov10+0.01)*(ndvi11-ndvi10), Ndvi); 
 # Ndvi = if (day ge daycov11, ndvi11 + (day - daycov11)/(daycov12-daycov11+0.01)*(ndvi12-ndvi11), Ndvi); 
 # Ndvi = if (day ge daycov11, ndvi12, Ndvi); 

  Ndvi = cover(Ndvi,0)*mask;
  NDVI = min(1.0,max(0,Ndvi));                         
  #make sure NDVI is between 0 and 1!
  
  #Cover = 1 -exp(-2*NDVI/(1.5 - NDVI));              # Converting NDVI values to Cover factor according to Van der Knijff et al., 1999
  a_ = 4.257;
  b_ = 100.719;
  c_ = -5.439;
  Cover = min(0.99,max(0.0,a_*NDVI*NDVI + b_*NDVI + c_)/100.0)*mask;
  report Cover.tss = timeoutput(nLU, areaaverage(Cover, nLU));    
 # report avg cover per land use type in a graph
  Litter = Cover;
  Covcum = Covcum + Cover;							   # cumulative cover for average

  N = 0.005 + 0.01*rr.map + 0.1*Cover;   # mannings n depends on cover
  N = if(chanmask eq 1,0.05,N);

  Cohesion = coh_cal*max(-1,Coh + Cover*CohVeg);
  Cohesion = if (chanmask eq 1,-1,Cohesion);

  #################################
  ######### Interception ##########
  #################################

  interception = interception + rain - ETp;                 # add rainfall and subtract evaporation from interception
  interception = min(Smax, interception);                   # fill up the interception with rain, make sure water has Smax = 0   
  interception = max(0, interception);                      # interception cannot be less than 0
  Pe = max(rain - interception, 0);                         # Estimating effective rainfall (mm)
  intccum = intccum + interception;                         # cumulating interception
  report intccum_tss = maptotal(intccum)/nrCells;    

  ### NOTE: Infiltration and runoff and flow and splash erosion are solved in 10 min intervals with a kinematic wave for Q and V
  
  ###############################
  ### Infiltration and Runoff ###
  ###############################

  count = 0;
  maxcount = 3600/dt*hoursrain; # 10 min timesteps for kin wave

  W = dx/3; # flow width m, water is concentrated in 1/3 of the gridcell
  Infilday = 0;  # to calc daily infiltration

  WHavg = 0*mask;
  Qavg = 0*mask;

  ###loop that does the fast processes per 10 min and uses a kin wave for routing
  repeat {

    Pinterval = Pe*exp(-0.693*(count+1)); #mm/dt
    # this gives a rainfall for each count that is Pe/2, Pe/4, Pe/8 etc. an exponential decline

    store = max(0.0, theta_s-theta)*SD1;                      # effective storage in mm
    Infilcap = infil_cal*hoursrain/maxcount*Ksat*mask;        # infil is a fraction of Ksat 
    Infilcap = min(store, Infilcap);                          # infiltration in mm, is smallest of storage or rainfall
    Infil = min(Infilcap, Pinterval);
    Infilday = Infilday + Infil;
  
    infilcum = infilcum + Infil;
    # infil cumulative in mm on points
    ### Runoff water height during hours of rain
    WH = cover(max(0, Pinterval-Infilcap)*0.001,0)*mask;  #waterheight in m avg per hour
  
    WH *= dx/W;  # assume the flow is concentrated over 1/3 of the gridcell

    # manning's velocity and discharge and alpha: A = alphaQ^0.6
    V = WH**(2/3)*sqrt(Grad)/N;
    Q = V * WH * W;
    alpha =(N/sqrt(Grad) * (W**(2.0/3.0)))**0.6;
    ts = ordinal(1);         # not sure what this is
    Qn = kinematic(LDD, Q, 0, alpha, 0.6, ts, dt, dx);

    # new Q, alpha and WH after kin wave
    Qn = cover(Qn,0)*mask;
    A = alpha*(Qn**0.6);
    WH = A/W;
    
    # new velocity for flow erosion
    V = WH**(2/3)*sqrt(Grad)/N;
    Qn = V*W*WH;
	
    # volume water needed for sed concentration
    Vol = WH*W*dx;
  
    WHavg = WHavg + WH;

    # runoff leaving the area
    #rocumout = maptotal(rocum * if(pit(LDD) ne 0,scalar(1),0));
    #rocumout = rocumout + maptotal(if(pit(LDD) ne 0,Qn*dt/(dx*dx)*1000,0));
    
    Qavg = Qavg + Qn;
 
    #############################
    ########### Erosion #########
    #############################
  
    #######################################################
    ### Kinetic energy of rainfall in J/m2 and SPLASH Ds ##
    #######################################################
   
    #Int = rain/hoursrain; # avg rainfall intensity in mm/h
    Intensity = Pinterval*3600/dt;
    # kinetic energy calculations for splash detachment
    # DT = direct rainfall on open part of the gridcell, 
    # LD is leaf drainage from covered part
    KE_DT = (1-Cover)*28.3*(1-0.52*exp(-0.042*Intensity))*Pinterval;          # kinetic energy of direct rainfall     
    KE_LD = (1-Litter)*Cover*Pinterval*max(0, 15.3*sqrt(PH) - 5.87);            # kinetic energy of LD in J/m2 proposed by Brandt (1990)
    
    KE = if(rain eq 0, 0, KE_DT + KE_LD);                     # total kinetic energy of rainfall
    splashfactor = (W+0.15*(dx-W))/dx;                        # splash from the dry part in a gridcell to the runoff flowing part
    Ds = splashfactor*Kfactor*KE*0.001*cellarea();  # detachment by raindrop g/J * J/m2 * 0.001 = kg/m2 * cell size = kg
    Ds = if(Kfactor eq -1, 0,Ds);
    #Ds = Ds * exp(-1.48*WH*1000);                   # no splash when the water is deep then a few mm
    Ds = if (WH gt 0, Ds, 0);                       # splash erosion when there is no water at the surface does not count 
    Sed = Sed + Ds;                                 # add splash sediment to flow and cal concentration for flow detachment  
    Conc = min(840, if (Vol gt 1e-6, Sed/Vol, 0));  # 840 kg/m3 is the highest conc before the flow is considered mudflow and mannings eq is no longer valid
  
    ##############################################
    ### Estimating soil particle detachment ######
    ##############################################
  
    ### transport capacity in kg/m3

	# set V to ero for channel to avoid erosion in channel  
    #V = if (chanmask eq 0, min(10,V),0); 

    rho_s = 2650; # specific density grains of sand kg/m3
    omega = Grad*max(0,V-0.004);  # unit stream power 
   
    TC = rho_s*Cg*(omega**Dg);  # transport capacity Govers et al. 1995
   
    Y = cover(if (Cohesion le 0, 0, 1/(2*Cohesion)),0)*mask; 
    # soil strength factor 0-1, from MMF
   
    # Y = cover(if (Cohesion le 0, 0, 0.79*exp(-0.85*Cohesion)),0)*mask; 
    # soil strength factor 0-1, from Eurosem
   
    Df  = max(TC-Conc,0)*Y*Vs*dt*W*dx;  # flow detachment in kg/celll:  kg/sec/m3 * m/s *m*m* sec = kg/cell  
    Df = if(chanmask > 1, 0, Df);

    # potential deposition in kg/cell (negative)   #kg/m3 * (-) * m3/cell  = kg/cell
    Dep = min(TC-Conc,0)*Vs*dt*W*dx;
    Dep = if (WH gt 0.0001, min(TC-Conc,0)*(1-exp(-dt*Vs/WH))*Vol, 0);  
   # Dep = min(0, max(Dep, -(Sed+Df)))*mask;     # not more deposition than there is sediment
    
    #exclude channel cells, these are not channel equations
    Df = if(chanmask eq 1,0,Df);   
    Dep = if(chanmask eq 1,0,Dep);
  
    # cumulative for display
    Dfcum = Dfcum + Df;    
    Depcum = Depcum + Dep;
    Dscum = Dscum + Ds;

    erosion = Df+Ds+Dep + erosion;         # cumulating soil loss  (detachment - deposition, is net loss)                       
  
    count = count + 1;

    
  } until count >= maxcount; #1 or more hours

  ### NOW THE REST OF THE WATER BALANCE ###

  #################################
  ### Actual Evapotranspiration ###
  #################################

  # actual evapotranspiration ETa linear with soil moisture content (mm)
  ETpoint = theta_wp + (theta_fc - theta_wp)*2/3;
  ETfactor = if (theta gt ETpoint, 1.0, 0.0);
  ETfactor = if (theta lt ETpoint and theta ge theta_wp,(theta-theta_wp)/(ETpoint-theta_wp), ETfactor);
  ETfactor = if (theta lt theta_wp, 0.0, ETfactor);
  #ETfactor = cover(ETfactor, 0.0001); #????
         
  Ta = ETp * ETfactor * Cover;                                # actual transpiration (mm)
  Ea = ETp * theta/theta_s * (1-Cover);                       # actual soil evaporation (mm)
  ETa = Ea + Ta;                                              # ETa sum of the Evap and Transp  
  #ETa = if(landuse eq RICE or chanmask eq 1, ETp, ETa);                        # ETa equals ETp in case of water body
  ETa = min(ETa, SoilMoisture);                               # cannot be more than soil moisture present
 
  # graphs with average and cumulative average ETa of all cells
  ETacum = ETacum + ETa;
  ETFactor = ETa/ETp;  
  #ETfcum = ETfCum + ETFactor;
                                       # recalculating true ETFactor based on final ETa/ETp 
  ETE0fac = areaaverage(ETfactor, nLU);

  #################################
  ########### Percolation #########
  #################################
  
  Perc = 24.0*Ksat*(theta/theta_s)**(2+3/lambda);  
  Perc = if(theta lt theta_wp, 0, Perc);
  Perc = min(SoilMoisture, Perc);
  perccum=perccum+Perc;

  #################################################
  ########### Computing new soil moisture #########
  #################################################

  SoilMoisture = SoilMoisture + (Infilday - ETa - Perc);
  SoilMoisture = min(SoilMoisture, theta_s*SD1);
  SoilMoisture = max(0, SoilMoisture);                       # soil moisture cannot be negative
  theta = SoilMoisture/SD1*mask;                      # findng soil moisture at surface
  se = theta/theta_s;
  report sepoints.tss = timeoutput(nLU,se);

  #### reporting cumulative graphs (spatial average) 

  report etpcum_tss = maptotal(ETpcum)/nrCells;
  report  eta_tss = maptotal(ETa)/nrCells;                    
  report etacum_tss = maptotal(ETacum)/nrCells;
  report ETfact_tss = maptotal(ETfactor)/nrCells;           # computing time series ETFactor

  #report WHavg = WHavg/maxcount;
  # spatial average runoff in mm cumulative of the whole area
  roavg = maptotal(WHavg*1000)/nrCells;
  rocum = rocum + roavg; 
  report runoffLU.tss=timeoutput(nLU,rocum);  

  report infil.tss= timeoutput(nLU, infilcum); # cumulative infil per landuse   

  #report outflow_mm.tss=rocumout; # outflowmout of the area, NOT the runoff that is on the surface 

  #report Qavg = Qavg/maxcount;
  
  #report rofract.tss = rocumout/((maptotal(Pcum)+1)/nrCells);
 
  report percolation.tss = maptotal(perccum)/nrCells;

  # converting from kg/gridcell to ton/ha kg: kg to ton factor 1000 
  # gridcell area to ha: # kg/gridcell /cellarea() *10000 = kg/ha /1000 = ton/ha
  # so (10000/1000)/cellarea()
  factor = 10/cellarea();
  report ds.map=Dscum*factor;
  report df.map=Dfcum*factor;
  report dep.map=Depcum*factor;

  # soil loss is the net positive detachment in ton/ha, expressed per gridcell. Note that a gridcell is less than a ha 
  soilloss = max(0, erosion*factor);            
        
  report Soilloss.map = soilloss;
  report soillosslu.tss = timeoutput(nLU, soilloss); # soilloss per land unit ton/ha

  # classify with FAO classification
  erosrate = scalar(
    if(soilloss ge 0 and soilloss lt 2, 0, 
    if(soilloss ge 2 and soilloss lt 5, 1, 
    if(soilloss ge 5 and soilloss lt 10, 2, 
    if(soilloss ge 10 and soilloss lt 20, 3, 
    if(soilloss ge 20 and soilloss lt 50, 4, 
    if(soilloss ge 50 and soilloss lt 100, 5, 6)))))));
  report erosrate.map = ordinal(erosrate);
  erosrate = scalar(
    if(soilloss ge 0 and soilloss lt 5, 0, 
    if(soilloss ge 2 and soilloss lt 10, 1, 
    if(soilloss ge 10 and soilloss lt 20, 2, 
    if(soilloss ge 20 and soilloss lt 30, 3, 
    if(soilloss ge 30 and soilloss lt 40, 4, 
    if(soilloss ge 40 and soilloss lt 60, 5, 6)))))));
  report erosrate1.map = ordinal(erosrate);
  # set everything to zero for the next day, 
  # we assume there is no runoff that lasts until the ext day
  Sed = 0;
  Conc = 0;
  Vol = 0;
  WH = 0;  

  day = day + 1;
