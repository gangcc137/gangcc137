library(raster)
library(sp)
library(ncdf4)
library(dplyr)
library(gridExtra)
library(ggplot2)
setwd("D:\\TRENDYv9\\DLEM")
source ("D:/LC_CMIP6/kons_general_help.R") 
dat=read.csv("dat_npp_rh.csv")
dat2=read.csv("dat_S2.csv")
# fixme:
# Your parameters will most likely differ but you can still use the
# distinctions between different sets of parameters. The aim is to make
# the different tasks in the code more obvious. In principal you could
# have a lot of overlapping sets and just have to keep them consistent. 
# 'namedtuples' can be used just as normal tuples by functions
# that are not aware of the names. They can still use the positions like 
# in the original code

# @Kostia and the 'R'tists: 
# It is not necessary to replicate the complete functionality of the #
# namedtuple classes. A simple 'R' approximation is a list with named entries 
# pa=list(C_leaf_0=2,...)

# This set is used by the functions that produce the 
# specific ingredients (functions) that will be run by
# mcmc alg.
UnEstimatedParameters = list(
  C_leaf_0=0,
  C_sp_0=0,#C_wood_0=0
  # C_fr_0=0,#C_root_0=0,
  # C_cr_0=0,#C_litter_above_0=0,
  C_AOM1_0=0,#C_litter_below_0=0,
  #C_AOM2_0=0,#C_fastsom_0=0,
  C_SB1_0=0,#C_slowsom_0=0,
  C_SR_0=0,#C_passsom_0=0,
  C_SB2_0=0,#C_slowsom_0=0,
  C_NM_0=0,#C_passsom_0=0,
  C_DM_0=0,
  #C_PSM_0=0,
  npp=0,
  number_of_months=0
  # tsl=0, # soil temperature - for dynamic turnover rates
  # mrso=0, # soil moisture - for dynamic turnover rates
  # ts=0 # air temperature - for dynamic turnover rates
)

# This set is used (as argument) by the functions that are called
# inside the mcmc
EstimatedParameters = list(
  beta_leaf=0,    #  1 (indices uses in original code) 
  beta_sp=0,#beta_wood=0,    #  2
  f_leaf2aom1=0, #f_leaflit2fastsom=0,  #  3
  f_sp2aom1=0, #f_leaflit2slowsom=0,#  4
  f_fr2aom1=0, #f_leaflit2passsom=0,#  5
  f_cr2aom1=0 #f_woodlit2fastsom=0,  #  6
  # f_woodlit2slowsom=0,#  7
  # f_woodlit2passsom=0,#  8
  # f_rootlit2fastsom=0,  #  9
  # f_rootlit2slowsom=0,#  10
  # f_rootlit2passsom=0,#  11
  # k_leaf=0,       #  12
  # k_wood=0,       #  13
  # k_root=0,       #  14
  # k_leaflit=0,	#  15
  # k_woodlit=0,	#  16
  # k_rootlit=0,	#  17
  # k_fastsom=0,	#  18
  # k_slowsom=0,	# 19
  # k_passsom=0,	# 20
  # C_leaflit_0=0,	# 21
  # T0=0,	# 22
  # E=0,	# 23
  # KM=0  # 24
)

# This is the set off all 
Parameters=c(EstimatedParameters,UnEstimatedParameters)

# This set defines the order of the c pools
# The order is crucial for the compatibility
# with the matrices (B and b) If you change it
# the matrix changes
StateVariables = list(
  C_leaf=0,
  C_sp=0,#C_wood=0,
  C_froot=0,#C_root=0,
  C_croot=0,#C_leaflit=0,
  C_AOM1=0,#C_woodlit=0,
  C_AOM2=0,#C_rootlit=0,
  C_SMB1=0,
  C_SR=0,
  C_SMB2=0,
  C_NOM=0,
  C_DOM=0,
  C_PSOM=0
  # C_fastsom=0,
  # C_slowsom=0,
  # C_passsom=0
)
Observables = list(
  CVeg=0,
  #C_sp=0,#C_wood=0,
  #C_froot=0,#C_froot=0,
  #C_croot=0,#C_croot=0,
  # C_AOM1=0,#C_litter_above=0,
  # C_AOM2=0,#C_litter_below=0,    
  # C_SB1=0,
  # C_SR=0,
  # C_SB2=0,
  # C_NM=0,
  # C_DM=0,
  # C_PSM=0,
  # C_fastsom=0,
  # C_slowsom=0,
  # C_passsom=0,
  rh=0,
  cLitter=0,#litter pool
  cSoil=0 #soil pool
)

# We define another set of parameters which describes
# the parameters of the matrices A,K and the vector b
# and drivers like npp (in form of arrays)
# but does not include start values and hyperparameters like the 'number_of_months'
# This distinction is helpful for the forward simulation where the
# distinction between estimated and constant is irrelevant.
ModelParameters = Parameters[!( names(Parameters) %in% list(
  'C_leaf_0',
  'C_sp_0',#'C_root_0',
  'C_froot_0',#'C_wood_0',
  'C_croot_0',#'C_litter_above_0',
  'C_AOM1_0',#'C_litter_below_0',
  'C_AOM2_0',#'C_fastsom_0',
  'C_SMB1_0',#'C_slowsom_0',
  'C_SR_0',#'C_passsom_0',
  'C_SMB2_0',
  'C_NOM_0',
  'C_DOM_0',
  'C_PSOM_0',
  'rh_0',
  #'C_leaflit_0',
  'number_of_months')                    
)] 

# to do
# 1.) make a namedtuple for the yycable data or use xarray to create a multifile dataset

#    # Read NetCDF data  ******************************************************************************************************************************

get_example_site_vars<-function(dataPath, lon, lat){
  
  # pick up 1 site
  point<-SpatialPoints(cbind(lon,lat), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
  
  # load all netCDF files in the folder
  files<-list.files(dataPath, pattern="..nc") # list all netCDF files in the folder
  (var_names<-substr(files, 0 ,regexpr("_",files)-1)) # derive variable names from file names
  dat<-data.frame(point)
  
  for (i in 1:length(var_names)){
    r<-stack(paste0(dataPath,"/",files[i]))
    r_dat<-extract(r,point)
    dat<-cbind(dat,r_dat[1,])
  }
  names(dat)<-c("lon","lat", var_names) # assign variable names
  
  #calculate wood pool from veg, leaf and root
  dat$cWood=dat$cVeg-dat$cLeaf-dat$cRoot
  # correct fluxes from per second to per day
  dat$gpp<-dat$gpp*86400
  dat$npp<-dat$npp*86400
  dat$nppLeaf<-dat$nppLeaf*86400
  dat$nppWood<-dat$nppWood*86400
  dat$nppRoot<-dat$nppRoot*86400
  dat$fVegLitter<-dat$fVegLitter*86400
  dat$fLitterSoil<-dat$fLitterSoil*86400
  dat$rh<-dat$rh*86400
  
  # explore the data
  print("Importing data:")
  names(dat)
  #save the data to a file
  
  output=data.frame(
    C_leaf=dat$cLeaf,
    C_wood=dat$cWood,
    C_root=dat$cRoot,
    C_litter_above=dat$cLitterAbove,
    C_litter_below=dat$cLitterBelow,
    C_fastsom=dat$cSoilFast,
    C_slowsom=dat$cSoilMedium,
    C_passsom=dat$cSoilSlow,
    npp=dat$npp,
    rh=dat$rh,
    f_veg2litter=dat$fVegLitter,
    f_litter2som=dat$fLitterSoil,
    # tsl=dat$tsl, # soil temperature - for dynamic turnover rates
    # mrso=dat$mrsos, # soil moisture - for dynamic turnover rates
    # ts=dat$ts # air temperature - for dynamic turnover rates
  )
  write.csv(output,paste0(dataPath,"/dat.csv"))
  return(output)
  
} 
datapath="D:/TRENDYv9/DLEM"
# read data from a csv file if previously saved  
get_data_from_file<-function(dataPath){
  dat<-read.csv(paste0(dataPath,"/dat.csv")) 
  #output=list(dat$npp,dat$nppLeaf,dat$nppWood,dat$nppRoot,dat$fVegLitter,dat$fLitterSoil,dat$fLitterSoil, dat$rh,
  #            dat$cLeaf, dat$cVeg, dat$cRoot, dat$cLitterAbove, dat$cLitterBelow, dat$cSoilFast, dat$cSoilMedium, dat$cSoilSlow, 
  #            dat$tsl, dat$mrso, dat$ts)
  return(dat)
}

# function to check if mcmc-generated parameters can be accepted
make_param_filter_func<-function(C_max, C_min) {
  
  isQualified<-function(c){
    # fixme
    #   this function is model specific: It discards parameter proposals
    #   where beta1 + beta3 are > 1 and 
    #   where sum of transfer coefficients from the same pool > 1
    paramNum = length(c)
    flag = T
    for (i in 1:paramNum){
      if(c[i] > C_max[i] || c[i] < C_min[i]) {flag = F; break}
      if(c[1] + c[2] > 1){flag = F; break}
      # if(c[3] + c[4] + c[5] > 1){flag = F; break}
      # if(c[6] + c[7] + c[8] > 1){flag = F; break}
      if(c[10] + c[11] + c[12] + c[13] + c[14] > 1){flag = F; break}
    }
    return (flag)
  }
  return (isQualified)
}

make_weighted_cost_func<-function(obs) {
  
  costfunction<-function(out_simu){
    # fixme
    #   as indicated by the fact that the function lives in this
    #   model-specific module it is not apropriate for (all) other models.
    #   There are model specific properties:
    #   1.) The weight for the different observation streams
    #
    
    # we assume the model output to be in the same shape and order
    # as the observation
    # this convention has to be honored by the forward_simulation as well
    # which in this instance already compresses the 3 different litter pools
    # to C_litter and the 3 different soil pools to one
    
    J_obj1 = mean (( out_simu_annual$cVeg - obs$cVeg)**2)/(2*var(obs$cVeg))
    # J_obj2 = mean (( out_simu$C_wood - obs$C_wood )**2)/(2*var(obs$C_wood))
    # J_obj3 = mean (( out_simu$C_root - obs$C_root )**2)/(2*var(obs$C_root))
    # J_obj4 = mean (( out_simu$C_litter_above - obs$C_litter_above )**2)/(2*var(obs$C_litter_above))
    # J_obj5 = mean (( out_simu$C_litter_below - obs$C_litter_below )**2)/(2*var(obs$C_litter_below))
    # J_obj6 = mean (( out_simu$C_fastsom - obs$C_fastsom )**2)/(2*var(obs$C_fastsom))
    # J_obj7 = mean (( out_simu$C_slowsom - obs$C_slowsom )**2)/(2*var(obs$C_slowsom))
    # J_obj8 = mean (( out_simu$C_passsom - obs$C_passsom )**2)/(2*var(obs$C_passsom))
    J_obj2 = mean (( out_simu_annual$rh - obs$rh_year )**2)/(2*var(obs$rh_year))
    J_obj3 = mean (( out_simu_annual$cLitter - obs$cLitter )**2)/(2*var(obs$cLitter))
    J_obj4 = mean (( out_simu_annual$csoil - obs$cSoil )**2)/(2*var(obs$cSoil))
    
    J_new= (J_obj1 + J_obj2 + J_obj3 + J_obj4)/200
    
    return (J_new)
  }
  return (costfunction)
}

make_param2res<-function(cpa){
  # """The returned function 'param2res' is used by the metropolis hastings algorithm
  #   to produce a result from a proposed parameter(tuple).
  #   'param2res' will be called only with the parameters to be estimated 'pe'.
  #   The constant parameters are the parameters of the present function 'pc' and are automatically 
  #   available in the returned 'param2res' (closure).
  # 
  #   In the case of this model 'param2res' has to perform the following
  #   tasks.
  #   -   use both sets of parameters 'pc' and 'pe' to build a forward model
  #   -   run the forward model to produce output for the times
  #       when observations are available.(Here monthly)
  #   -   project the output of the forward model to the observed variables.
  #       (here summation of all different soil-pools to compare to the single
  #       observed soil-pool and the same for the litter pools)
  # 
  #   In this version of the function all tasks are performed at once as in
  #   the original script.
  #   See the alternative implementation to see how all the tasks can be
  #   delegated to smaller functions.
  #   """
  # fixme
  # This function is model-specific in several ways:
  # 0. Which are the fixed and which are the estimated variables
  # 1. The matrices for the forward simulation
  # 2. the driver (here monthly npp)
  # 3. the projection of the simulated pool values to observable variables
  #    (here summation of
  
  param2res<-function(epa){
    # pa is a numpy array when pa comes from the predictor
    # so we transform it to be able to use names instead of positions
    #epa=EstimatedParameters
    days = c(31,28,31,30,31,30,31,31,30,31,30,31)# Construct b vector
    
    # leaf, root, wood
    # beta1=epa$beta_leaf; beta2=epa$beta_wood 
    # beta3 = 1- beta1- beta2
    # B = c(beta1, beta2, beta3, 0, 0, 0, 0,0,0)   
    beta1=epa$beta_leaf; beta2=epa$beta_sp 
    beta3=0.6*(1- beta1- beta2)
    beta4=0.4*(1- beta1- beta2)
    B = c(beta1, beta2, beta3, beta4, 0, 0, 0,0,0,0,0,0) # allocation
    # transfer coefficients
    # f41=1; f52=1; f63=1
    # f74 = epa$f_leaflit2fastsom; f84=epa$f_leaflit2slowsom;  f94=epa$f_leaflit2passsom;
    # f75=epa$f_woodlit2fastsom; f85=epa$f_woodlit2slowsom;  f95=epa$f_woodlit2passsom;
    # f76=epa$f_rootlit2fastsom; f86=epa$f_rootlit2slowsom; f96=epa$f_rootlit2passsom
    f51 = epa$f_leaf2aom1;   f52 = epa$f_sp2aom1;        f53 = epa$f_fr2aom1; f54 = epa$f_cr2aom1; 
    f61 = 1-epa$f_leaf2aom1; f62 = epa$f_sp2aom2; f63 = 1-epa$f_fr2aom1; f64 = epa$f_cr2aom2;
    f75=0.04368; f76=0.0624; f78=0.8; f710=0.368; f711=0.5; f712=0.4;
    f89=0.5;
    f95=0.09632; f96=0.1376;
    f105=0.21; f107=0.6435;f1011=0.3
    f115=0.013; f116=0.64;f1110=0.00015; 
    f127=0.0065; f1210=0.032;
    A = c(-1,  0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  #leaf
          0,  -1,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  #sw
          0,   0,  -1,   0,   0,   0,   0,   0,   0,  0,  0,  0,  #fr
          0,   0,   0,  -1,   0,   0,   0,   0,   0,  0,  0,  0,  #cr
          f51, f52,f53, f54, -1,   0,   0,   0,   0,  0,  0,  0,  #AOM1
          f61, f62,f63, f64,  0,  -1,   0,   0,   0,  0,  0,  0,  #AOM2
          0,   0,   0,   0,  f75, f76,  -1, f78,  0,f710,f711,f712,  #SMB1
          0,   0,   0,   0,   0,   0,   0,  -1,  f89, 0,  0,  0,  #SMR
          0,   0,   0,   0,   f95, f96, 0,   0,  -1,  0,  0,  0,  #SMB2
          0,   0,   0,   0,   f105,0,  f107, 0,   0, -1,f1011,0,  #NOM
          0,   0,   0,   0,   f115,f116,0,   0,   0,f1110,-1, 0,  #DOM
          0,   0,   0,   0,   0,   0,  f127, 0,   0,f1210,0, -1)/12 #PSOM
    
    A = matrix(A, nrow = 12, byrow = TRUE)
    
    #turnover rate per day of pools:
    #temp = c(epa$k_leaf,epa$k_wood,epa$k_root, epa$k_leaflit, epa$k_woodlit, epa$k_rootlit, epa$k_fastsom, epa$k_slowsom, epa$k_passsom)
    clay=0.2028
    temp = c(0.001826,0.000137,0.002740, 0.0000685,0.000608, 0.000148, 0.000458/clay, 
             0.000375, 0.000458,0.010959,0.958904,0.002192)
    # K matrix
    # K = rep(0,81)
    # K = matrix(K, nrow = 9, byrow = TRUE)
    # for (i in 1:9) { K[i,i] = temp[i] }
    K = rep(0,144)
    K = matrix(K, nrow = 12, byrow = TRUE)
    for (i in 1:12) { K[i,i] = temp[i] }
    
    # X vector
    x=rep(0,cpa$number_of_months)
    x_fin=data.frame(x,x,x,x,x,x,x,x,x,x,x,x); 
    #names(x_fin)=c("leaf","wood","root","leaflit","woodlit","rootlit","fastsom","slowsom","passsom")
    names(x_fin)=c("leaf","sp","froot","croot","AOM1","AOM2","SMB1","SR","SMB2","NOM","DOM","PSOM")
    
    rh_fin=rep(0,cpa$number_of_months);   
    f_veg_lit_fin=rep(0,cpa$number_of_months);
    f_lit_soil_fin=rep(0,cpa$number_of_months)
    
    # "leaf","wood","root","leaf_litter","wood_litter","root_litter","soil_fast","soil_slow","soil_passive"
    x_init = c(epa$C_leaf_0, # leaf
               epa$C_sp_0,   #sapwood
               cpa$C_froot_0,   #froot
               cpa$C_croot_0,   #croot
               epa$C_AOM1_0, #litter-AOM1
               cpa$C_AOM2_0, #litter-AOM2
               epa$C_SMB1_0,  #soil-SMB1
               epa$C_SR_0,   #soil-SR
               epa$C_SMB2_0,  #soil-SMB2
               epa$C_NOM_0,   #soil-NOM
               epa$C_DOM_0,   #soil-DOM
               cpa$C_PSOM_0)  #soil-PSOM
               # cpa$C_wood_0, # stem
               # cpa$C_root_0, # root
               # epa$C_leaflit_0, # leaf litter - unknown
               # cpa$C_litter_above_0[1]-epa$C_leaflit_0, # wood litter
               # cpa$C_litter_below_0, # root litter
               # cpa$C_fastsom_0, # soil active
               # cpa$C_slowsom_0, # soil intermediate
               # cpa$C_passsom_0)  # soil passive  
    X=x_init   # initialize carbon pools 
    
    # initialize first respiration value
    #co2_rh=cpa.rh_0
    # fixme:
    # slight change to the original
    # I would like to start the solution with the initial values
    # m=0 means after 0 moths = in the initial step
    #B=A@K
    #pa=Parameters.from_EstimatedParametersAndUnEstimatedParameters(epa,cpa)
    #B=make_compartmental_matrix_func(pa)(0,X)
    
    ############################ Additional functions specific for VISITe ##############
    
    # # modifier for leaf pool turnover rate for deciduous forest depending on air temp
    # leaf_fall<-function (TA)
    # {
    #   lf=0
    #   if (TA<281) {lf<-1}
    #   return (lf)
    # }
    # # modifier for soil pools turnover rate to implement time-dependent soil respiration
    # Rh_calc<-function (TS, M, T0, E, KM)
    # {
    #   TS<-TS-273.15
    #   if (TS>T0) {rh_out=exp(E*(1/(10-T0)-1/(TS-T0))) *M/(KM+M)}  else {rh_out=0}
    #   return(rh_out)
    # }
    ######################################################################################
    
    for (m in 1:cpa$number_of_months){
      npp_in = cpa$npp[m] 
      co2_rh = 0; 
      f_veg_lit = 0;
      f_lit_soil = 0
      
      # environmental factor ksi: modifies leaf fall and turnover rates / respiration 
      # depending on monthly surface t "ts", soil t "tsl" and soil moisture data "mrso"
      # rh_modifier - VISIT-specific - remove if not relevant for your model
      # rh_modifier=Rh_calc(cpa$tsl[m], cpa$mrso[m], epa$T0, epa$E, epa$KM)
      # # ksi vector (remove if no environmental modifiers)
      # ksi=c(leaf_fall(cpa$ts[m]), # leaf pool
      #       1, # wood pool
      #       1, # root pool
      #       rh_modifier, # leaf liter
      #       rh_modifier, # wood litter
      #       rh_modifier, # root litter
      #       rh_modifier, # fast soil
      #       rh_modifier, # slow soil
      #       rh_modifier) # passive soil
      
      for (d in 1:days[m%%12+1]) {
        
        # matrix equation with ksi (remove ksi if no environmental modifiers)
        X = X + B * npp_in + A %*% K %*% X #* ksi
        # deriving rh from each litter and soil pool as 1 - sum of transfer coefficients 
        # co2_rate = c(0,0,0, 
        #              (1-f74-f84-f94)*K[4,4]*ksi[4],
        #              (1-f75-f85-f95)*K[5,5]*ksi[5],
        #              (1-f76-f86-f96)*K[6,6]*ksi[6], 
        #              K[7,7]*ksi[7], 
        #              K[8,8]*ksi[8], 
        #              K[9,9]*ksi[9])
        co2_rate = c(0,0,0,0,
                     (1-f75-f95-f105-f115)*K[5,5],
                     (1-f76-f96-f116)*K[6,6],
                     (1-f107-f127)*K[7,7],
                     (1-f78)*K[8,8], 
                     (1-f89)*K[9,9],
                     (1-f710-f1110-f1210)*K[10,10],
                     (1-f711)*K[11,11],
                     (1-f712)*K[12,12])
        co2=sum(co2_rate*X)
        co2_rh = co2_rh + co2/days[m%%12+1]   # monthly average rh
        # deriving litterfall
        # litterfall_rate = c(f41*K[1,1]*ksi[1],f52*K[2,2]*ksi[2],f63*K[3,3]*ksi[3],0,0,0,0,0,0)
        litterfall_rate = c((f51+f61)*K[1,1],(f52+f62)*K[2,2],(f53+f63)*K[3,3],(f54+f64)*K[4,4],0,0,0,0,0,0,0,0)
        litterfall=sum(litterfall_rate*X)
        f_veg_lit=f_veg_lit+litterfall/days[m%%12+1]
        # deriving humus formation
        # litter_to_soil_rate = c(0,0,0,f74*K[4,4]*ksi[4]+f84*K[4,4]*ksi[4]+f84*K[4,4]*ksi[4],
        #                         f75*K[5,5]*ksi[5]+ f85*K[5,5]*ksi[5]+ f95*K[5,5]*ksi[5],
        #                         f76*K[6,6]*ksi[6]+ f86*K[6,6]*ksi[6]+ f96*K[6,6]*ksi[6],
        #                         0,0,0)
        litter_to_soil_rate = c(0,0,0,0,
                              (f75+f95+f105+f115)*K[5,5],
                              (f76+f96+116)*K[6,6],
                              (f107+f127)*K[7,7],
                              f78*K[8,8],
                              f89*K[9,9],
                              (f710+f1110+f1210)*K[10,10],
                              (f711+f1011)*K[11,11],
                              f712*K[10,10])
        litter_to_soil=sum(litter_to_soil_rate*X)
        f_lit_soil=f_lit_soil+litter_to_soil/days[m%%12+1]
      }
      x_fin[m,]=X
      rh_fin[m]=co2_rh
      f_veg_lit_fin[m]=f_veg_lit
      f_lit_soil_fin[m]=f_lit_soil
    }
    
    # We create an output that has the same shape
    # as the observations to make the costfunctions
    # easier.
    # To this end we project our 10 output variables of the matrix simulation
    # onto the 6 data streams by summing up the 3 litter pools into one
    # and also the 3 soil pools into one
    # out_simu = list(
    #   C_leaf=x_fin$leaf,
    #   C_wood=x_fin$wood,
    #   C_root=x_fin$root,
    #   C_litter_above=x_fin$leaflit+x_fin$woodlit,
    #   C_litter_below=x_fin$rootlit,
    #   C_fastsom=x_fin$fastsom,
    #   C_slowsom=x_fin$slowsom,
    #   C_passsom=x_fin$passsom,
    #   rh=rh_fin,
    #   f_veg2litter=f_veg_lit_fin,
    #   f_litter2som=f_lit_soil_fin)
    #"leaf","sp","froot","croot","AOM1","AOM2","SMB1","SR","SMB2","NOM","DOM","PSOM"
      month_pool = list(
        # C_leaf=x_fin$leaf,
        # C_sp=x_fin$sp,
        # C_froot=x_fin$froot,
        # C_croot=x_fin$croot,
        # C_AOM1=x_fin$AOM1,
        # C_AOM2=x_fin$AOM2,
        # C_SMB1=x_fin$SMB1,
        # C_SR=x_fin$SR,
        # C_SMB2=x_fin$SMB2,
        # C_NOM=x_fin$NOM,
        # C_DOM=x_fin$DOM,
        # C_PSOM=x_fin$PSOM,
        C_veg=x_fin$leaf+x_fin$sp+x_fin$froot+x_fin$croot,
        C_litter=x_fin$AOM1+x_fin$AOM2,
        C_soil=x_fin$SMB1+x_fin$SR+x_fin$SMB2+x_fin$NOM+x_fin$DOM+x_fin$PSOM,
        rh=rh_fin,
        f_veg2litter=f_veg_lit_fin,
        f_litter2soc=f_lit_soil_fin
    )
      #annual Veg pool
      C_veg_annual=NULL
      for (i in 0:319) {
        numveg=sum(month_pool$C_veg[(12*i+1):(12*i+12)])
        C_veg_annual=c(C_veg_annual,numveg)
      }
      #annual litter pool
      C_litter_annual=NULL
      for (i in 0:319) {
        numlit=sum(month_pool$C_litter[(12*i+1):(12*i+12)])
        C_litter_annual=c(C_litter_annual,numlit)
      }
      #annual soil pool
      C_soil_annual=NULL
      for (i in 0:319) {
        numsoil=sum(month_pool$C_soil[(12*i+1):(12*i+12)])
        C_soil_annual=c( C_soil_annual,numsoil)
      }

      #annual rh
      rh_annual=NULL
      for (i in 0:319) {
        numrh=sum(month_pool$rh[(12*i+1):(12*i+12)])
        rh_annual=c(rh_annual,numrh)
      }

      out_simu=list(
        cVeg=C_veg_annual,
        rh=rh_annual,
        cLitter=C_litter_annual,
        csoil=C_soil_annual)
      
    return (out_simu)
  }
  return (param2res)
}
################################################################################
# # alternative implementation not yet translated to R
