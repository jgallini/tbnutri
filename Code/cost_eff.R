
# This function runs cost effectiveness analyses for a nutritional intervention 
# in TB patients in India. The inputs are iterations (any integer), 
# and a set of logical (TRUE/FALSE) inputs for whether or not to hold a 
# particular parameter constant. These options should be used for 
# deterministic sensitivity analyses. The default for all of these inputs is 
# FALSE, so each parameter will be varied according to the appropriate 
# distribution if no option is specified.

cost_eff<-function(iterations,tb_death_c=FALSE,
                   oth_death_c=FALSE, failed_c=FALSE, ltf_c=FALSE,
                   ltf_relap_c=FALSE,
                   earpt_relap_c=FALSE,ltpt_relap_c=FALSE,
                   food_cost_c=FALSE,tb_cost_c=FALSE,tb_daly_c=FALSE,
                   posttb_daly_c=FALSE,oth_death_mult_c=FALSE){
  
  #specifying required libraries
  library(lhs)
  library(mc2d)
  
  #setting seed and number of parameters that will vary
  set.seed(84937)
  params<-12 
  
  #Setting up using with Latin hypercube sampling to get the matrix of control 
  #probabilities.
  #The iteration we're on will correspond to the row of this matrix
  A <- randomLHS(iterations,params)
  B <- matrix(nrow = iterations, ncol = params) 
  
  #control parameters
  #TB death
  p1<-0.128 #6 month probability from Pranay (1.29.25 email)
  r1<-(-log(1-p1)/6) #6 month rate
  mu1<-1-exp(-r1) #1 month probability
  var1<-0.01**2 #just using variance of 6 month distribution to 
  #allow for more variability than the 1 month would
  alpha1 <- ((1 - mu1) / var1 - 1 / mu1) * mu1 ^ 2
  beta1 <- alpha1 * (1 / mu1 - 1)
  percs1<-qbeta(c(0.025, 0.975),alpha1,beta1)
  B[,1] <- qbeta(A[,1], shape1 = alpha1, shape2= beta1) 
  #Other (background) death
  p3<-0.023669 #5 year (60 month) probability of background death from life 
  #table ages 45-49
  r2<-(-log(1-p3)/60) #60 month rate
  p4<-1-exp(-r2) #1 month probability: 0.000399
  B[,2] <- qpert(A[,2],mean=0.000399,max=1.2*0.000399,min=0.8*0.000399) 
  #taking point estimate plus/minus 20%
  percs2<-qpert(c(0.025,0.975),mean=p4,max=1.2*p4,min=0.8*p4)
  #LTF
  p2<-0.052 #6 month probability from Pranay (1.29.25 email)
  r2<-(-log(1-p2)/6) #6 month rate
  mu2<-1-exp(-r2) #1 month probability
  var2<-0.003**2 #using 6 month variance to increase variability
  alpha2 <- ((1 - mu2) / var2 - 1 / mu2) * mu2 ^ 2
  beta2 <- alpha2 * (1 / mu2 - 1)
  percs3<-qbeta(c(0.025, 0.975),shape1=alpha2,shape2=beta2)
  B[,3] <- qbeta(A[,3], shape1=alpha2,shape2=beta2) 
  #Food cost
  B[,4]<- qgamma(A[,4],shape=23,scale=1) 
  percs4<-qgamma(c(0.025,0.975),shape=23,scale=1)
  #changed after meeting 8.19.24 to
  #mean of 23 and var of 23 from mean of 43 based on Pranay's calculations
  #TB cost
  B[,5]<- qgamma(A[,5],shape=25,scale=0.8)
  percs5<-qgamma(c(0.025, 0.975),shape=25,scale=0.8)
  #TB DALY value
  B[,6]<- qpert(A[,6],min=(0.224/12),max=(0.454/12),mean=(0.333/12))
  percs6<-qpert(c(0.025, 0.975),min=(0.224/12),max=(0.454/12),mean=(0.333/12))
  #Long term increased DALY for former TB patients
  B[,7]<- qgamma(A[,7],shape=100,scale=0.000044)
  percs7<-qgamma(c(0.025, 0.975),shape=100,scale=0.000044)
  #failed: distribution from Pranay email 1/29/25
  p8<-0.004 #6 month probability
  r8<-(-log(1-p8)/6) #6 month rate
  mu8<-1-exp(-r8) #1 month probability
  var8<-0.001**2 #using 6 month variance to increase variability
  alpha8 <- ((1 - mu8) / var8 - 1 / mu8) * mu8 ^ 2
  beta8 <- alpha8 * (1 / mu8 - 1)
  percs8<-qbeta(c(0.025, 0.975),shape1 = alpha8, shape2= beta8)
  B[,8] <- qbeta(A[,8], shape1 = alpha8, shape2= beta8)
  #ltf_relap (from LTF to relapse)
  p9<-((158/(158+1210))*1.61) #12 month beta distribution from Pranay: 
  #(2.1.24 email)
  r9<-(-log(1-p9)/12) #12 month rate
  mu9<-1-exp(-r9) #1 month probability
  var9<-(158*1210)/(((158+1210)**2)*(158+1210+1))/(36*(1.61**2)) 
  #converting variance of 6 month distribution to 
  #1 month by dividing by 36 and also 1.61 multiplier squared
  alpha9 <- ((1 - mu9) / var9 - 1 / mu9) * mu9 ^ 2
  beta9 <- alpha9 * (1 / mu9 - 1)
  percs9<-qbeta(c(0.025, 0.975),shape1 = alpha9, shape2= beta9)
  B[,9] <- qbeta(A[,9], shape1 = alpha9, shape2= beta9)
  #earpt_relap (from early post treatment to relapse)
  p10<-(158/(158+1210)) #12 month beta distribution from Pranay (1.30.24 email)
  r10<-(-log(1-p10)/12) #12 month rate
  mu10<-1-exp(-r10) #1 month probability
  var10<-(158*1210)/(((158+1210)**2)*(158+1210+1))/(36) 
  #converting variance of 6 month distribution to 
  #1 month by dividing by 36
  alpha10 <- ((1 - mu10) / var10 - 1 / mu10) * mu10 ^ 2
  beta10 <- alpha10 * (1 / mu10 - 1)
  percs10<-qbeta(c(0.025, 0.975),shape1 = alpha10, shape2= beta10)
  B[,10] <- qbeta(A[,10], shape1 = alpha10, shape2= beta10)
  #ltpt_relap (from late post treatment to relapse)
  B[,11]<-qpert(A[,11],min=((169/100000)/12),max=((231/100000)/12),
                mean=((199/100000)/12))
  percs11<-qpert(c(0.025, 0.975),min=((169/100000)/12),max=((231/100000)/12),
                 mean=((199/100000)/12))
  #other death multiplier (Pranay email 10.25.24)
  B[,12]<-exp(qpert(A[,12],min=log(1.02),max=log(1.34),mean=log(1.14)))
  
  #naming columns
  colnames(B)<-c("tb_death","oth_death","ltf",
                 "food_cost","tb_cost","tb_daly","posttb_daly",
                 "failed","ltf_relap","earpt_relap","ltpt_relap",
                 "oth_death_mult")
  
  #This section takes into account the user inputs for which variables should be 
  #held constant. Parameters are overwritten with their theoretical mean for all 
  #values if it's specified that parameter be held constant. If LTF is specified 
  #constant all months are constant.
  
  if (tb_death_c==TRUE) {B[,"tb_death"]<-mu1}
  if (oth_death_c==TRUE) {B[,"oth_death"]<-0.000399}
  if (ltf_c==TRUE) {B[,"ltf"]<-mu2}
  if (food_cost_c==TRUE) {B[,"food_cost"]<-23}
  if (tb_cost_c==TRUE) {B[,"tb_cost"]<-25*0.8}
  if (tb_daly_c==TRUE) {B[,"tb_daly"]<-0.333/12}
  if (posttb_daly_c==TRUE) {B[,"posttb_daly"]<-(100*0.000044)}
  if (failed_c==TRUE) {B[,"failed"]<-mu8}
  if (ltf_relap_c==TRUE) {B[,"ltf_relap"]<-mu9}
  if (earpt_relap_c==TRUE) {B[,"earpt_relap"]<-mu10}
  if (ltpt_relap_c==TRUE) {B[,"ltpt_relap"]<-(199/100000)/12}
  if (oth_death_mult_c==TRUE) {B[,"oth_death_mult"]<-1.14}
  
  
  #Now using Latin hypercube sampling to create a matrix of Treatment Inputs.
  
  #Treatment
  C<-randomLHS(iterations,params)
  D<-matrix(nrow = iterations, ncol = params) 
  
  #TB death: from Pranay email 10.30.24
  #draw random risk ratio from beta distribution to keep less than 1
  rr<-qbeta(C[,1],shape1=4,shape2=2.25)
  #rrpercs1<-qbeta(c(0.025,0.975),shape1=4,shape2=2.25)
  #multiply risk ratios by control parameters to get probability 
  #of TB death for treatment group
  D[,1]<-rr*B[,1] 
  #Other death: same as control other death
  D[,2] <- B[,2]
  #LTF, now using risk ratio from Pranay email 10.15.24, multiplying by
  #control parameter
  D[,3] <- exp(qpert(C[,3],min=log(0.28),max=log(0.77),mean=log(0.45)))*B[,3] 
  #rrpercs2<-exp(qpert(c(0.025, 0.975),min=log(0.28),max=log(0.77),
  #mean=log(0.45)))
  #Food cost (same as control)
  D[,4]<- B[,4]
  #TB cost (same as control)
  D[,5]<- B[,5]
  #TB DALY value (same as control)
  D[,6]<- B[,6]
  #Long term increased DALY for former TB patients (same as control)
  D[,7]<- B[,7]
  #failed, now using risk ratio from email Pranay 10.15.24, multiplying by
  #control group
  D[,8]<-exp(qpert(C[,8],min=log(0.03),max=log(0.4),mean=log(0.1)))*B[,8]
  #rrpercs3<-exp(qpert(c(0.025, 0.975),min=log(0.03),max=log(0.4),mean=log(0.1)))
  #ltf_relap (same as control group)
  D[,9] <- B[,9]
  #earpt_relap (same as control group)
  D[,10] <- B[,10]
  #ltpt_relap (same as control group)
  D[,11]<-B[,11]
  #oth_death_mult (same as control group)
  D[,12]<-B[,12]
  
  #naming columns
  colnames(D)<-c("tb_death","oth_death","ltf",
                 "food_cost","tb_cost","tb_daly","posttb_daly",
                 "failed","ltf_relap","earpt_relap","ltpt_relap",
                 "oth_death_mult")
  
  
  #This section takes into account the user inputs for which variables should be 
  #held constant. Parameters are overwritten with their theoretical mean for all 
  #values if it's specified that parameter be held constant. 
  
  if (tb_death_c==TRUE) {D[,"tb_death"]<-mu1*0.64}
  if (oth_death_c==TRUE) {D[,"oth_death"]<-0.000399}
  if (ltf_c==TRUE) {D[,"ltf"]<-mu2*0.45}
  if (food_cost_c==TRUE) {D[,"food_cost"]<-23}
  if (tb_cost_c==TRUE) {D[,"tb_cost"]<-25*0.8}
  if (tb_daly_c==TRUE) {D[,"tb_daly"]<-0.333/12}
  if (posttb_daly_c==TRUE) {D[,"posttb_daly"]<-(100*0.000044)}
  if (failed_c==TRUE) {D[,"failed"]<-mu8*0.1}
  if (ltf_relap_c==TRUE) {D[,"ltf_relap"]<-mu9}
  if (earpt_relap_c==TRUE) {D[,"earpt_relap"]<-mu10}
  if (ltpt_relap_c==TRUE) {D[,"ltpt_relap"]<-(199/100000)/12}
  if (oth_death_mult_c==TRUE) {D[,"oth_death_mult"]<-1.14}
  
  #storing the final inputs in a matrix
  Control_Inputs<-B
  Treatment_Inputs<-D
  
  #initiating vectors needed in loop
  cost_per_daly<-vector()
  Control_DALYS<-vector()
  Treatment_DALYS<-vector()
  ctrl_csum<-vector()
  trt_csum<-vector()
  cost_dif<-vector()
  daly_av2<-vector()

  Subjects_trt<-list()
  Subjects_cont<-list()
  
  # This is a function that will iterate through a certain number of times
  for(i in 1:iterations) {
      
    # now calling custom function that creates transition matrices for treatment
    # and control groups
    source(here("Code/Transition_function_wide.R"))
    Subjects_trt[[i]]<-TB_state_wide(Treatment_Inputs[i,1],
                                     Treatment_Inputs[i,2],
                                Treatment_Inputs[i,3],Treatment_Inputs[i,3],
                                Treatment_Inputs[i,3],Treatment_Inputs[i,3],
                                Treatment_Inputs[i,3],Treatment_Inputs[i,3], 
                                #same LTF value for all 6 months of treatment
                                Treatment_Inputs[i,8],Treatment_Inputs[i,9],
                                Treatment_Inputs[i,10],Treatment_Inputs[i,11],
                                Treatment_Inputs[i,12])
    
    Subjects_cont[[i]]<-TB_state_wide(Control_Inputs[i,1],
                                      Control_Inputs[i,2],
                                 Control_Inputs[i,3],Control_Inputs[i,3],
                                 Control_Inputs[i,3],Control_Inputs[i,3],
                                 Control_Inputs[i,3],Control_Inputs[i,3],
                                 #same LTF value for all 6 months of treatment
                                 Control_Inputs[i,8],Control_Inputs[i,9],
                                 Control_Inputs[i,10],Control_Inputs[i,11],
                                 Control_Inputs[i,12])
    
    #subsetting into relevant states for control active TB DALY calculation
    Control_d<-subset(Subjects_cont[[i]],select=-c(die_other_cause_tb,
                                                   die_other_cause_posttb,
                                                   die_tb,ltfu,
                                                   early_posttrt,late_posttrt))
    
    #subsetting into relevant states for control post TB DALY calculation
    Control_postd<-subset(Subjects_cont[[i]],select=c(early_posttrt,
                                                      late_posttrt,ltfu))
    
    #summing control TB DALY's by month
    s1<-Control_Inputs[i,"tb_daly"]*Control_d #DALY's from having active TB
    s2<-Control_Inputs[i,"posttb_daly"]*Control_postd #DALY's post TB
    sum_d<-as.vector(rowSums(s1))+as.vector(rowSums(s2))
    
    #control background deaths, selecting out states that are relevant
    #ignoring background deaths that occurred while they had TB
    Control_backded<-subset(Subjects_cont[[i]],select=c(die_other_cause_posttb))
    
    #control tb deaths, selecting out states that are relevant
    Control_tbded<-subset(Subjects_cont[[i]],select=c(die_tb))
    
    #summing control tb deaths 
    s12tb<-Control_tbded*(1/12)
    #summing control post tb deaths scaled by 0.14 since that's the excess
    #death beyond what we'd usually expect
    s12post<-Control_backded*(1/12)*(0.14)
    ded2<-as.vector(rowSums(s12tb))+as.vector(rowSums(s12post))
    
    #setting discounting rate
    rate<-0.03
    
    #summing DALY's by month
    dvec<-ded2+sum_d
    
    #splits months into groups by year
    dcheck1<-split(dvec,ceiling(seq_along(dvec)/12))
    
    #summing DALY's by year
    l1<-list()
    for (j in 1:length(dcheck1)){
      l1[j]<-sum(dcheck1[[j]])
    }
    u1<-unlist(l1)
    
    #applying discounting
    u2<-list()
    for (j in 1:length(u1)){
      u2[j]<-u1[j]/((1+rate)**j)
    }
    Control_DALYS[i]<-sum(unlist(u2))
    
    
    ####### treatment DALY's
    
    #subsetting into relevant states for treatment DALY calculation for active 
    #TB
    Trt_d<-subset(Subjects_trt[[i]],select=-c(die_other_cause_tb,
                                              die_other_cause_posttb,
                                              die_tb,ltfu,early_posttrt,
                                              late_posttrt))
    
    #subsetting into relevant states for post TB DALY's
    Trt_dp<-subset(Subjects_trt[[i]],select=c(early_posttrt,late_posttrt,ltfu))
    
    #summing treatment TB DALY's by month
    s1t<-Treatment_Inputs[i,"tb_daly"]*Trt_d #active TB DALY's
    s2t<-Treatment_Inputs[i,"posttb_daly"]*Trt_dp #post TB DALY's
    sum_dt<-as.vector(rowSums(s1t))+as.vector(rowSums(s2t))
    
    #treatment background deaths, selecting out states that are relevant
    #ignoring background deaths that occur during TB states
    Treatment_dedpost<-subset(Subjects_trt[[i]],
                              select=c(die_other_cause_posttb))
    
    #treatment tb deaths, selecting states that are relevant
    Treatment_dedtb<-subset(Subjects_trt[[i]],select=c(die_tb))
    
    #summing treatment deaths
    s12ttb<-(1/12)*Treatment_dedtb
    s12tpost<-(1/12)*Treatment_dedpost*(0.14)
    #again scaling by excess death post tb
    sum_d2t<-as.vector(rowSums(s12ttb))+as.vector(rowSums(s12tpost))
    
    #setting discounting rate
    rate<-0.03
    
    #summing DALY's by month
    dvect<-sum_d2t+sum_dt
    
    #splits months into groups by year
    dcheck1t<-split(dvect,ceiling(seq_along(dvect)/12))
    
    #summing DALY's by year
    l1<-list()
    for (j in 1:length(dcheck1t)){
      l1[j]<-sum(dcheck1t[[j]])
    }
    u1t<-unlist(l1)
    
    #applying discounting
    u2t<-list()
    for (j in 1:length(u1t)){
      u2t[j]<-u1t[j]/((1+rate)**j)
    }
    Treatment_DALYS[i]<-sum(unlist(u2t))
    
    #Now calculating total DALY's averted.
    #total DALY's averted
    daly_av2[i]<-Control_DALYS[i]-Treatment_DALYS[i]
    #daly_av3<-sort(daly_av2)
    
    # Next the total costs for the control group and the intervention group are 
    # calculated. The cost of TB treatment was summed over patients in active 
    # treatment states over the 60 year period. This was done for both control 
    # and intervention groups. Additionally for the intervention group, costs of 
    # food were summed over patients in active treatment states over the 60 year
    # period. Discounting was applied at a 3% rate for both intervention and 
    # control. The costs of treatment and food were summed for intervention 
    # patients. Then, the difference in costs was subtracted between the 
    # intervention and control groups.
    
    #control costs
    #selecting out states that are relevant
    Control_c<-subset(Subjects_cont[[i]],select=-c(die_other_cause_tb,
                                                   die_other_cause_posttb,
                                                   die_tb,ltfu,failed,
                                                   early_posttrt,late_posttrt,
                                              relapse))
    
    #summing control costs by month
    s1tc<-Control_Inputs[i,"tb_cost"]*Control_c #active TB DALY's
    sum_s1tc<-as.vector(rowSums(s1tc))
    
    #splits months into groups by year
    dcheck1c<-split(sum_s1tc,ceiling(seq_along(sum_s1tc)/12))
    
    #summing costs by year
    l1<-list()
    for (j in 1:length(dcheck1c)){
      l1[j]<-sum(dcheck1c[[j]])
    }
    u1c<-unlist(l1)
    
    #applying discounting
    u2c<-list()
    for (j in 1:length(u1c)){
      u2c[j]<-u1c[j]/((1+rate)**j)
    }
    ctrl_csum[i]<-sum(unlist(u2c))
    
    #treatment costs
    #selecting out states that are relevant
    Intervention_c<-subset(Subjects_trt[[i]],select = -c(die_other_cause_tb,
                                                         die_other_cause_posttb,
                                                         die_tb,ltfu,failed,
                                                    early_posttrt,late_posttrt,
                                                    relapse))
    
    total_cost<-Treatment_Inputs[i,"food_cost"]+Treatment_Inputs[i,"tb_cost"]
    
    #summing treatment costs by month
    s1tt<-total_cost*Intervention_c 
    sum_dtt<-as.vector(rowSums(s1tt))
    
    #splits months into groups by year
    dcheck1tt<-split(sum_dtt,ceiling(seq_along(sum_dtt)/12))
    
    #summing costs by year
    l1t<-list()
    for (j in 1:length(dcheck1tt)){
      l1t[j]<-sum(dcheck1tt[[j]])
    }
    u1tt<-unlist(l1t)
    
    #applying discounting
    u2tt<-list()
    for (j in 1:length(u1tt)){
      u2tt[j]<-u1tt[j]/((1+rate)**j)
    }
    trt_csum[i]<-sum(unlist(u2tt))
    
    #cost difference
    cost_dif[i]<-trt_csum[i]-ctrl_csum[i]

    #cost per DALY
    cost_per_daly[i]<-cost_dif[i]/daly_av2[i]
    
      #print(paste(i))
  }

  #outputting data we need to use in results file
  return(list(cost_per_daly=cost_per_daly,Control_DALYS=Control_DALYS,
              Treatment_DALYS=Treatment_DALYS,ctrl_csum=ctrl_csum,
              trt_csum=trt_csum,Control_Inputs=Control_Inputs,
              Treatment_Inputs=Treatment_Inputs,
              Subjects_cont=Subjects_cont,Subjects_trt=Subjects_trt))

}

