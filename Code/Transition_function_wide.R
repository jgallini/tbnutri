

# Creating function for Multi-State Markov Model for TB nutritional 
# supplementation

# This function allows inputs for probability of: 
# LTF from TB tx Months 1-6, death from TB, 
# death from other causes, treatment failure, relapsing from LTF state, 
# relapsing from early post-treatment state, and relapsing from late 
# post-treatment state

# It defaults to 10,000 subjects all entering at the beginning and 
# runs for 60 years (720 months)

# This version of the function creates a wide 720x20 data frame with months 
# as the rows and 
# all 20 states as the columns

# Coder: Julia Gallini. PI: Pranay Sinha.

TB_state_wide<-function(tb_death,oth_death,
                        ltf1,ltf2,ltf3,ltf4,ltf5,ltf6,
                        failed,ltf_relap,earpt_relap,ltpt_relap,
                        oth_death_mult){
  
  #forming columns
  On_TB_tx_M1<-c(rep(0,20))
  On_TB_tx_M2<-c(-1,rep(0,19))
  On_TB_tx_M3<-c(0,-1,rep(0,18))
  On_TB_tx_M4<-c(0,0,-1,rep(0,17))
  On_TB_tx_M5<-c(0,0,0,-1,rep(0,16))
  #adding in increased risk of death for cured patients with other death 
  #multiplier
  On_TB_tx_M6<-c(0,0,0,0,-1,rep(0,15)) 
                     0,rep(oth_death,5),0)
  die_other_cause_posttb<-c(rep(0,7),1,0,0,0,(oth_death_mult*oth_death),
                            (oth_death_mult*oth_death),
                            rep(0,7))
  die_tb<-c(rep(tb_death,6),0,0,1,0,0,0,0,0,rep(tb_death,5),0)
  ltfu<-c(ltf1,ltf2,ltf3,ltf4,ltf5,ltf6,0,0,0,-1,rep(0,10))
  failed<-c(rep(0,5),failed,rep(0,14))
  early_posttrt<-c(rep(0,5),-1,0,0,0,0,0,-1,rep(0,7),1) 
  late_posttrt<-c(rep(0,11),-1,-1,rep(0,7))
  #cured<-c(rep(0,5),-1,0,0,0,0,-1,rep(0,6),1)
  relapse<-c(rep(0,9),ltf_relap,0,earpt_relap,ltpt_relap,rep(0,7))
  repeat_tb_tx_M1<-c(rep(0,10),1,0,0,1,rep(0,6))
  repeat_tb_tx_M2<-c(rep(0,14,),-1,rep(0,5))
  repeat_tb_tx_M3<-c(rep(0,15,),-1,rep(0,4))
  repeat_tb_tx_M4<-c(rep(0,16,),-1,rep(0,3))
  repeat_tb_tx_M5<-c(rep(0,17,),-1,0,0)
  repeat_tb_tx_M6<-c(rep(0,18,),-1,0)
  
  #creating matrix
  mat<-cbind(On_TB_tx_M1,On_TB_tx_M2,On_TB_tx_M3,On_TB_tx_M4,On_TB_tx_M5,
             On_TB_tx_M6,die_other_cause_tb,die_other_cause_posttb,
             die_tb,ltfu,failed,early_posttrt,late_posttrt,relapse,
             repeat_tb_tx_M1,repeat_tb_tx_M2,repeat_tb_tx_M3,repeat_tb_tx_M4,
             repeat_tb_tx_M5,repeat_tb_tx_M6)
  
  #solving for values that are based on other values
  mat[1,2]<- 1-sum(mat[1,-2])
  mat[2,3]<- 1-sum(mat[2,-3])
  mat[3,4]<- 1-sum(mat[3,-4])
  mat[4,5]<- 1-sum(mat[4,-5])
  mat[5,6]<- 1-sum(mat[5,-6])
  mat[10,10]<- 1-sum(mat[10,-10])
  mat[6,12]<- 1-sum(mat[6,-12])
  mat[12,12]<- (11/12)*(1-sum(mat[12,c(-12,-13)]))
  mat[12,13]<- (1/12)*(1-sum(mat[12,c(-12,-13)]))
  mat[13,13]<- 1-sum(mat[13,-13])
  mat[15,16]<- 1-sum(mat[15,-16])
  mat[16,17]<- 1-sum(mat[16,-17])
  mat[17,18]<- 1-sum(mat[17,-18])
  mat[18,19]<- 1-sum(mat[18,-19])
  mat[19,20]<- 1-sum(mat[19,-20])
  
  
  #initiating matrix of 10,000 people to start in TB treatment month 1
  Subjects_M1<-t(as.vector(c(10000,rep(0,19))))
  
  #looping through any number of months of time, creates "Subjects" table
  Subjects<-data.frame(Subjects_M1)
  rownames(Subjects)[1]<-"Month_1"
  for (i in 1:719){
    ind<-ceiling(i/60)
    mat2<-mat
    
    #different oth_death option for each age range
    mat2[,"die_other_cause_tb"]<-(1.50**(ind-1))*mat2[,"die_other_cause_tb"]
    mat2[,"die_other_cause_posttb"]<-
      (1.50**(ind-1))*mat2[,"die_other_cause_posttb"]
    #increasing risk of background death by 50% every 5 years 
    #since we're starting at age 45
    
    mat2[,"die_other_cause_tb"]<-ifelse(mat2[,"die_other_cause_tb"]>1,
                                        1,mat2[,"die_other_cause_tb"])
    mat2[,"die_other_cause_posttb"]<-ifelse(mat2[,"die_other_cause_posttb"]>1,
                                            1,mat2[,"die_other_cause_posttb"])
    #capping transition probabilities at 1
    
    #now we need to adjust all of the other probabilities that depend 
    # on other death so that the total probabilities add to 1
    mat2[1,2]<- 1-sum(mat2[1,-2])
    mat2[2,3]<- 1-sum(mat2[2,-3])
    mat2[3,4]<- 1-sum(mat2[3,-4])
    mat2[4,5]<- 1-sum(mat2[4,-5])
    mat2[5,6]<- 1-sum(mat2[5,-6])
    mat2[10,10]<- 1-sum(mat2[10,-10])
    mat2[6,12]<- 1-sum(mat2[6,-12])
    mat2[12,12]<- (11/12)*(1-sum(mat2[12,c(-12,-13)]))
    mat2[12,13]<- (1/12)*(1-sum(mat2[12,c(-12,-13)]))
    mat2[13,13]<- 1-sum(mat2[13,-13])
    mat2[15,16]<- 1-sum(mat2[15,-16])
    mat2[16,17]<- 1-sum(mat2[16,-17])
    mat2[17,18]<- 1-sum(mat2[17,-18])
    mat2[18,19]<- 1-sum(mat2[18,-19])
    mat2[19,20]<- 1-sum(mat2[19,-20])
    
    #now doing matrix calculation
    Subjects[i+1,]<-data.matrix(Subjects[i,])%*%data.matrix(mat2)

    rownames(Subjects)[i+1]<-paste0("Month","_",i+1)
    colnames(Subjects)<-colnames(mat)
    }  
  
  #outputting resulting matrix
  Subjects
  
}









