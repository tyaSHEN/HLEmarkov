load("output/dem.Rdata")
rm("Raw")
library(tidyverse)
library(doParallel)
library(nnet)

#### DDMSLT ####
co = c("26","36") 
# CUT = 0 # 5 or 0
iteration = 500
for (cut in c(0,3,5)) {
  
  for (CO in co) {
    Age = YAW %>% filter(FY == paste0("19",CO)) %>% pull(FA)
    a = Age+5
    # a = Age+cut
    
    registerDoParallel(min(detectCores(),14))
    BASELINE = foreach(iter=c(1:iteration),.packages= c("tidyverse"),.combine = rbind) %dopar% {
      
      imp = sample(c(1:10),1)
      load(paste0("output/Sex/ALL",CO," (duration",imp,").Rdata"))
      
      
      C_CON.D_cut = C_CON.D
      # #constrain the sample for all
      # C_CON.D_cut$unknown = ifelse(C_CON.D_cut$unknown <CUT,NA,C_CON.D_cut$unknown)
      #
      C_CON.D_cut$duration = ifelse(C_CON.D_cut$duration >=cut,cut,C_CON.D_cut$duration)
      C_CON.D_cut$duration = ifelse(!is.na(C_CON.D_cut$unknown) & C_CON.D_cut$unknown >=cut,cut,C_CON.D_cut$duration)
      C_CON.D_cut$ADL.D = paste0(C_CON.D_cut$ADL.D,C_CON.D_cut$duration)
      C_CON.D_cut$ADL.D = ifelse(C_CON.D_cut$ADL.D =="HNA","H",C_CON.D_cut$ADL.D)
      C_CON.D_cut$ADL.D = ifelse(C_CON.D_cut$ADL.D %in% c("ANA","LNA","NANA"),NA,C_CON.D_cut$ADL.D)
      C_CON.D_cut$pre_ADL.D = lag(C_CON.D_cut$ADL.D)
      
      if(iter==1){
        tem = data.frame(rahhidpn = unique(C_CON.D_cut$rahhidpn))
      }else{
        tem = data.frame(rahhidpn = sample(unique(C_CON.D_cut$rahhidpn),length(unique(C_CON.D_cut$rahhidpn)),replace = T))
      }
      C_CON = left_join(tem,C_CON.D_cut)
      
      C_CON$pre_ADL.D = factor(C_CON$pre_ADL.D,levels = c("A0","A1","A2","A3","A4","A5","A6","A7","A8","A9","L0","L1","L2","L3","L4","L5","L6","L7","L8","L9"))
      C_CON$pre_ADL.D = factor(C_CON$pre_ADL.D,levels = unique(sort(C_CON$pre_ADL.D)))
      
      C_CON$ADL.D = factor(C_CON$ADL.D,levels = c("A0","A1","A2","A3","A4","A5","A6","A7","A8","A9","L0","L1","L2","L3","L4","L5","L6","L7","L8","L9","H"))
      C_CON$ADL.D = factor(C_CON$ADL.D,levels = unique(sort(C_CON$ADL.D)))
      
      tab = C_CON %>% filter(!is.na(state),!is.na(pre_state),!is.na(pre_ADL.D)) %>% mutate(x= paste0(substring(pre_ADL.D,0,1),substring(ADL.D,0,1)))
      nrow(tab)
      table(tab$x)
      mean(tab$age)
      table(tab$Sex)/sum(table(tab$Sex))*100
      table(tab$race)/sum(table(tab$race))*100
      table(tab$edu)/sum(table(tab$race))*100
      
      
      BL = C_CON %>% filter(age %in% c((a-5):(a+4)),wt >0,!is.na(ADL.D)) %>% group_by(rahhidpn) %>% arrange(age) %>% slice_head(n=1)
      Baseline = left_join(tem,BL) %>% filter(!is.na(wt),ADL.D != "H")
      Baseline %>% group_by(Sex) %>% summarise(tot = sum(wt)) %>% filter(!is.na(Sex)) %>% mutate(pro = tot/sum(tot))
      Baseline %>% group_by(race) %>% summarise(tot = sum(wt)) %>% filter(!is.na(race)) %>% mutate(pro = tot/sum(tot))
      Baseline %>% group_by(edu) %>% summarise(tot = sum(wt)) %>% filter(!is.na(edu)) %>% mutate(pro = tot/sum(tot))
      tab = BL %>% filter(!is.na(wt),ADL.D != "H")
      table(tab$Sex)/sum(table(tab$Sex))*100
      table(tab$race)/sum(table(tab$race))*100
      table(tab$edu)/sum(table(tab$edu))*100
      table(tab$ADL.D)
      sum(table(tab$ADL.D))
      
      # Baseline = left_join(tem,BL) %>% filter(!is.na(wt),ADL.D != "H")
      Baseline = Baseline %>% group_by(ADL.D,Sex) %>% summarise(tot = sum(wt)) %>% filter(!is.na(Sex)&!is.na(ADL.D)) %>% mutate(ADL.D = as.character(ADL.D))
      Baseline = Baseline %>% ungroup() %>% mutate(pro = tot/sum(tot)) %>% select(-tot)
      Baseline$iter = iter
      Baseline
    }
    stopImplicitCluster()
    write_csv(BASELINE,paste0("output/Sex/Res/BASELINE.",CO,".",cut,".csv"))
    
    registerDoParallel(min(detectCores(),14))
    PROB = foreach(iter=c(1:iteration),.packages= c("tidyverse","nnet"),.combine = rbind) %dopar% {
      imp = sample(c(1:10),1)
      load(paste0("output/Sex/ALL",CO," (duration",imp,").Rdata"))
      
      C_CON.D_cut = C_CON.D
      #constrain the sample for all
      # C_CON.D_cut$unknown = ifelse(C_CON.D_cut$unknown <CUT,NA,C_CON.D_cut$unknown)
      #
      C_CON.D_cut$duration = ifelse(C_CON.D_cut$duration >=cut,cut,C_CON.D_cut$duration)
      C_CON.D_cut$duration = ifelse(!is.na(C_CON.D_cut$unknown) & C_CON.D_cut$unknown >=cut,cut,C_CON.D_cut$duration)
      C_CON.D_cut$ADL.D = paste0(C_CON.D_cut$ADL.D,C_CON.D_cut$duration)
      C_CON.D_cut$ADL.D = ifelse(C_CON.D_cut$ADL.D =="HNA","H",C_CON.D_cut$ADL.D)
      C_CON.D_cut$ADL.D = ifelse(C_CON.D_cut$ADL.D %in% c("ANA","LNA","NANA"),NA,C_CON.D_cut$ADL.D)
      C_CON.D_cut$pre_ADL.D = lag(C_CON.D_cut$ADL.D)
      
      if(iter==1){
        tem = data.frame(rahhidpn = unique(C_CON.D_cut$rahhidpn))
      }else{
        tem = data.frame(rahhidpn = sample(unique(C_CON.D_cut$rahhidpn),length(unique(C_CON.D_cut$rahhidpn)),replace = T))
      }
      C_CON = left_join(tem,C_CON.D_cut)
      
      C_CON$pre_ADL.D = factor(C_CON$pre_ADL.D,levels = c("A0","A1","A2","A3","A4","A5","A6","A7","A8","A9","L0","L1","L2","L3","L4","L5","L6","L7","L8","L9"))
      C_CON$pre_ADL.D = factor(C_CON$pre_ADL.D,levels = unique(sort(C_CON$pre_ADL.D)))
      
      C_CON$ADL.D = factor(C_CON$ADL.D,levels = c("A0","A1","A2","A3","A4","A5","A6","A7","A8","A9","L0","L1","L2","L3","L4","L5","L6","L7","L8","L9","H"))
      C_CON$ADL.D = factor(C_CON$ADL.D,levels = unique(sort(C_CON$ADL.D)))
      
      tmsa = multinom(formula = ADL.D ~ pre_ADL.D + age + agesq + as.factor(Sex) + as.factor(Sex):age,
                      weights = wt*w, data = C_CON %>% mutate(agesq = age^2) %>% filter(age>=50),na.action=na.omit,maxit = 1000)
      
      dwrite <- crossing(pre_ADL.D = unique(sort(C_CON$pre_ADL.D)), Sex = c("1","2"),age=c(55:105))
      dwrite$agesq = dwrite$age ^ 2
      
      pp.s <- cbind(dwrite, predict(tmsa, newdata = dwrite, type = "probs", se = TRUE))
      lpp.s <- pivot_longer(pp.s,c(5:ncol(pp.s)))
      lpp.s$value = ifelse(lpp.s$value<0.0005,NA,lpp.s$value)
      lpp.s = lpp.s %>% group_by(pre_ADL.D,Sex,age) %>% mutate(value = value/sum(value,na.rm = T)) %>% filter(!is.na(value))
      lpp.s$iter = iter
      lpp.s
    }
    stopImplicitCluster()
    write_csv(PROB,paste0("output/Sex/Res/PROB",CO,".",cut,".csv"))
    
    registerDoParallel(min(detectCores(),14))
    PROB = foreach(iter=c(1:iteration),.packages= c("tidyverse","nnet"),.combine = rbind) %dopar% {
      imp = sample(c(1:10),1)
      load(paste0("output/Sex/ALL",CO," (duration",imp,").Rdata"))
      
      C_CON.D_cut = C_CON.D
      #constrain the sample for all
      # C_CON.D_cut$unknown = ifelse(C_CON.D_cut$unknown <CUT,NA,C_CON.D_cut$unknown)
      #
      C_CON.D_cut$duration = ifelse(C_CON.D_cut$duration >=cut,cut,C_CON.D_cut$duration)
      C_CON.D_cut$duration = ifelse(!is.na(C_CON.D_cut$unknown) & C_CON.D_cut$unknown >=cut,cut,C_CON.D_cut$duration)
      C_CON.D_cut$ADL.D = paste0(C_CON.D_cut$ADL.D,C_CON.D_cut$duration)
      C_CON.D_cut$ADL.D = ifelse(C_CON.D_cut$ADL.D =="HNA","H",C_CON.D_cut$ADL.D)
      C_CON.D_cut$ADL.D = ifelse(C_CON.D_cut$ADL.D %in% c("ANA","LNA","NANA"),NA,C_CON.D_cut$ADL.D)
      C_CON.D_cut$pre_ADL.D = lag(C_CON.D_cut$ADL.D)
      
      if(iter==1){
        tem = data.frame(rahhidpn = unique(C_CON.D_cut$rahhidpn))
      }else{
        tem = data.frame(rahhidpn = sample(unique(C_CON.D_cut$rahhidpn),length(unique(C_CON.D_cut$rahhidpn)),replace = T))
      }
      C_CON = left_join(tem,C_CON.D_cut)
      
      C_CON$pre_ADL.D = factor(C_CON$pre_ADL.D,levels = c("A0","A1","A2","A3","A4","A5","A6","A7","A8","A9","L0","L1","L2","L3","L4","L5","L6","L7","L8","L9"))
      C_CON$pre_ADL.D = factor(C_CON$pre_ADL.D,levels = unique(sort(C_CON$pre_ADL.D)))
      
      C_CON$ADL.D = factor(C_CON$ADL.D,levels = c("A0","A1","A2","A3","A4","A5","A6","A7","A8","A9","L0","L1","L2","L3","L4","L5","L6","L7","L8","L9","H"))
      C_CON$ADL.D = factor(C_CON$ADL.D,levels = unique(sort(C_CON$ADL.D)))
      
      tmsa = multinom(formula = state ~ pre_state + age + agesq + as.factor(Sex) + as.factor(Sex):age,
                      weights = wt*w, data = C_CON %>% mutate(agesq = age^2) %>% filter(age>=50,!is.na(pre_ADL.D)),na.action=na.omit,maxit = 1000)
      
      dwrite <- crossing(pre_state = unique(sort(C_CON$pre_state)), Sex = c("1","2"),age=c(55:105))
      dwrite$agesq = dwrite$age ^ 2
      
      pp.s <- cbind(dwrite, predict(tmsa, newdata = dwrite, type = "probs", se = TRUE))
      lpp.s <- pivot_longer(pp.s,c(5:ncol(pp.s)))
      lpp.s$value = ifelse(lpp.s$value<0.0005,NA,lpp.s$value)
      lpp.s = lpp.s %>% group_by(pre_state,Sex,age) %>% mutate(value = value/sum(value,na.rm = T)) %>% filter(!is.na(value))
      lpp.s$iter = iter
      colnames(lpp.s)[1] = "pre_ADL.D"
      lpp.s
    }
    stopImplicitCluster()
    write_csv(PROB,paste0("output/Sex/Res/PROB",CO,".",cut," (Markov).csv"))
  }
  
}


#### semi-markov ####
co = c("26","36") 
# cut = 2
iteration = 500
for (CO in co) {
  Age = YAW %>% filter(FY == paste0("19",CO)) %>% pull(FA)
  a = Age+5
  
  registerDoParallel(min(detectCores(),14))
  PROB = foreach(iter=c(1:iteration),.packages= c("tidyverse"),.combine = rbind) %dopar% {
    imp = sample(c(1:10),1)
    load(paste0("output/Sex/ALL",CO," (duration",imp,").Rdata"))
    C_CON.D = C_CON.D %>% group_by(rahhidpn) %>% mutate(pre_ADL.D = lag(ADL.D),DUR = lag(duration))
    C_CON.D = C_CON.D %>% select(-ragender,-raracem,-rahispan)
    # 
    if(iter==1){
      tem = data.frame(rahhidpn = unique(C_CON.D$rahhidpn))
    }else{
      tem = data.frame(rahhidpn = sample(unique(C_CON.D$rahhidpn),length(unique(C_CON.D$rahhidpn)),replace = T))
    }
    C_CON = left_join(tem,C_CON.D)
    
    C_CON$pre_ADL.D = factor(C_CON$pre_ADL.D,levels = c("A","L"))
    C_CON$pre_ADL.D = factor(C_CON$pre_ADL.D,levels = unique(sort(C_CON$pre_ADL.D)))
    # 
    C_CON$ADL.D = factor(C_CON$ADL.D,levels = c("A","L","H"))
    C_CON$ADL.D = factor(C_CON$ADL.D,levels = unique(sort(C_CON$ADL.D)))
    # 
    BL = C_CON %>% filter(age %in% c((a-5):(a+4)),wt >0, !is.na(duration)) %>% group_by(rahhidpn) %>% arrange(age) %>% slice_head(n=1)
    Baseline = left_join(tem,BL) %>% filter(!is.na(wt),ADL.D != "H")
    tab = BL %>% filter(!is.na(wt),ADL.D != "H")
    table(tab$Sex)/sum(table(tab$Sex))*100
    table(tab$race)/sum(table(tab$race))*100
    table(tab$edu)/sum(table(tab$edu))*100
    table(tab$ADL.D,tab$duration)
    sum(table(tab$ADL.D))
    
    Baseline = Baseline %>% group_by(ADL.D,duration,Sex) %>% summarise(tot = sum(wt))%>% filter(!is.na(Sex)&!is.na(ADL.D)) %>% mutate(ADL.D = as.character(ADL.D))
    Baseline = Baseline %>% ungroup() %>% mutate(pro = tot/sum(tot)) %>% select(-tot)
    Baseline$iter = iter
    Baseline
    
  }
  stopImplicitCluster()
  write_csv(PROB,paste0("output/Sex/Res/BASELINE.",CO,".SEM.csv"))
  
  
  registerDoParallel(min(detectCores(),14))
  PROB = foreach(iter=c(1:iteration),.packages= c("tidyverse","nnet"),.combine = rbind) %dopar% {
    imp = sample(c(1:10),1)
    load(paste0("output/Sex/ALL",CO," (duration",imp,").Rdata"))
    C_CON.D = C_CON.D %>% group_by(rahhidpn) %>% mutate(pre_ADL.D = lag(ADL.D),DUR = lag(duration))
    C_CON.D = C_CON.D %>% select(-ragender,-raracem,-rahispan)
    # 
    if(iter==1){
      tem = data.frame(rahhidpn = unique(C_CON.D$rahhidpn))
    }else{
      tem = data.frame(rahhidpn = sample(unique(C_CON.D$rahhidpn),length(unique(C_CON.D$rahhidpn)),replace = T))
    }
    C_CON = left_join(tem,C_CON.D)
    
    C_CON$pre_ADL.D = factor(C_CON$pre_ADL.D,levels = c("A","L"))
    C_CON$pre_ADL.D = factor(C_CON$pre_ADL.D,levels = unique(sort(C_CON$pre_ADL.D)))
    # 
    C_CON$ADL.D = factor(C_CON$ADL.D,levels = c("A","L","H"))
    C_CON$ADL.D = factor(C_CON$ADL.D,levels = unique(sort(C_CON$ADL.D)))
    # 
    tab = C_CON %>% filter(!is.na(state),!is.na(pre_state),!is.na(DUR)) %>% mutate(x= paste0(pre_ADL.D,ADL.D))
    nrow(tab)
    table(tab$x)
    mean(tab$age)
    table(tab$Sex)/sum(table(tab$Sex))
    table(tab$race)/sum(table(tab$race))*100
    table(tab$edu)/sum(table(tab$race))*100
    
    tmsa = multinom(formula = ADL.D ~ pre_ADL.D + age + agesq +DUR:age + pre_ADL.D:DUR + DUR + DURsq + as.factor(Sex) + as.factor(Sex):DUR+ as.factor(Sex):age,
                    weights = wt*w, data = C_CON %>% filter(age>=50,!is.na(DUR)) %>% mutate(age = age-DUR,agesq = age^2,DURsq = DUR^2) ,na.action=na.omit,maxit = 1000)
    
    dwrite <- crossing(pre_ADL.D = unique(sort(C_CON$pre_ADL.D)), Sex = c("1","2"), DUR = c(0:15),age=c(55:105))
    dwrite$agesq = dwrite$age ^ 2
    dwrite$DURsq = dwrite$DUR ^ 2
    
    # pp.s <- cbind(dwrite, predict(tmsa, newdata = dwrite, "response"))
    pp.s <- cbind(dwrite, predict(tmsa, newdata = dwrite, type = "probs", se = TRUE))
    lpp.s <- pivot_longer(pp.s,c(7:ncol(pp.s)))
    lpp.s = lpp.s %>% group_by(pre_ADL.D,DUR,Sex,age) %>% mutate(value = value/sum(value,na.rm = T)) %>% filter(!is.na(value))
    lpp.s$iter = iter
    lpp.s
  }
  stopImplicitCluster()
  PROB$duration = ifelse(PROB$pre_ADL.D==PROB$name,PROB$DUR+1,0)
  PROB$age = PROB$age+PROB$DUR
  PROB$State1 = paste0(PROB$pre_ADL.D,PROB$DUR)
  PROB$State2 = paste0(PROB$name,PROB$duration)
  PROB$State2 = ifelse(PROB$State2=="H0","H",PROB$State2)
  
  write_csv(PROB,paste0("output/Sex/Res/PROB",CO,".SEM (new).csv"))
  
  registerDoParallel(min(detectCores(),14))
  PROB = foreach(iter=c(1:iteration),.packages= c("tidyverse","nnet"),.combine = rbind) %dopar% {
    imp = sample(c(1:10),1)
    load(paste0("output/Sex/ALL",CO," (duration",imp,").Rdata"))
    C_CON.D = C_CON.D %>% group_by(rahhidpn) %>% mutate(pre_ADL.D = lag(ADL.D),DUR = lag(duration))
    C_CON.D = C_CON.D %>% select(-ragender,-raracem,-rahispan)
    # 
    if(iter==1){
      tem = data.frame(rahhidpn = unique(C_CON.D$rahhidpn))
    }else{
      tem = data.frame(rahhidpn = sample(unique(C_CON.D$rahhidpn),length(unique(C_CON.D$rahhidpn)),replace = T))
    }
    C_CON = left_join(tem,C_CON.D)
    
    C_CON$pre_ADL.D = factor(C_CON$pre_ADL.D,levels = c("A","L"))
    C_CON$pre_ADL.D = factor(C_CON$pre_ADL.D,levels = unique(sort(C_CON$pre_ADL.D)))
    # 
    C_CON$ADL.D = factor(C_CON$ADL.D,levels = c("A","L","H"))
    C_CON$ADL.D = factor(C_CON$ADL.D,levels = unique(sort(C_CON$ADL.D)))
    # 
    tab = C_CON %>% filter(!is.na(state),!is.na(pre_state),!is.na(DUR)) %>% mutate(x= paste0(pre_ADL.D,ADL.D))
    nrow(tab)
    table(tab$x)
    mean(tab$age)
    table(tab$Sex)/sum(table(tab$Sex))
    table(tab$race)/sum(table(tab$race))*100
    table(tab$edu)/sum(table(tab$race))*100

    tmsa = multinom(formula = ADL.D ~ pre_ADL.D + age + agesq +as.factor(Sex) + as.factor(Sex):age,
                    weights = wt*w, data = C_CON %>% filter(age>=50,!is.na(DUR)) %>% mutate(agesq = age^2) ,na.action=na.omit,maxit = 1000)
    
    
    dwrite <- crossing(pre_ADL.D = unique(sort(C_CON$pre_ADL.D)), Sex = c("1","2"),age=c(55:105))
    dwrite$agesq = dwrite$age ^ 2
    
    # pp.s <- cbind(dwrite, predict(tmsa, newdata = dwrite, "response"))
    pp.s <- cbind(dwrite, predict(tmsa, newdata = dwrite, type = "probs", se = TRUE))
    lpp.s <- pivot_longer(pp.s,c(5:ncol(pp.s)))
    lpp.s = lpp.s %>% group_by(pre_ADL.D,Sex,age) %>% mutate(value = value/sum(value,na.rm = T)) %>% filter(!is.na(value))
    lpp.s$iter = iter
    lpp.s
  }
  stopImplicitCluster()
  # PROB$duration = ifelse(PROB$pre_ADL.D==PROB$name,PROB$DUR+1,0)
  PROB$State1 = paste0(PROB$pre_ADL.D,"0")
  PROB$State2 = paste0(PROB$name,"0")
  PROB$State2 = ifelse(PROB$State2=="H0","H",PROB$State2)
  
  write_csv(PROB,paste0("output/Sex/Res/PROB",CO,".SEM (Markov).csv"))
}


Full = C_CON %>% filter(!is.na(pre_ADL.D),!is.na(ADL.D))
Semi = C_CON %>% filter(!is.na(DUR),!is.na(duration))
nrow(Semi)/nrow(Full)

# Semi-M is older
mean(Semi$age)
mean(Full$age)
# Semi-M has more women
Full %>% filter(!duplicated(rahhidpn)) %>% group_by(sex) %>% summarise(c=n()) %>% mutate(c/sum(c))
Semi %>% filter(!duplicated(rahhidpn)) %>% group_by(sex) %>% summarise(c=n()) %>% mutate(c/sum(c))
# Semi-M has slightly lower education (likely because of more turbulence in life)
Full %>% filter(!duplicated(rahhidpn)) %>% group_by(edu) %>% summarise(c=n()) %>% mutate(c/sum(c))
Semi %>% filter(!duplicated(rahhidpn)) %>% group_by(edu) %>% summarise(c=n()) %>% mutate(c/sum(c))


#### DDMSLT (history) ####
co = c("26","36") 
# CUT = 0
iteration = 500
for (cut in c(3,5)) {
  
  for (CO in co) {
    Age = YAW %>% filter(FY == paste0("19",CO)) %>% pull(FA)
    a = Age+5
    # a = Age+cut
    
    registerDoParallel(min(detectCores(),14))
    BASELINE = foreach(iter=c(1:iteration),.packages= c("tidyverse"),.combine = rbind) %dopar% {
      
      imp = sample(c(1:10),1)
      load(paste0("output/Sex/ALL",CO," (duration",imp,").Rdata"))
      
      
      C_CON.D_cut = C_CON.D
      #constrain the sample for all
      # C_CON.D_cut$unknown = ifelse(C_CON.D_cut$unknown <CUT,NA,C_CON.D_cut$unknown)
      #
      C_CON.D_cut$duration = ifelse(C_CON.D_cut$duration >=cut,cut,C_CON.D_cut$duration)
      C_CON.D_cut$duration = ifelse(!is.na(C_CON.D_cut$unknown) & C_CON.D_cut$unknown >=cut,"X",C_CON.D_cut$duration)
      C_CON.D_cut$duration = ifelse(C_CON.D$ADL.D == "L" & C_CON.D_cut$duration == "X", cut , C_CON.D_cut$duration)
      C_CON.D_cut$ADL.D = paste0(C_CON.D_cut$ADL.D,C_CON.D_cut$duration)
      C_CON.D_cut$ADL.D = ifelse(C_CON.D_cut$ADL.D =="HNA","H",C_CON.D_cut$ADL.D)
      C_CON.D_cut$ADL.D = ifelse(C_CON.D_cut$ADL.D %in% c("ANA","LNA","NANA"),NA,C_CON.D_cut$ADL.D)
      C_CON.D_cut$pre_ADL.D = lag(C_CON.D_cut$ADL.D)
      
      if(iter==1){
        tem = data.frame(rahhidpn = unique(C_CON.D_cut$rahhidpn))
      }else{
        tem = data.frame(rahhidpn = sample(unique(C_CON.D_cut$rahhidpn),length(unique(C_CON.D_cut$rahhidpn)),replace = T))
      }
      C_CON = left_join(tem,C_CON.D_cut)
      
      C_CON$pre_ADL.D = factor(C_CON$pre_ADL.D,levels = c("A0","A1","A2","A3","A4","A5","A6","A7","A8","A9","AX","L0","L1","L2","L3","L4","L5","L6","L7","L8","L9"))
      C_CON$pre_ADL.D = factor(C_CON$pre_ADL.D,levels = unique(sort(C_CON$pre_ADL.D)))
      
      C_CON$ADL.D = factor(C_CON$ADL.D,levels = c("A0","A1","A2","A3","A4","A5","A6","A7","A8","A9","AX","L0","L1","L2","L3","L4","L5","L6","L7","L8","L9","H"))
      C_CON$ADL.D = factor(C_CON$ADL.D,levels = unique(sort(C_CON$ADL.D)))
      
      BL = C_CON %>% filter(age %in% c((a-5):(a+4)),wt >0,!is.na(ADL.D)) %>% group_by(rahhidpn) %>% arrange(age) %>% slice_head(n=1)
      
      Baseline = left_join(tem,BL) %>% filter(!is.na(wt),ADL.D != "H")
      Baseline %>% group_by(Sex) %>% summarise(tot = sum(wt)) %>% filter(!is.na(Sex)) %>% mutate(pro = tot/sum(tot))
      Baseline %>% group_by(race) %>% summarise(tot = sum(wt)) %>% filter(!is.na(race)) %>% mutate(pro = tot/sum(tot))
      Baseline %>% group_by(edu) %>% summarise(tot = sum(wt)) %>% filter(!is.na(edu)) %>% mutate(pro = tot/sum(tot))
      tab = BL %>% filter(!is.na(wt),ADL.D != "H")
      table(tab$Sex)/sum(table(tab$Sex))*100
      table(tab$race)/sum(table(tab$race))*100
      table(tab$edu)/sum(table(tab$edu))*100
      table(tab$ADL.D)
      sum(table(tab$ADL.D))
      
      # Baseline = left_join(tem,BL) %>% filter(!is.na(wt),ADL.D != "H")
      Baseline = Baseline %>% group_by(ADL.D,Sex) %>% summarise(tot = sum(wt))%>% filter(!is.na(Sex)&!is.na(ADL.D)) %>% mutate(ADL.D = as.character(ADL.D))
      Baseline = Baseline %>% ungroup() %>% mutate(pro = tot/sum(tot)) %>% select(-tot)
      Baseline$iter = iter
      Baseline
    }
    stopImplicitCluster()
    write_csv(BASELINE,paste0("output/Sex/Res/BASELINE.",CO,".",cut," (his).csv"))
    
    registerDoParallel(min(detectCores(),14))
    PROB = foreach(iter=c(1:iteration),.packages= c("tidyverse","nnet"),.combine = rbind) %dopar% {
      imp = sample(c(1:10),1)
      load(paste0("output/Sex/ALL",CO," (duration",imp,").Rdata"))
      
      C_CON.D_cut = C_CON.D
      #constrain the sample for all
      # C_CON.D_cut$unknown = ifelse(C_CON.D_cut$unknown <CUT,NA,C_CON.D_cut$unknown)
      #
      C_CON.D_cut$duration = ifelse(C_CON.D_cut$duration >=cut,cut,C_CON.D_cut$duration)
      C_CON.D_cut$duration = ifelse(!is.na(C_CON.D_cut$unknown) & C_CON.D_cut$unknown >=cut,"X",C_CON.D_cut$duration)
      C_CON.D_cut$duration = ifelse(C_CON.D$ADL.D == "L" & C_CON.D_cut$duration == "X", cut , C_CON.D_cut$duration)
      C_CON.D_cut$ADL.D = paste0(C_CON.D_cut$ADL.D,C_CON.D_cut$duration)
      C_CON.D_cut$ADL.D = ifelse(C_CON.D_cut$ADL.D =="HNA","H",C_CON.D_cut$ADL.D)
      C_CON.D_cut$ADL.D = ifelse(C_CON.D_cut$ADL.D %in% c("ANA","LNA","NANA"),NA,C_CON.D_cut$ADL.D)
      C_CON.D_cut$pre_ADL.D = lag(C_CON.D_cut$ADL.D)
      
      if(iter==1){
        tem = data.frame(rahhidpn = unique(C_CON.D_cut$rahhidpn))
      }else{
        tem = data.frame(rahhidpn = sample(unique(C_CON.D_cut$rahhidpn),length(unique(C_CON.D_cut$rahhidpn)),replace = T))
      }
      C_CON = left_join(tem,C_CON.D_cut)
      
      C_CON$pre_ADL.D = factor(C_CON$pre_ADL.D,levels = c("A0","A1","A2","A3","A4","A5","A6","A7","A8","A9","AX","L0","L1","L2","L3","L4","L5","L6","L7","L8","L9"))
      C_CON$pre_ADL.D = factor(C_CON$pre_ADL.D,levels = unique(sort(C_CON$pre_ADL.D)))
      
      C_CON$ADL.D = factor(C_CON$ADL.D,levels = c("A0","A1","A2","A3","A4","A5","A6","A7","A8","A9","AX","L0","L1","L2","L3","L4","L5","L6","L7","L8","L9","H"))
      C_CON$ADL.D = factor(C_CON$ADL.D,levels = unique(sort(C_CON$ADL.D)))
      
      tmsa = multinom(formula = ADL.D ~ pre_ADL.D + age + agesq + as.factor(Sex) + as.factor(Sex):age,
                      weights = wt*w, data = C_CON %>% mutate(agesq = age^2) %>% filter(age>=50),na.action=na.omit,maxit = 1000)
      
      dwrite <- crossing(pre_ADL.D = unique(sort(C_CON$pre_ADL.D)), Sex = c("1","2"),age=c(55:105))
      dwrite$agesq = dwrite$age ^ 2
      
      pp.s <- cbind(dwrite, predict(tmsa, newdata = dwrite, type = "probs", se = TRUE))
      lpp.s <- pivot_longer(pp.s,c(5:ncol(pp.s)))
      lpp.s$value = ifelse(lpp.s$value<0.0005,NA,lpp.s$value)
      lpp.s = lpp.s %>% group_by(pre_ADL.D,Sex,age) %>% mutate(value = value/sum(value,na.rm = T)) %>% filter(!is.na(value))
      lpp.s$iter = iter
      lpp.s
    }
    stopImplicitCluster()
    write_csv(PROB,paste0("output/Sex/Res/PROB",CO,".",cut," (his).csv"))
    
  }
  
}