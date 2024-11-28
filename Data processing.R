rm(list=setdiff(ls(), "Raw"))
library(tidyverse)
library(nnet)
library(doParallel)
library(data.table)
# registerDoParallel(detectCores())

load("HRS2020.RData")
rm(list=setdiff(ls(), "Raw"))
ADL = Raw %>% select("rahhidpn",
                     # paste0("r",c(2:14),"walk1a"),
                     # paste0("r",c(2:14),"stoopa"),
                     # paste0("r",c(2:14),"sita"),
                     # paste0("r",c(2:14),"chaira"),
                     # paste0("r",c(2:14),"clim1a"),
                     # paste0("r",c(2:14),"lifta"),
                     # paste0("r",c(2:14),"dimea"),
                     # paste0("r",c(2:14),"armsa"),
                     # paste0("r",c(2:14),"pusha"),
                     # paste0("r",c(2:14),"iadlza"),
                     # paste0("r",c(2:14),"adlwa"),
                     paste0("r",c(2:15),"adl5a"))
# the functional item are set to no issue if missing
# ADL[,c(2:(13*9+1))][is.na(ADL[,c(2:(13*9+1))])] <- 0
ADL = pivot_longer(ADL,c(2:ncol(ADL)))
ADL$wave= parse_number(ADL$name)
ADL$name = gsub("[0-9]","",ADL$name)
ADL = pivot_wider(ADL,names_from=name,values_from = value)
# ADL = ADL %>% mutate(radl = rwalka+rstoopa+rsita+rchaira+rclima+rlifta+
#                        rdimea+rarmsa+rpusha+4*radlwa)
# ADL = ADL[,c(1,2,14)]
colnames(ADL)[3]="radl"

SW = Raw %>% select("rahhidpn",paste0("r",c(2:14),"wtcrnh"))
SW = pivot_longer(SW,c(2:14))
SW$wave= parse_number(SW$name)
SW$name = gsub("[0-9]","",SW$name)
SW = pivot_wider(SW,names_from=name,values_from = value)

dem = Raw %>% select("rahhidpn","raeduc","ragender","raracem","rahispan")
dem$raeduc = haven::as_factor(dem$raeduc)
dem$raeduc = as.character(dem$raeduc)
dem$raeduc = substring(dem$raeduc,0,1)
dem$edu = ifelse(dem$raeduc %in% c("2","3"),3,ifelse(dem$raeduc %in% c("4","5"),5,ifelse(dem$raeduc %in% c("1"),1,NA)))
dem$race = dem$raracem
dem$race = ifelse(dem$race == 3, 4 ,dem$race)
dem$race = ifelse(dem$rahispan == 1, 3 ,dem$race)
# dem = dem %>% filter(race != 4)

YAW = data.frame(FY = c(1936,1926),
                 FA = c(60,70),
                 WA = c(5,5))

save.image(file = "output/dem.Rdata")

for (CO in unique(YAW$FY)) {
  Cohort = Raw %>% select("rahhidpn","rabyear","radyear",paste0("r",c(5:15),"agey_e")) %>% 
    filter(rabyear %in% c(CO:(CO+9))) %>%
    mutate(WD = as.integer(radyear/2-995),WD = ifelse(is.na(WD),15,WD))
  Cohort = pivot_longer(Cohort,c(4:14),values_to = "age")
  Cohort$wave= parse_number(Cohort$name)
  Cohort = Cohort %>% select(-name)
  
  Cohort = left_join(Cohort,ADL)
  Cohort$radl=ifelse(Cohort$wave>Cohort$WD,999,Cohort$radl)
  
  Cohort$state=ifelse(Cohort$radl <1,"A",ifelse(Cohort$radl<8,"L",ifelse(Cohort$radl<100,"L","H")))
  
  unique(Cohort$state)
  Cohort$state = factor(Cohort$state,level= c("A","L","H"))
  
  ### remove empty waves and fill age
  Cohort[which(Cohort$state=="H"),"age"] = Cohort$radyear[which(Cohort$state=="H")]-Cohort$rabyear[which(Cohort$state=="H")]
  c = Cohort
  Mis = Cohort[which(is.na(Cohort$state)),]
  for (i in 1:nrow(Mis)) {
    I = Mis$rahhidpn[i]
    W = Mis$wave[i]
    if(length(is.na(Cohort[which(Cohort$rahhidpn == I & Cohort$wave==W+1),"state"]))==0|
       length(is.na(Cohort[which(Cohort$rahhidpn == I & Cohort$wave==W-1),"state"]))==0){next()}
    if(is.na(Cohort[which(Cohort$rahhidpn == I & Cohort$wave==W+1),"state"])|
       is.na(Cohort[which(Cohort$rahhidpn == I & Cohort$wave==W-1),"state"])){next()}
    if(Cohort[which(Cohort$rahhidpn == I & Cohort$wave==W+1),"state"]=="H"){
      c[which(c$rahhidpn==I&c$wave==W),"state"]=Cohort[which(Cohort$rahhidpn==I&Cohort$wave==W-1),"state"]
      c[which(c$rahhidpn==I&c$wave==W),"age"] = ifelse((Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W-1)]+1)==(Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W+1)]-1),(Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W+1)]-1),
                                                       sample(c((Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W-1)]+1):(Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W+1)]-1)),1))
    }else{
      n=rnorm(1)
      c[which(c$rahhidpn==I&c$wave==W),"state"]=ifelse(n>0,Cohort[which(Cohort$rahhidpn==I&Cohort$wave==W-1),"state"],Cohort[which(Cohort$rahhidpn==I&Cohort$wave==W+1),"state"])
      c[which(c$rahhidpn==I&c$wave==W),"age"] = ifelse((Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W-1)]+1)==(Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W+1)]-1),(Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W+1)]-1),
                                                       sample(c((Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W-1)]+1):(Cohort$age[which(Cohort$rahhidpn==I&Cohort$wave==W+1)]-1)),1))
    }
  }
  
  c_aw = left_join(c,dem)
  c_aw$pre_state = c_aw$state
  c_aw$state = lead(c_aw$state)
  c_aw$aw = ifelse(!is.na(c_aw$pre_state) & is.na(c_aw$state),1,0)
  c_aw$aw = ifelse(c_aw$pre_state=="H",NA,c_aw$aw)
  c_aw$aw = ifelse(c_aw$wave==15,NA,c_aw$aw)
  denom.fit = glm(aw ~ as.factor(ragender) + as.factor(race) + as.factor(edu) + age + as.factor(pre_state), 
                  family = binomial(), data = c_aw, na.action = na.exclude)
  # pd <- predict(denom.fit, type = "response")
  
  numer.fit <- glm(aw~1, family = binomial(), data = c_aw, na.action = na.exclude)
  # pn <- predict(numer.fit, type = "response")
  
  Mis = c[which(is.na(c$state)),]
  c_hh=c()
  for (id in unique(c$rahhidpn)) {
    tem = c %>% filter(rahhidpn==id) 
    if(length(which(tem$state=="H")[-1])>0){
      tem = tem[-which(tem$state=="H")[-1],]
    }
    if(nrow(tem) >1){
      c_hh = rbind(c_hh,tem)
    }
  }
  
  ### first imputation (randomly selected from the observation before and after)
  c_con = c()
  for (id in unique(c_hh$rahhidpn)) {
    tem = c_hh %>% filter(rahhidpn==id)
    tem$wave = ifelse(tem$state=="H",NA,tem$wave)
    if("TRUE" %in% duplicated(tem$age)){
      if(!"TRUE" %in% is.na(tem[which(tem$age == tem$age[which(duplicated(tem$age))]),"wave"]))
        tem = tem[-which(duplicated(tem$age)),]
    }
    tem = full_join(tem,data.frame(age= c(50:105)),by = "age")
    tem = arrange(tem,age)
    tem$rahhidpn=id
    Mis = tem[which(is.na(tem$state)),]
    for (i in 1:nrow(Mis)) {
      A = Mis$age[i]
      if(length(is.na(tem[which(tem$age==A-1),"state"]))==0|
         length(is.na(tem[which(tem$age==A+1),"state"]))==0|
         length(is.na(tem[which(tem$age==A+2),"state"]))==0){next()}
      if(is.na(tem[which(tem$age==A+1),"state"])[1] & 
         is.na(tem[which(tem$age==A+2),"state"])[1]){next()}
      if(is.na(tem[which(tem$age==A+1),"state"])[1] & 
         !is.na(tem[which(tem$age==A+2),"state"])[1]){
        tem[which(tem$age==A),"state"]=tem[which(tem$age==A-1),"state"]
      }else if("TRUE" %in% (tem[which(tem$age==A+1),"state"]=="H")){
        tem[which(tem$age==A),"state"]=tem[which(tem$age==A-1),"state"]
      }else if(!is.na(tem[which(tem$age==A-1),"state"])&
               !is.na(tem[which(tem$age==A+1),"state"])){
        tem[which(tem$age==A),"state"]=sample(c(tem[which(tem$age==A-1),"state"],tem[which(tem$age==A+1),"state"]),1)
      }
    }
    c_con = rbind(c_con,tem)
  }
  c_con = c_con %>% filter(!is.na(age))
  
  c_con = c_con %>% group_by(rahhidpn) %>% mutate(pre_state = lag(state))
  c_con$dis = substring(c_con$state,0,1)
  c_con$adl = substring(c_con$state,2,3)
  c_con$pre_dis = substring(c_con$pre_state,0,1)
  c_con$pre_adl = substring(c_con$pre_state,2,3)
  # c_con = c_con %>% select(-raeduc,-ragender,-raracem,-edu,-rahispan,-race)
  
  
  c_con = left_join(c_con,dem)
  
  ### sample weight
  c_con = left_join(c_con,SW)
  c_con = c_con %>% group_by(rahhidpn) %>% mutate(lag_wt = lag(rwtcrnh)) %>% 
    mutate(wt = coalesce(rwtcrnh,lag_wt)) %>% mutate(lag_wt = lag(wt)) %>% 
    mutate(wt = coalesce(wt,lag_wt)) %>% mutate(lag_wt = lag(wt)) %>% 
    mutate(wt = coalesce(wt,lag_wt)) %>% select(-lag_wt,-rwtcrnh)
  
  c_con$pre_state = as.character(c_con$pre_state)
  c_con$pre_state = ifelse(c_con$pre_state =="H",NA,c_con$pre_state)
  c_con$pre_state = factor(c_con$pre_state,levels = c("A","L"))
  pd <- predict(denom.fit, newdata = c_con, type = "response")
  pn <- predict(numer.fit, newdata = c_con, type = "response")
  c_con$w <- (1-pn)/(1-pd)
  
  save(list = c("Cohort","c_con"),file = paste0("output/C",substr(CO,3,4)," (10).Rdata"))
  
  
} 

# second imputation (allow interwave different from the neighboring observations)
co = c("26","36")
load("output/dem.Rdata")
rm("Raw")

for (CO in co) {
  for(iter in c(1:10)){
  yaw = YAW %>% filter(FY==paste0("19",CO))
  load(paste0("output/C",CO," (10).RData"))
  c_con$sex = ifelse(c_con$ragender==1,1,ifelse(c_con$ragender == 2,2, NA))
  Age = (yaw[1,"FA"]-6):(yaw[1,"FA"]+25)
  iteration = 20
  L = length(Age)
  registerDoParallel(min(detectCores(),14))
    #TP = c()
    ALL = foreach(iter=1:iteration,.packages= c("tidyverse","nnet")) %dopar% {
      
      if(iteration == 1){
        C_CON = c_con
      }else{
        tem = data.frame(rahhidpn = sample(unique(c_con$rahhidpn),length(unique(c_con$rahhidpn)),replace = T))
        C_CON = left_join(tem,c_con)
      }
      
      BL = C_CON %>% filter(age%in% c((yaw[1,"FA"]-5):(yaw[1,"FA"]+4))) %>% filter(wave == yaw[1,"WA"])%>% filter(wt >0)
      BL = BL %>% group_by(state,sex) %>% summarise(tot = sum(wt))
      # BL$sex = as.numeric(as.character(BL$ragender))
      BL = BL %>% filter(!is.na(sex)&!is.na(state))
      BL = BL %>% mutate(state = as.character(state))
      # BL = BL %>% ungroup() %>% mutate(pro = round(tot/sum(tot)*1000000,0)) %>% filter(pro>0)
      BL$state = recode(BL$state,"A"=1,"L"=2)
      BL$tot = 25000
      
      bl = BL %>% mutate_at(c("sex"),as.numeric)
      bl = as.matrix(bl[,c(1,2,3)])
      bl[,3] = as.numeric(bl[,3])
      
      tms = multinom(formula = state ~ pre_state + age + agesq + as.factor(sex) +
                       as.factor(sex)*age, 
                     weights = wt*w, data = C_CON %>% mutate(agesq = age^2),na.action=na.omit,maxit = 1000)
      
      dwrite <- crossing(pre_state = c("A","L"), sex = c("1","2"),age=c(50:95))
      dwrite$agesq = dwrite$age ^ 2
      pp.writes <- cbind(dwrite, predict(tms, newdata = dwrite, type = "probs", se = TRUE))
      
      lpp.writes <- pp.writes %>% pivot_longer(c(5:7),names_to = "state",values_to = "prob")
      lpp.writes$pre_state= recode(lpp.writes$pre_state,"A"=1,"L"=2)
      lpp.writes$state= recode(lpp.writes$state,"A"=1,"L"=2,"H"=3)
      lpp.writes = lpp.writes %>% group_by(pre_state,sex,age) %>% mutate(c = cumsum(prob))
      lpp.writes = lpp.writes %>% mutate(c = ifelse(state=="3", 1, c))
      lpp.writes = lpp.writes %>% select(-agesq,-prob)
      ## add in dead to dead prob
      H = crossing(pre_state= 3,sex= unique(lpp.writes$sex),
                   age = unique(lpp.writes$age),state = c(1,2,3),c = 0)
      H = H %>% mutate(c = ifelse(pre_state %in% c(3) & state %in% c(3), 1, c))
      lpp.writes = bind_rows(lpp.writes,H)
      lpp.writes = xtabs(c ~ pre_state+state+age+sex,data = lpp.writes)
      
      # write_csv(pp.writes,paste0("tmp/pp.",(z-1)*iteration+iter,".csv"))
      # write_csv(bl,paste0("tmp/bl.",(z-1)*iteration+iter,".csv"))
      
      # microsimulation
      Res = c()
      for (n in 1:nrow(bl)) {
        tem = matrix(c(rep(bl[n,2], bl[n,3]),rep(bl[n,1], bl[n,3])),nrow = bl[n,3])
        Res = rbind(Res,tem)
      }
      
      for(age in Age){
        Y = apply(Res, 1, function(x) lpp.writes[as.character(x[age-(Age[1]-2)]),,as.character(age),as.character(x[1])])
        Y = t(Y)
        
        NEX = apply(Y,1,function(y){
          RAM = runif(1)
          s = ifelse(RAM<y[1],1,ifelse(RAM<y[2],2,3))
          s
        })
        Res= cbind(Res,NEX)
      }
      Res
    }
    stopImplicitCluster() 
    # assign(paste0("ALL",CO,yaw[1,"FA"],".",z),ALL)
    # save(list = c(paste0("ALL",CO,yaw[1,"FA"],".",z)),file = paste0("output/Sex/ALL",CO,yaw[1,"FA"],".",z,".RData"))
    # rm(list = paste0("ALL",CO,yaw[1,"FA"],".",z))
    
    Age = (yaw[1,"FA"]-6):(yaw[1,"FA"]+26)
    All = do.call(rbind,ALL) %>% as.data.frame()
    colnames(All) = c("Sex",Age) 
    All$id = 1:nrow(All)
    All = pivot_longer(All,2:(ncol(All)-1),names_to = "Age",values_to = "State")
    All$State = recode(All$State,`1`="A",`2`="L",`3`="H")
    All$Age = as.numeric(All$Age)
    All = setDT(All)
    setkey(All, Sex, Age, State)
    
    load(paste0("output/C",CO," (10).RData"))
    Cohort = left_join(Cohort,dem)
    Cohort$sex = ifelse(Cohort$ragender==1,1,ifelse(Cohort$ragender == 2,2, NA))
  
    registerDoParallel(min(detectCores(),14))
    gap_filled = foreach(id = unique(Cohort$rahhidpn),.packages= c("tidyverse","data.table"),.combine = "rbind") %dopar% {
      tem = filter(Cohort,rahhidpn==id) 
      x = tem %>% filter(state %in% c("A","L"))
      if(nrow(x)<2 & nrow(x)>0){
        x = tem %>% filter(state %in% c("H"))
        if(nrow(x)<1){
          return(NULL)
        }
      }else if(nrow(x)<1){
        return(NULL)
      }
      age = tem %>% filter(state == "H") %>% slice(1) %>% pull(age) 
      tem = setDT(tem)
      setkey(tem, sex, age, state)
      # Perform the left join
      tem <- All[tem, on = .(Sex = sex, Age = age, State = state), nomatch = 0]
      # tem = left_join(tem, All ,by = c("sex"="Sex","age"="Age","state"="State"))
      if(length(age)>0){
        x2=0
        repeat {
          # Calculate the table of 'id' frequencies and get the maximum frequency
          id_freq <- tem[, .N, by = id]
          max_freq <- max(id_freq$N)
          
          # Get the names of 'id' values with the maximum frequency
          max_freq_ids <- id_freq[N == max_freq, id]
          
          # Randomly select one 'id' value from the ones with the maximum frequency
          x <- sample(max_freq_ids, 1)
          tem2 <- All[id == x & Age == age - 1 & State != "H"]
          if(nrow(tem2)>0|x2==x){
            break
          }
          x2 = x
        }
      }else{
        # Calculate the table of 'id' frequencies and get the maximum frequency
        id_freq <- tem[, .N, by = id]
        max_freq <- max(id_freq$N)
        
        # Get the names of 'id' values with the maximum frequency
        max_freq_ids <- id_freq[N == max_freq, id]
        
        # Randomly select one 'id' value from the ones with the maximum frequency
        x <- sample(max_freq_ids, 1)
      }
      tem2 <- All[id == x,] %>% select(-id)
      
      filled = filter(c_con,rahhidpn==id)
      filled = left_join(filled,tem2,by = c("age"="Age"))
      filled
    }
    stopImplicitCluster() 
    Sys.time()
    
    gap_filled$State = ifelse(is.na(gap_filled$state),NA, gap_filled$State)
    gap_filled$State = ifelse(!is.na(gap_filled$radl),as.character(gap_filled$state),gap_filled$State)
    length(unique(gap_filled$rahhidpn))
    
    save(list = "gap_filled",file = paste0("output/Sex/ALL",CO," (imp",iter,").RData"))
    
  }
}


#### generate duration column ####

co = c("26","36")
for (iter in 1:10) {
  for (CO in co) {
    load(paste0("output/Sex/ALL",CO," (imp",iter,").Rdata"))
    
    gap_filled$ADL.D = gap_filled$State
    C_CON.D <- list()
    for (id in unique(gap_filled$rahhidpn)) {
      tem = gap_filled %>% filter(rahhidpn == id)
      
      tem$duration = NA
      tem$unknown = NA
      # tem = tem %>% mutate(pre_adl1 = lag(ADL.D,1))
      tem$ADL.D = as.character(tem$ADL.D)
      tem$ADL.D[which(is.na(tem$ADL.D))] = "X"
      
      for (n in 2:nrow(tem)) {
        if(tem$ADL.D[n] == "X"){
          next()
        }else if(tem$ADL.D[n] == "H"){
          break()
        }else if(tem$ADL.D[n-1]=="X"){
          tem$unknown[n] = 0
        }else if(tem$ADL.D[n-1]==tem$ADL.D[n]){
          tem$unknown[n] = tem$unknown[n-1]+1
          tem$duration[n] = tem$duration[n-1]+1
        }else if(tem$ADL.D[n-1]!=tem$ADL.D[n]){
          tem$duration[n] = 0
        }
        
      }
      
      C_CON.D <- c(C_CON.D, list(tem))
      
    }
    C_CON.D <- do.call(rbind, C_CON.D)
    
    C_CON.D$ADL.D = ifelse(C_CON.D$ADL.D=="X",NA,C_CON.D$ADL.D)
    C_CON.D$ADL.D = factor(C_CON.D$ADL.D,levels = c("A","L","H"))
    save(C_CON.D,file = paste0("output/Sex/ALL",CO," (duration",iter,").Rdata"))
  }
}


