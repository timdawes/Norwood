# "Predicting clinical trajectory after the Norwood procedure: focus on aortic and pulmonary re-interventions"
# Authors: Haapanen H, Dawes TJW, Brown K, Giardini A, Dedieu N, Shetty P, Tsang V, Kostolny M.


# SetUp Code



# Load patient data
      cat("\n Loading the Excel spreadsheets...\n")
      col_types = as.character(read.csv("Data/col_types.csv")[1:357,1])
      d<- read_excel("Data/Norwood.xlsx", sheet = "Mastersheet", col_types=col_types)
      
      col_types9 = as.character(read.csv("Data/col_types.csv")[1:61,9])
      d9<- read_excel("Data/Norwood.xlsx", sheet = "Ao_PA_reop_dates", col_types=col_types9)

# Exclude the hybrid cases
      
      # Find MRNs for cases which are not a hybrid
          MRNs_included<- d$Hosp_no[which(d$Suc_hyb == 0)]
          no.cases<- length(MRNs_included)
          
      # Include only non-hybrid cases
          d<- d[match(MRNs_included, d$Hosp_no),]
          d9<- d9[match(MRNs_included, d9$Hosp_no),]
          
        cat("\n ",no.cases,"cases included for analysis...\n")
          

cat("\n Extracting the PA re-operation data...\n")

# PA reoperation
      # Construct long-database for analysis
            d9.long<- NULL              
            
            PA.reop.cols<- match(c("PA_1st_reop_date","PA_2nd_reop_date","PA_3rd_reop_date","PA_4th_reop_date", "PA_5th_reop_date","PA_6th_reop_date"), colnames(d9))
            Ao.reop.cols<- match(c("Ao_1st_reop_date","Ao_2nd_reop_date","Ao_3rd_reop_date","Ao_4th_reop_date"), colnames(d9))
            
            sv.operation.cols<- match(c("Norwood_date","Glenn_date", "TCPC_date"), colnames(d9))
            FU.cols<- match(c("DOD_or_transplant_date", "Last_FU_date"), colnames(d9))
            DorT.cols<- match("Dead_or_transplant", colnames(d9))


# Transfer matrix
      
                  
              # All re-operations: complex transition matrix
                  
                  # Initialize a 14x14 matrix filled with NA
                          mat.Allc <- matrix(NA, nrow = 14, ncol = 14)
                  
                          mat.Allc[1,]<- c(NA,1,NA,2,3,NA,NA,NA,NA,NA,NA,NA,NA,4)
                          mat.Allc[2,]<- c(NA,NA,5,6,7,NA,NA,NA,NA,NA,NA,NA,NA,8)
                          mat.Allc[3,]<- c(NA,NA,NA,9,10,NA,NA,NA,NA,NA,NA,NA,NA,11)
                          mat.Allc[4,]<- c(NA,NA,NA,NA,NA,12,NA,NA,13,14,NA,NA,NA,15)
                          mat.Allc[5,]<- c(NA,NA,NA,NA,NA,16,NA,NA,17,18,NA,NA,NA,19)
                          mat.Allc[6,]<- c(NA,NA,NA,NA,NA,NA,20,NA,21,22,NA,NA,NA,23)
                          mat.Allc[7,]<- c(NA,NA,NA,NA,NA,NA,NA,24,25,26,NA,NA,NA,27)
                          mat.Allc[8,]<- c(NA,NA,NA,NA,NA,NA,NA,NA,28,29,NA,NA,NA,30)
                          mat.Allc[9,]<- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,31,NA,NA,32)
                          mat.Allc[10,]<- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,33,NA,NA,34)
                          mat.Allc[11,]<- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,35,NA,36)
                          mat.Allc[12,]<- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,37,38)
                          mat.Allc[13,]<- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,39)
                          mat.Allc[14,]<- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
                          
                  # Format to a transfer matrix
                          tmat.Allc<- mstate::transMat(tapply(mat.Allc, rep(1:nrow(mat.Allc), ncol(mat.Allc)), function(i) which(is.na(i)==FALSE)),
                                                       names = c("N","R1","R2","Glenn","Glenn_R","R3","R4","R5","TCPC","TCPC_R","R6","R7","R8","DorT"))
                          
                        
                  
      # Multi-state modelling matrix wide versions
            

                    # Define the covariates we're interested in 
                                covs<- c("Weight", "Preterm", "era","Ndays","SysV","DGcat","SystV_f","AVVRcat","SanovsBT","Interdigitation","N1ECMO")
                                covs.original.names<- c("N1weight", "Preterm", "era","Ndays","SysV", "DG", "SystV_f", "AVVR", "SanovsBT", "N1adigit")
            

                   
                   # All re-interventions
                         msm1.Allc<- msm2.Allc<- NULL
                         
                         msm2.Allc.wide<- matrix(0, nrow=no.cases, ncol=41)
                         paths.Allc<- paths.Allc.DorT<- rep(NA, length=no.cases)
                         
                         rownames(msm2.Allc.wide)<- 1:no.cases
                         colnames(msm2.Allc.wide)<- c("id",
                                                     sapply(2:14, function(x) {paste(c("state","censor"),x, sep="")}),
                                                     c("N1Age","Weight","Preterm","era","Ndays","SysV","DG","DGcat","SystV_f","AVVR","AVVRcat","SanovsBT","Interdigitation","N1ECMO"))
                         
                         msm2.Allc.wide[,1]<- 1:no.cases
                         msm2.Allc.wide[,seq(2,26,2)]<- Inf
                         
                         did.this.transition.happen.Allc<- matrix("-", nrow=14, ncol=14)
                         rownames(did.this.transition.happen.Allc)<- rownames(tmat.Allc)
                         colnames(did.this.transition.happen.Allc)<- rownames(tmat.Allc)
                         
                         FirstReinterventionDate<- matrix(NA, nrow=no.cases, ncol=1)
                         d$FUy_to_firstreint<- NA
                         d$Reint<- 0
                   
                   # Multi-stage modelling for Ao interventions
                   
                   cat("\n Extracting the joint Ao/PA re-operation data for complex reoperation model...\n")
                   
                   for (i in 1:nrow(d9))
                   {
                      # Patients 9,32,51,62,66,68 had both Ao and PA interventions
                     
                     # Update progress bar
                          setTxtProgressBar(pb,i)
                     
                     # Extract dates for SV pathway operations and FU  
                         Ao.reop.dates<- data.frame(d9[i,Ao.reop.cols[!is.na(d9[i,Ao.reop.cols])]])
                         PA.reop.dates<- data.frame(d9[i,PA.reop.cols[!is.na(d9[i,PA.reop.cols])]])
                         
                     # Amalgamate PA and Ao operations
                         All.reop.dates<- cbind(Ao.reop.dates, PA.reop.dates)
                         unique.reop.dates<- which(!duplicated(as.character(All.reop.dates)))
                         All.reop.dates<- All.reop.dates %>% select(all_of(unique.reop.dates))
                         
                    # Find when the first reintervention was done
                         if (ncol(All.reop.dates) == 0) {FirstReinterventionDate[i]<- NA} else {
                           FirstReinterventionDate[i,]<- as.character(All.reop.dates[,order(unlist(matrix(All.reop.dates[1,])))[1]])
                           d$FUy_to_firstreint[i]<- (ymd(FirstReinterventionDate[i,]) - ymd(d$DOB[i])) / 365
                           d$Reint[i]<- 1
                           }
                         
                         op.names<- colnames(All.reop.dates)
                         op.names<- gsub("PA_|Ao_", "", op.names)
                         op.names<- gsub("1st_|2nd_|3rd_|4th_|5th_|6th_", "", op.names)
                         colnames(All.reop.dates)<- op.names
                         
                     
                     sv.dates<- d9[i,sv.operation.cols]
                     which.are.not.nas<- which(!is.na(sv.dates)[1,])
                     sv.dates<- sv.dates[,which.are.not.nas]
                     FU.dates<- t(data.frame(d9[i,FU.cols]))[,1]
                     FU<- FU.dates[2]
                     DorT<- FU.dates[1]
                     
                     DorT<- t(data.frame(d9[i,DorT.cols]))[,1]
                     
                     all.op.dates<- cbind(sv.dates, All.reop.dates, FU)
                     
                     o<- order(unlist(matrix(all.op.dates[1,])))
                     
                     all.op.dates<- all.op.dates[,o]
                     
                     unique.values<- which(!duplicated(as.character(all.op.dates)))
                     not.unique.values<- which(duplicated(as.character(all.op.dates)))
                     
                     z<- colnames(all.op.dates)
                     
                     # Store all the paths observed in this dataset
                           path1<- paste(" ->", names(all.op.dates)[unique.values])
                           path2<- paste(" AND", names(all.op.dates)[not.unique.values])
                           path3<- paste(c(path1,path2)[order(c(unique.values, not.unique.values))], collapse = "")
                           path4<- gsub("_date", "", path3)
                           paths.Allc[i]<- path4
                           
                     
                     # Any re-operation during SV operation
                     
                     for (date.counter in unique.values[-length(unique.values)])
                     {
                       date.counter.plus.one<- unique.values[match(date.counter, unique.values) + 1]
                       date.counter.plus.two<- unique.values[match(date.counter, unique.values) + 2]
                       
                       # Calculate which row of matrix 'm' you should jump to next
                       op.from<- match(z[date.counter], c("Norwood_date","Glenn_date","TCPC_date",
                                                          "reop_date","reop_date.1","reop_date.2","reop_date.3","reop_date.4","reop_date.5","reop_date.6","reop_date.7","reop_date.8","FU")) 
                       
                       # Define how many re-operations were done before each SV stage
                       # Before Glenn
                           glenn.number<- match("Glenn_date", z)
                           if (is.na(glenn.number) == TRUE) {glenn.number<- match("FU", z)}
                           
                       # Before TCPC  
                           tcpc.number<- match("TCPC_date", z)
                           if (is.na(tcpc.number) == TRUE) {tcpc.number<- match("FU", z)}
                           
                          reop.numbers<- grep("reop", z)
                       
                       if (op.from>3)
                       {
                         
                         if (!is.na(match("Glenn_date", z))) {after.glenn<- all.op.dates[date.counter] > all.op.dates[match("Glenn_date", z)]}
                         if (!is.na(match("Glenn_date", z))) {with.glenn<- all.op.dates[date.counter] == all.op.dates[match("Glenn_date", z)]}
                         
                         if (!is.na(match("TCPC_date", z))) {after.tcpc<- all.op.dates[date.counter] > all.op.dates[match("TCPC_date", z)]}
                         if (!is.na(match("TCPC_date", z))) {with.tcpc<- all.op.dates[date.counter] == all.op.dates[match("TCPC_date", z)]}
                         
                         if (length(reop.numbers) > 0) {
                           reop.numbers.before.tcpc<- reop.numbers[sapply(reop.numbers, function(x) {all.op.dates[x] <= all.op.dates[tcpc.number]})]
                           reop.numbers.after.glenn<- reop.numbers[sapply(reop.numbers, function(x) {all.op.dates[x] > all.op.dates[glenn.number]})]
                           reop.numbers.before.glenn<- reop.numbers[sapply(reop.numbers, function(x) {all.op.dates[x] <= all.op.dates[glenn.number]})]
                           
                           interstage2<- length(intersect(reop.numbers.after.glenn, reop.numbers.before.tcpc))
                           interstage1<- length(reop.numbers.before.glenn)
                           
                         } else {interstage2<- interstage1<- 0}
                         
                         
                         if (after.tcpc==TRUE) {matrix.row<- 7 + op.from - interstage2 - interstage1} # reop after TCPC
                         if (with.tcpc==TRUE) {matrix.row<- 10} # Reop after TCPC
                         
                         if (after.glenn==TRUE && after.tcpc==FALSE) {matrix.row<- 2 + op.from - interstage1} # reop after Glenn
                         if (with.glenn==TRUE) {matrix.row<- 5}
                         
                         if (after.glenn==FALSE && with.glenn==FALSE) {matrix.row<- op.from - 2}
                         
                         
                       }  else {
                         
                         if (!is.na(match("Glenn_date", z))) {reop.with.glenn<- all.op.dates["Glenn_date"] == all.op.dates[match("Glenn_date", z) + 1]}
                         if (!is.na(match("TCPC_date", z))) {reop.with.tcpc<- all.op.dates["TCPC_date"] == all.op.dates[match("TCPC_date", z) + 1]}
                         
                         if (op.from == 1) {matrix.row<- 1}
                         
                         if (op.from == 2 && reop.with.glenn == FALSE) {matrix.row<- 4}
                         if (op.from == 2 && reop.with.glenn == TRUE) {matrix.row<- 5}
                         
                         if (op.from == 3 && reop.with.tcpc == FALSE) {matrix.row<- 9}
                         if (op.from == 3 && reop.with.tcpc == TRUE) {matrix.row<- 10}
                         
                         
                       }
                       
                       
                       # Calculate which column of matrix 'mat.Allc' the patient had next
                               matrix.col<- NA
                               op.to<- match(z[date.counter.plus.one], c("Norwood_date","Glenn_date","TCPC_date","reop_date","reop_date.1","reop_date.2","reop_date.3","reop_date.4","reop_date.5","reop_date.6","reop_date.7","reop_date.8","FU")) 
                               
                       # Which date is the Glenn?
                               glenn.number<- match("Glenn_date", z)
                               if (is.na(glenn.number) == TRUE) {glenn.number<- match("FU", z)}
                               
                       # Which date is the TCPC?
                               tcpc.number<- match("TCPC_date", z)
                               if (is.na(tcpc.number) == TRUE) {tcpc.number<- match("FU", z)}
                               
                       # What are the reoperation numbers?
                              reop.numbers<- grep("reop", z)
                       
                       if (op.to>3)
                       {
                         if (!is.na(match("Glenn_date", z))) {after.glenn<- all.op.dates[date.counter.plus.one] > all.op.dates[match("Glenn_date", z)]} else {after.glenn<- FALSE}
                         if (!is.na(match("Glenn_date", z))) {with.glenn<- all.op.dates[date.counter.plus.one] == all.op.dates[match("Glenn_date", z)]} else {after.glenn<- FALSE}
                         
                         if (!is.na(match("TCPC_date", z))) {after.tcpc<- all.op.dates[date.counter.plus.one] > all.op.dates[match("TCPC_date", z)]} else {after.tcpc<- FALSE}
                         if (!is.na(match("TCPC_date", z))) {with.tcpc<- all.op.dates[date.counter.plus.one] == all.op.dates[match("TCPC_date", z)]} else {with.tcpc<- FALSE}
                         
                         if (length(reop.numbers) > 0) {
                           
                           reop.numbers.before.tcpc<- reop.numbers[sapply(reop.numbers, function(x) {all.op.dates[x] <= all.op.dates[tcpc.number]})]
                           reop.numbers.after.glenn<- reop.numbers[sapply(reop.numbers, function(x) {all.op.dates[x] > all.op.dates[glenn.number]})]
                           reop.numbers.before.glenn<- reop.numbers[sapply(reop.numbers, function(x) {all.op.dates[x] <= all.op.dates[glenn.number]})]
                           
                           interstage2<- length(intersect(reop.numbers.after.glenn, reop.numbers.before.tcpc))
                           interstage1<- length(reop.numbers.before.glenn)
                           
                         } else {interstage2<- interstage1<- 0}
                         
                         
                         if (after.tcpc==TRUE) {matrix.col<- 7 + op.to - interstage2 - interstage1} # reop after TCPC
                         if (with.tcpc==TRUE) {matrix.col<- 10} # Ao reop after TCPC
                         
                         if (after.glenn==TRUE && after.tcpc==FALSE) {matrix.col<- 2 + op.to - interstage1} # Ao reop after Glenn
                         if (with.glenn==TRUE) {matrix.col<- 5}
                         
                         if (after.glenn==FALSE && with.glenn==FALSE) {matrix.col<- op.to - 2}
                         
                         if (op.to == 13 && DorT == 1) {matrix.col<- 14}
                         if (op.to == 13 && DorT == 0) {matrix.col<- 999}
                         
                       } else {
                         
                         if (!is.na(match("Glenn_date", z))) {reop.with.glenn<- all.op.dates["Glenn_date"] == all.op.dates[date.counter.plus.one + 1]}
                         if (!is.na(match("TCPC_date", z))) {reop.with.tcpc<- all.op.dates["TCPC_date"] == all.op.dates[date.counter.plus.one + 1]}
                         
                         if (op.to == 1) {matrix.col<- 1}
                         
                         if (op.to == 2 && reop.with.glenn == FALSE) {matrix.col<- 4}
                         if (op.to == 2 && reop.with.glenn == TRUE) {matrix.col<- 5}
                         
                         if (op.to == 3 && reop.with.tcpc == FALSE) {matrix.col<- 9}
                         if (op.to == 3 && reop.with.tcpc == TRUE) {matrix.col<- 10}
                         
                       }
                       
                       
                       
                       # Make a multi-state matrix for THIS patient only (msm.for.one.patient)
                             list.of.tos<- which(!is.na(mat.Allc[matrix.row,]))
                             l<- length(list.of.tos)
                             
                             msm.for.one.patient<- data.frame(matrix(ncol=22, nrow=l))
                             colnames(msm.for.one.patient)<- c("ID", "from","to", "trans", "Tstart", "Tstop", "time",
                                                               "status","N1Age","Preterm","Weight","era","Ndays","SysV",
                                                               "DG","DGcat","SystV_f","AVVR","AVVRcat","SanovsBT","Interdigitation","N1ECMO")
                                                              
                       
                       
                       # Construct matrix for this patient (j) column by column
                       # ID
                          id<- i
                       
                       # From 
                          from<- matrix.row
                       
                       # Selected covariates
                             N1Age<- d$N1age[i]
                             Preterm<- d$Preterm[i]
                             Weight<- d$N1weight[i]
                             era<- match(1, c(d[i,]$era1, d[i,]$era2, d[i,]$era3))
                             if (is.na(d$N1_date[i]) == TRUE) {Ndays<- ymd(d$Hyb_d[i]) - ymd("2000-01-01")} else {Ndays<- ymd(d$N1_date[i]) - ymd("2000-01-01")}
                             SysV<- d$SysV[i]
                             DG<- d$DG[i]
                             DGcat<- as.numeric(DG<5) # All the diagnoses 1-4 are subtypes of HLHS
                             SystV_f<- d$SystV_f[i]
                             AVVR<- d$AVVR[i]
                             AVVRcat<- as.numeric(AVVR<3) # AVVR 0=none, 1=trivial, 2=mild, 3=mod, 4=severe
                             
                       # Surgical covariates
                             SanovsBT<- d$SanovsBT[i]
                             Interdigitation<- d$N1adigit[i]
                             N1ECMO<- as.numeric((d$N1ECMOin[i] == 1) | (d$N1ECMOpo[i] == 1))  # AVVR 0=none, 1=trivial, 2=mild, 3=mod, 4=severe
                             
                       
                       # To, trans, Tstart, Tstop, time, status
                             for (k in 1:l)
                             {
                               to<- list.of.tos[k]
                               trans<- mat.Allc[matrix.row,to]
                               Tstart<- difftime(all.op.dates[,date.counter], all.op.dates[,1], units = "days")
                               Tstop<- difftime(all.op.dates[,date.counter.plus.one], all.op.dates[,1], units = "days")
                               time<- Tstop - Tstart
                               if (to == matrix.col) {status<- 1} else {status<- 0}
                               msm.for.one.patient[k,]<- c(id, from, to, trans, Tstart, Tstop, time, status, N1Age, Preterm, Weight, era, Ndays, SysV, DG, DGcat, SystV_f, AVVR, AVVRcat, SanovsBT, Interdigitation, N1ECMO)
                             }
                             
                       # Register that this transition happened, at least in this patient  
                            if (matrix.col != 999) {did.this.transition.happen.Allc[matrix.row, matrix.col]<- "Yes"}
                       
                       
                       # Construct 'wide' matrix for this patient
                       # id, state1 -> 13, time & censoring for each
                       
                             # Define the column mappings for simplicity
                                   columns <- c("state" = "to", "censor" = "status", "N1Age" = "N1Age", "Preterm" = "Preterm", 
                                                "Weight" = "Weight", "era" = "era", "Ndays" = "Ndays", "SysV" = "SysV", 
                                                "DG" = "DG", "DGcat" = "DGcat", "SystV_f" = "SystV_f", "AVVR" = "AVVR", 
                                                "AVVRcat" = "AVVRcat", "SanovsBT" = "SanovsBT", "Interdigitation" = "Interdigitation", 
                                                "N1ECMO" = "N1ECMO")
                                   
                             for (k in 1:l) {
                               time <- msm.for.one.patient[k,]$time
                               
                               # Update mapped columns dynamically
                               for (col_name in names(columns)) {
                                 value <- msm.for.one.patient[k,][[columns[col_name]]]
                                 col_idx <- which(colnames(msm2.Allc.wide) == col_name)
                                 msm2.Allc.wide[i, col_idx] <- value
                               }
                               
                               # Additional specific updates
                               state_col <- which(colnames(msm2.Allc.wide) == paste("state", msm.for.one.patient[k,]$to, sep = ""))
                               censor_col <- which(colnames(msm2.Allc.wide) == paste("censor", msm.for.one.patient[k,]$to, sep = ""))
                               msm2.Allc.wide[i, state_col] <- Tstop
                               msm2.Allc.wide[i, censor_col] <- msm.for.one.patient[k,]$status
                               msm2.Allc.wide[i, which(msm2.Allc.wide[i,] == Inf)] <- max(Tstop)
                             }
                             
                       
                       msm1.Allc<- rbind(msm1.Allc, msm.for.one.patient)
                     }
                     
                   }
                   
                   cat(" done.")
                   close(pb)
                   
                   
                   # Illegal transitions
                           illegal.transitions<- which(apply(cbind(c(did.this.transition.happen.Allc), c(tmat.Allc)), 1, function(x) {x[1]=="Yes" && is.na(x[2])==TRUE}) == TRUE)
                           if (length(illegal.transitions!=0)) {cat("\n Transitions that happened, but which are NOT allowed:", illegal.transitions)}
                           
                   # Transitions that did NOT happen, but which are allowed:
                           missing.transitions.Allc<- c(tmat.Allc)[which(apply(cbind(c(did.this.transition.happen.Allc), c(tmat.Allc)), 1, function(x) {x[1]=="-" && is.na(x[2])==FALSE}) == TRUE)]
                           #for (i in 1:14) {for (j in 1:14) {if (tmat.Allc[i,j] %in% missing.transitions.Allc) {cat("\n", paste(rownames(tmat.Allc)[i], "->", colnames(tmat.Allc)[j]))}}}
                           
                      # Impute and assign the correct attributes to allows redrank.TD to analyse it
                           imputation.df<- data.frame(
                                   Weight = msm2.Allc.wide[,29],  SanovsBT = as.factor(msm2.Allc.wide[,39]), Interdigitation = as.factor(msm2.Allc.wide[,40]), 
                                   N1Age = msm2.Allc.wide[,28], 
                                   era = msm2.Allc.wide[,31],
                                   Ndays = msm2.Allc.wide[,32], SysV = as.factor(msm2.Allc.wide[,33]),
                                   DG = as.factor(msm2.Allc.wide[,34]),
                                   DGcat = as.factor(msm2.Allc.wide[,35]),
                                   SystV_f = as.factor(msm2.Allc.wide[,36]),
                                   AVVR = as.factor(msm2.Allc.wide[,37]),
                                   AVVRcat = as.factor(msm2.Allc.wide[,38]),
                                   N1ECMO = as.factor(msm2.Allc.wide[,41]),
                                   Preterm = msm2.Allc.wide[,30], 
                                   
                             # Extra variables which improve imputation
                                   Weight_N2 = d$N2weight, Weight_N3 = d$N3weight, ShuntSize = as.factor(d$s_size),
                                   AA_MA = d$AA_MA, AA_MS = d$AA_MS, AS_MA = d$AS_MA, AS_MS = d$AS_MS, d$N1coarc,
                                   d$N1surg, Sex = d$Sex)
                           
                           msm2.Allc.wide.i<- cbind(msm2.Allc.wide[,1:28], missForest(imputation.df)$ximp[,1:15])
                           #msm2.Allc.wide.i[,which(colnames(msm2.Allc.wide.i) == "SystV_f")]<- round(msm2.Allc.wide.i[,which(colnames(msm2.Allc.wide.i) == "SystV_f")])
                           
                           
                      # Multi-stage modelling matrix made by mstate package
                            msm2.Allc.i  <- msprep(data = msm2.Allc.wide.i,
                                              trans = tmat.Allc,
                                              time = c(NA, paste("state",2:14,sep="")),
                                              status = c(NA, paste("censor",2:14,sep="")),
                                              id = 1:no.cases,
                                              keep = c("N1Age","Preterm","Weight","era","Ndays","SysV","DG","DGcat","SystV_f","AVVR","AVVRcat","SanovsBT","Interdigitation", "N1ECMO"))
        
                            class(msm2.Allc.i) <- c("msdata", "data.frame")
                            attr(msm2.Allc.i, "trans")<- tmat.Allc
                           
                      # Convert the times into years
                           msm2.Allc.i[,c("Tstart","Tstop","time")]<- msm2.Allc.i[,c("Tstart","Tstop","time")] / 365.25
                           
                      # Convert 'FU' to 'DorT' for some of the cases
                           for (i in 1:no.cases) {if (d$DorT[i] == 1) {paths.Allc.DorT[i]<- gsub("FU", "DorT", paths.Allc.DorT[i])}}
                
            
      # Convert the covariates into transition-specific one-hot variables
            msm2.Allc.covs<- expand.covs.msdata.TD(msm2.Allc.i, c("N1Age", "Preterm", "Weight", "Ndays","SysV","DG","DGcat","SystV_f","AVVR","AVVRcat","SanovsBT","Interdigitation","N1ECMO"), longnames = FALSE)
            
            
            
            
            
            
            
            
            
            
            
      