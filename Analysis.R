# "Predicting clinical trajectory after the Norwood procedure: focus on aortic and pulmonary re-interventions"
# Authors: Haapanen H, Dawes TJW, Brown K, Giardini A, Dedieu N, Shetty P, Tsang V, Kostolny M.


# Set ... to Desktop folder
      working.dir<- setwd("...")
      setwd(working.dir)
      
# Make a folder on the Desktop called 'Norwood' with multiple subfolders called: Code, Data

# Load libraries and input data from Excel sheet
      source("Norwood/Code/Functions.R")
      source("Norwood/Code/SetUp.R")
      

    #  Non-parametric models
                # These ignore the influence of covariates
                # Aim is to calculate transition intensities for the model (= hazard rate for multi-state models)
                # i.e. describes the rate of transitioning from one state to another state, so can be >1
                # Do this using Cox PH models, with separate baseline hazards for each transition
                # When using discrete models, the term 'transition hazard' is often used instead of 'transition intensity'

                c0.Allc <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = msm2.Allc.i, method = "breslow")
                
                # Estimate the CUMULATIVE TRANSITION HAZARDS for each possible transition in the msm
                # Also associated (co)variances - NB variances are the diagonal of the matrix
                        par(mfrow=c(1,1))
                        msf0.Allc <- msfit(object = c0.Allc, vartype = "greenwood", trans = tmat.Allc)
                        par(mfrow=c(1,1), mar=c(4,4,2,2))
                        plot(msf0.Allc, lwd=4)


                # Estimate the TRANSITION PROBABILITIES 
                        pt0.Allc <- probtrans(msf0.Allc, predt = 0, method = "greenwood")

                # ...and plot the transition probabilities  
                        par(mfrow=c(1,1))
                        plot(pt0.Allc, xlab = "Years since Norwood", las = 1, type = "filled", col = rainbow(20))


      # Reduced rank models and simulation   
            
            
              # Reduced rank regression
            
              # Iterative fitting, so repeated over 's' different seeds
            
              # Set up the matrices
                    s<- 1000
                    loglik.mat<- matrix(0, nrow=s, ncol=2)
                    
                    cox.coef1<- matrix(0, nrow=s, ncol=45, dimnames = list(paste("iter",1:s,sep=""), 
                                                 paste(rep(c("coef","se","se2","Chisq","p"), each = 9),
                                                       rep(c("Weight","Ndays","DGcat","SysV","SystV_f1","AVVRcat","SanovsBT","Interdigitation","N1ECMO"), 5), sep = "_")))
                    
                    cox.coef2<- matrix(0, nrow=s, ncol=195, dimnames = list(paste("iter",1:s,sep=""), 
                                                 paste(rep(c("coef","se(coef)","se2","Chisq","p"), each = 39),
                                                      rep(paste("AlphaX",1:39,sep=""), 5), sep="_")))
                    
                    folder.rr<- c("Norwood/Results")
                    write.table(0, file=paste(folder.rr, "rr_iteration_number.csv", sep="/"), col.names=FALSE, row.names=FALSE)
                    pb <- txtProgressBar(min = 0, max = s, style = 3, char = "=")
                    seeds<- 1:s
                    b<- 0
             
              # Iterate through the seeds
                    repeat {
                              # Allow easy exit
                                  stop.check<- read.table(paste(folder.rr, "stop_check.csv", sep="/"), col.names=FALSE)[1,]
                                  if (stop.check != "go") {break}
                                  
                            # Read in iteration number and set seed
                                  b<- read.table(file=paste(folder.rr, "rr_iteration_number.csv", sep="/"), col.names=FALSE)[1,]
                                  b<- b + 1
                                  setTxtProgressBar(pb, b)
                                  set.seed(seeds[b])
                                  rr1<- NULL
                                  
                            
                            
                            
                             while(is.null(rr1)){
                                 try(rr1 <- redrank.TD(Surv(Tstart, Tstop, status) ~ Weight + Ndays + DGcat + SysV + SystV_f + AVVRcat + SanovsBT + Interdigitation + N1ECMO,
                                                       data = msm2.Allc.covs, R = 1, max.iter = 200, eps=1e-8, print.level = 2, remove.NAs = FALSE, scale = TRUE))}
                            
                            # Store the regression on alpha for this iteration
                                    cox.coef1[b,]<- c(unlist(summary(rr1$cox.itr1)$coefficients[,1:4]), p.adjust(summary(rr1$cox.itr1)$coefficients[,6], method = "fdr"))
                                    
                            # Store the regression on gamma for this iteration
                                    cox.coef2[b,]<- c(unlist(summary(rr1$cox.itr2)$coefficients[,1:4]), p.adjust(summary(rr1$cox.itr2)$coefficients[,6], method = "fdr"))
                                                  
                                   loglik.mat[b,]<- c(seeds[b],rr1$loglik)
                                   
                            # Save some pictures
                               
                               if (b>1) {
                                 # Log likelihood variation
                                       pdf(file = paste(working.dir,"/Norwood/Results/Loglik.pdf",sep=""), height = 5, width = 10)
                                       par(mar=c(8,8,2,2))
                                       plot(loglik.mat[1:b,], pch=19, lwd=4, col="blue", bty='n', axes=F, xlab="", ylab="", ylim=10*c(floor(min(loglik.mat[1:b,2])/10), ceiling(max(loglik.mat[1:b,2])/10)))
                                       axis(side = 1, at = seq(0,1000,200), lwd = 4, cex.axis = 1.5)
                                       axis(side = 2, at = seq(-1540,-1500,10), las = 1, lwd = 4, cex.axis = 1.5, line = 0)
                                       mtext("Log likelihood", side = 2, cex=2, line = 5)
                                       mtext("Iteration", side = 1, cex=2, line = 3)
                                       dev.off()
                                       
                                 # Chi-squared values for covariates and for transitions
                                           par(mfrow=c(1,1))
                                           
                                       # For alpha values (covariate coefficients)
                                               p.values.alpha.fdr.corrected<- t(apply(cox.coef1[1:b, 37:45], 1, function(x) {p.adjust(x, method="fdr")}))
                                               p.values.alpha.fdr.corrected.minus.log10<- -log10(p.values.alpha.fdr.corrected)
                                               chi.sq.alpha<- round(apply(cox.coef1[1:b, 28:36], 2, median),1)
                                               
                                               # Draw the boxplot
                                                     pdf(file = paste(working.dir,"/Norwood/Results/coxcoef1.pdf",sep=""), height = 5, width = 10)
                                                         par(mar=c(15,4,2,2))
                                                         boxplot(p.values.alpha.fdr.corrected.minus.log10, las=2, frame = F, axes=F, col = "lightblue", lwd=3, ylim=c(0,10))
                                                         axis(side = 1, at = 1:9, padj = -0.5, lwd = 4, cex.axis = 1, line = 0, labels = c("Weight", "", "Diagnosis", "", "Ventricular function", "", "Shunt type", "", "ECMO"))
                                                         axis(side = 1, at = c(2,4,6,8), padj = 1, lwd = 0, cex.axis = 1, labels = c("Era", "Systemic ventricle", "AVVR","Interdigitation"))
                                                         axis(side = 2, at = seq(0,10,2), las = 1, lwd = 4, cex.axis = 1.5, line = -1)
                                                         mtext(text=expression("-log"[10]*" p-values"), side = 2, cex=2, line = 1.5)
                                                         mtext(text="Covariate", side = 1, cex=2, line = 4)
                                                         abline(h=-log10(0.05), lwd=4, lty=5, col="red")
                                                     dev.off()         
                                           
                                      # For gamma values (transition coefficients)
                                               p.values.gamma.median<- apply(cox.coef2[1:b, 157:195], 2, function(x) {median(x)})
                                               p.values.gamma.fdr.corrected<- p.adjust(p.values.gamma.median, method="fdr")
                                               p.values.gamma.fdr.corrected.minus.log10<- -log10(p.values.gamma.fdr.corrected)
                                               p.values.gamma.fdr.corrected.minus.log10.rescaled<- rescale.yaxis(p.values.gamma.fdr.corrected.minus.log10)
                                               chi.sq.gamma<- round(apply(cox.coef2[1:b, 118:156], 2, median),3)
                                               
                                                     pdf(file = paste(working.dir,"/Norwood/Results/coxcoef2.pdf",sep=""), height = 5, width = 10)
                                                         par(mar=c(6,6,4,2))
                                                         boxplot(p.values.gamma.fdr.corrected.minus.log10.rescaled, las=2, frame = F, axes=F, col="red", lwd=1, ylim=c(0,2), cex=0.3)
                                                         axis(side = 1, at = seq(1,39,2), labels = seq(1,39,2), padj = -0.3, las = 1, lwd = 4, cex.axis = 1, line = 0)
                                                         axis(side = 1, at = seq(2,38,2), labels = seq(2,38,2), padj = 0.3, las = 1, lwd = 0, cex.axis = 1)
                                                         mtext(text=expression("-log"[10]*" p-values"), side = 2, cex=2, line = 3)
                                                         mtext("Transition Number", side = 1, cex=2, line = 3)
                                                         abline(h=rescale.yaxis(-log10(0.05)), lwd=4, lty=5, col="red")
                                                         
                                                         # Y-axis to be re-scaled
                                                             y.labels<- c(seq(0,1,0.1), 2, 3, seq(4,10,1))
                                                             y.labels.trans<- rescale.yaxis(y.labels)
                                                             axis(side = 2, at = y.labels.trans, labels = y.labels, las = 1, lwd = 4, cex.axis = 1.5, line = -1)
                                                    dev.off()  
                                 }
                               
                               # Save data and iteration number
                                     write.table(cox.coef1, file=paste(folder.rr, "coxcoef1.csv", sep="/"), col.names=FALSE, row.names=FALSE)
                                     write.table(cox.coef2, file=paste(folder.rr, "coxcoef2.csv", sep="/"), col.names=FALSE, row.names=FALSE)
                                     write.table(loglik.mat, file=paste(folder.rr, "loglik.csv", sep="/"), col.names=FALSE, row.names=FALSE)
                                     write.table(b, file=paste(folder.rr, "rr_iteration_number.csv", sep="/"), col.names=FALSE, row.names=FALSE)
                                     
                                     if (b>1) {
                                               write.table(round(chi.sq.alpha,1), file = paste(folder.rr, "coxcoef1chisq.csv", sep="/"))
                                               write.table(apply(p.values.alpha.fdr.corrected.minus.log10, 2, median), file = paste(folder.rr, "coxcoef1pvalsfdrcorrected.csv", sep="/"))
                                               write.table(round(chi.sq.gamma,1), file = paste(folder.rr, "coxcoef2chisq.csv", sep="/"))
                                               write.table(apply(p.values.gamma.fdr.corrected.minus.log10, 2, median), file = paste(folder.rr, "coxcoef2pvalsfdrcorrected.csv", sep="/"))
                                               }
                                     
                                     
                              if (b == s) {break}              
                    }
                     
            
                    # Chi-squared and p-values for Table 3 for paper
                        # Chi-squared values
                              round(apply(cox.coef1[1:1000,28:36],2,mean),1)
                                       
                        # P-values
                              round(p.adjust(sapply(apply(cox.coef1[1:1000,28:36],2,mean), function(x) {pchisq(x, 1, lower.tail = FALSE)}), method="fdr"),5)
                      
                              
                              
                              
                              
                              
      # ------------- Demonstrate the direction of effects --------------------
                              
                              
      # Negative alpha coefficients correspond to better survival if you have more of the independent variable
      
            interpretation.matrix<- matrix(c("Preterm1","If the baby is term","If the baby is preterm",
                                             "Weight","If the baby has low birthweight", "If the baby has higher birthweight",
                                             "Ndays","If the operation is nearer 2000", "If the operation was done more recently",
                                             "era2","If the operation is from era 1 compared to era 2","If the operation was from era 2 when compared to era 1",
                                             "era3","If the operation is from era 1 compared to era 3","If the operation was from era 3 when compared to era 1",
                                             "DGcat1", "If the diagnosis is not HLHS", "If the diagnosis is HLHS",
                                             "SysV1","If the systemic ventricle is a left ventricle","If the systemic ventricle is a right ventricle",
                                             "SystV_f","If the systemic ventricle has good function","If the systemic ventricle has poor function",
                                             "AVVRcat1","If there is less AVVR","If there is severe AVVR",
                                             "SanovsBT1","If there is a Sano shunt","If there is a BT shunt",
                                             "Interdigitation1", "If surgical technique does not use interdigitation", "If surgical technique uses interdigitation",
                                             "N1ECMO1","If the baby does not go on ECMO","If the baby receives ECMO"),
                                           dimnames=list(paste("Covariate",1:12,sep=""), c("Covariate","Interpretation if covariate is zero", "Interpretation if covariate is one")),
                                           ncol=3, byrow=T)
            

      # Correct for multiple testing
            
            # Which covariates are significant?
                  significant.covariates.multiple.testing.corrected<- p.adjust(summary(rr1$cox.itr1)$coefficients[,6], method = "fdr")
                  significant.covariates<- which(significant.covariates.multiple.testing.corrected < 0.05)
            
            # Which transitions are significant?
                  significant.transitions.multiple.testing.corrected<- p.adjust(summary(rr1$cox.itr2)$coefficients[,6], method = "fdr")
                  significant.transitions<- which(significant.transitions.multiple.testing.corrected < 0.05)
                  
            
            interp.rows<- match(rownames(rr1$Alpha)[significant.covariates], interpretation.matrix[,1])
            interp.cols<- 2.5-sign(rr1$Alpha[significant.covariates])/2
            
      # Establish which way round the effects go
            # Step 1: Define poor survival and good survival cases
                  clinical.risk<- matrix(c(2.5,3.5,
                                           677,7699,
                                           1,0,
                                           1,0,
                                           3,0,
                                           3,0,
                                           1,0,
                                           0,1,
                                           1,0), ncol=9, dimnames = list(c("Poor survival", "Good survival"), c("Weight","Ndays","DGcat","SysV","Syst_f","AVVRcat","SanovsBT","Interdigitation","N1ECMO")))
            
                  clinical.risk.ce<- t(apply(clinical.risk,1, function(x) x - rr1$center))
                  clinical.risk.ce.sc<- t(apply(clinical.risk.ce, 1, function(x) {x / rr1$scale}))
                  
            # Step 2: Calculate the prognostic score (ALPHA[T] . Z) for each
                  prognostic.scores<- clinical.risk.ce.sc %*% rr1$Alpha
                  prognostic.scores
            
            # Step 3: Define the effect direction
                  direction<- c("worse","better")[(1 + as.numeric(prognostic.scores[1,1] > prognostic.scores[2,1]))]
                  

      # Alpha: determines the effects of individual covariates 
            cat(paste("\n The baby has a",direction,"prognosis if: \n ----------------------------------"))
            for (i in 1:length(interp.rows)) {cat("\n", interpretation.matrix[interp.rows[i], interp.cols[i]],
                                                  "   p ~",
                                                  significant.covariates.multiple.testing.corrected[significant.covariates[i]])}


      # Gamma: Determines the size of the effect of the prognostic score on transition k
            t(rr1$Gamma) # How much do the coefficients influence these transitions?
            summary(rr1$cox.itr2)$coefficients[significant.transitions,]
            
          
                

 
            
            # -------------- Define direct, direct+, indirect and dead-before-TCPC pathways ------------------
                  
                  # Sample paths to estimate how many get through a SV pathway with no intervention
                          tvec<- seq(0, 20, 1/48)
                          summary.paths <- mssample.TD(Haz = msf.Allc$Haz, trans = tmat.Allc, clock = "reset", output = "path", tvec = tvec, M = 500)
                          paths.tmat.Allc<- paths(tmat.Allc)
                          
                  
                  # Direct: N, Gl, T 
                          direct.vector.list<- list(c(1,rep(NA,11)),
                                                    c(1,4,rep(NA,10)),
                                                    c(1,4,9,rep(NA,9)),
                                                    c(1,4,9,14,rep(NA,8)))
                          
                          path.direct<- which(apply(paths(tmat.Allc), 1, function(x) {any(sapply(direct.vector.list, FUN = identical, x))}) == TRUE)
                          col.direct<- match(paste("ppath",path.direct,sep=""), colnames(summary.paths))
                  
                  # Directplus: N, Gl/Gl-R, N/N-R      
                          directplus.vector.list<- list(c(1,5,rep(NA,10)),   # N + Gl-R
                                                        c(1,5,9,rep(NA,9)),  # N + Gl-R + TCPC
                                                        c(1,5,9,14,rep(NA,8)),  # N + Gl-R + TCPC + died
                                                        c(1,5,10,rep(NA,9)),  # N + Gl-R + TCPC + died
                                                        c(1,5,10,14,rep(NA,8)),  # N + Gl-R + TCPC + died
                                                        c(1,4,10,rep(NA,9)),    # N + Gl + TCPC-R
                                                        c(1,4,10,14,rep(NA,8))) # N + Gl + TCPC-R + died
                          path.directplus<- which(apply(paths(tmat.Allc), 1, function(x) {any(sapply(directplus.vector.list, FUN = identical, x))}) == TRUE)
                          col.directplus<- match(paste("ppath",path.directplus,sep=""), colnames(summary.paths))
                  
                  # Indirect: Additional re-intervention at any stage
                  
                          # Which paths have a reintervention?
                                  path.reintervention<- which(apply(paths.tmat.Allc, 1, function(x) {any(c(2,3,6,7,8,11,12,13) %in% x)}) == TRUE)
                          
                          # Now need to pick out two groups:
                          #     1. those who had (a) a reintervention (b) didn't die
                          #     2. those who had (a) a reintervention (b) TCPC first (c) then died
                          
                                  # What is the index of dying for all cases?
                                          path.dead.index<- apply(paths(tmat.Allc),1, function(x) {match(14, x)})
                                  
                                  #  What is the index of the TCPC for all cases?
                                          path.TCPC.index<- apply(paths.tmat.Allc, 1, function(x) {match(c(9,10), x)})
                                          for (i in 1:438) {if (is.na(path.TCPC.index[1,i]) == TRUE) {path.TCPC.index[1,i]<- path.TCPC.index[2,i]}}
                                          path.TCPC.index<- path.TCPC.index[1,]
                                          
                                  # Group 1
                                          path.didnt.die<- which(apply(paths(tmat.Allc),1, function(x) {c(14) %in% x}) == FALSE)
                                          path.indirect.group1<- intersect(path.reintervention, path.didnt.die)
                                          
                                  # Group 2
                                          path.TCPC.before.dead<- which(path.TCPC.index < path.dead.index)
                                          path.indirect.group2<- intersect(path.reintervention, path.TCPC.before.dead)
                                  
                          path.indirect<- sort(c(path.indirect.group1, path.indirect.group2))
                          col.indirect<- match(paste("ppath",path.indirect,sep=""), colnames(summary.paths))
                          
                  
                  # Died before TCPC
                  
                          # Who died before TCPC?
                                  path.did.die<- which(apply(paths(tmat.Allc),1, function(x) {c(14) %in% x}) == TRUE)
                                  path.did.not.have.TCPC<- which(apply(paths.tmat.Allc, 1, function(x) {any(c(9,10) %in% x)}) == FALSE)
                                  path.died.before.TCPC<- intersect(path.did.die, path.did.not.have.TCPC)
                                  col.died.before.TCPC<- match(paste("ppath",path.died.before.TCPC,sep=""), colnames(summary.paths))
                                  
                  
                  
                  # Collect together direct, directplus, indirect and died-before-TCPC pathways into one dataframe
                          variables.mean<- data.frame(direct = apply(data.matrix(summary.paths[,col.direct]),1,sum),
                                                      directplus = apply(data.matrix(summary.paths[,col.directplus]),1,sum),
                                                      indirect = apply(data.matrix(summary.paths[,col.indirect]),1,sum),
                                                      incomplete = apply(data.matrix(summary.paths[,col.died.before.TCPC]),1,sum))
                          
                          pathway.cols<- list(direct = col.direct, directplus = col.directplus, indirect = col.indirect, died.before.TCPC = col.died.before.TCPC)
                          
                          
                  # Check all the paths are accounted for
                          missing.paths<- setdiff(1:438, c(path.direct, path.directplus, path.indirect, path.died.before.TCPC))
                          missing.paths
                          par(mfrow=c(1,1))
                          plot(1:961, apply(variables.mean,1,sum), ylim = c(0,1), xlab="Timepoints", ylab="Fraction of patients accounted for by the 4 pathways")
                          
                  
            
           

            # ----------------------------------------------------------

                          
                          
                          
                          
            # Finally, use simulation to calculate transition probabilities not assuming Markov
                              
            
                              bootstrap.reps<- 1000
                              
                    # Reset variables for outcome estimation
                              time<- seq(0,20,length.out=961)
                              surv.times<- c(1/12, 6/12, 1, 5, 10, 20)
                              msfsample.Allc.path<- list()
                              
                              survival.overall<- reintfreesurvival.overall<- pathwaysurvival.overall<- matrix(0, nrow=961, ncol=12, dimnames=list(1:961, paste(rep(c("All","Era1","Era2","Era3"), each=3), c("Mu","UpperCI","LowerCI"), sep="_")))
                              survival.output<- reintfreesurvival.output<- pathwaysurvival.output<- matrix(0, nrow=24, ncol=5, dimnames = list(paste(rep("era",24), rep(0:3, each = 6),sep=""), c("Era","Time","MeanSurv","LowerCI","UpperCI")))
                              
                              sampled.outcomes<- matrix(0, nrow=bootstrap.reps, ncol=16, dimnames = list(1:bootstrap.reps, c(paste(rep(paste("era",0:3,sep="_"), each = 4), c("Direct","Direct+","Indirect","Died before TCPC"), sep="."))))
                              pathwaysampled.outcomes<- matrix(0, nrow=4, ncol=16, dimnames = list(1:4, c(paste(rep(paste("era",0:3,sep="_"), each = 4), c("Direct","Direct+","Indirect","Died.before.TCPC"), sep="."))))
                              
                              pathwaysurvival.by.bootstrap.reps.direct<- pathwaysurvival.by.bootstrap.reps.directplus<- pathwaysurvival.by.bootstrap.reps.indirect<- pathwaysurvival.by.bootstrap.reps.died.before.TCPC<- matrix(0, nrow=3844, ncol=1000, dimnames = list(paste(rep("era",3844),rep(0:3, each=961),sep="_"), paste("b",1:1000,sep="")))
                              stacked_values.mean.by.bootstrap.reps<- stacked_values.lowerCI.by.bootstrap.reps<- stacked_values.upperCI.by.bootstrap.reps<- list()
                              b<- 0
                              
                              survival.by.bootstrap.reps<- reintfreesurv.by.bootstrap.reps<- matrix(0, nrow=961*4, ncol=bootstrap.reps)
                              boxcox.vector<- rep(1,961)
                              
                    # Define folders to save the result in
                              
                              folder.survival<-paste("/Users/", tims.name, "Data/Era_bootstrap", sep="")
                              folder.reintfreesurv<-paste("/Users/", tims.name, "Data/Reintervention_bootstrap", sep="")
                              folder.pathwaysurv<-paste("/Users/", tims.name, "Data/Pathway_bootstrap", sep="")
                              
                              
                    # Initiate some blank variables in both folders
                              save.list<- list(reintfreesurvival.overall = reintfreesurvival.overall, reintfreesurvival.output = reintfreesurvival.output, survival.by.bootstrap.reps = survival.by.bootstrap.reps, b = b)
                              save(save.list, file=paste(folder.reintfreesurv, "Reintsurvloopvariables.RData", sep="/"))
                              
                              save.list<- list(survival.overall = survival.overall, survival.output = survival.output, reintfreesurv.by.bootstrap.reps = reintfreesurv.by.bootstrap.reps, b = b)
                              save(save.list, file=paste(folder.survival, "Survivalloopvariables.RData", sep="/"))
                             
                              save.list<- list(pathwaysurvival.overall = pathwaysurvival.overall, pathwaysurvival.output = pathwaysurvival.output, b = b)
                              for (route in 1:4) {save(save.list, file=paste(folder.pathwaysurv, "/Pathwaysurvloopvariables_",route,".RData", sep=""))}
                              
                              write.table(0, file=paste(folder.survival, "iteration_number.csv", sep="/"), col.names=FALSE, row.names=FALSE)
                              write.table(0, file=paste(folder.reintfreesurv, "iteration_number.csv", sep="/"), col.names=FALSE, row.names=FALSE)
                              write.table(0, file=paste(folder.pathwaysurv, "iteration_number.csv", sep="/"), col.names=FALSE, row.names=FALSE)
                              
                              
                        # Progress bar
                              pb <- txtProgressBar(min = 0, max = bootstrap.reps, style = 3, char = "=")
                              bw<- 1
                              
                        # Define the outcomes for the 'survival' and 'survival without reintervention' endpoints                              
                            # Survival outcome (14)
                                  path.number.ends.in.death<- which(apply(paths(tmat.Allc), 1, function(x) {c(14) %in% x}) == TRUE)
                                  path.number.does.not.end.in.death<- setdiff(1:438, path.number.ends.in.death)
                                  col.number.does.not.end.in.death<- path.number.does.not.end.in.death + 1
                                  col.number.ends.in.death<- path.number.ends.in.death + 1

                            # Intervention-free survival (2,3,5,6,7,8,10,11,12,13) or death (14)
                                  
                                  path.number.includes.reintervention<- which(apply(paths(tmat.Allc), 1, function(x) {any(c(2,3,5,6,7,8,10,11,12,13) %in% x)}) == TRUE)
                                  path.number.does.not.include.reintervention<- setdiff(1:438, path.number.includes.reintervention)
                                  col.number.does.not.include.reintervention<- path.number.does.not.include.reintervention + 1
                                  col.number.includes.reintervention<- path.number.includes.reintervention + 1
                                  



              
                 
          # Define the cumulative transition hazard models () 
              for (era in 0:3){
                        cat("...",era)
                
                    # Find the IDs which correspond to the era in question
                        if (era == 0) {ids<- 1:no.cases} else {ids<- unique(msm2.Allc.covs$id[which(msm2.Allc.covs$era == era)])}
                
                    # Find all the rows from these patients
                        era.rows<- which(msm2.Allc.covs$id %in% ids)
                                
                    # Model the baseline hazard from these era of patients
                        coxmod<- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = msm2.Allc.covs[era.rows,], method = "breslow")
                                
                    # Compute the cumulative transition hazard for each transition
                        newpts<- data.frame(msm2.Allc.covs[era.rows,])
                        colnames(newpts)[which(colnames(newpts) == "trans")]<- "strata(trans)"
                        cum_haz<- msfit(coxmod, newdata = newpts, variance = FALSE, vartype = "aalen", trans = tmat.Allc)$Haz
                        save(cum_haz, file=paste("Norwood/Data/Cum_haz_era_",era,".RData",sep=""))
                        
              }
                                
                      
                      

      repeat {
                      b<- read.table(file=paste(folder.survival, "iteration_number.csv", sep="/"), col.names=FALSE)[1,]
                      b<- b + 1
                      setTxtProgressBar(pb, b)
                      reintervention.dataframe<- survival.dataframe<- NULL
                    
                                              
                        # Simulate occupancy probabilities from predicted cumulative hazards
                                      for (era in 0:3)
                                      {
                                        load(file=paste(working.dir, "/Norwood/Data/Cum_haz_era_",era,".RData",sep=""))
                                        msfsample.Allc.path[[era+1]]<- mssample.TD(Haz = cum_haz, trans = tmat.Allc, clock = "forward", output = "path", tvec = seq(0,20,1/48), M = 100)
                                      }
                          
                      
                      
                    # Survival outcomes
                      
                                      for (era in 0:3)
                                      {
                                          
                                                # Load the survival data
                                                      load(file=paste(folder.survival, "Survivalloopvariables.RData",sep="/"))
                                                      columns.by.era<- (1+era*3):(3+era*3)
                                                   
                                            # # Approach 2: Update the mean and SD by storing all the data
                                                      
                                                      start<- 1 + era*961
                                                      stop<- start + 960
                                                      survival.by.bootstrap.reps[start:stop,b]<- apply(msfsample.Allc.path[[era+1]][,col.number.does.not.end.in.death],1,sum)
                                                      
                                                      # Measure median and 95% CIs from bootstraps
                                                          if ((b/100) == round(b/100)) {boxcox.vector<- apply(survival.by.bootstrap.reps[start:stop,1:b], 1, function(x) {BoxCox.lambda(x)})}
                                                          for (k in 1:961){survival.overall[k,columns.by.era]<- mean.sd.after.boxcox(survival.by.bootstrap.reps[(start:stop)[k], 1:b], boxcox.vector[k])}
                                                      
                                                          new.rows<- data.frame(time = time, 
                                                                                mean = 100 * survival.overall[,columns.by.era[1]],
                                                                                ymin = 100 * survival.overall[,columns.by.era[2]],
                                                                                ymax = 100 * survival.overall[,columns.by.era[3]], 
                                                                                Era = as.factor(c("All", "2001 - 2008", "2008 - 2014", "2014 - 2021")[era+1]), era = as.factor(era))
                                                        
                                                        # Smooth the results and cap at 0 and 100
                                                                new.rows<- smoother(new.rows, bw)
                                                                survival.dataframe<- rbind(survival.dataframe, new.rows) 
                                                                
                                                # Update the survival.output matrix         
                                                      for (i in surv.times)
                                                           {
                                                                m<- which.min(abs(time - i))
                                                                surv.time.index<- match(i,surv.times)
                                                                output.row<- era*6 + surv.time.index
                                                                survival.output[output.row,]<- c(era, i, round(100* survival.overall[m, columns.by.era],1))
                                                            }
                                                      
                                                 # Save progress
                                                        save.list<- list(survival.overall = survival.overall, survival.output = survival.output, survival.by.bootstrap.reps, b = b)
                                                        save(save.list, file=paste(folder.survival, "Survivalloopvariables.RData", sep="/"))
                                          }                              
                                        
                                        # Save some individual files
                                                write.csv(survival.output, file=paste(folder.survival, "survival_output.csv", sep="/"))
                                                write.csv(survival.dataframe, file=paste(folder.survival, "survival_dataframe.csv", sep="/"))
                                                write.table(b, file=paste(folder.survival, "iteration_number.csv", sep="/"), col.names=FALSE, row.names=FALSE)
                                        
                                        # Draw the survival curves
                                                rows.by.era<- which(survival.dataframe$era != 0)
                                                
                                                g1<- ggplot(data=survival.dataframe[rows.by.era,], aes(x=time, y=mean_smoothed, color=Era, fill=Era)) + 
                                                  labs(x = "Time (years)", y = "Transplant-Free Survival (%)") +
                                                  geom_ribbon(show.legend = TRUE, aes(ymin=ymin_smoothed, ymax=ymax_smoothed, fill = Era), linetype = "dashed", alpha=.3) +
                                                  geom_line(show.legend = TRUE, size = 2, alpha=1) + guides(fill = "none") + ylim(0, 101) +
                                                  scale_fill_manual(name = 'Era',
                                                                    values = c("All" = "black", "2001 - 2008" = "#009E73","2008 - 2014" = "#2c45ff", "2014 - 2021" = "#c10012"),
                                                                    drop = FALSE) +
                                                  scale_color_manual(name = "Era",
                                                                     values = c("All" = "black", "2001 - 2008" = "#126600","2008 - 2014" = "#0a2472", "2014 - 2021" = "#c10012"),
                                                                     drop = FALSE) +
                                                  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                     axis.line = element_line(colour = "black"), legend.position="bottom", axis.text.x = element_text(size = 14),
                                                                     axis.text.y = element_text(size = 14))
                                                
                                                g3<- ggplot(data=survival.dataframe[-rows.by.era,], aes(x=time, y=mean_smoothed, color=factor(Era, levels = c("All","2001 - 2008", "2008 - 2014", "2014 - 2021")),
                                                                                                                                 fill=factor(Era, levels = c("All","2001 - 2008", "2008 - 2014", "2014 - 2021")))) + 
                                                  labs(x = "Time (years)", y = "Transplant-Free Survival (%)") +
                                                  geom_ribbon(show.legend = TRUE, aes(ymin=ymin_smoothed, ymax=ymax_smoothed, fill = Era), linetype = "dashed", alpha=.3) +
                                                  geom_line(show.legend = TRUE, size = 2, alpha=1) + guides(fill = "none") + ylim(0,101) +
                                                  scale_fill_manual(name = 'Era',
                                                                    values = c("All" = "black", "2001 - 2008" = "#009E73","2008 - 2014" = "#2c45ff", "2014 - 2021" = "#c10012"),
                                                                    drop = FALSE) +
                                                  scale_color_manual(name = 'Era',
                                                                     values = c("All" = "black", "2001 - 2008" = "#126600","2008 - 2014" = "#0a2472", "2014 - 2021" = "#c10012"),
                                                                     drop = FALSE) +
                                                  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                     axis.line = element_line(colour = "black"), legend.position="bottom", axis.text.x = element_text(size = 14),
                                                                     axis.text.y = element_text(size = 14))
                                                
                              # Re-intervention-free survival outcome
                              
                                            for (era in 0:3)
                                            {   
                                                        # Load the survival data
                                                            load(file=paste(folder.reintfreesurv, "Reintsurvloopvariables.RData",sep="/"))
                                                            columns.by.era<- (1+era*3):(3+era*3)
                                              
                                                         # Approach 2: Update the mean and SD by storing all the data
                                                            start<- 1 + era*961
                                                            stop<- start + 960
                                                            reintfreesurv.by.bootstrap.reps[start:stop,b]<- apply(msfsample.Allc.path[[era+1]][,col.number.does.not.include.reintervention],1,sum)
                                                            
                                                            # Measure median and 95% CIs from bootstraps
                                                                for (k in 1:961){reintfreesurvival.overall[k,columns.by.era]<- mean.sd.after.boxcox(reintfreesurv.by.bootstrap.reps[(start:stop)[k], 1:b], boxcox.vector[k])}
                                                                
                                                          # Calculate new rows, smooth and rbind to reintervention.dataframe
                                                                new.rows<- data.frame(time = time, 
                                                                                      mean = 100 * reintfreesurvival.overall[,columns.by.era[1]],
                                                                                      ymin = 100 * reintfreesurvival.overall[,columns.by.era[2]],
                                                                                      ymax = 100 * reintfreesurvival.overall[,columns.by.era[3]], 
                                                                                      Era = as.factor(c("All", "2001 - 2008", "2008 - 2014", "2014 - 2021")[era+1]), era = as.factor(era))
                                                                
                                                                new.rows<- smoother(new.rows, bw)
                                                                reintervention.dataframe<- rbind(reintervention.dataframe, new.rows) 
                                                                
                                                        # Add discrete survival figures
                                                            for (i in surv.times)
                                                            {
                                                               m<- which.min(abs(time - i))
                                                               surv.time.index<- match(i,surv.times)
                                                               output.row<- era*6 + surv.time.index
                                                               reintfreesurvival.output[output.row,]<- c(era, i, round(100 * reintfreesurvival.overall[m,columns.by.era],1))
                                                             }
                                                            
                                                            # Save progress
                                                            save.list<- list(reintsurvival.overall = reintfreesurvival.overall, 
                                                                             reintfreesurvival.output = reintfreesurvival.output,
                                                                             reintfreesurv.by.bootstrap.reps = reintfreesurv.by.bootstrap.reps, b = b)
                                                            
                                                            save(save.list, file=paste(folder.reintfreesurv, "Reintsurvloopvariables.RData", sep="/"))
                                                            
                                                                          
                                            }
                                            
                                            # Save some individual files
                                                    write.csv(reintfreesurvival.output, file=paste(folder.reintfreesurv, "reintfreesurvival_output.csv", sep="/"))
                                                    write.csv(reintervention.dataframe, file=paste(folder.reintfreesurv, "reintervention_dataframe.csv", sep="/"))
                                                    write.table(b, file=paste(folder.reintfreesurv, "iteration_number.csv", sep="/"), col.names=FALSE, row.names=FALSE)
                                                    
                                            rows.by.era<- which(reintervention.dataframe$era != 0)
                                                    
                                            g2<- ggplot(data=reintervention.dataframe[rows.by.era,], aes(x=time, y=mean_smoothed, color=factor(Era, levels = c("All","2001 - 2008", "2008 - 2014", "2014 - 2021")))) + 
                                                labs(x = "Time (years)", y = "Reintervention-Free (%)") +
                                                geom_ribbon(show.legend = TRUE, aes(ymin=ymin_smoothed, ymax=ymax_smoothed, fill = Era), linetype = "dashed", alpha=.3) +
                                                geom_line(show.legend = TRUE, size = 2, alpha=1) + guides(fill = "none") + ylim(0,101) +
                                                scale_fill_manual(name = 'Era',
                                                                  values = c("All" = "black", "2001 - 2008" = "#009E73","2008 - 2014" = "#2c45ff", "2014 - 2021" = "#c10012"),
                                                                  drop = FALSE) +
                                                scale_color_manual(name = 'Era',
                                                                  values = c("All" = "black", "2001 - 2008" = "#126600","2008 - 2014" = "#0a2472", "2014 - 2021" = "#c10012"),
                                                                  drop = FALSE) +
                                                theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 axis.line = element_line(colour = "black"), legend.position="bottom", axis.text.x = element_text(size = 14),
                                                                 axis.text.y = element_text(size = 14))
                                                
                                            g4<- ggplot(data=reintervention.dataframe[-rows.by.era,], aes(x=time, y=mean_smoothed, color=factor(Era, levels = c("All","2001 - 2008", "2008 - 2014", "2014 - 2021")))) + 
                                                labs(x = "Time (years)", y = "Reintervention-Free (%)") +
                                                geom_ribbon(show.legend = TRUE, aes(ymin=ymin_smoothed, ymax=ymax_smoothed, fill = Era), linetype = "dashed", alpha=.3) +
                                                geom_line(show.legend = TRUE, size = 2, alpha=1) + guides(fill = "none") + ylim(0,101) +
                                                scale_fill_manual(name = "Era",
                                                                  values = c("All" = "black", "2001 - 2008" = "#009E73","2008 - 2014" = "#2c45ff", "2014 - 2021" = "#c10012"),
                                                                  drop = FALSE) +
                                                scale_color_manual(name = "Era",
                                                                   values = c("All" = "black", "2001 - 2008" = "#009E73","2008 - 2014" = "#2c45ff", "2014 - 2021" = "#c10012"),
                                                                   drop = FALSE) +
                                                theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 axis.line = element_line(colour = "black"), legend.position="bottom", axis.text.x = element_text(size = 14),
                                                                 axis.text.y = element_text(size = 14))
                                                
                      # Plot joint-pane survival curves
                            suppressMessages(ggarrange(g3, g1, ncol=2, common.legend = TRUE, legend = "bottom", labels = c("A","B"), hjust = -3.5, font.label = list(size = 24)) %>% 
                                               ggexport(filename = paste("/Users/", tims.name, "Figure3AB_Survival Overall.pdf", sep="")))
                            
                      # Plot individual survival curves
                            suppressMessages(ggexport(g3, filename = paste("/Users/", tims.name, "Figure3A_Survival Overall.pdf", sep="")))
                            suppressMessages(ggexport(g1, filename = paste("/Users/", tims.name, "Figure3B_Survival Overall.pdf", sep="")))
                      
                      # Plot joint-pane reintervention-free survival curves
                            suppressMessages(ggarrange(g4, g2, ncol=2, common.legend = TRUE, legend = "bottom", labels = c("A","B"), hjust = -3.5, font.label = list(size = 24)) %>% 
                                               ggexport(filename = paste("/Users/", tims.name, "Figure4AB_Survival by Era.pdf",sep="")))
                            
                      
                      
                      # Use the sampled paths to estimate how many get through a SV pathway with no intervention
                      
                                                      
                                        era<- 0
                                        pathwaysurvival.dataframe.direct<- pathwaysurvival.dataframe.directplus<- pathwaysurvival.dataframe.indirect<- pathwaysurvival.dataframe.died.before.TCPC<- NULL
                                        pathwaysurvival.overall.direct<- pathwaysurvival.overall.directplus<- pathwaysurvival.overall.indirect<- pathwaysurvival.overall.died.before.TCPC<- matrix(0, nrow=961, ncol=12, dimnames = list(paste("time",1:961,sep=""), paste(rep(paste(rep("Era",4),0:3,sep="_"),each=3), rep(c("Mean","LowerCI","UpperCI"),4), sep="_")))
                                        pathwaysurvival.output.direct<- pathwaysurvival.output.directplus<- pathwaysurvival.output.indirect<- pathwaysurvival.output.died.before.TCPC<- matrix(0, nrow=24, ncol=5, dimnames = list(paste(rep("time",24), round(rep(surv.times,4),1)), c("Era","Survivaltime","Mean","UpperCI","LowerCI")))
                                        
                      
                                        for (era in 0:3)
                                              {
                                                                  columns.by.era<- (1+era*3):(3+era*3)
                                                                  start<- 1 + era*961
                                                                  stop<- start + 960
                                                                  
                                                            # Route 1: Direct
                                                                        pathwaysurvival.by.bootstrap.reps.direct[start:stop,b]<- apply(data.matrix(msfsample.Allc.path[[era+1]][,pathway.cols$direct]), 1, sum)
                                                                  
                                                                  # Measure median and 95% CIs from bootstraps
                                                                        pathwaysurvival.overall.direct[,columns.by.era]<- t(apply(data.matrix(pathwaysurvival.by.bootstrap.reps.direct[start:stop,1:b]), 1, function(x) {quantile(x, probs = c(0.5, 0.025, 0.975))}))
                                                                        
                                                                        new.rows.direct<- data.frame(time = time, 
                                                                                              mean = 100 * pathwaysurvival.overall.direct[,columns.by.era[1]],
                                                                                              ymin = 100 * pathwaysurvival.overall.direct[,columns.by.era[2]],
                                                                                              ymax = 100 * pathwaysurvival.overall.direct[,columns.by.era[3]], 
                                                                                              Era = as.factor(c("All", "2001 - 2008", "2008 - 2014", "2014 - 2021")[era+1]), era = as.factor(era))
                                                                        
                                                            # Route 2: Directplus
                                                                        pathwaysurvival.by.bootstrap.reps.directplus[start:stop, b]<- apply(msfsample.Allc.path[[era+1]][,pathway.cols$directplus], 1, sum)
                                                                  
                                                                  # Measure median and 95% CIs from bootstraps
                                                                        pathwaysurvival.overall.directplus[,columns.by.era]<- t(apply(data.matrix(pathwaysurvival.by.bootstrap.reps.directplus[start:stop,1:b]), 1, function(x) {quantile(x, probs = c(0.5, 0.025, 0.975))}))
                                                                        
                                                                        new.rows.directplus<- data.frame(time = time, 
                                                                                                     mean = 100 * pathwaysurvival.overall.directplus[,columns.by.era[1]],
                                                                                                     ymin = 100 * pathwaysurvival.overall.directplus[,columns.by.era[2]],
                                                                                                     ymax = 100 * pathwaysurvival.overall.directplus[,columns.by.era[3]], 
                                                                                                     Era = as.factor(c("All", "2001 - 2008", "2008 - 2014", "2014 - 2021")[era+1]), era = as.factor(era))
                                                                  
                                                                
                                                            # Route 3: Indirect
                                                                        pathwaysurvival.by.bootstrap.reps.indirect[start:stop, b]<- apply(msfsample.Allc.path[[era+1]][,pathway.cols$indirect], 1, sum)
                                                                  
                                                                  # Measure median and 95% CIs from bootstraps
                                                                        pathwaysurvival.overall.indirect[,columns.by.era]<- t(apply(data.matrix(pathwaysurvival.by.bootstrap.reps.indirect[start:stop,1:b]), 1, function(x) {quantile(x, probs = c(0.5, 0.025, 0.975))}))
                                                                  
                                                                        new.rows.indirect<- data.frame(time = time, 
                                                                                                         mean = 100 * pathwaysurvival.overall.indirect[,columns.by.era[1]],
                                                                                                         ymin = 100 * pathwaysurvival.overall.indirect[,columns.by.era[2]],
                                                                                                         ymax = 100 * pathwaysurvival.overall.indirect[,columns.by.era[3]], 
                                                                                                         Era = as.factor(c("All", "2001 - 2008", "2008 - 2014", "2014 - 2021")[era+1]), era = as.factor(era))
                                                                  
                                                                  
                                                            # Route 4: Died.before.TCPC
                                                                        pathwaysurvival.by.bootstrap.reps.died.before.TCPC[start:stop, b]<- apply(msfsample.Allc.path[[era+1]][,pathway.cols$died.before.TCPC], 1, sum)
                                                                  
                                                                  # Measure median and 95% CIs from bootstraps
                                                                        pathwaysurvival.overall.died.before.TCPC[,columns.by.era]<- t(apply(data.matrix(pathwaysurvival.by.bootstrap.reps.died.before.TCPC[start:stop,1:b]), 1, function(x) {quantile(x, probs = c(0.5, 0.025, 0.975))}))
                                                                        
                                                                        new.rows.died.before.TCPC<- data.frame(time = time, 
                                                                                                         mean = 100 * pathwaysurvival.overall.died.before.TCPC[,columns.by.era[1]],
                                                                                                         ymin = 100 * pathwaysurvival.overall.died.before.TCPC[,columns.by.era[2]],
                                                                                                         ymax = 100 * pathwaysurvival.overall.died.before.TCPC[,columns.by.era[3]], 
                                                                                                         Era = as.factor(c("All", "2001 - 2008", "2008 - 2014", "2014 - 2021")[era+1]), era = as.factor(era))
                                                                        
                                                                  
                                                                        
                                                                  # Rescale totals to account for median =/= mean and totals must always add up to 100%
                                                                        rescale.factor<- apply(cbind(new.rows.direct$mean, new.rows.directplus$mean, new.rows.indirect$mean, new.rows.died.before.TCPC$mean),1,sum) / 100
                                                                        new.rows.direct<- cbind(time = new.rows.direct$time, sweep(new.rows.direct[,c("mean", "ymin", "ymax")], 1, rescale.factor, "/"), Era = new.rows.direct$Era, era = new.rows.direct$era)
                                                                        new.rows.directplus<- cbind(time = new.rows.directplus$time, sweep(new.rows.directplus[,c("mean", "ymin", "ymax")], 1, rescale.factor, "/"), Era = new.rows.directplus$Era, era = new.rows.directplus$era)
                                                                        new.rows.indirect<- cbind(time = new.rows.indirect$time, sweep(new.rows.indirect[,c("mean", "ymin", "ymax")], 1, rescale.factor, "/"), Era = new.rows.indirect$Era, era = new.rows.indirect$era)
                                                                        new.rows.died.before.TCPC<- cbind(time = new.rows.indirect$time, sweep(new.rows.died.before.TCPC[,c("mean", "ymin", "ymax")], 1, rescale.factor, "/"), Era = new.rows.died.before.TCPC$Era, era = new.rows.died.before.TCPC$era)
                                                                        
                                                                        
                                                                    
                                                                  pathwaysurvival.dataframe.direct<- rbind(pathwaysurvival.dataframe.direct, smoother(new.rows.direct, bw))
                                                                  pathwaysurvival.dataframe.directplus<- rbind(pathwaysurvival.dataframe.directplus, smoother(new.rows.directplus, bw))
                                                                  pathwaysurvival.dataframe.indirect<- rbind(pathwaysurvival.dataframe.indirect, smoother(new.rows.indirect, bw))
                                                                  pathwaysurvival.dataframe.died.before.TCPC<- rbind(pathwaysurvival.dataframe.died.before.TCPC, smoother(new.rows.died.before.TCPC, bw))
                                                                  
                                                            
                                                            # Add discrete survival figures
                                                                  for (i in surv.times)
                                                                  {
                                                                    m<- which.min(abs(time - i))
                                                                    surv.time.index<- match(i,surv.times)
                                                                    output.row<- era*6 + surv.time.index
                                                                    pathwaysurvival.output.direct[output.row,]<- c(era, i, round(100 * pathwaysurvival.overall.direct[m,columns.by.era],1))
                                                                    pathwaysurvival.output.directplus[output.row,]<- c(era, i, round(100 * pathwaysurvival.overall.directplus[m,columns.by.era],1))
                                                                    pathwaysurvival.output.indirect[output.row,]<- c(era, i, round(100 * pathwaysurvival.overall.indirect[m,columns.by.era],1))
                                                                    pathwaysurvival.output.died.before.TCPC[output.row,]<- c(era, i, round(100 * pathwaysurvival.overall.died.before.TCPC[m,columns.by.era],1))
                                                                  }
                                                            
                                                          
                                                          }                              
                      
                      # Save progress
                            save.list.direct<- list(pathwaysurvival.overall.direct, pathwaysurvival.output.direct, pathwaysurvival.by.bootstrap.reps.direct, b = b)
                            save.list.directplus<- list(pathwaysurvival.overall.directplus, pathwaysurvival.output.directplus, pathwaysurvival.by.bootstrap.reps.directplus, b = b)
                            save.list.indirect<- list(pathwaysurvival.overall.indirect, pathwaysurvival.output.indirect, pathwaysurvival.by.bootstrap.reps.indirect, b = b)
                            save.list.died.before.TCPC<- list(pathwaysurvival.overall.died.before.TCPC, pathwaysurvival.output.died.before.TCPC, pathwaysurvival.by.bootstrap.reps.died.before.TCPC, b = b)
                            
                            save(save.list.direct, file = paste(folder.pathwaysurv, "Pathwaysurvivalloopvariables_direct.RData", sep="/"))
                            save(save.list.directplus, file = paste(folder.pathwaysurv, "Pathwaysurvivalloopvariables_directplus.RData", sep="/"))
                            save(save.list.indirect, file = paste(folder.pathwaysurv, "Pathwaysurvivalloopvariables_indirect.RData", sep="/"))
                            save(save.list.died.before.TCPC, file = paste(folder.pathwaysurv, "Pathwaysurvivalloopvariables_died_before_TCPC_RData", sep="/"))
                            
                      # Save some individual files
                            write.csv(pathwaysurvival.output.direct, file=paste(folder.pathwaysurv, "/pathwaysurvival_output_direct.csv", sep=""))
                            write.csv(pathwaysurvival.output.directplus, file=paste(folder.pathwaysurv, "/pathwaysurvival_output_directplus.csv", sep=""))
                            write.csv(pathwaysurvival.output.indirect, file=paste(folder.pathwaysurv, "/pathwaysurvival_output_indirect.csv", sep=""))
                            write.csv(pathwaysurvival.output.died.before.TCPC, file=paste(folder.pathwaysurv, "/pathwaysurvival_output_died_before_TCPC.csv", sep=""))
                            
                            write.csv(pathwaysurvival.dataframe.direct, file=paste(folder.pathwaysurv, "/pathwaysurvival_dataframe_direct.csv", sep=""))
                            write.csv(pathwaysurvival.dataframe.directplus, file=paste(folder.pathwaysurv, "/pathwaysurvival_dataframe_directplus.csv", sep=""))
                            write.csv(pathwaysurvival.dataframe.indirect, file=paste(folder.pathwaysurv, "/pathwaysurvival_dataframe_indirect.csv", sep=""))
                            write.csv(pathwaysurvival.dataframe.died.before.TCPC, file=paste(folder.pathwaysurv, "/pathwaysurvival_dataframe_died_before_TCPC.csv", sep=""))
                            
                                        
                      
                
                      # Draw Figure 5 (Central Figure) real-time
                                    variables.mean<- variables.lowerCI<- variables.upperCI<- matrix(0, nrow=961, ncol=4, dimnames = list(1:961, c("Direct","Directplus","Indirect","Died before TCPC")))
                                    
                                    for (i in 1:4)
                                    {
                                      route<- c("direct","directplus","indirect","died_before_TCPC")[i]
                                      filename<- paste(working.dir, "Data/pathwaysurvival_dataframe_",route,".csv",sep="")
                                      t<- read.csv(file = filename)[1:961,]
                                      
                                      variables.mean[,i]<- t[,3]
                                      variables.lowerCI[,i]<- t[,5]
                                      variables.upperCI[,i]<- t[,4]
                                    }
                                    
                                    # Create stacked data for areas
                                          stacked_values <- t(apply(variables.mean, 1, cumsum))
                                          stacked_values.lowerCI <- t(apply(variables.lowerCI, 1, cumsum))
                                          stacked_values.upperCI <- t(apply(variables.upperCI, 1, cumsum))
                                          
                                          stacked_values.mean.by.bootstrap.reps[[b]]<- stacked_values
                                          stacked_values.lowerCI.by.bootstrap.reps[[b]]<- stacked_values.lowerCI
                                          stacked_values.upperCI.by.bootstrap.reps[[b]]<- stacked_values.upperCI
                                          
                                    # Base plot
                                          par(mfrow=c(1,1))
                                          filename<- paste(working.dir, "Figure x Transition Plot (also CENTRAL).pdf",sep="")
                                          pdf(file = filename, width = 8, height = 8) 
                                          
                                                plot(NA, xlim=range(tvec), ylim=c(0, 100), xlab="Time (years)", ylab="Probability of being on the pathway", type="n", main = "")
                                                
                                                # Add stacked areas
                                                      for (i in seq_len(ncol(variables.mean))) {
                                                        if (i == 1) {
                                                          polygon(c(tvec, rev(tvec)), c(rep(0, length(tvec)), rev(stacked_values[, i])),
                                                                  col=rainbow(5)[i], border=NA)
                                                        } else {
                                                          polygon(c(tvec, rev(tvec)), c(stacked_values[, i - 1], rev(stacked_values[, i])),
                                                                  col=rainbow(5)[i], border=NA)
                                                        }
                                                      }
                                            
                                            
                                                # Add lines on top of the stacked areas
                                                      cols<- sapply(1:5, function(x) {colorRampPalette(c(rainbow(5)[x], "black"))(100)[80]})
                                                      matlines(tvec, stacked_values, lty=1, lwd=4, col=cols)
                                                      
                                                # Labels (NB: Direct = Treatment, Direct+ = Treatment+, Indirect = Off Treatment)
                                                      stckd.timepoint<- stacked_values[which(tvec == 10),] 
                                                      mean.y.position<- 0.5*(c(0, stckd.timepoint) + c(stckd.timepoint,0))[1:4]
                                                      mean.y.position<- mean.y.position - 2
                                                      for (i in 1:4) {text(10.5,mean.y.position[i], labels = c("Treatment","Treatment Plus","Off Treatment","Died before TCPC")[i], cex=2, adj=c(0.5, 0))}
                                                      
                                          dev.off()
                                          
                                          
                                    # Save stacked values by bootstrap reps
                                          save.list.pathways<- list(stacked_values.mean.by.bootstrap.reps, stacked_values.lowerCI.by.bootstrap.reps, stacked_values.upperCI.by.bootstrap.reps, b = b)
                                          save(save.list.pathways, file = paste(folder.pathwaysurv, "Stacked_values.RData", sep="/"))
                                          
                      write.table(variables.mean, file=paste(folder.pathwaysurv, "variables_mean.csv", sep="/"), col.names=TRUE, row.names=TRUE)
                      write.table(variables.upperCI, file=paste(folder.pathwaysurv, "variables_upperCI.csv", sep="/"), col.names=TRUE, row.names=TRUE)
                      write.table(variables.lowerCI, file=paste(folder.pathwaysurv, "variables_lowerCI.csv", sep="/"), col.names=TRUE, row.names=TRUE)
                      write.table(b, file=paste(folder.pathwaysurv, "iteration_number.csv", sep="/"), col.names=FALSE, row.names=FALSE)
                                                            
                                
                      if (b == bootstrap.reps) {break}              
      }
                
                
                   
                                  
                                  
                                  
                                  
                                  
             
        

                  