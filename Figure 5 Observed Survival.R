# "Predicting clinical trajectory after the Norwood procedure: focus on aortic and pulmonary re-interventions"
# Authors: Haapanen H, Dawes TJW, Brown K, Giardini A, Dedieu N, Shetty P, Tsang V, Kostolny M.


# Observed Clinical Trajectories: treatment, treatment plus, off-treatment and incomplete pathway patients
      path.type<- rep("", no.cases)
      paths.Allc.DorT<- paths.Allc
      for (i in 1:no.cases) {if (d$DorT[i] == 1) {paths.Allc.DorT[i]<- gsub("FU", "DorT", paths.Allc.DorT[i])}}


# Vectors which help define the three pathways

# Treatment ("Direct")
      direct.path.list<- list(c(" -> Norwood -> Glenn -> TCPC -> FU"),
                              c(" -> Norwood -> Glenn -> TCPC -> DorT"),
                              c(" -> Norwood -> Glenn -> FU"),
                              c(" -> Norwood -> FU"))
      
      path.direct<- sapply(paths.Allc.DorT, function(x) {any(sapply(direct.path.list, function(y) {identical(x,y)}))})


# Treatment+ ("Direct+")
      directplus.path.list<- list(c(" -> Norwood -> Glenn AND reop -> FU"),
                                  c(" -> Norwood -> Glenn AND reop -> TCPC -> FU"),
                                  c(" -> Norwood -> Glenn AND reop -> TCPC -> DorT"),
                                  c(" -> Norwood -> Glenn AND reop -> TCPC AND reop -> FU"),
                                  c(" -> Norwood -> Glenn AND reop -> TCPC AND reop -> DorT"),
                                  c(" -> Norwood -> Glenn -> TCPC AND reop -> FU"),
                                  c(" -> Norwood -> Glenn -> TCPC AND reop -> DorT"))
      
      path.directplus<- sapply(paths.Allc.DorT, function(x) {any(sapply(directplus.path.list, function(y) {identical(x,y)}))})

# Off-treatment: ("Indirect")
      indirect<- c("-> reop ->", "-> reop.1 ->", "-> reop.2 ->", "-> reop.3 ->", "-> reop.4 ->", "-> reop.5 ->")
# Which paths have a reintervention?
      path.reintervention<- sapply(paths.Allc.DorT, function(path) {any(sapply(indirect, function(x) {grepl(x, path)}))})

# Now need to pick out two groups:
#     1. those who had (a) a reintervention (b) didn't die
#     2. those who had (a) a reintervention (b) TCPC (c) then died

# What is the index of dying for all cases?
      path.dead.index<- sapply(paths.Allc.DorT, function(x) {gregexpr("DorT", x)[[1]]})
      path.dead.index[path.dead.index < 0]<- NA

#  What is the index of the TCPC for all cases?
      path.TCPC.index<- sapply(paths.Allc.DorT, function(x) {gregexpr("TCPC", x)[[1]]})
      path.TCPC.index[path.TCPC.index < 0]<- NA

# Group 1
      path.didnt.die<- sapply(paths.Allc.DorT, function(x) {!grepl("DorT", x)})
      path.indirect.group1<- which(path.reintervention & path.didnt.die)

# Group 2
      path.TCPC.before.dead<- path.TCPC.index < path.dead.index
      path.indirect.group2<- which(path.reintervention & path.TCPC.before.dead)
      paths.Allc.DorT[path.indirect.group2]

      path.indirect<- sort(c(path.indirect.group1, path.indirect.group2))

# Died before TCPC

# Who died before TCPC?
      path.did.die<- sapply(paths.Allc.DorT, function(x) {grepl("DorT", x)})
      path.did.not.have.TCPC<- sapply(paths.Allc.DorT, function(x) {!grepl("TCPC", x)})
      path.died.before.TCPC<- which(path.did.die & path.did.not.have.TCPC)


      path.type<- rep(NA, no.cases)
      path.type[path.direct]<- "Treatment"
      path.type[path.directplus]<- "Treatmentplus"
      path.type[path.indirect]<- "Offtreatment"
      path.type[path.died.before.TCPC]<- "DiedBeforeTCPC"
      table(path.type)





# Weight, N1Age, Ndays
# Categorical: N1ECMO, Interdigitation, SanovsBT, AVVRcat, AVVR, SystV_f, DGcat, SysV, era

      df<- NULL
      for (var in c("N1ECMO", "Interdigitation", "SanovsBT", "AVVRcat", "AVVR", "SystV_f", "DGcat", "SysV", "era"))
      {
        t<- table(path.type, vars[[var]])
        t<- t[-which(rownames(t) == "Incomplete_not_dead"),]
        t<- 100 * round(sweep(t, 1, apply(t,1,sum), FUN = "/"), 3)
        Freq = c(t)
        Label = rep(rownames(t), length(na.omit(unique(vars[[var]]))))
        Label2 = as.factor(rep(na.omit(unique(vars[[var]])), each = nrow(t)))
        
        df<- rbind(df, data.frame(Freq = Freq, Pathway = Label, Category = Label2, Group = var))
      }
      
      p<- ggplot(df, aes(x = Pathway, y = Freq)) + geom_bar(
        aes(fill = Category), stat = "identity", color = "white",
        position = "stack") +
        facet_wrap(~Group) + fill_palette("jco")
      p
      

