## Figure
# If we want to use this to make a figure, we can pass a sequence
# of values to x1. Remember, x2 must always be 1 row. Let's 
# visualize the RSS for forage vs. mean forage.
x1 <- data.frame(forage = seq(0, 800, length.out = 100), 
                 predator = 0, 
                 sl_ = 50, log_sl_ = log(50),
                 cos_ta_ = 0)

x2 <- data.frame(forage = mean(values(hab$forage)), 
                 predator = 0, 
                 sl_ = 50, log_sl_ = log(50),
                 cos_ta_ = 0)

logRSS <- log_rss(m1, x1, x2, ci = "se", ci_level = 0.95)

# We have a plot method for 'log_rss' objects to make a very basic figure.
plot(logRSS)

# But if we want more control, we can use ggplot with the 'df' data.frame.
ggplot(logRSS$df, aes(x = forage_x1, y = exp(log_rss), 
                      ymin = exp(lwr), ymax = exp(upr))) +
  geom_ribbon(color = "black", fill = "gray80", linetype = "dashed") +
  geom_line() +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  xlab("Forage at x1") +
  ylab("RSS vs. Mean Forage") +
  theme_bw()


#################################################################
Efxplot2 <- function(ModelList,
                    Sig = TRUE, StarLoc = NULL,
                    Alpha1 = 1, Alpha2 = 1,
                    PointOutline = T,
                    ModelNames = NULL,
                    VarNames = NULL, VarOrder = NULL,
                    Intercept = TRUE, PointSize = 1,
                    BarWidth = 0.1,
                    tips = 0.2){
  
  require(dplyr); require(ggplot2); require(INLA); require(MCMCglmm); require(MASS)
  
  Graphlist <- list()
  
  if(!class(ModelList) == "list"){
    
    ModelList <- list(ModelList)
    
  }
  
  if(is.null(ModelNames) & !is.null(names(ModelList))){
    
    ModelNames <- names(ModelList)
    
  }
  
  for(i in 1:length(ModelList)){
    
    model <- ModelList[[i]]
    
    if(class(model) == "inla"){
      
      Graph <- as.data.frame(summary(model)$fixed)
      colnames(Graph)[which(colnames(Graph)%in%c("0.025quant","0.975quant"))] <- c("Lower","Upper")
      colnames(Graph)[which(colnames(Graph)%in%c("0.05quant","0.95quant"))] <- c("Lower","Upper")
      colnames(Graph)[which(colnames(Graph)%in%c("mean"))] <- c("Estimate")
      
      rownames(Graph) %<>% str_replace_all(":", "_")
      
    }
    
    if(class(model) == "MCMCglmm"){
      Graph <- as.data.frame(summary(model)$solutions)
      colnames(Graph)[1:3] <- c("Estimate","Lower","Upper")
    }
    
    if(any(class(model) %>% str_detect("bam|gam"))){
      
      Graph  <-
        summary(model)[1:4] %>%
        map(~.x[1:length(summary(model)[[1]])]) %>%
        bind_cols() %>%
        rename(Estimate = p.coeff, P = p.pv) %>%
        # mutate(Lower = Estimate - se, Upper = Estimate + se) %>%
        as.data.frame
      
      Coefs <- coef(model)
      VC <- vcov(model)
      
      sim <- mvrnorm(1000, mu = Coefs, Sigma = VC)
      
      MCMCEffects <- sim %>% as.data.frame %>%
        map(~.x %>% as.mcmc %>% HPDinterval)
      
      rownames(Graph) <- summary(model)[[1]] %>% attr("names")
      
      Graph <- MCMCEffects[1:nrow(Graph)] %>%
        map(~tibble(Lower = .x[1],
                    Upper = .x[2])) %>%
        # bind_rows(.id = "Name") %>%
        bind_rows() %>%
        bind_cols(Graph, .)
      
    }
    
    Graph$Model <- i
    Graph$Factor <- rownames(Graph)
    
    Graphlist[[i]] <- Graph
    
  }
  
  Graph <- bind_rows(Graphlist)
  
  Graph$Sig <- with(Graph, ifelse(Lower*Upper>0, "*", ""))
  
  Graph$Model <- as.factor(Graph$Model)
  
  if(!is.null(ModelNames)){
    levels(Graph$Model) <- ModelNames
  }
  
  position <- ifelse(length(unique(Graph$Model))  ==  1, "none", "right")
  
  if(is.null(VarOrder)) VarOrder <- rev(unique(Graph$Factor))
  if(is.null(VarNames)){
    
    VarNames <- c(VarOrder)
    names(VarNames) <- c(VarOrder)
    
  }
  
  Graph$Factor %<>% str_replace_all("(Intercept)", "Intercept")
  
  Graph$Factor <- factor(Graph$Factor, levels = VarOrder)
  
  levels(Graph$Factor) %<>% str_replace_all(VarNames)
  
  Graph %<>% as.data.frame %>% filter(!is.na(Factor))
  
  if(!Intercept){
    
    VarNames <- VarNames[!str_detect(VarNames, "ntercept")]
    
    Graph <- Graph %>% filter(Factor %in% VarNames)
    
  }
  
  Graph$starloc <- NA
  
  min <- min(Graph$Lower,na.rm = T)
  max <- max(Graph$Upper,na.rm = T)
  
  if(Sig == TRUE){
    
    Graph$starloc <- max + (max - min)/10
    
  }
  
  if(!is.null(StarLoc)){
    
    Graph$starloc <- StarLoc
    
  }
  
  Graph$Alpha <- with(Graph, ifelse(Lower*Upper>0, Alpha1, Alpha2))
  
  Graph %>%
    mutate(SigAlpha = factor(as.numeric(Lower*Upper > 0),
                             levels = c(1, 0))) ->
    
    Graph
  
  if(PointOutline){
    
    PointOutlineAlpha <- Alpha1
    
  }else{
    
    PointOutlineAlpha <- 0
    
  }
  return(Graph)
}
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")  
  
  ggplot(Graph,
         aes(x = as.factor(Factor),
             y = Estimate,
             alpha = SigAlpha)) +
    geom_point(position = position_dodge(w = 0.5), size = 1) +
    geom_errorbar(position = position_dodge(w = 0.5),
                  aes(ymin = Lower, ymax = Upper), size = 0.3,
                  width = tips) +
    geom_hline(aes(yintercept = 0),lty = 2) + labs(x = NULL) + coord_flip() +
    theme(legend.position = position) +
    geom_text(aes(label = Sig, y = starloc),
              position = position_dodge(w = 0.5),
              show.legend = F) +
    scale_alpha_manual(values = c(1, 1)) +
    guides(alpha = "none") +
    geom_point(position = position_dodge(w = 0.5), size = PointSize*(4/3),
               alpha = PointOutlineAlpha) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
                  width = 0.1,
                  position = position_dodge(w = 0.5)) +
    geom_point(position = position_dodge(w = 0.5), size = PointSize) +
  scale_colour_manual(values=cbbPalette)


Efxplot2(r.inla.global)
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# To use for line and point colors, add
scale_colour_manual(values=cbPalette)


GraphGlobal <- GraphGlobal[-10,]
ggplot(GraphGlobal,
       aes(x = as.factor(Factor),
           y = Estimate,
           alpha = SigAlpha)) +
  geom_point(position = position_dodge(w = 0.5), size = 2) +
  geom_hline(aes(yintercept = 0),lty = 2, colour="#56B4E9") + labs(x = NULL) + coord_flip() +
  geom_text(aes(label = Sig, y = starloc),
            position = position_dodge(w = 0.5),colour="#009E73",
            show.legend = F, size=8) +
  scale_alpha_manual(values = c(1, 1)) +
  guides(alpha = "none") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, colour="#E69F00"),
                width = 0,show.legend=F,
                position = position_dodge(w = 0.5)) +
  geom_point(position = position_dodge(w = 0.5), size = 1) +
  scale_x_discrete(labels = c("Log Step Length x Rolling Avg",'Log Step Length','Cosine Turning Angle','Step Length','Distance to Trail','TRI','Wetland','Herbaceous','Shrub','Developed')) +
  ggtitle("Posterior Beta Estimates for Global Model") +
  theme(plot.title = element_text(hjust = 0.5),axis.text=element_text(size=12),panel.background = element_blank())



#Day and Night
DayNightList <- list(Day=r.inla.global.day, Night=r.inla.global.night)
DayNightGraph <- Efxplot2(DayNightList)
DayNightGraph <- DayNightGraph[c(-10,-21),]
position <- ifelse(length(unique(DayNightGraph$Model))  ==  1, "none", "right")

ggplot(DayNightGraph,
       aes(x = as.factor(Factor),
           y = Estimate,
           group = Model,
           colour = Model,
           alpha = SigAlpha)) +
  geom_point(position = position_dodge(w = 0.5), size = 2) +
  geom_errorbar(position = position_dodge(w = 0.5),
                aes(ymin = Lower, ymax = Upper), width = 0, show.legend=T) +
  geom_hline(aes(yintercept = 0),lty = 2) + labs(x = NULL) + coord_flip() +
  theme(legend.position = position) +
  geom_text(aes(label = Sig, y = starloc),
            position = position_dodge(w = 0.5),
            show.legend = F, size=8) +
  scale_alpha_manual(values = c(1, 1)) +
  guides(alpha = "none") +
  geom_point(position = position_dodge(w = 0.5), size = 2,
             alpha = 1) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                width = 0,
                position = position_dodge(w = 0.5)) +
  geom_point(position = position_dodge(w = 0.5), size = 1) +
  #scale_x_discrete(labels = c("Log Step Length x Rolling Avg",'Log Step Length','Cosine Turning Angle','Step Length','Distance to Trail','TRI','Wetland','Herbaceous','Shrub','Developed')) +
  ggtitle("Posterior Beta Estimates for Day and Night Models") +
  labs(color= "Legend") +
  scale_colour_manual(values= cbbPalette) +
  theme(plot.title = element_text(hjust = 0.5),axis.text=element_text(size=12),panel.background = element_blank())
###################################################################################
Day2 <- Day[,c(2,4,5,7,30,37,38)]
seq <- rep(g$round, times=9)
#sl_
x <- rep(g$`RndSteps5$x`, times=9)
newobs <- data.frame(sl_=0,log_sl_=0,cos_ta_=0,x=x,UnscaleDist=seq,Stratum=seq(424:792),case=NA,ANIMAL_ID1=rep(1:9, each=41))
newsobs <- newobs[order(newobs$ANIMAL_ID1,newobs$x),]

#log_sl_
range(Night$sl_)
range(Night$log_sl_)
newobs <- data.frame(sl_=, log_sl_=,cos_ta=mean(Night$cos_ta_),
                     x=mean(Night$x),RA=mean(Night$RA),TRI=mean(Night$TRI),
                     lc=c(41,22,52,71,90))
                     
                      
newobs$developed<-ifelse(RndSteps4$lc %in% c(21,22,23,24), 1, 0) #This is call dummy code 
newobs$forest<-ifelse(RndSteps4$lc %in% c(41,42), 1, 0)
newobs$shrub<-ifelse(RndSteps4$lc %in% c(52), 1, 0)
newobs$herb <- ifelse(RndSteps4$lc %in% c(71,81,82), 1, 0)
newobs$wetlands <- ifelse(RndSteps4$lc %in% c(90,95), 1, 0)