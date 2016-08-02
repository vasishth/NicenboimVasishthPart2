## ----setup,include=FALSE,cache=FALSE-------------------------------------
library(knitr)
library(coda)
library(Cairo)
# set global chunk options, put figures into folder
options(replace.assign=TRUE,show.signif.stars=FALSE)
opts_chunk$set(cache=TRUE,fig.path='figures2/figure-', fig.align='center', fig.show='hold', cache.path='cache/graphics-')
options(replace.assign=TRUE,width=75)
opts_chunk$set(dev='postscript')
#opts_chunk$set(dev='pdf')


options(digits = 2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(9991)



# knit_hooks$set(source = function(x, options) {
#     paste("\\singlespacing\n\\begin{listing}\n\\begin{Verbatim}[numbers=left,frame=single,fontfamily=courier,fontsize=\\footnotesize]\n
# ", x, 
#         "\n\\end{Verbatim}\n\\end{listing}\n\\doublespacing\n", sep = "")
# })

knit_hooks$set(source = function(x, options) {
    paste("\\begin{Verbatim}[numbers=left,fontfamily=courier,fontsize=\\footnotesize, firstnumber=last]\n", x, 
        "\n \\end{Verbatim}\n", sep = "")
}, output = function(x, options) {
    paste("\\begin{Verbatim}[fontfamily=courier,fontsize=\\footnotesize, firstnumber=last]\n", x, 
        " \\end{Verbatim}\n", sep = "")
}    )

#to reset the counter for the code, works?
# knit_hooks$set(reset = function(before, options, envir){
# if(before){
#     return("\\setcounter{lstnumber}{1}")
# }
# })

# , output = function(x, options) {
#     paste("\\begin{lstlisting}[basicstyle={\\ttfamily}]\n", x, 
#         "\\end{lstlisting}\n", sep = "")
# }
# , warning = hook_lst_bf, message = hook_lst_bf, error = hook_lst_bf)

# options$capT
 # "\\caption{", options$capT, "}\n",
# \singlespacing
# \begin{listing}
# \begin{Verbatim}[numbers=left,frame=single,fontfamily=courier,fontsize=\footnotesize]
# \end{Verbatim}
# \caption{Code for the fixed effects Model 1.}\label{fig:Model1code}
# \end{listing}
# \doublespacing




## ----fig01posteriorbinomial,echo=FALSE,fig.height=4,fig.width=5----------
library(ggplot2)
library(Cairo)
BT <- function(A,B,Z,N){
    l <- 1000
    prior <- data.frame(P= seq(0,1,length=l),density= dbeta(seq(0,1,length=l), A, B),Distribution="prior")  

    likelihood <- data.frame(P= seq(0,1,length=l),density= dbeta(seq(0,1,length=l), Z, N-Z),Distribution="likelihood")  

    posterior <- data.frame(P= seq(0,1,length=l),density= dbeta(seq(0,1,length=l), Z+A, N -Z + B),Distribution="posterior") 
    rbind(prior,likelihood,posterior)
}

BT1 <- BT(1,1,80,100)
BT1$title <- "Flat prior; N=100"
BT2 <- BT(1,1,8,10)
BT2$title <- "Flat prior; N=10"
BT3 <- BT(5,8,80,100)
BT3$title <- "Weakly informative prior; N=100"
BT4 <- BT(5,8,8,10)
BT4$title <- "Weakly informative prior; N=10"
BT <- rbind(BT1,BT2,BT3,BT4)
BT$Distribution <-  factor(BT$Distribution, levels=c("prior","likelihood","posterior"))

# plot <- ggplot(data=BT,aes(x=P,y=density,linetype=Distribution)) +geom_line()+
#   geom_ribbon(data=BT,aes(ymax=density,fill=Distribution),ymin=0,
#               colour=NA,alpha=1)+ #0.5 
#                theme_bw()+xlab("Estimate")+ylab("Density")+xlim(c(0,1))+  theme(legend.position="bottom")+facet_wrap(~title,nrow=2) #+scale_fill_manual(values=c("gray","white","black"))

plot <- ggplot(data=BT,aes(x=P,y=density,linetype=Distribution)) +geom_line()+
               theme_bw()+xlab("Estimate")+ylab("Density")+xlim(c(0,1))+  theme(legend.position="bottom")+facet_wrap(~title,nrow=2)+
               scale_linetype_manual(values=c("longdash","dotted","solid"))

                #+scale_fill_manual(values=c("gray","white","black"))


suppressWarnings(print(plot))

## ----model0,include=TRUE, message=F, warning=F,hide=T,background="white"----
library(rstanarm)
dgw <- read.table("gibsonwu2012data.txt")
dgw_hn <- subset(dgw, subset = region == "headnoun")
dgw_hn$cond <- ifelse(dgw_hn$type == "obj-ext", 1, -1)
m1 <- stan_lmer(formula = log(rt) ~ cond + (cond | subj) + (cond | item), 
                prior_intercept = normal(0, 10), 
                prior = normal(0, 1),            
                prior_covariance = decov(regularization = 2), 
                data = dgw_hn, 
                chains = 4, 
                iter = 2000, 
                cores = 4)
#summary(m1) # Very long summary with all the parameters in the model

## ----include=FALSE-------------------------------------------------------
save(m1,file="gibsonwu_m1.Rda")

## ----include=TRUE,background="white"-------------------------------------
samples_m1 <- as.data.frame(m1) # It saves all the samples from the model.
posterior_condition <- samples_m1$cond
options(digits = 4) 
mean(posterior_condition) 
median(posterior_condition)  

## ----include=FALSE-------------------------------------------------------
options(digits = 2)
posterior_Intercept <- samples_m1$'(Intercept)'
Pbelow0 <- mean(posterior_condition < 0 )
#difference between conditions in ms
diffms <- 20
expeffect <- round(-log(diffms/2 + exp(mean(posterior_Intercept))) +mean(posterior_Intercept),3)
Pbelow_expeffect <- mean(posterior_condition<expeffect)
condition_dens <- with(density(posterior_condition,adjust=2),data.frame(estimate=x,density=y,Distribution="posterior"))
prior_dens <- with(density(rnorm(4000,0,1),adjust=2),data.frame(estimate=x,density=y,Distribution="prior"))
both_dens <- rbind(condition_dens,prior_dens)

posterior_mean <- mean(posterior_condition) 



## ----include=TRUE,background="white"-------------------------------------
mean(posterior_condition < 0) # Posterior probability that lies below zero. 

## ----fig02ORSRdiff,echo=FALSE,fig.height=4,fig.width=5-------------------
library(ggplot2)
both_dens$Distribution <- factor(both_dens$Distribution,levels=c("prior","posterior"))
plot <- ggplot(data=both_dens,aes(x=estimate,y=density,linetype=Distribution)) +geom_line()+
 # geom_ribbon(data=both_dens,aes(ymax=density,fill=Distribution),ymin=0,
              #colour=NA,alpha=0.5)+ 
               theme_bw()+xlab("Estimate")+ylab("Density")+xlim(c(-2,2))+
                              scale_linetype_manual(values=c("longdash","solid"))

               #+scale_fill_manual(values=c("#d95f02","#1b9e77","#7570b3"))
suppressWarnings(print(plot))


## ----fig03posteriorbelows,echo=FALSE,fig.height=2.5,fig.width=5, fig.show='hold'----
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

plot1<-ggplot(data=condition_dens,aes(x=estimate,y=density)) +geom_line()+
  geom_ribbon(data=subset(condition_dens, estimate<0),aes(ymax=density),ymin=0,
              fill="gray",colour=NA,alpha=1)+ theme_bw()+xlab("Estimate")+ylab("Posterior Density")+
   geom_vline(xintercept = 0,linetype="dashed")+scale_fill_manual(values=c("#1b9e77","#d95f02","#7570b3")) +ggtitle("(a)")

plot2<-ggplot(data=condition_dens,aes(x=estimate,y=density)) + geom_line()+
  geom_ribbon(data=subset(condition_dens, estimate < expeffect),aes(ymax=density),ymin=0,
              fill="gray",colour=NA,alpha=1)+ theme_bw()+xlab("Estimate")+ylab("Posterior Density")+
   geom_vline(xintercept = expeffect,linetype="dashed") +scale_fill_manual(values=c("#1b9e77","#d95f02","#7570b3"))  +ggtitle("(b)")

multiplot(plot1,plot2,cols=2)

## ----include=TRUE,background="white"-------------------------------------
options(digits = 4) 
posterior_interval(m1, par = "cond", prob = 0.95) # 95% Percentile Interval
library(SPIn) # For calculating the HPDI
bootSPIn(posterior_condition)$spin # 95% HPDI  

## ----include=FALSE-------------------------------------------------------
options(digits = 2)

## ----fig04CrI,echo=FALSE,fig.height=4,fig.width=5------------------------
CrI <- as.numeric(posterior_interval(m1, par = "cond", prob = 0.95)) # 95% Percentile Interval

ggplot(data=condition_dens,aes(x=estimate,y=density)) +geom_line()+
  geom_ribbon(data=subset(condition_dens, estimate>CrI[1] & estimate<CrI[2]),aes(ymax=density),ymin=0,
              fill="gray",colour=NA,alpha=1)+ theme_bw()+xlab("Estimate")+ylab("Posterior Density")+
   geom_vline(xintercept = CrI,linetype="dashed")
   #+scale_fill_manual(values=c("#1b9e77","#d95f02","#7570b3"))

## ----include=FALSE-------------------------------------------------------

#for the weakly informative priors
realistic_effms <- 200 
realistic_eff <- round(-log(realistic_effms/2 + 
                exp(mean(log(dgw_hn$rt)))) +mean(log(dgw_hn$rt)),3)


realistic_eff2 <- round(-log(realistic_effms/4 + 
                exp(mean(log(dgw_hn$rt)))) +mean(log(dgw_hn$rt)),3)


#for the informative priors
expected_effms <- 40 

expected_eff <- round(-log(expected_effms/2 +
                exp(mean(log(dgw_hn$rt)))) +mean(log(dgw_hn$rt)),3)



priors <- list(very_weakly_informative = list(beta_mu = 0, beta_sd = 1),
                weakly_informative = list(beta_mu = 0, beta_sd = abs(realistic_eff)), 
                      informative = list(beta_mu = 0, beta_sd =abs(realistic_eff2)), 
                  off_informative = list(beta_mu = expected_eff*4, beta_sd =.02),
             opposite_informative = list(beta_mu = -expected_eff, beta_sd =.02),
             too_informative = list(beta_mu = expected_eff, beta_sd =.02)
             )


## ----include=FALSE,eval=TRUE---------------------------------------------

 
    diff_models <- plyr::llply(priors, function(pr){
       
        beta_mu <<- pr$beta_mu #Ugly way to make it work in knitr
        beta_sd <<- pr$beta_sd 

        model <- stan_lmer(formula = log(rt) ~ cond + (cond | subj) + (cond | item), 
                prior_intercept = normal(0, 10), # Prior for the intercept,
                                                # or grand mean since the
                                                # contrast is orthogonal.
                prior = normal(beta_mu, beta_sd), # Prior for the estimates, in this
                                      # case only condition.           
                prior_covariance = decov(regularization = 2), # Prior for
                                                              # covariance
                                                              # matrix.
                data = dgw_hn, 
                chains = 4, # Default number of chains.
                iter = 2000, # The number of iterations can be reduced to  
                             # speed up the computation to test the model, 
                             # but (at least) 2000 are recommended.
                cores = 4   # This argument is optional and on a multicore
                            # system the user may want to set it to speed   
                            # up the computation.
                )


        pr$model <- model
        return(pr)
    })

save(diff_models,file="gibsonwu_diff_models.Rda")

## ----include=FALSE-------------------------------------------------------

estimates <- plyr::ldply(diff_models, function(s){
  samples <- as.data.frame(s$model)
  mean <-round(mean(samples$cond),2)
  Cr <-round(quantile(samples$cond,c(0.025,0.975)),2)
  Pbelow0 <-round( mean(samples$cond < 0),2)
  prior <- paste("Normal(",round(s$beta_mu,2),", ",round(s$beta_sd,2),")",sep="")
  data.frame(prior = prior ,Cr_l = Cr[1],Cr_h= Cr[2],Pbelow0 = Pbelow0 ,mean = mean)
})

table <- paste(paste(paste("$",estimates$prior,"$"),
                      paste("$",estimates$Cr_l,"$"),
                      paste("$",estimates$Cr_h,"$"),
                      paste("$",estimates$Pbelow0,"$"),
                      paste("$",estimates$mean,"$"),
                      sep=" & "),collapse=" \\\\ ")


## ----fig05diffpriors,echo=FALSE,fig.height=4,fig.width=5-----------------

densities <- plyr::ldply(diff_models, function(s){
  model <- paste("Normal(",s$beta_mu,", ",s$beta_sd,")",sep="")
  samples <- as.data.frame(s$model)
  posterior_beta <- samples$cond
  condition_dens<- with(density(posterior_beta,adjust=2),data.frame(beta=x,density=y,Distribution="posterior",model=model))
  prior_dens <- with(density(rnorm(4000,s$beta_mu,s$beta_sd),adjust=2),data.frame(beta=x,density=y,Distribution="prior",model=model))
  rbind(condition_dens,prior_dens)
})
densities$Distribution <- factor(densities$Distribution,levels=c("prior","posterior"))
plot <- ggplot(data=densities,aes(x=beta,y=density,linetype=Distribution)) +
facet_wrap(~model,nrow=2) +geom_line() +
# geom_ribbon(data=densities,aes(ymax=density,fill=Distribution),ymin=0,
# colour=NA,alpha=1)+
 xlim(c(-0.25,0.25))+
theme_bw()+  theme(legend.position="bottom") +xlab("Estimate")+ylab("Density") + geom_vline(xintercept =
0,linetype="dashed")  + scale_linetype_manual(values=c("longdash","solid"))

#+scale_fill_manual(values=c("#d95f02","#1b9e77","#7570b3"))

suppressWarnings(print(plot))

## ----include=FALSE-------------------------------------------------------
options(digits = 4)

## ----include=TRUE,background="white"-------------------------------------
library(polspline)
fit_posterior <- logspline(posterior_condition)
posterior <- dlogspline(0, fit_posterior) # Height of the posterior at 0
prior     <- dnorm(0, 0, 1) # Height of the prior at 0
(BF01 <- posterior/prior) #BF01 shows clear support for H0 

## ----include=FALSE-------------------------------------------------------
options(digits = 2)

## ----include=FALSE-------------------------------------------------------
library(polspline)
BF01 <- plyr::ldply(diff_models, function(s){
  samples <- as.data.frame(s$model)
  beta_posterior <- samples$cond
  fit_posterior <- logspline(beta_posterior)
  posterior <- dlogspline(0, fit_posterior) # this gives the pdf at point delta = 0
  prior     <- dnorm(0,s$beta_mu,s$beta_sd)
  prior_d <- paste("Normal(",s$beta_mu,", ",s$beta_sd,")",sep="")

  data.frame( H0= ifelse(posterior>prior,round(posterior/prior,2),1) ,H1=ifelse(posterior<prior,round(prior/posterior,2),1),prior=prior_d)
})

BF01


BF01_1<- BF01[c(2,3),]
BF01_2<- BF01[c(4,5,6),]
table1 <- paste(paste(BF01_1$H0,BF01_1$H1,BF01_1$prior,sep=" & "),collapse=" \\\\ ")
table2 <- paste(paste(BF01_2$H0,BF01_2$H1,BF01_2$prior,sep=" & "),collapse=" \\\\ ")

#BF_H0/BF_{H_opp} / (BF_H0/BF_{H_inf})  = BF_{H_inf}/BF_{H_opp}
BFHinfHopp  <- BF01[BF01$'.id'=="opposite_informative",]$H0/BF01[BF01$'.id'=="too_informative",]$H0

