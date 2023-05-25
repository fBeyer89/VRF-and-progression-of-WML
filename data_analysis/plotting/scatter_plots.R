library(ggeffects)
library(ggplot2)
library(cowplot)
library(dplyr)
comp_plot <- function(model, predictor){
    #function to plot effects for one model and predictor
    wd="/data/pt_life_whm/Results/VRF_cSVD/LME/results/Conf/"

    
    #load imputed dataset
    miceadds::load.Rdata(objname = "imp", "/data/pt_life_whm/Results/Tables/imputed_data_08.5.23/imputed_data_08.5.23.Rdata")
    
    #Generate values to plot for individual predictors
    #load lme4 results
    print(paste0(wd,"workspace_",model,"_freq_imp_plot.RData"))
    load(paste0(wd,"workspace_",model,"_freq_imp_plot.RData"))
    
    #Laurenz way: does not work because data is not saved in the analyses 
    #effects:effect(term="asinh_wml_base", mod=plot$analyses[[1]]) 
    
    #ggpredict returns predicted values and confidence intervals
    #range = seq(min(imp$data$asinh_wml_change),max(imp$data$asinh_wml_change),length.out=10)
    #ggpredict(plot$analyses[[1]], "asinh_wml_change [range]", interval="prediction")
    
    #Use this tutorial to provide pooled estimation for all imputations:
    #https://solomonkurz.netlify.app/blog/2021-10-21-if-you-fit-a-model-with-multiply-imputed-data-you-can-still-plot-the-line/

    #Prepare imputed data for plotting effects
    #Effects to be plotted from model M1   
    fitted_lines <-
      tibble(.imp = 1:imp$m) %>% 
      mutate(p = map(.imp, ~ ggpredict(plot$analyses[[.]], terms=paste0(predictor, " [n=10]")) %>% 
                       data.frame())
      )
    
    fitted_lines <- fitted_lines %>% 
      unnest(p) 

    fitted_lines=fitted_lines %>% 
      group_by(x) %>% #group by predictor values
      summarise(fit_bar = mean(predicted),
                v_w     = mean(std.error^2),
                v_b     = sum((predicted - fit_bar)^2) / (imp$m - 1),
                v_p     = v_w + v_b * (1 + (1 / imp$m)),
                se_p    = sqrt(v_p))%>% 
      # use the _p suffix to indicate these are pooled
      mutate(lwr_p = fit_bar - se_p * 1.96,
             upr_p = fit_bar + se_p * 1.96) 
    
    #Significances based on frequentist results
    #load lmerMod results
    load(paste0(wd,"workspace_",model,"_freq_imp_res.RData"))
    est=summary(mice::pool(res))
    asterisiks="not significant"
    sign <- c(est[est$term==paste0(predictor),"p.value"])
    if(sign < 0.0001){
      asterisiks <- "****"
    } else if(sign < 0.001){
      asterisiks <- "***"
    } else if(sign < 0.01){
      asterisiks <- "**"
    } else if(sign < 0.05){
      asterisiks <- "*"
    } 
    
    #Bayes factors
    bf <- read.csv(paste0(wd, model,'_model_res_bayes.csv'))
    if (model=="M1_VRF"){
      tmp <- bf %>% mutate(full = mean_one_sided_bf/(mean_one_sided_bf + 1)*100,
                    null = 1/(mean_one_sided_bf+1)*100)
      tmp$pred=c("age_change:WHR_base", "age_change:DBP_base", "WHR_change", "DBP_change")
      df1 <- tmp %>% pivot_longer(cols=full:null)
      colnames(df1)[7]="model"
    } else {
    df1 <- data.frame(model = c("full", "null"), 
                      value= c(bf[bf$pred==paste0(predictor),"mean_one_sided_bf"]/(bf[bf$pred==paste0(predictor),"mean_one_sided_bf"]+1)*100,
                               1/(bf[bf$pred==paste0(predictor),"mean_one_sided_bf"]+1)*100))
    }
    
  
    #Actual PLOTTING
    if (predictor=="asinh_wml_change"){
      x_axis_label="Change of WML volume (asinh mm³)"
      col2pred="#fdae61" 
    } else if (predictor=="asinh_wml_base"){
      x_axis_label="Baseline of WML volume (asinh mm³)"
      col2pred="#fdae61" 
    }else if (predictor=="DBP_base"){
      x_axis_label="Baseline DBP (mmHg)"
      col2pred="#fdae61" 
    } else if (predictor=="DBP_change"){
      x_axis_label="Change of DBP (mmHg)"
      col2pred="#fdae61" 
    } else if (predictor=="WHR_base"){
      x_axis_label="Baseline WHR"
      col2pred="#fdae61" 
    } else if (predictor=="age_change:DBP_base"){
        x_axis_label="Baseline DBP (mmHg)"
        col2pred="#fdae61"
    } else if (predictor=="age_change:WHR_base"){
      x_axis_label="Baseline WHR"
      col2pred="#fdae61"
    } else if (predictor=="WHR_change"){
      x_axis_label="Change of WHR"
      col2pred="#fdae61" 
    }
    
    if (model=="M2_exfunct"){
      outcome="exfunct"
      y_axis_label="executive function (Z-scored"
    } else if(model == "M3_globalcog"){
      outcome="globalcog"
      y_axis_label="global cognitive function (Z-scored)"
    } else if(model == "M1_VRF"){
      outcome="asinh_wml"
      y_axis_label="WML volume (asinh-transformed)"
    }
    
    xmin=min(fitted_lines$x)
    xmax=max(fitted_lines$x)
    ymax=max(max(imp$data[,paste0(outcome)],na.rm=T))
    
    scatter = 
      fitted_lines %>% 
      ggplot(aes(x = x)) +
      geom_ribbon(aes(ymin = lwr_p, ymax = upr_p),
                  alpha = 1/2) +
      geom_line(aes(y = fit_bar), 
                size = 1/2, col=col2pred) +
      # add the observed data for good measure
      geom_point(data = imp$data,
                 aes_string(x=predictor, y = outcome), col = "gray80")+
      #add astericks depending on significance
      geom_text(data = data.frame(x = xmin, y = ymax), aes(x = x, y = y), label = asterisiks, cex = 8, hjust = 0)  +
      theme_bw() +
      #theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
      labs(x=x_axis_label, y=y_axis_label)+
      xlim(xmin, xmax) 
    
      p1 <- ggdraw(scatter)

      if (model=="M1_VRF" & (predictor != "DBP_base" & predictor != "WHR_base")){
      pie1 <- ggplot(df1[df1$pred==predictor,], aes(x="", y=value, fill=model)) +
        theme_void() + theme(legend.position="none") +
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0) +
        scale_fill_manual(values= c("#39e75f", "#000000"))
      
      p1 = p1 + draw_plot(pie1,x=0.8, y=0.8, width=.2, height=.2)
      return(p1)
      } else if (model == "M2_exfunct" | model == "M3_globalcog"){
      pie1 <- ggplot(df1, aes(x="", y=value, fill=model)) +
          theme_void() + theme(legend.position="none") +
          geom_bar(stat="identity", width=1) +
          coord_polar("y", start=0) +
          scale_fill_manual(values= c("#39e75f", "#000000"))
        
        p1 = p1 + draw_plot(pie1,x=0.8, y=0.8, width=.2, height=.2)
      } else
      {print("no pieplot created")
       return(p1)
        }
      
       
      
}
 
