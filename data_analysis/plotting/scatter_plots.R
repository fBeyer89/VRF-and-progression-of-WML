library(ggeffects)
library(ggplot2)
library(cowplot)
library(dplyr)
comp_plot <- function(model, predictor, dir, impdata){
    #function to plot effects for one model and predictor
    wd=dir

    
    #load imputed dataset
    miceadds::load.Rdata(objname = "imp", paste0(impdata))

    
    
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
    
    #replace missing by mean of imputed values for plotting
    comp_imp=mice::complete(imp, "long")
    comp_sum=comp_imp %>% select(subj,time,.imp, DBP_base, WHR_base) %>% 
      mutate(as.factor(.imp)) %>%
      group_by(subj,time) %>% summarize_at(c("DBP_base", "WHR_base"),
                                           mean)
    imp_data_s=imp$data %>% select(-DBP_base,-WHR_base) 
    imp_data=merge(imp_data_s, comp_sum, by=(c("subj", "time")), all.x=T)  
    DBP_cutpoints = quantile(imp_data$DBP_base, probs = c(0,.33,.66,1),  na.rm = TRUE)
    WHR_cutpoints = quantile(imp_data$WHR_base, probs = c(0,.33,.66,1),  na.rm = TRUE)
    WMH_cutpoints = quantile(imp_data$asinh_wml_base, probs = c(0,.33,.66,1),  na.rm = TRUE)
    
    imp_data= imp_data %>%    
      mutate(group_DBP = cut(DBP_base, breaks = DBP_cutpoints))%>%    
      mutate(group_WHR = cut(WHR_base, breaks = WHR_cutpoints))%>%    
      mutate(group_WMH = cut(asinh_wml_base, breaks = WMH_cutpoints))
    
    imp_data[is.na(imp_data$group_DBP),"group_DBP"]="(54,72]"
    imp_data[is.na(imp_data$group_WHR),"group_WHR"]="(0.728,0.902]"
    imp_data[is.na(imp_data$group_WMH),"group_WMH"]="(0.0296,0.507]"
    levels(imp_data$group_WHR)=c("low", "middle", "high")
    levels(imp_data$group_DBP)=c("low", "middle", "high")
    levels(imp_data$group_WMH)=c("low", "middle", "high")
    
    ## Plotting dependent on predictor
    if (predictor=="age_change:DBP_base"){
      print("DBP time interaction")
     
    #select group variable
    imp_data$group=imp_data$group_DBP
    
    terms=c("age_change",
            paste0("DBP_base [",
                   ((DBP_cutpoints[1]+DBP_cutpoints[2])/2),",",
                   (DBP_cutpoints[2]+DBP_cutpoints[3])/2,",", 
                   (DBP_cutpoints[3]+DBP_cutpoints[4])/2,
                   "]"))
    print(terms)
    
    fitted_lines <-
      tibble(.imp = 1:imp$m) %>% 
      mutate(p = map(.imp, ~ ggpredict(plot$analyses[[.]], 
                                       terms=terms) %>% 
                       data.frame()))
    fitted_lines <- fitted_lines %>% 
      unnest(p) 
    
     levels(imp_data$group)=c("low", "middle", "high")
     levels(fitted_lines$group)=c("low", "middle", "high")
     fitted_lines=fitted_lines %>% 
       group_by(x,group) %>% #group by predictor values
       summarise(fit_bar = mean(predicted),
                 v_w     = mean(std.error^2),
                 v_b     = sum((predicted - fit_bar)^2) / (imp$m - 1),
                 v_p     = v_w + v_b * (1 + (1 / imp$m)),
                 se_p    = sqrt(v_p))%>% 
       # use the _p suffix to indicate these are pooled
       mutate(lwr_p = fit_bar - se_p * 1.96,
              upr_p = fit_bar + se_p * 1.96) 
      
    } 
    else if(predictor=="asinh_wml_base"){
      print("asinh_wml_base time interaction")
      
      #select group variable
      imp_data$group=imp_data$group_WMH
      
      
      terms=c("age_change",
              paste0("asinh_wml_base [",
                     ((WMH_cutpoints[1]+WMH_cutpoints[2])/2),",",
                     (WMH_cutpoints[2]+WMH_cutpoints[3])/2,",", 
                     (WMH_cutpoints[3]+WMH_cutpoints[4])/2,
                     "]"))
      print(terms)
      
      fitted_lines <-
        tibble(.imp = 1:imp$m) %>% 
        mutate(p = map(.imp, ~ ggpredict(plot$analyses[[.]], 
                                         terms=terms) %>% 
                         data.frame()))
      fitted_lines <- fitted_lines %>% 
        unnest(p) 
      
      levels(imp_data$group)=c("low", "middle", "high")
      levels(fitted_lines$group)=c("low", "middle", "high")
      fitted_lines=fitted_lines %>% 
        group_by(x,group) %>% #group by predictor values
        summarise(fit_bar = mean(predicted),
                  v_w     = mean(std.error^2),
                  v_b     = sum((predicted - fit_bar)^2) / (imp$m - 1),
                  v_p     = v_w + v_b * (1 + (1 / imp$m)),
                  se_p    = sqrt(v_p))%>% 
        # use the _p suffix to indicate these are pooled
        mutate(lwr_p = fit_bar - se_p * 1.96,
               upr_p = fit_bar + se_p * 1.96) 
      
    
    } 
    else if(predictor=="age_change:WHR_base"){
      print("WHR time interaction")
      
      #select group variable
      imp_data$group=imp_data$group_WHR
      
      
      terms=c("age_change",
              paste0("WHR_base [",
                     ((WHR_cutpoints[1]+WHR_cutpoints[2])/2),",",
                     (WHR_cutpoints[2]+WHR_cutpoints[3])/2,",", 
                     (WHR_cutpoints[3]+WHR_cutpoints[4])/2,
                     "]"))
      print(terms)
      
      fitted_lines <-
        tibble(.imp = 1:imp$m) %>% 
        mutate(p = map(.imp, ~ ggpredict(plot$analyses[[.]], 
                                         terms=terms) %>% 
                         data.frame()))
      fitted_lines <- fitted_lines %>% 
        unnest(p) 
      
      levels(imp_data$group)=c("low", "middle", "high")
      levels(fitted_lines$group)=c("low", "middle", "high")
      fitted_lines=fitted_lines %>% 
        group_by(x,group) %>% #group by predictor values
        summarise(fit_bar = mean(predicted),
                  v_w     = mean(std.error^2),
                  v_b     = sum((predicted - fit_bar)^2) / (imp$m - 1),
                  v_p     = v_w + v_b * (1 + (1 / imp$m)),
                  se_p    = sqrt(v_p))%>% 
        # use the _p suffix to indicate these are pooled
        mutate(lwr_p = fit_bar - se_p * 1.96,
               upr_p = fit_bar + se_p * 1.96) 
      
      
    }
    else {
      
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
    }

    #Significances based on frequentist results
    #load lmerMod results
    load(paste0(wd,"workspace_",model,"_freq_imp_res.RData"))
    est=summary(mice::pool(res))
    asterisiks=""
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
    #Pizza plot explanation: 17 times more likely -> 17 x 1 of area should be colored
    #https://jasp-stats.org/wp-content/uploads/2020/03/The_JASP_Data_Library_Version_3.pdf
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
      x_axis_label="Change of WMH (asinh mm³)"
      col2pred="#fdae61" 
      if (model=="M2_exfunct"){
        outcome="exfunct"
        y_axis_label=""
      } else if(model == "M3_globalcog"){
        outcome="globalcog"
        y_axis_label=""}
    } else if (predictor=="asinh_wml_base"){
      x_axis_label="Time since baseline (y)"
      lablabel="Baseline WMH \n volume"
      if (model=="M2_exfunct"){
        outcome="exfunct"
        y_axis_label="executive function\n (Z-scored)"
      } else if(model == "M3_globalcog"){
        outcome="globalcog"
        y_axis_label="global cognitive function\n (Z-scored)"}
    }else if (predictor=="DBP_base"){
      x_axis_label="Baseline DBP (mmHg)"
      col2pred="#fdae61" 
      y_axis_label="WMH volume\n (asinh-transformed)"
    } else if (predictor=="DBP_change"){
      x_axis_label="Change of DBP (mmHg)"
      col2pred="#fdae61" 
      xpos=0.8
      outcome="asinh_wml"
      y_axis_label=""
    } else if (predictor=="WHR_base"){
      x_axis_label="Baseline WHR"
      col2pred="#fdae61" 
      outcome="asinh_wml"
      y_axis_label="WMH volume\n (asinh-transformed)"
    } else if (predictor=="age_change:DBP_base"){
        x_axis_label="Time since baseline (y)"
        lablabel="Baseline DBP"
        xpos=0.8
        outcome="asinh_wml"
        y_axis_label="WMH volume\n (asinh)"
    } else if (predictor=="age_change:WHR_base"){
      x_axis_label="Time since baseline (y)"
      lablabel="Baseline WHR"
      xpos=0.8
      outcome="asinh_wml"
      y_axis_label="WMH volume\n (asinh)"
    } else if (predictor=="WHR_change"){
      x_axis_label="Change of WHR"
      col2pred="#fdae61" 
      xpos=0.8
      outcome="asinh_wml"
      y_axis_label=""
      
    }
    
    

    
    xmin=min(fitted_lines$x)
    xmax=max(fitted_lines$x)
    ymax=max(max(imp_data[,paste0(outcome)],na.rm=T))
    
    if (predictor=="age_change:WHR_base"|predictor=="age_change:DBP_base") {
      scatter = 
        ggplot() +
        geom_blank() +
        geom_ribbon(aes(x = x, ymin = lwr_p, ymax = upr_p,group=group,fill=group),
                    alpha = 0.2, data=fitted_lines) +
        geom_line(aes(x=x, y = fit_bar, group=group, color=group), 
                  size = 1/2, data=fitted_lines) +
        # add the observed data for good measure
        geom_point(data = imp_data,
                   aes_string(x="age_change", y = outcome, group="group", color="group"), 
                   alpha=0.5)+
        #add astericks depending on significance
        geom_text(data = data.frame(x = xmin, y = ymax), aes(x = x, y = y), 
                  label = asterisiks, cex = 4, hjust = 0)  +
        theme_bw(base_size = 12) +
        theme(legend.position="left") +
        labs(x=x_axis_label, y=y_axis_label)+
        xlim(xmin, xmax) +
        scale_color_manual(name=lablabel, labels = c("low", "middle", "high"), values = hue_pal()(3)) +
        scale_fill_manual(name=lablabel, labels = c("low", "middle", "high"), values = hue_pal()(3)) 
      
      p1 <- ggdraw(scatter)
    } 
    else if (predictor=="WHR_change"){
      scatter = 
      fitted_lines %>% 
      ggplot(aes(x = x)) +
      # add the observed data for good measure
      geom_point(data = imp_data,
                 aes_string(x=predictor, y = outcome, color="group_WHR"), alpha=0.5)+
        geom_ribbon(aes(ymin = lwr_p, ymax = upr_p),
                    alpha = 1/2) +
        geom_line(aes(y = fit_bar), 
                  size = 1/2, col=col2pred) +
      #add astericks depending on significance
      geom_text(data = data.frame(x = xmin, y = ymax), aes(x = x, y = y), label = asterisiks, cex = 4, hjust = 0)  +
      theme_bw(base_size = 12) +
        theme(legend.position="none") +
      labs(x=x_axis_label, y=y_axis_label)+
      xlim(xmin, xmax) + scale_colour_discrete(name = "Baseline WHR")
    
    p1 <- ggdraw(scatter)
    } 
    else if (predictor=="DBP_change"){
      scatter = 
      fitted_lines %>% 
      ggplot(aes(x = x)) +
      # add the observed data for good measure
      geom_point(data = imp_data,
                 aes_string(x=predictor, y = outcome, color="group_DBP"),
                 alpha=0.5)+
      geom_ribbon(aes(ymin = lwr_p, ymax = upr_p),
                  alpha = 1/2) +
      geom_line(aes(y = fit_bar), 
                size = 1/2, col=col2pred) +
      #add astericks depending on significance
      geom_text(data = data.frame(x = xmin, y = ymax), aes(x = x, y = y), 
                label = asterisiks, cex= 5, hjust = 0)  +
      theme_bw(base_size = 12) +
      theme(legend.position="none") +
      labs(x=x_axis_label, y=y_axis_label)+
      xlim(xmin, xmax)+ scale_colour_discrete(name = "Baseline DBP")
    
    p1 <- ggdraw(scatter)}
    else if (predictor=="asinh_wml_change"){
      scatter = 
        fitted_lines %>% 
        ggplot(aes(x = x)) +
        # add the observed data for good measure
        geom_point(data = imp_data,
                   aes_string(x=predictor, y = outcome, color="group_WMH"), alpha=0.5)+
        geom_ribbon(aes(ymin = lwr_p, ymax = upr_p),
                    alpha = 1/2) +
        geom_line(aes(y = fit_bar), 
                  size = 1/2, col=col2pred) +
        #add astericks depending on significance
        geom_text(data = data.frame(x = xmin, y = ymax), aes(x = x, y = y), 
                  label = asterisiks, cex= 5, hjust = 0)  +
        theme_bw(base_size = 12) +
        theme(legend.position="none") +
        labs(x=x_axis_label, y=y_axis_label)+
        xlim(xmin, xmax)
      
      p1 <- ggdraw(scatter)}
    else if (predictor=="asinh_wml_base"){
      scatter = 
        ggplot() +
        geom_blank() +
        geom_ribbon(aes(x = x, ymin = lwr_p, ymax = upr_p,group=group,fill=group),
                    alpha = 0.2, data=fitted_lines) +
        geom_line(aes(x=x, y = fit_bar, group=group, color=group), 
                  size = 1/2, data=fitted_lines) +
        # add the observed data for good measure
        geom_point(data = imp_data,
                   aes_string(x="age_change", y = outcome, group="group", color="group"), 
                   alpha=0.5)+
        #add astericks depending on significance
        geom_text(data = data.frame(x = xmin, y = ymax), aes(x = x, y = y), 
                  label = asterisiks, cex = 4, hjust = 0)  +
        theme_bw(base_size = 12) +
        theme(legend.position="bottom") +
        labs(x=x_axis_label, y=y_axis_label)+
        xlim(xmin, xmax) +
        scale_color_manual(name=lablabel, labels = c("low", "middle", "high"), values = hue_pal()(3)) +
        scale_fill_manual(name=lablabel, labels = c("low", "middle", "high"), values = hue_pal()(3)) 

        p1 <- scatter}

      if (model=="M1_VRF" & (predictor != "DBP_base" & predictor != "WHR_base")){
      pie1 <- ggplot(df1[df1$pred==predictor,], aes(x="", y=value, fill=model)) +
        theme_void() + theme(legend.position="none") +
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0) +
        scale_fill_manual(values= c("#39e75f", "#000000"))
      
      p1p = p1 + draw_plot(pie1,x=xpos, y=0.8, width=.2, height=.2)
      return(list(p1p, p1, pie))
      } else if (model == "M2_exfunct" | model == "M3_globalcog"){
      pie1 <- ggplot(df1, aes(x="", y=value, fill=model)) +
          theme_void() + theme(legend.position="none") +
          geom_bar(stat="identity", width=1) +
          coord_polar("y", start=0) +
          scale_fill_manual(values= c("#39e75f", "#000000"))
        
        p1p = ggdraw(p1) + draw_plot(pie1,x=0.8, y=0.8, width=.2, height=.2)
        return(list(p1p, p1, pie1))
      } else
      {print("no pieplot created")
       return(p1)
        }
      
       
      
}
 
