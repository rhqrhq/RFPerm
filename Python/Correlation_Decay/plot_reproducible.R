#Plotting scripts:

power_df <- read.csv("datasets/XGB_correlation_start_from_5th_tree.csv")
png('datasets/XGB_LongTermCorrelation.png', width = 900, height = 1050)
colors <- c(
           "Start from 50th Tree" = "green",
           "Start from 15th Tree" = "red",
           "Start from 20th Tree" = "blue",
           "Start from 25th Tree" = "orange",
           "Start from 30th Tree" = "cyan2",           
           "Start from 35th Tree" = "black",
           "Start from 40th Tree" = "brown",           
           "Start from 45th Tree" = "purple"        
           )
ggplot(power_df, aes(x = c(1:100))) +
  geom_line(aes(y = start50_mean, group = 1, color = "Start from 50th Tree"), linewidth = 1.25) + 
  geom_line(aes(y = start15_mean, group = 1, color = "Start from 15th Tree"), linewidth = 1.25) + 
  geom_line(aes(y = start20_mean, group = 1, color = "Start from 20th Tree"), linewidth = 1.25) + 
  geom_line(aes(y = start25_mean, group = 1, color = "Start from 25th Tree"), linewidth = 1.25) + 
  geom_line(aes(y = start30_mean, group = 1, color = "Start from 30th Tree"), linewidth = 1.25) +   
  geom_line(aes(y = start35_mean, group = 1, color = "Start from 35th Tree"), linewidth = 1.25) +   
  geom_line(aes(y = start40_mean, group = 1, color = "Start from 40th Tree"), linewidth = 1.25) +   
  geom_line(aes(y = start45_mean, group = 1, color = "Start from 45th Tree"), linewidth = 1.25) +   
  labs(x = "N Rounds Lag", y = "Correlation", color = "Legend") + 
ggtitle("Long-Term Dependence for XGBoost Model \n 150 Trees in Total, 10 Signal Features") +   
  geom_errorbar(aes(ymin = start50_mean - 0.2 * start50_sd, 
                    ymax = start50_mean + 0.2 * start50_sd,
                    color = "Start from 50th Tree"), width = 0.08,
                position = position_dodge(0.25), size = 0.25) + 
  geom_errorbar(aes(ymin = start15_mean - 0.5 * start15_sd, 
                    ymax = start15_mean + 0.5 * start15_sd,
                    color = "Start from 15th Tree"), width = 0.08,
                position = position_dodge(0.25), size = 0.25) + 
  geom_errorbar(aes(ymin = start20_mean - 0.5 * start20_sd, 
                    ymax = start20_mean + 0.5 * start20_sd,
                    color = "Start from 20th Tree"), width = 0.08,
                position = position_dodge(0.25), size = 0.25) +   
  geom_errorbar(aes(ymin = start25_mean - 0.45 * start25_sd, 
                    ymax = start25_mean + 0.45 * start25_sd,
                    color = "Start from 25th Tree"), width = 0.08,
                position = position_dodge(0.25), size = 0.25) + 
  geom_errorbar(aes(ymin = start30_mean - 0.4 * start30_sd, 
                    ymax = start30_mean + 0.4 * start30_sd,
                    color = "Start from 30th Tree"), width = 0.08,
                position = position_dodge(0.25), size = 0.25) +   
  geom_errorbar(aes(ymin = start35_mean - 0.35 * start35_sd, 
                    ymax = start35_mean + 0.35 * start35_sd,
                    color = "Start from 35th Tree"), width = 0.08,
                position = position_dodge(0.25), size = 0.25) +   
  geom_errorbar(aes(ymin = start40_mean - 0.3 * start40_sd, 
                    ymax = start40_mean + 0.3 * start40_sd,
                    color = "Start from 40th Tree"), width = 0.08,
                position = position_dodge(0.25), size = 0.25) +   
  geom_errorbar(aes(ymin = start45_mean - 0.25 * start45_sd, 
                    ymax = start45_mean + 0.25 * start45_sd,
                    color = "Start from 45th Tree"), width = 0.08,
                position = position_dodge(0.25), size = 0.25) +                                                     
    theme(legend.title = element_blank()) + theme_bw() + 
     guides(colour = guide_legend(nrow = 4, byrow = TRUE)) +        
    scale_color_manual(values = colors) + 
     theme(plot.title = element_text(size = 30)) + 
   theme(axis.text = element_text(size = 30)) + 
   theme(axis.title = element_text(size = 30)) + 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(legend.position='bottom') + 
  theme(legend.title = element_text(size = 32), 
        legend.text = element_text(size = 32)) + 
  theme(legend.title = element_blank()) +   
  theme(legend.key.size = unit(1, 'cm'))   
dev.off()

power_df <- read.csv("datasets/xgb_longterm_correlation_results_startfrom55.csv")
png('datasets/XGB_LongTermCorrelation_55trees.png', width = 900, height = 1050)
colors <- c(
           "Start from 55th Tree" = "red",
           "Start from 60th Tree" = "red",
           "Start from 65th Tree" = "blue",
           "Start from 70th Tree" = "orange",
           "Start from 75th Tree" = "orange",           
           "Start from 80th Tree" = "black",
           "Start from 85th Tree" = "green",           
           "Start from 90th Tree" = "purple"        
           )
ggplot(power_df, aes(x = c(1:58))) +
  geom_line(aes(y = start55_mean, group = 1, color = "Start from 55th Tree"), linewidth = 1.25) + 
  geom_line(aes(y = start60_mean, group = 1, color = "Start from 60th Tree"), linewidth = 1.25) + 
  geom_line(aes(y = start65_mean, group = 1, color = "Start from 65th Tree"), linewidth = 1.25) + 
  geom_line(aes(y = start70_mean, group = 1, color = "Start from 70th Tree"), linewidth = 1.25) + 
  geom_line(aes(y = start75_mean, group = 1, color = "Start from 75th Tree"), linewidth = 1.25) +   
  geom_line(aes(y = start80_mean, group = 1, color = "Start from 80th Tree"), linewidth = 1.25) +   
  geom_line(aes(y = start85_mean, group = 1, color = "Start from 85th Tree"), linewidth = 1.25) +   
  geom_line(aes(y = start90_mean, group = 1, color = "Start from 90th Tree"), linewidth = 1.25) +   
  labs(x = "N Rounds Lag", y = "Correlation", color = "Legend") + 
ggtitle("Long-Term Dependence for XGBoost Model \n 150 Trees in Total, 10 Signal Features") +   
  geom_errorbar(aes(ymin = start55_mean - 0.2 * start55_sd, 
                    ymax = start55_mean + 0.2 * start55_sd,
                    color = "Start from 55th Tree"), width = 0.08,
                position = position_dodge(0.25), size = 0.25) + 
  geom_errorbar(aes(ymin = start60_mean - 0.5 * start60_sd, 
                    ymax = start60_mean + 0.5 * start60_sd,
                    color = "Start from 60th Tree"), width = 0.08,
                position = position_dodge(0.25), size = 0.25) + 
  geom_errorbar(aes(ymin = start65_mean - 0.5 * start65_sd, 
                    ymax = start65_mean + 0.5 * start65_sd,
                    color = "Start from 65th Tree"), width = 0.08,
                position = position_dodge(0.25), size = 0.25) +   
  geom_errorbar(aes(ymin = start70_mean - 0.45 * start70_sd, 
                    ymax = start70_mean + 0.45 * start70_sd,
                    color = "Start from 70th Tree"), width = 0.08,
                position = position_dodge(0.25), size = 0.25) + 
  geom_errorbar(aes(ymin = start75_mean - 0.4 * start75_sd, 
                    ymax = start75_mean + 0.4 * start75_sd,
                    color = "Start from 75th Tree"), width = 0.08,
                position = position_dodge(0.25), size = 0.25) +   
  geom_errorbar(aes(ymin = start80_mean - 0.35 * start80_sd, 
                    ymax = start80_mean + 0.35 * start80_sd,
                    color = "Start from 80th Tree"), width = 0.08,
                position = position_dodge(0.25), size = 0.25) +   
  geom_errorbar(aes(ymin = start85_mean - 0.3 * start85_sd, 
                    ymax = start85_mean + 0.3 * start85_sd,
                    color = "Start from 85th Tree"), width = 0.08,
                position = position_dodge(0.25), size = 0.25) +   
  geom_errorbar(aes(ymin = start90_mean - 0.25 * start90_sd, 
                    ymax = start90_mean + 0.25 * start90_sd,
                    color = "Start from 90th Tree"), width = 0.08,
                position = position_dodge(0.25), size = 0.25) +                                                     
    theme(legend.title = element_blank()) + theme_bw() + 
     guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +        
    scale_color_manual(values = colors) + 
     theme(plot.title = element_text(size = 30)) + 
   theme(axis.text = element_text(size = 30)) + 
   theme(axis.title = element_text(size = 30)) + 
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(legend.position='bottom') + 
  theme(legend.title = element_text(size = 32), 
        legend.text = element_text(size = 32)) + 
  theme(legend.title = element_blank()) +   
  theme(legend.key.size = unit(1, 'cm'))   
dev.off()



