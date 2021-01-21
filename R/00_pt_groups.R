library(dplyr)
library(ggplot2)
library(ggrepel)
library(forcats)
library(scales)

CDE <- read.csv(file = '/Users/senosam/Documents/Massion_lab/CyTOF_summary/CDE_TMA36_2020FEB25_SA.csv')

pie_chart <- function(dt, col_name, plot_title){
    x <- data.frame(table(dt[col_name]))
    colnames(x) <- c('Group', 'value')
    x <- x %>% 
        arrange(desc(value)) %>%
        mutate(prop = percent(value / sum(value)))

    print(ggplot(x, aes(x = "", y = value, fill = fct_inorder(Group))) +
        geom_bar(width = 1, stat = "identity") +
        coord_polar("y", start = 0) +
        geom_label_repel(aes(label = value), size=5, show.legend = F, nudge_x = 1) +
        #geom_label_repel(aes(label = prop), size=5, show.legend = F, nudge_x = 1) + #%
        guides(fill = guide_legend(title = "Group")) +
        labs(title=plot_title) +
        theme(plot.title = element_text(hjust = 0.5, size=22)))
}


pie_chart(CDE, 'CANARY', 'CANARY')
pie_chart(CDE, 'Smoking_Status', 'Smoking Status')
pie_chart(CDE, 'Stages_simplified', 'Stage')
pie_chart(CDE, 'Hist_predominant', 'Histology')
pie_chart(CDE, 'Death_st', 'Death Status')

flag = c()

for (i in 1:nrow(CDE)){
    # if SILA >=0.5 stage shouldn't be 0 (or 1), histology shouldn't be lepidic
    if (CDE$SILA[i]>=0.5) {
        if( CDE$Stages_simplified[i] %in% c('Stage 0', 'Stage I') || 
            CDE$Hist_predominant[i] %in% c('lepidic')) {
            cat('flag!: ', CDE$pt_ID[i], '\n')
            flag = c(flag, i)
        }
    }
    
    if (CDE$SILA[i]<0.5) {
        if ( CDE$Stages_simplified[i] %in% c('Stage III', 'Stage IV') || 
            CDE$Hist_predominant[i] %in% c('solid', 'micropapillary') || 
            CDE$Death_st[i] == 'Yes' || 
            CDE$Recurrence_st[i] == 'Yes' || 
            CDE$Progression_st[i] == 'Yes' ||
            CDE$CANARY[i] == 'P' ){
            cat('flag!: ', CDE$pt_ID[i], '\n')
            flag = c(flag, i)
        }

    }
}



fl_df <- CDE[flag, c('pt_ID','SILA', 'CANARY', 'Stages_simplified', 'X8th_ed_path_stage', 
    'Hist_predominant', 'Hist_other_patterns', 'Death_st', 'Recurrence_st', 
    'Progression_st')]
fl_df[which(fl_df$CANARY == 'G'),]
fl_df[which(fl_df$CANARY == 'I'),]
fl_df[which(fl_df$CANARY == 'P'),]