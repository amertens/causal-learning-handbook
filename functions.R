
forest_plot <- function(tab, arm_strat=F){
  
  tab$label_f[tab$label_f == "Total Carbohydrate"] <- "Total Carbohydrates"
  tab$sigFDR <- factor(tab$sigFDR, levels=c(0,1))
  tab$sig <- factor(tab$sig, levels=c(0,1))
  tab$study <- factor(tab$study, levels=c("Misame", "Vital","Elicit"))
  tab$visit <- factor(tab$visit, levels=c( "14-21 days", "1-2 mo.","3-4 mo.", "1.5 mo.","2 mo.","1 mo.","5 mo."))
  tab=tab %>% mutate(sigcat=case_when(
    sigFDR == 1 & sig == 1 ~ "Sig",
    sigFDR == 0 & sig == 1 ~ "Sig before FDR",
    sigFDR == 0 & sig == 0 ~ "Not Significant"
  ), sigcat=factor(sigcat, levels=c("Not Significant", "Sig before FDR", "Sig")))
  
  unique(tab$visit)
  unique(tab$study)
  
  if(arm_strat){
    ggplot(tab, aes(x = reorder(label_f, -est), y = est, color=contrast, shape=sigcat, group=contrast)) +
      geom_point(position = position_dodge(width = 0.5), size = 2) +
      geom_errorbar(aes(ymin = cil, ymax = ciu), width = 0.2, position = position_dodge(width = 0.5)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      coord_flip() +
      scale_shape_manual(values = c(1, 19, 17), 
                         labels = c("Not Significant","Sig before FDR", "Sig")) +
      scale_color_manual(values = c(tableau10)) +
      facet_wrap(study~visit) +
      labs(title = "Forest Plot of Estimates",
           x = "Variable",
           y = "Estimate") +
      theme_minimal() + theme(legend.position = "right")
  }else{
    ggplot(tab, aes(x = reorder(label_f, -est), y = est, color=sigFDR, shape=sig, group=contrast)) +
      geom_point(position = position_dodge(width = 0.5), size = 2) +
      geom_errorbar(aes(ymin = cil, ymax = ciu), width = 0.2, position = position_dodge(width = 0.5)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      coord_flip() +
      scale_shape_manual(values = c(1, 19), 
                         labels = c("Not Significant", "Significant")) +
      scale_color_manual(values = c("grey60", tableau10[2]), 
                         labels = c("Not Significant", "Significant")) +
      facet_wrap(study~visit) +
      labs(title = "Forest Plot of Estimates",
           x = "Variable",
           y = "Estimate") +
      theme_minimal() + theme(legend.position = "none")
  }
  
  
}


clean_tab <- function(tab){
  rownames(tab)=NULL
  
  tab <- tab %>% 
    mutate(ATE = paste0(round(est, 2), " (", round(cil, 2), ", ", round(ciu, 2), ")"),
           pval_cat=case_when(
             pval < 0.001 ~ "***",
             pval < 0.01 ~ "**",
             pval < 0.05 ~ "*",
             TRUE ~ ""
           ),
           pval_adj_cat=case_when(
             pval_adj < 0.001 ~ "***",
             pval_adj < 0.01 ~ "**",
             pval_adj < 0.05 ~ "*",
             TRUE ~ ""
           ),
           pval=paste0(round(pval, 3),pval_cat), pval_adj = paste0(round(pval_adj, 3),pval_adj_cat),
    ) %>%
    subset(., select = -c( cil, ciu)) %>%
    select(study, visit, contrast, label_f, ATE, pval, pval_adj)
  
  tab = DT::datatable(tab) %>%
    formatStyle(columns = colnames(tab), fontSize = '11px')
  
  return(tab)
  
}
