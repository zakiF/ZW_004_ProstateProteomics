x <- c(10,20,15,16,28,12,30)
y <- c(1,2,1,6,8,2,3)


t.test(y,x, alternative = "less")
uu <- unique(dat_NPX_anno2$Panel)
d1 <- dat_NPX_anno_all %>% dplyr::filter(Panel == uu[2])

dat_tt <- 
  olink_ttest(df = dat_NPX_anno_all,
              variable = 'c_a_vs_b')
as.data.frame(dat_tt) %>% head(n=10) %>% select(Assay, p.value, Adjusted_pval, Threshold)
#dplyr::filter(dat_tt, Assay == "HPGDS") %>% select(Assay, p.value, Adjusted_pval, Threshold)

dat_tt_2 <- 
  olink_ttest_2(df = dat_NPX_anno_all,
                variable = 'c_a_vs_b')
as.data.frame(dat_tt_2) %>% head(n=10) %>% select(Assay, p.value, Adjusted_pval, Threshold)


dd <- d1 %>% filter(Assay == "INSL3")
ggplot(dd, aes(x=c_a_vs_b, y=NPX, fill=c_a_vs_b)) +
  geom_boxplot(outlier.alpha = 0,  width = 0.2, alpha=0.5, #, coef=0,
               position = position_dodge(width = 0.9)) +
  geom_quasirandom(dodge.width=0.9, colour="grey10", alpha=0.5) +
  scale_fill_manual(values=unique(pal.tmp)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"),
        axis.text.x = element_text(angle = 0, hjust = 1)) +
  ylab("NPX") +
  xlab("") +
  NULL



dd <- dat_NPX_anno_all %>% filter(Assay == "GSTA1")
ggplot(dd, aes(x=cohort_simple, y=NPX, fill=cohort_simple)) +
  geom_boxplot(outlier.alpha = 0,  width = 0.2, alpha=0.5, #, coef=0,
               position = position_dodge(width = 0.9)) +
  geom_quasirandom(dodge.width=0.9, colour="grey10", alpha=0.5) +
  scale_fill_manual(values=unique(pal.cohort_main)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"),
        axis.text.x = element_text(angle = 0, hjust = 1)) +
  ylab(paste(unique(dd$Assay))) +
  facet_wrap(~Assay, ncol=5, scale="free_y") +
  xlab("") +
  NULL

dd <- dat_NPX_anno_all %>% filter(OlinkID %in% head(row_4))
dd <- dat_NPX_anno_all %>% filter(Assay == "INSL3") %>%
  filter(!chaarted_volume == "N/A")
ggplot(dd, aes(x=cohort_simple, y=NPX, fill=chaarted_volume)) +
  geom_boxplot(outlier.alpha = 0,  width = 0.2, alpha=0.5, #, coef=0,
               position = position_dodge(width = 0.9)) +
  geom_quasirandom(dodge.width=0.9, colour="grey10", alpha=0.5) +
  #scale_fill_manual(values=c("darkgreen", "grey90")) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"),
        axis.text.x = element_text(angle = 0, hjust = 1)) +
  ylab(paste(unique(dd$Assay))) +
  xlab("") +
  facet_wrap(~Assay, ncol=5, scale="free_y") +
  NULL



dplyr::filter(df_group_de, Cat == "c_b_vs_a", Assay == "CDH6")





dat_NPX_anno_all <- dat_NPX_anno2 %>%
  dplyr::mutate(
    #c_a = case_when(cohort_simple == "A" ~ "Y", TRUE ~ "N"),
    #c_a_vs_b = case_when(cohort_simple == "A" ~ "Y", cohort_simple == "B" ~ "N"),
    #c_a_vs_c = case_when(cohort_simple == "A" ~ "Y", cohort_simple == "C" ~ "N"),
    c_a_vs_d = case_when(cohort_simple == "A" ~ "Y", cohort_simple == "D" ~ "N"),
    
    #c_b = case_when(cohort_simple == "B" ~ "Y", TRUE ~ "N"),
    #c_b_vs_a = case_when(cohort_simple == "B" ~ "Y", cohort_simple == "A" ~ "N"),
    #c_b_vs_c = case_when(cohort_simple == "B" ~ "Y", cohort_simple == "C" ~ "N"),
    c_b_vs_d = case_when(cohort_simple == "B" ~ "Y", cohort_simple == "D" ~ "N"),
    
    #c_c = case_when(cohort_simple == "C" ~ "Y", TRUE ~ "N"),
    #c_c_vs_a = case_when(cohort_simple == "C" ~ "Y", cohort_simple == "A" ~ "N"),
    #c_c_vs_b = case_when(cohort_simple == "C" ~ "Y", cohort_simple == "B" ~ "N"),
    c_c_vs_d = case_when(cohort_simple == "C" ~ "Y", cohort_simple == "D" ~ "N"),
    
    #c_d = case_when(cohort_simple == "D" ~ "Y", TRUE ~ "N"),
    #c_d_vs_a = case_when(cohort_simple == "D" ~ "Y", cohort_simple == "A" ~ "N"),
    #c_d_vs_b = case_when(cohort_simple == "D" ~ "Y", cohort_simple == "B" ~ "N"),
    #c_d_vs_c = case_when(cohort_simple == "D" ~ "Y", cohort_simple == "C" ~ "N"),
  )



dplyr::filter(df_group_de, Cat == "c_b_vs_a")
dplyr::filter(combined_o, Assay == "INSL3")
