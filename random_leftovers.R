#Unused / old code i might want


# PLot a histogram to see the frequency of when sampling takes place
ccp_blood_data_long %>% 
  filter(onset_to_lab >= 0) %>% 
  select(onset_to_lab, assay) %>% 
  filter(assay %in% c("lymphocyte_count", "neutrophil_count", "platelet_count")) %>% 
  ggplot(aes(onset_to_lab), colour = assay) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(0, 50)) +
  facet_wrap(vars(assay)) +
  scale_colour_manual(values = c("#e31a1c", "#ff7f00", "#c51b7d",
                                 "#023858", "#3690c0", "#f768a1")) +
  labs(x = "Onset to sample (days)", y = "Number of samples", colour = "Assay") + 
  theme(legend.position = "top")
