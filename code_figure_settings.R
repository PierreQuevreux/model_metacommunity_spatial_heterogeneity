# THEMES # ----
theme<-theme_gray()+
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(colour='grey'),
        panel.grid.major.y = element_line(colour='grey'),
        strip.background = element_rect(fill=NULL),
        text = element_text(family="serif",size=20),
        axis.text = element_text(size=20),
        axis.line = element_line(),
        legend.key=element_blank(),
        plot.title = element_text(hjust = 0.5))

theme_raster<-theme_gray()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family="serif",size=20),
        axis.text = element_text(size=20),
        axis.line = element_line(),
        legend.key=element_blank(),
        legend.text.align=0,
        plot.title = element_text(hjust = 0.5))

# SCALES # ----
corr_colour_TL_2<-scale_colour_manual(values=c("dodgerblue3","chocolate1"),
                                      labels=c("1","2"),
                                      guide = guide_legend(reverse = TRUE),
                                      name='trophic\nlevel')
corr_colour_TL_3<-scale_colour_manual(values=c("dodgerblue3","chocolate1","chartreuse4"),
                                      labels=c("1","2","3"),
                                      guide = guide_legend(reverse = TRUE),
                                      name='trophic\nlevel')
corr_colour_TL_4<-scale_colour_manual(values=c("dodgerblue3","chocolate1","chartreuse4","red"),
                                      labels=c("1","2","3","4"),
                                      guide = guide_legend(reverse = TRUE),
                                      name='trophic\nlevel')
corr_colour_TL_4_fill<-scale_fill_manual(values=c("dodgerblue3","chocolate1","chartreuse4","red"),
                                         labels=c("1","2","3","4"),
                                         guide = guide_legend(reverse = TRUE),
                                         name='trophic\nlevel')
corr_colour_TL_2_TS<-scale_colour_manual(values=c("dodgerblue1","chocolate1","dodgerblue4","chocolate4"),
                                         guide = guide_legend(reverse = TRUE),
                                         name='trophic\nlevel')
patch_line<-scale_linetype_manual(values=c("solid","22"),
                                  name='patch',
                                  labels=c("fast","slow"))
corr_line_TL_2_TS<-scale_linetype_manual(values=c("solid","solid","22","22"),
                                         guide = guide_legend(reverse = TRUE),
                                         name='Trophic\nlevel')
CV_colour_grad<-scale_fill_gradient2(low = "red",
                                     mid = "white",
                                     high = "blue",
                                     midpoint = 0,
                                     name="covariance")
corr_colour_grad<-scale_fill_gradient2(low = "red",
                                       mid = "white",
                                       high = "blue",
                                       midpoint = 0,
                                       limits = c(-1,1),
                                       name="correlation\ncoefficient")
corr_diff_colour_grad<-scale_fill_gradient2(low = "red",
                                            mid = "white",
                                            high = "blue",
                                            na.value = "grey92",
                                            midpoint = 0,
                                            name="correlation\ndifference")
biomass_diff_colour_grad<-scale_fill_gradient2(low = "red",
                                               mid = "white",
                                               high = "blue",
                                               na.value = "grey92",
                                               midpoint = 0,
                                               name="biomass\ndifference")
scaled_biomass_diff_colour_grad<-scale_fill_gradient2(low = "red",
                                               mid = "white",
                                               high = "blue",
                                               na.value = "grey92",
                                               midpoint = 0,
                                               name="scaled\nbiomass\ndifference")
CV_diff_colour_grad<-scale_fill_gradient2(low = "red",
                                          mid = "white",
                                          high = "blue",
                                          na.value = "grey92",
                                          #trans = "log10",
                                          midpoint = 0,
                                          name="CV\ndifference")
jacobian_colour_grad<-scale_fill_gradient2(low = "red",
                                           mid = "white",
                                           high = "blue",
                                           na.value = "grey92",
                                           midpoint = 0,
                                           name="Jacobian\ncoef.")
CV_poptot_colour<-scale_colour_manual(values=c("darkseagreen3","darkorchid3"),
                                      labels=c(expression("CV"["pop"]),expression("CV"["tot"])))
perturbation_line<-scale_linetype_manual(values=c("solid","22","dotted"),
                                         labels=c("fast patch","slow patch","fast and slow patches"),
                                         name='perturbation in')
perturbation_pred_line<-scale_linetype_manual(values=c("solid","22","dotted"),
                                              labels=c("fast patch","slow patch","fast and slow patches"),
                                              name='perturbation\nof pred in the')
perturbation_prey_line<-scale_linetype_manual(values=c("solid","22","dotted"),
                                              labels=c("fast patch","slow patch","fast and slow patches"),
                                              name='perturbation\nof prey in the')
perturbation_prey_colour<-scale_colour_manual(values=c("cyan4","darkgoldenrod1"),
                                              labels=c("prey in patch #1","prey in patch #2"),
                                              name='perturbation of')
pert_colour_TL_2<-scale_colour_manual(values=c("dodgerblue3","chocolate1"),
                                      labels=c("1","2"),
                                      guide = guide_legend(reverse = TRUE),
                                      name='perturbed\nspecies')
fill_resilience<-scale_fill_viridis(name="Asymptotic\nresilience",
                                    trans = "log10")
corr_alpha<-scale_alpha_discrete(name="correlation of\nperturbations",
                                 range=c(0.5,1))
biomass_colour_grad<-scale_fill_viridis(name="Biomass",
                                        trans = "log10")

# AXIS # ----
x_axis_gamma<-scale_x_continuous(breaks=seq(1,10,3))
x_axis_log10_short<-scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
y_axis_log10_short<-scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
x_axis_log10<-scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))
y_axis_log10<-scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))
# LABELS # ----
# ea and ma factors
ma<-expression(italic(ma))
ma001<-expression(italic(ma)*"=0.01")
ma01<-expression(italic(ma)*"=0.1")
ma1<-expression(italic(ma)*"=1")
ma10<-expression(italic(ma)*"=10")
ma100<-expression(italic(ma)*"=100")
ea<-expression(italic("\u03B5")*italic(a))
ea001<-expression(italic("\u03B5")*italic(a)*"=0.01")
ea01<-expression(italic("\u03B5")*italic(a)*"=0.1")
ea1<-expression(italic("\u03B5")*italic(a)*"=1")
ea10<-expression(italic("\u03B5")*italic(a)*"=10")
ea100<-expression(italic("\u03B5")*italic(a)*"=100")

# labels
label_dispersal<-expression("Scaled dispersal rate "*italic(d["i"]))
label_correlation<-"Correlation between the two patches"
label_correlation_intra<-"Correlation between predator and prey"
label_CV<-"Coefficient of variation (CV)"
label_CV_pop<-"Biomass CV at population scale"
label_CV_metapop<-"Biomass CV at metapopulation scale"
label_CV_metacom<-"Biomass CV at metacommunity scale"
label_M1<-"Relative importance of dispersal"
label_lambda<-expression("Strength of top-down control "*italic("\u03BB"))
label_max_gamma<-expression("Maximum asymmetry coefficient "*italic("\u03B3"))
label_gamma<-expression("Asymmetry of interaction strength "*italic("\u03B3"))
label_omega<-expression("Asymmetry of resource supply "*italic("\u03C9"))
label_asymmetry<-expression("Asymmetry ("*italic("\u03B3 , \u03C9")*")")
label_asymmetry_short<-expression(italic("\u03B3 , \u03C9"))
label_contribution<-"Relative contribution to the lead eigenvector"
label_contribution_2lines<-expression(atop("Relative contribution","to the lead eigenvector"))
label_resilience<-"Asymptotic resilience"
label_jacobian<-"Direct effect"

# perturbations
pert_11<-"prey perturbed in fast patch"
pert_12<-"prey perturbed in slow patch"
pert_21<-"pred. perturbed in fast patch"
pert_22<-"pred. perturbed in slow patch"

# patches
patch_1<-"fast patch"
patch_2<-"slow patch"
patch_labels<-c('1'=patch_1, '2'=patch_2)
