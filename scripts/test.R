regression_att <- nls(formula = slope ~ slope[1] * exp(-(DEPTH-DEPTH[1])/star), data = iPOC, start = list(star=50))
star = summary(regression_att)$parameters[1]

z = c(0:800)
regression_data <- data.frame(z, slope = iPOC$slope[1] * exp(-(z-iPOC$DEPTH[1])/star))




regression_att <- nls(formula = slope ~ slope[1] * (DEPTH+1) ^(-b), data = iPOC, start = list(b=1))
b = summary(regression_att)$parameters[1]

z = c(0:800)
regression_data <- data.frame(z, slope = iPOC$slope[1] * (z) ^(-b))




if(i==1){
  Att_plots <- ggplot(data = iPOC, aes(x = DEPTH, y = slope)) +
                            geom_ribbon(aes(ymin = slope_min, ymax = slope_max), width = 18, col="darkgrey", alpha = 0.1) +
                            geom_point(col="black",fill = colors) +
                            # geom_point(shape=21, fill = colors, size=2) +
                            # geom_function(fun = f, args = list(a=a, b=b), color = "darkred", size = 1) +
                            geom_line(data = regression_data, aes(x = z, y = slope), color ="darkred", size = 0.5) +
                            scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(800,0)) +
                            scale_y_continuous(name = expression(paste(Seasonal~accumulation,"  [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0,0),limits = c(0,90)) +
                            coord_flip() +
                            # geom_line(data = regression_data, aes(x=z, y=pred_slope), color = 'red')+
                            theme_bw() +
                            theme(axis.text = element_text(size = 11, colour = "black"),
                                  axis.title = element_text(size = 11, colour = "black"),
                                  axis.text.x=element_text(angle=0),
                                  plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
                            annotation_custom(grobTree(textGrob(paste("b = ",format(b, scientific = F, digits = 3)), x=0.56,  y=0.75, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            # annotation_custom(grobTree(textGrob(paste("c = ",format(c, scientific = F, digits = 3)), x=0.56,  y=0.75, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            # annotation_custom(grobTree(textGrob(paste("l = ",format(b, scientific = F, digits = 4)), x=0.56,  y=0.70, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            # annotation_custom(grobTree(textGrob(paste("C = ",format(c, scientific = F, digits = 3)), x=0.56,  y=0.65, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            # annotation_custom(grobTree(textGrob(paste("max = ",format(slopes$slope[1], scientific = F, digits = 3)), x=0.56,  y=0.60, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            # annotation_custom(grobTree(textGrob(paste("min = ",format(slopes$slope[6], scientific = F, digits = 3)), x=0.56,  y=0.55, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            
                            ggtitle(names[i])
}

if(i!=5){
  Att_plots <- Att_plots +
                            geom_ribbon(aes(ymin = slope_min, ymax = slope_max), width = 18, col="darkgrey", alpha = 0.1) +
                            geom_point(col="black",fill = colors) +
                            # geom_point(shape=21, fill = colors, size=2) +
                            # geom_function(fun = f, args = list(a=a, b=b), color = "darkred", size = 1) +
                            geom_line(data = regression_data, aes(x = z, y = slope), color ="darkred", size = 0.5) +
                            scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(800,0)) +
                            scale_y_continuous(name = expression(paste(Seasonal~accumulation,"  [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0,0),limits = c(0,90)) +
                            coord_flip() +
                            # geom_line(data = regression_data, aes(x=z, y=pred_slope), color = 'red')+
                            theme_bw() +
                            theme(axis.text = element_text(size = 11, colour = "black"),
                                  axis.title = element_text(size = 11, colour = "black"),
                                  axis.text.x=element_text(angle=0),
                                  plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
                            annotation_custom(grobTree(textGrob(paste("b = ",format(b, scientific = F, digits = 3)), x=0.56,  y=0.75, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            # annotation_custom(grobTree(textGrob(paste("c = ",format(c, scientific = F, digits = 3)), x=0.56,  y=0.75, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            # annotation_custom(grobTree(textGrob(paste("l = ",format(b, scientific = F, digits = 4)), x=0.56,  y=0.70, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            # annotation_custom(grobTree(textGrob(paste("C = ",format(c, scientific = F, digits = 3)), x=0.56,  y=0.65, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            # annotation_custom(grobTree(textGrob(paste("max = ",format(slopes$slope[1], scientific = F, digits = 3)), x=0.56,  y=0.60, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            # annotation_custom(grobTree(textGrob(paste("min = ",format(slopes$slope[6], scientific = F, digits = 3)), x=0.56,  y=0.55, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            
                            ggtitle(names[i])
}

if(i==5){
  Att_plots <- Att_plots +
                            geom_col(col="black",fill = colors[-(7:8)]) +
                            geom_errorbar(aes(ymin = slope_min, ymax = slope_max), width = 18, col="darkgrey") +
                            # geom_point(shape=21, fill = colors[-(7:8)], size=2) +
                            # xlim(c(0,40)) +
                            # geom_function(fun = f, args = list(a=a, b=b), color = "darkred",size = 1) +
                            geom_line(data = regression_data[(1:551),], aes(x = z, y = flux), color ="darkred", size = 1) +
                            scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(800,0)) +
                            scale_y_continuous(name = expression(paste(Accumulation~saisonnière,"  [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0,0),limits = c(0,90)) +
                            coord_flip() +
                            # geom_line(data = regression_data, aes(x=z, y=pred_slope), color = 'red')+
                            theme_bw() +
                            theme(axis.text = element_text(size = 11, colour = "black"),
                                  axis.title = element_text(size = 11, colour = "black"),
                                  axis.text.x=element_text(angle=0),
                                  plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
                            annotation_custom(grobTree(textGrob(paste("b = ",format(b, scientific = F, digits = 3)), x=0.56,  y=0.75, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            # annotation_custom(grobTree(textGrob(paste("c = ",format(c, scientific = F, digits = 3)), x=0.56,  y=0.75, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            # annotation_custom(grobTree(textGrob(paste("l = ",format(b, scientific = F, digits = 4)), x=0.56,  y=0.70, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            # annotation_custom(grobTree(textGrob(paste("C = ",format(c, scientific = F, digits = 3)), x=0.56,  y=0.65, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            # annotation_custom(grobTree(textGrob(paste("max = ",format(slopes$slope[1], scientific = F, digits = 3)), x=0.56,  y=0.60, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            # annotation_custom(grobTree(textGrob(paste("min = ",format(slopes$slope[6], scientific = F, digits = 3)), x=0.56,  y=0.55, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) +
                            
                            ggtitle(names[i])
}
