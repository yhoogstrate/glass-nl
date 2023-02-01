#!/usr/bin/env R

# requires to have Arial: - https://github.com/wch/extrafont/issues/32
#
# https://superuser.com/questions/1153990/anyone-know-how-to-install-arial-fonts-on-centos-7
# sudo yum install curl cabextract xorg-x11-font-utils fontconfig
# wget http://www.itzgeek.com/msttcore-fonts-2.0-3.noarch.rpm
# rpm -Uvh msttcore-fonts-2.0-3.noarch.rpm
# 
# 1. update remotes / Rttf2pt1:
# install.packages('remotes')
# remotes::install_version("Rttf2pt1", version = "1.3.8")
# 
# 2. install and update fonts: - https://fromthebottomoftheheap.net/2013/09/09/preparing-figures-for-plos-one-with-r/
# install.packages("extrafont")
# library("extrafont")
# font_import(')
# font_import(pattern='Arial')
# 
# 3. check if Arial exists:
# loadfonts()
# sort(fonts())

theme_cellpress <- theme_bw() +
  theme(
  text =          element_text(size = 7, family = "Arial", face = "plain"),
  axis.text =     element_text(size = 7, family = "Arial", face = "plain"),
  axis.title.x =  element_text(size = 7, family = "Arial", face = "plain"),
  axis.title.y =  element_text(size = 7, family = "Arial", face = "plain"),
  strip.text =    element_text(size = 7, family = "Arial", face = "plain"),
  
  legend.title =  element_text(size = 7, family = "Arial", face = "plain"),
  legend.text =   element_text(size = 7, family = "Arial", face = "plain"),
  
  plot.title =    element_text(size = 7, family = "Arial", face = "plain"), # `title` covers both title and subtitle
  plot.subtitle = element_text(size = 7, family = "Arial", face = "italic"),
  plot.caption =  element_text(size = 7, family = "Arial", face = "italic"),
  
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank()
)




youri_gg_theme <- theme(
  text = element_text(family = 'Helvetica'),
  #axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = 'bottom',
  plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
  panel.background = element_rect(fill = 'white', colour = 'white'),
  axis.title = element_text(face = "bold",size = rel(1)),
  axis.title.y = element_text(angle=90,vjust =2),
  axis.title.x = element_text(vjust = -0.2),
  axis.text = element_text(),
  axis.line = element_line(colour="black"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_line(colour = 'grey20', linetype = 'dotted'),
  panel.grid.minor.y = element_line(colour = 'grey50', linetype = 'dotted')
) + theme(
  panel.grid.major.x = element_line(colour = 'grey20', linetype = 'dotted',size=0.25),
  panel.grid.minor.x = element_line(colour = 'grey50', linetype = 'dotted')
) 
