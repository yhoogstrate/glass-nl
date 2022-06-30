#!/usr/bin/env

library(ggplot2)
library(pathwork)
library(gghalves) #devtools::install_github('erocoar/gghalves')


expr <- data.frame(expr = c(
  13, 11,
  4, 3,
  12, 8,
  7, 6),
  pat = c('A','A','B','B','C','C','D','D'),
  cond = c('primary','recurrent','primary','recurrent','primary','recurrent','primary','recurrent')
) %>% 
  dplyr::mutate(uid = 1:n()) %>% 
  dplyr::mutate(x = (as.numeric(as.factor(pat))-1) * 4 + as.numeric(cond == "recurrent"))


# expr2 <- rbind(expr, expr %>%  dplyr::mutate(expr=NULL))  


p1 <- ggplot(expr, aes(x=x, y=expr, col=pat,shape=cond)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=expr),col="black") +
  geom_point(size=4) +
  theme_bw() + 
  labs(y="BRAF expression", x=NULL)
  
p2 <- ggplot(expr, aes(x=cond, y=expr, fill=cond)) +
  geom_violin(draw_quantiles = c(0.5), alpha=0.5) +
  ggbeeswarm::geom_beeswarm(size=4,side = 1L) +
  theme_bw() + 
  labs(y="BRAF expression", x=NULL)

p1 + p2

# ggsave("/tmp/vis1.png",width=12,height=4)



expr.norm <- expr %>%
  dplyr::group_by(pat) %>%
  dplyr::mutate(pat.mean = mean(expr)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(expr = expr - pat.mean + 8.5)

p1 <- ggplot(expr.norm, aes(x=x, y=expr, col=pat,shape=cond)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=expr),col="black") +
  geom_point(size=4) +
  theme_bw() + 
  labs(y="BRAF expression", x=NULL)

p2 <- ggplot(expr.norm, aes(x=cond, y=expr, fill=cond)) +
  geom_violin(draw_quantiles = c(0.5), alpha=0.5) +
  ggbeeswarm::geom_beeswarm(size=4,side = 1L) +
  theme_bw() + 
  labs(y="BRAF expression", x=NULL)

p1 + p2

#ggsave("/tmp/vis2.png",width=12,height=4)




