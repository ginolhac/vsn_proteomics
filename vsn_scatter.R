library("dplyr")
library("readr")
library("ggplot2")
library("reshape2")

mac <- read_tsv("fullmelt_24_2_1.25_drop_norm_arcsinh.tsv", col_names = TRUE) %>%
  mutate(platform = "mac")
deb <- read_tsv("deb_fullmelt_24_2_1.25_drop_norm_arcsinh.tsv", col_names = TRUE) %>%
  mutate(platform = "debian")
df <- inner_join(mac, deb, by = c("Accession", "condition", "variable"))

df %>%
  filter(variable == "Q.Value") %>%
  ggplot()+
  geom_point(aes(x = value.x, y = value.y))+
  scale_x_log10()+
  scale_y_log10()+
  theme_light(16)+
  xlab("log10(qvalue) mac")+
  ylab("log10(qvalue) debian")

df %>%
  filter(variable == "Q.Value") %>%
  mutate(diff = value.x - value.y) %>%
  ggplot()+
  geom_density(aes(x = diff), fill = "orange", alpha = 0.6)+
  theme_light(16)+
  xlim(-1e-3, 4e-3)

mac %>%
  dcast(Accession + platform + condition ~ variable) %>%
  tbl_df %>%
  filter(Q.Value <= 0.15)
# 722
deb %>%
  dcast(Accession + platform + condition ~ variable) %>%
  tbl_df %>%
  filter(Q.Value <= 0.15)
# 743
