library(ggplot2)
p2data <- "https://raw.githubusercontent.com/datavizpyr/data/master/palmer_penguin_species.tsv"
tbl.penguins <- read.table(p2data, sep="\t", header=TRUE)
head(tbl.penguins)
dim(tbl.penguins)           # 344 7
lapply(tbl.penguins, class)
  #           species:  "character"
  #            island:  "character"
  #  culmen_length_mm:  "numeric"
  #   culmen_depth_mm:  "numeric"
  # flipper_length_mm:  "integer"
  #       body_mass_g:  "integer"
  #              sex:  "character"

table(tbl.penguins$species)  #    Adelie Chinstrap    Gentoo
                             #       152        68       124

table(tbl.penguins$sex)      # FEMALE   MALE
                             #    165    168

tbl.penguins <- subset(tbl.penguins, !is.na(sex))
ggplot(tbl.penguins,
       aes(x = species,
           y = flipper_length_mm,
           fill = species)) +
    geom_violin()


quartz()
ggplot(tbl.penguins,
       aes(x = species,
           y = flipper_length_mm,
           fill = sex)) +
    geom_violin()
