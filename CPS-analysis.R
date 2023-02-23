#Title: Analysis of Combinatorial CPS States
#Author: Jason Saba
#Output: numerous descriptive files  

##Instructions: provide a folder with .xls extension files.
#Folder must also include file "AllCPScombos.xlsx"

library(readxl) #version 1.4.1
library(writexl) #version 1.4.0
library(ggridges) #version 0.5.4
library(gplots) #version 3.1.3
library(tidyverse) ###****tidyr version 1.1.4 is absolutely critical
library(remotes) #version 2.4.2.Use this to install versions specified here

##For installing necessary versions, use lines such as:
#install_version("tidyr", version = "1.1.4") 

#read in combos table for later analysis
all_combos <- readxl::read_excel(("AllCPScombos.xlsx"), 
                         sheet = "Sheet1", 
                         range = cell_cols("A")
)

##read data and define dataframe
fileNames <- Sys.glob("*.xls") ## Important!! this is specifically .xls; NOT .xlsx

for (fileName in fileNames) {

sampleID <- gsub(".xls*", "", fileName)

df <- readxl::read_excel(paste(sampleID, ".xls", sep=""), 
                         sheet = sampleID, 
                         range = cell_cols("A:O")
                         )

df <- df[,c(1,3,2,5,4,7,6,9,8,11,10,13,12,15,14)] #Swap column order to get PSX-ON then PSX-OFF
    
colnames(df) <- sub("...1", "Barcode", colnames(df))

#Separate 'ON' and 'OFF' Columns
df_on <- df %>% select(1, 2, 4, 6, 8, 10, 12, 14)
df_off <- df %>% select(1, 3, 5, 7, 9, 11, 13, 15)

##Append columns 
#First remove '_ON' and '_OFF' strings to have same header
colnames(df_on) <- sub('_ON', '', colnames(df_on))
colnames(df_off) <- sub('_OFF', '', colnames(df_off))

#Append rows
df_operon <- rbind(df_on, df_off)

##Determine basic info
barcount <- nrow(df) #total number of barcodes/cells

operon_totalRC  <- colSums(df_operon[,2:8]) #total reads per operon

operon_totalRC <- data.frame(as.list(operon_totalRC)) #convert to usable list for next step

operon_totalRC$TotalReads <- rowSums(operon_totalRC) #determine total reads and append column with name 'TotalReads'

operon_totalRC <- operon_totalRC %>% pivot_longer(everything()) #pivot longer for future transformations

operon_totalRC <- mutate(operon_totalRC, PercentTotalReads = 100*value/colSums(operon_totalRC[-8,2])) #determine operon % of total reads in library and append as column

#Make ridgeplot for each individual operon for reads greater than 0  
piv_oper <- df_operon %>% select(2:8) %>% pivot_longer(everything()) #Transpose table for the operon dfs 

oper_ridge <- filter(piv_oper, value > 0) %>% 
  mutate(name = fct_reorder(name, value, median, .desc = TRUE)) %>% 
  ggplot(aes(y = name, x = value, fill = name, group = name)) + 
  geom_density_ridges(
    scale = 1, 
    show.legend = FALSE, 
    quantile_lines = TRUE, 
    quantiles = 2, 
    from = 0, 
    to = 100
    ) +
  theme_bw() + labs(y = "", x = "read counts")

## Data Filtering
df2 <- df #replicate data to protect original 

df2$TotalReads <- rowSums(df2[2:15]) #determine total reads and append column with name 'TotalReads'

#Remove barcode if locus total reads < 1% of barcode total reads
find_less_one_percent_barcodes <- function(f) {
  f %>%
    # go from wide to long format
    pivot_longer(PSA_ON:PSH_OFF) %>%
    # split out the operon and position into two columns
    separate(name, into = c("operon", "position"), sep = "_") %>%
    group_by(Barcode, operon) %>%
    # determine locus reads percent of total reads per barcode
    summarize(percent_of_barcode_reads = sum(value)/TotalReads, .groups = "drop") %>%
    # look for barcodes that have a locus reads <1% of total reads per barcode 
    filter(percent_of_barcode_reads < 0.01)
}

filter_out_less_one_percent_barcodes <- function(g){
  # get the data for the bad barcodes
  bad_barcodes <- find_less_one_percent_barcodes(g)
  g %>% 
    filter(
      !(Barcode %in% 
          # pull the barcode column out of bad barcodes and 
          # only consider the unique ones
          (bad_barcodes %>% 
             pull(Barcode) %>% 
             unique()
          )
      )
    )
}

filter_out_greater_one_percent_barcodes <- function(g){
  # get the data for the bad barcodes
  bad_barcodes <- find_less_one_percent_barcodes(g)
  g %>% 
    filter(
      (Barcode %in% 
          # pull the barcode column out of bad barcodes and 
          # only consider the unique ones
          (bad_barcodes %>% 
             pull(Barcode) %>% 
             unique()
          )
      )
    )
}

df3 <- filter_out_less_one_percent_barcodes(df2)

oper_less_than_0.01_perBC <- nrow(df2)-nrow(df3) #counts the number of barcodes excluded by this step

#Call values for each orientation less than 1% of total 0, else = value 
df4 <- replace(df3, df3 < 1, 0)

df3 <- df3 %>% column_to_rownames(., var = 'Barcode')

df4 <- replace(df3, df3/df3$TotalReads < 0.01, 0)

df4 <- df4 %>% rownames_to_column()

df4 <- df4 %>% rename(Barcode = rowname) %>% tibble()

## Graph ridge plot after filtering 
#Separate 'ON' and 'OFF' Columns
f14_on <- df4 %>% select(1, 2, 4, 6, 8, 10, 12, 14)
f14_off <- df4 %>% select(1, 3, 5, 7, 9, 11, 13, 15)

#Remove '_ON' and '_OFF' strings to have same header
colnames(f14_on) <- sub('_ON', '', colnames(f14_on))
colnames(f14_off) <- sub('_OFF', '', colnames(f14_off))

#Append columns
f14_operon <- rbind(f14_on, f14_off)

#Transpose table for the operon dataframes
piv_f14 <- f14_operon %>% select(2:8) %>% pivot_longer(everything())

## Discretize data
#Replace values >0 with 1 
df5 <- df4[,1:15]

df5[,2:15][df5[,2:15] > 0] <- 1

#This step removes barcodes having both orientations ON or both OFF for any operon
find_bad_barcodes <- function(h) {
  h %>%
    # go from wide to long format
    pivot_longer(PSA_ON:PSH_OFF) %>%
    # split out the operon and position into two columns
    separate(name, into = c("operon", "position"), sep = "_") %>%
    group_by(Barcode, operon) %>%
    # get the total counts per operon and position
    summarize(total = sum(value), .groups = "drop") %>%
    # look for counts that aren't exactly 1
    filter(total != 1)
}

filter_out_bad_barcodes <- function(j){
  # get the data for the bad barcodes
  bad_barcodes <- find_bad_barcodes(j)
  j %>% 
    filter(
      !(Barcode %in% 
          # pull the barcode column out of bad barcodes and 
          # only consider the unique ones
          (bad_barcodes %>% 
             pull(Barcode) %>% 
             unique()
          )
      )
    )
}

filter_out_good_barcodes <- function(k){
  # get the data for the bad barcodes
  bad_barcodes <- find_bad_barcodes(k)
  k %>% 
    filter(
      (Barcode %in% 
          # pull the barcode column out of bad barcodes and 
          # only consider the unique ones
          (bad_barcodes %>% 
             pull(Barcode) %>% 
             unique()
          )
      )
    )
}

df6 <- filter_out_bad_barcodes(df5)

#Summarize tossed values
tossed <- c("Less than 1% reads" = oper_less_than_0.01_perBC, 
            "Both ON or OFF" = nrow(df5) - nrow(df6),
            "Sum" = oper_less_than_0.01_perBC + nrow(df5) - nrow(df6)
            )

tossed <- tossed %>% as.data.frame() %>% rownames_to_column() 

colnames(tossed) <- c("Removed Barcodes", "Count")

tossed_less_than_1pc <- filter_out_greater_one_percent_barcodes(df2)

tossed_bothONorOFF <- filter_out_good_barcodes(df5)

#Make ridgeplot after ALL filtering
un_discretized_final_data <- filter(df, df$Barcode %in% df6$Barcode)

undisc_operon_for_ridge <- filter(df_operon, df_operon$Barcode %in% df6$Barcode)

piv_undisc <- undisc_operon_for_ridge[, 2:8] %>% pivot_longer(everything())

post_filt_oper_ridge <- filter(piv_undisc, value > 0) %>% 
  mutate(name = fct_reorder(name, value, median, .desc = TRUE)) %>% 
  ggplot(aes(y = name, x = value, fill = name, group = name)) + 
  geom_density_ridges(scale = 1, show.legend = FALSE, quantile_lines = TRUE, quantiles = 2, from = 0, to = 100) +
  theme_bw() + labs(y = "", x = "read counts")

post_filt_oper_ridge 

## Analysis
#Isolate ON columns for analysis
on <- df6 %>% select(2, 4, 6, 8, 10, 12, 14)

#Change dataframe to character so can convert value of 1 to letters for next steps
on_char <- on
on_char <- data.frame(lapply(on_char, as.character), stringsAsFactors = FALSE)

#Duplicate for making table for Cytoscape network plot (using dashes in place of '')
on_char_net <- on_char

#Set value of '1' = 3rd letter in column name of numeric dataframe "on"; column name is 1-dimensional in
  #"on_char" dataframe after switching it to character 

for(i in 1:ncol(on_char)) {
  for(j in 1:nrow(on_char)) {
    if(on_char[j, i] == 1) {
      on_char[j,i] <- str_sub(colnames(on[,i]), 3, 3)}
    else { on_char[j,i] <- ""}
  }}

for(i in 1:ncol(on_char_net)) {
  for(j in 1:nrow(on_char_net)) {
    if(on_char_net[j, i] == 1) {
      on_char_net[j,i] <- str_sub(colnames(on[,i]), 3, 3)}
    else { on_char_net[j,i] <- "-"}
  }}

#Concatenate strings and convert to tibble for processing
combos <- unite(on_char, combo, sep = '') %>% tibble() 
combos_for_network <- unite(on_char_net, combo, sep = '') %>% tibble() 


#Create column with counts of each combo
combos_count <- combos %>% group_by(combo) %>% mutate(count = n())
combos_count_piv <- combos_count %>% pivot_wider()

#Same for combos network... 
combos_countNET <- combos_for_network %>% group_by(combo) %>% mutate(count = n())
combos_countNET_piv <- combos_countNET %>% pivot_wider()

#Add "OFF" label for empty values... 
for(i in 1:nrow(combos_count)){
  if(combos_count[i, 1] == '') {
    combos_count[i, 1] <- 'ALL OFF'
  }}

#Plot combo vs proportion 
combo_plot <- ggplot(combos_count, aes(x = reorder(combo, count))) +
  geom_bar() +
  theme_classic() + # removes all gridlines 
  theme(
    axis.ticks.y = element_blank(), # removes flipped y-axis ticks
    panel.grid.major.x = element_line(color = "gray", size = 0.1, linetype = "solid")) +
  scale_x_discrete(name = "Combo") +
  scale_y_continuous(
    trans=scales::pseudo_log_trans(base = 10),
    breaks = c(1:10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000),
    labels = c(1,'','','','','','','','', 10, '','','','','','','','', 100, '', '', '', '','','','','', 1000, '',''),
    expand = c(0,0)) + #forces y axis intercept to 0,0
  coord_flip()

##Determine number of promoters simultaneously on plot
#Create column which adds '1's from the ON table (on)
on_sums <-  on #first replicate 

on_sums$sums <- apply(on_sums, 1, sum)

#Plot histogram promoters simultaneously ON
simul_on_plot <- ggplot(on_sums, aes(x = sums)) +
  geom_bar() +
  theme_classic() +
  scale_x_continuous(
    name = "Promoters Simultaneously Oriented ON", 
    breaks = c(0, 1, 2, 3, 4, 5, 6, 7),
    labels = c(0, 1, 2, 3, 4, 5, 6, 7)) +
  scale_y_continuous(expand = c(0,0))  #forces y axis intercept to 0,0)

simul_on_plot

#Make table of count for promoters simultaneously ON
SimultONCount <- on_sums %>% group_by(sums) %>% mutate(count = n()) 
SimultONCount <- tibble(SimultONCount[,8:9]) %>% pivot_wider() 
SimultONCount <- arrange(SimultONCount, sums)

#norm by proportion of population
SimultONCount$PercentOfPop <- 100* SimultONCount$count/sum(SimultONCount$count)

##Compute Levenshtein distance between strings and make barplots
library(stringdist)
#"stringdistmatrix" counts character distance between combos, and "table" counts the elements

lv <- stringdistmatrix(combos_countNET_piv$combo, combos_countNET_piv$combo)
lvt <- table(lv)
lvt_df <- as.data.frame(lvt)

colnames(lv) <- combos_countNET_piv$combo
rownames(lv) <- combos_countNET_piv$combo

#have to divide by 2 b/c matrix doubles counts
lvt_df[-1,2] <- lvt_df[-1,2]/2

LVT_plot_distance <- ggplot(lvt_df[-1,], aes(lv, Freq)) +
  geom_col() +
  theme_classic() + 
  scale_x_discrete(
    name = "Hamming Distance (HD)"
  ) +
  scale_y_continuous(
    expand = c(0,0),
    name = "Subpopulation Pairs")  #forces y axis intercept to 0,0)

if(nrow(lvt_df) == 1) {
  LVT_plot_distance <- as.data.frame("NA")
} 

#set values not equal to 1 to 0, then sum values and append sum to original cleaned matrix
lv_df <- as.data.frame(lv)
lv_df_ones_only <- replace(lv_df, lv_df != 1, 0)
lv_df$one_sum = rowSums(lv_df_ones_only)

#Determine subpopulation's proportion in population (%) versus # of immediate neighbors
count_vs_onesum <- tibble(
  COMBO = rownames(lv_df), 
  ChkSameCombo = combos_countNET_piv$combo, 
  count = combos_countNET_piv$count,
  count_frac = combos_countNET_piv$count/sum(combos_countNET_piv$count), 
  onesum = lv_df$one_sum
  )

ggplot(count_vs_onesum, aes(x=onesum, y=count_frac)) + 
  geom_point() + 
  theme_classic() + 
  scale_y_log10(limits = c(0.0001,1), 
                labels = scales::percent,
                breaks = c(0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0)) +
  expand_limits(y = 0.01) +
  scale_x_discrete(limit = c('1', '2', '3', '4', '5', '6', '7')) +
  xlab("# of immediate neighbors (HD=1)") + 
  ylab("Subpopulation Size (%)")

count_vs_onesum$onesum <- as.factor(count_vs_onesum$onesum)

box_count_vs_onesum <- ggplot(count_vs_onesum, aes(x=onesum, y=count_frac)) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_y_log10(limits = c(0.0001,1), 
                labels = scales::percent,
                breaks = c(0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0)) +
  expand_limits(y = 0.01) +
  scale_x_discrete(limit = c('1', '2', '3', '4', '5', '6', '7')) +
  xlab("# of connected nodes") + 
  ylab("combo's proportion in population")

dot_count_vs_onesum <- box_count_vs_onesum + geom_dotplot(binaxis = 'y', dotsize = 0.5, stackdir = 'center')
dot_count_vs_onesum

#Determine count vs # of connected nodes (population barplot)
lvt_onesum <- table(lv_df$one_sum) %>% as.data.frame()

LVT_plot_node_number <- ggplot(lvt_onesum, aes(Var1, Freq)) +
  geom_col() +
  theme_classic() + 
  scale_x_discrete(
    name = "# of neighbors with HD = 1"
  ) +
  scale_y_continuous(
    expand = c(0,0),
    name = "Count")  #forces y axis intercept to 0,0)

if(nrow(lvt_onesum) == 1) {
  LVT_plot_node_number <- as.data.frame("NA")
}

##Plot HD from most abundant state vs. proportion in population
#Determine most abundant state
ind = which(count_vs_onesum[,'count'] == max(count_vs_onesum$count))
most_abundant = count_vs_onesum$COMBO[ind]

#Calculate Levenshtein distance from most abundant state and append value as column, then plot
count_vs_onesum$lv_from_max <- stringdist(most_abundant, count_vs_onesum$COMBO)

ggplot(count_vs_onesum, aes(x=lv_from_max, y=count_frac)) + 
  geom_point() + 
  theme_classic() + 
  scale_y_log10(limits = c(0.0001,1), 
                labels = scales::percent,
                breaks = c(0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0)) +
  expand_limits(y = 0.01) +
  scale_x_discrete(limit = c('1', '2', '3', '4', '5', '6', '7')) +
  xlab("HD from Most Abundant Subpopulation") + 
  ylab("Subpopulation Size (%)")

count_vs_onesum$lv_from_max <- as.factor(count_vs_onesum$lv_from_max)

box_count_vs_onesum <- ggplot(count_vs_onesum, aes(x=lv_from_max, y=count_frac)) + 
  geom_boxplot() + 
  theme_classic() + 
  scale_y_log10(limits = c(0.0001,1), 
                labels = scales::percent,
                breaks = c(0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0)) +
  expand_limits(y = 0.01) +
  scale_x_discrete(limit = c('1', '2', '3', '4', '5', '6', '7')) +
  xlab("HD from Most Abundant Subpopulation") + 
  ylab("Subpopulation Size (%)")

dot_count_vs_lv_from_max <- box_count_vs_onesum + geom_dotplot(binaxis = 'y', dotsize = 0.5, stackdir = 'center')
dot_count_vs_lv_from_max

no_zeros_lv_from_max <- count_vs_onesum %>% subset(count_vs_onesum$lv_from_max != "0")

#Determine Spearman correlation for proportion in population vs. HD from Most Abundant Subpopulation
if(nrow(no_zeros_lv_from_max) <= 1) {
  SpearCor_df_fix <- as.data.frame("NA")
  } else{
    SpearCor <- cor.test(
      no_zeros_lv_from_max$count_frac,
      as.numeric(no_zeros_lv_from_max$lv_from_max),
      method = "spearman",
      exact = FALSE);
    SpearCor_df <- data.frame(t(sapply(SpearCor,c)));
    SpearCor_df_fix <- data.frame(SpearCor_df$estimate, SpearCor_df$p.value)
  }

#column1 = state, col2 = lv distance from most abundant state, col3=prop in pop 

##Count number of barcodes before and after filtering
barcount_after_filter <- nrow(df6)

BarCount <- c("Before Filtering" = barcount, 
              "After Filtering" = barcount_after_filter, 
              "Excluded"=barcount - barcount_after_filter
              )

BarCount <- as.data.frame(BarCount) %>% rownames_to_column()

#Determine capsule as proportion of population irrespective of others
a <- 100*colSums(on[,1])/barcount_after_filter
b <- 100*colSums(on[,2])/barcount_after_filter
d <- 100*colSums(on[,3])/barcount_after_filter
e <- 100*colSums(on[,4])/barcount_after_filter
f <- 100*colSums(on[,5])/barcount_after_filter
g <- 100*colSums(on[,6])/barcount_after_filter
h <- 100*colSums(on[,7])/barcount_after_filter

prop_of_pop <- setNames(c(a, b, d, e, f, g, h), c("A", "B", "D", "E", "F", "G", "H"))

#Make heatmap for combos
all_combos_fracs <- all_combos %>% full_join(count_vs_onesum, by = c("Combo" = "COMBO"), keep = TRUE)
all_combos_fracs_clean <- cbind(all_combos_fracs[, 1], all_combos_fracs[,5])
all_combos_fracs_clean[is.na(all_combos_fracs_clean)] <- 0

as.matrix(cbind(all_combos_fracs_clean[,2],all_combos_fracs_clean[,2]))

my_palette2 <- colorRampPalette(c("white", "light green", "green", "dark green"))(n = 1000)

heatmap.2(
  as.matrix(cbind(all_combos_fracs_clean[1:64,2],all_combos_fracs_clean[1:64,2])),
  dendrogram = "none",
  Colv = NULL,
  Rowv = NULL,
  labRow = all_combos_fracs_clean$Combo,
  trace = "n",
  col = my_palette2,
  labCol = "",
  margins = c(5,10),
  keysize = 1,
  offsetRow = 1,
  cexRow = 0.5,
  cexCol = 0.1,
  key.title = NA
)

##Export file 
write_xlsx(list(
  BarCount = BarCount, 
  TossedValues = tossed,
  Discretized_ALL = df6, 
  Discretized_ON = on,
  Filtered_Data = un_discretized_final_data,
  CombosTableClean = combos_count_piv, 
  CombosTableForNetwork = combos_countNET_piv, 
  CapsuleProps = tibble(prop_of_pop), 
  CountSimultON = SimultONCount, 
  LvFromMax = count_vs_onesum, 
  LvFromMaxSpear = SpearCor_df_fix, 
  LvPairDistanceBarPlot = lvt_df, 
  LvOneSum = lvt_onesum,
  Tossed_less1percent = tossed_less_than_1pc,
  Tossed_bothONorOFF = tossed_bothONorOFF),
  paste(sampleID,"_Summary.xlsx", sep="")
  )

write_xlsx(list(
  LvFromMax = count_vs_onesum, 
  Spearman = SpearCor_df_fix), 
  paste(sampleID,"_LvFromMax.xlsx", sep="")
  )

##Make a heatmap
heat <- as.matrix(on)
col <- colorRampPalette(c("light gray","forest green"))

#Colv and Rowv turn off clustering. revC puts results in correct order top-down. col(2) makes binary color map. 
png(file= paste(sampleID, "_Heatmap.png", sep="")) 
heatmap.2(
  heat[1:128,1:7],
  Colv = NA, 
  dendrogram = "none",
  trace = "none",
  col = col(2), 
  scale = "none",
  labRow = NA,
  labCol = paste(
    round(prop_of_pop, digits = 0),
    "%",
    " ",
    str_sub(colnames(on), 3, 3)
  ),
  main = sampleID,
  ylab = "Single Cells",
  key = FALSE,
  lhei = c(5, 30),
  lwid = c(2, 5),
  reorderfun=function(d, w) reorder(d, w, agglo.FUN = max),
  margins = c(8,2)
)
dev.off()

#Other figures compiled in PDF
require(lattice)
pdf(file= paste(sampleID, "_Figures.pdf", sep=""))
print(oper_ridge)
print(post_filt_oper_ridge)
print(combo_plot)
print(simul_on_plot)
print(LVT_plot_distance)
print(LVT_plot_node_number)
print(dot_count_vs_onesum)
print(dot_count_vs_lv_from_max)

dev.off()

}

#Make files containing combo counts for network analyses
fileNames2 <- Sys.glob("*_Summary.xlsx")

for (fileName2 in fileNames2) {

sampleID2 <- gsub("_Summary.xlsx", "", fileName2)

df_network <- readxl::read_excel(paste(
  sampleID2, "_Summary.xlsx", sep=""),
  sheet = "CombosTableForNetwork",
  range = cell_cols("A:B")
  )

write_xlsx(list(CombosTableForNetwork = df_network), paste(sampleID2,"_CombosForNetwork.xlsx", sep=""))
}
