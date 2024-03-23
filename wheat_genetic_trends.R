# Genetic trends analysis ----

# Objective ----
# - Organize and check phenotypic data
# - Run single trial BLUP models to calculate reliability
# - Run single trial BLUE models to estimate genotype's blues
# - Calculate the genetic correlations among traits
# - Trend analysis (genetic, non-genetic, and total phenotypic)

rm(list=objects()) # clean workspace

# Packages ----
library(tidyverse) # R packages for data science
library(janitor) # Simple Tools for Examining and Cleaning Dirty Data
library(asreml) # ASReml-R package
library(asremlPlus) # Augments 'ASReml-R'
library(corrplot) # A visualization of a correlation matrix
library(ggtext) # Improved text rendering support for ggplot2
asreml.options(maxit=250) # Maximum number of iterations
asreml.options(trace=F) # Report convergence monitoring

# Phenotypic data ----

## Raw data ----
# Read, clean, and organize data from Breedbase

dat1w <- #raw data wide format
  read.csv('data/gen_trend_dat_2022nov10.csv', skip = 3) |> 
  remove_empty(which = c('cols')) |>#remove columns entirely empty
  clean_names() |> # clean names
  # select variables to be used
  dplyr::select(observation_unit_name,study_name,study_year,location_name,
                block_number, row_number, col_number, germplasm_name,
                heading_time_julian_date_jd_co_321_0001233,
                plant_height_cm_co_321_0001301,
                grain_yield_kg_ha_co_321_0001218,
                grain_test_weight_g_l_co_321_0001210) |>
  # rename variables
  rename(id=observation_unit_name,
         year=study_year,
         loc=location_name,
         env=study_name,
         blk=block_number,
         row=row_number,
         col=col_number,
         gen=germplasm_name,
         heading_time=heading_time_julian_date_jd_co_321_0001233,
         plant_height=plant_height_cm_co_321_0001301,
         grain_yield=grain_yield_kg_ha_co_321_0001218,
         test_weight=grain_test_weight_g_l_co_321_0001210) |>
  arrange(year,loc,gen,blk,row,col) |> #arrange data
  # replace location name
  mutate(loc = str_replace(loc,'St. Jacob Township, IL','St. Jacob, IL')) |>
  # convert to factors
  mutate_at(vars(id:gen),as.factor) |>
  # convert to numeric
  mutate_at(vars(heading_time:test_weight),as.numeric) |>
  # replace 0 by NA for response variables
  mutate_at(vars(heading_time:test_weight), ~ifelse(.==0,NA,.)) |>
  # filter out AdvHY trials (Advanced High Yield trial)
  filter(!grepl('AdvHY',env)) |>
  # drop unused levels
  droplevels() |>
  glimpse()

# Raw data 1 long format

dat1l <- dat1w |> #raw data long format
  pivot_longer(cols = heading_time:test_weight,
               names_to = 'var', #variables
               values_to = 'val', #values
               values_drop_na = T) |>
  mutate(var=as.factor(var)) |>
  glimpse()

## Data summary ----

# Number of observations per trait
# Number of years, locations, environments, genotypes, and total observations per trait
dat1l |>
  group_by(var) |>
  mutate(tot=length(id)) |>
  summarise(nyear=length(unique(year)),
            nloc=length(unique(loc)),
            nenv=length(unique(env)),
            ngen=length(unique(gen)),
            tot=unique(tot)) |>
  mutate(var=c('Grain yield', 'Heading time', 'Plant height', 'Test weight'))

# Total number of observations
nrow(dat1l)

### Fig. traits by location and year ----
dat1l |>
  group_by(year,loc,var) |>
  summarise(no_gen=length(unique(gen))) |>
  pivot_wider(names_from = loc,
              values_from = no_gen) |>
  pivot_longer(names_to = 'loc',
               values_to = 'no_gen',
               cols = 3:8) |>
  ggplot(aes(year,loc)) +
  geom_tile(aes(fill=no_gen), color='black', na.rm = T) +
  facet_wrap(~var, 
             labeller=labeller(var=c('grain_yield'='Grain yield',
                                     'heading_time'='Heading time',
                                     'plant_height'='Plant height',
                                     'test_weight'='Test weight'))) +
  scale_fill_gradient(name = 'Number of \ngenotypes',
                      high = '#FF552E', low = '#13294B',
                      na.value='white') +
  xlab('Year') + ylab('Location') +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.5, vjust=0.5, size = 8),
        plot.margin = margin(10,25,0,0),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14), panel.grid = element_blank())
ggsave("figures/traits_loc.png", width = 6.5, height = 4, units = "in", dpi=320)

## Single trial analysis ----
asreml.options(trace=F) # Report convergence monitoring

### BLUP model ----
# BLUP models to calculate blup for each genotype within trial, blups will be used to calculate reliability

# Model
mod.blups <-
  function(dat){
    yr <- unique(as.character(dat$year)) # year being analyzed
    if('2021' %in% yr){
      # Model for 2021 (including row and column)
      mod <- asreml(fixed = val~1+blk,
                    random = ~gen+row+col,
                    data = dat, na.action = na.method(y='include', x='include'))
    }else{
      # Model for other years (2001 to 2020)
      mod <- asreml(fixed = val~1+blk,
                    random = ~gen,
                    data = dat, na.action = na.method(y='include', x='include'))
    }
    # Update model
    mod <- update(mod)
    # Blups
    blups<- predict(mod, classify='gen', ignore=c('(Intercept)','blk'))$pvals
    # Prediction error variance
    pev<- blups[,'std.error']^2
    # Genetic variance
    Vg<- summary(mod)$varcomp['gen','component']
    # Reliability for each gen
    rel<- 1-(pev/Vg)
    # Mean reliability
    m_rel<- mean(rel)
    # Result
    res<- data.frame(m_rel=m_rel)
    # Result that will be printed
    return(res)
  }

# Run model
out_blups <- dat1l|> # raw data long format
  # group by env and var
  group_by(env,var) |> 
  # number of blocks within env and var
  mutate(no_blk=length(unique(blk))) |> 
  # remove the ones with only one block
  filter(no_blk>1) |> 
  # remove no_blk column
  dplyr::select(-no_blk) |>
  nest() |>
  # run the model for each combination of env and var using the function
  mutate(blups=map(data,
                   ~mod.blups(.x))) |>
  # unnest reliability
  unnest_wider(blups) |>
  ungroup() |>
  # round reliability
  mutate_if(is.numeric,~round(.,3)) |>
  glimpse()

#### Reliability ----
# Summary
out_blups |>
  dplyr::select(-c(data)) |>
  group_by(var) |>
  summarise(min=min(m_rel),
            max=max(m_rel),
            mean=mean(m_rel),
            median=(median(m_rel))) |>
  mutate_if(is.numeric, ~round(.,2))

# Boxplot
out_blups  |>
  dplyr::select(-c(data)) |>
  group_by(var) |>
  mutate(mrel=format(round(mean(m_rel),2),nsmall=2)) |>
  ggplot(aes(x=var,y=m_rel)) + 
  geom_boxplot(fill='#FF552E',color='#13294B',width=0.5) + 
  geom_point(aes(y=as.numeric(mrel)),shape=4, size=2, color='#13294B') +
  ylab("Reliability") +
  scale_y_continuous(breaks = seq(0,1,by=0.1),limits=c(0,1),expand=c(0.01,0)) + 
  xlab(NULL) + 
  scale_x_discrete(labels=c('Grain yield','Heading time',
                            'Plant height','Test weight')) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust=0.5, vjust=0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = 'none', panel.grid = element_blank())
ggsave("figures/r2.png", width = 6, height = 4, units = "in", dpi=320)

### BLUE model ----
# Estimate genotype blues and weights within each trial

# Model
mod.blues <- # Function name
  function(dat){
    yr <- unique(as.character(dat$year)) # year being analyzed
    
    if('2021' %in% yr){
      # Model for 2021 (including row and column)
      mod <- asreml(fixed = val~1+gen+blk,
                    random = ~row+col,
                    data = dat,
                    na.action = na.method(y='include', x='include'))
    }else{
      # Model for other years (2001 to 2020)
      mod <- asreml(fixed = val~1+gen+blk,
                    data = dat,
                    na.action = na.method(y='include', x='include'))
    }
    # Update model
    mod <- update(mod)
    # Blues
    blues<- predict(mod, classify='gen')$pvals
    # Result
    return(blues)
  }

# Run model and organize output
out_blues <- # create df with blues
  out_blups |> # previous model output
  # unnest draw data used in previous model
  unnest(data) |>
  # group by env and var
  group_by(env,var) |>
  # nest by each combination of env and var
  nest() |>
  # run the model for each combination of env and var using the function
  mutate(blues=map(data,
                   ~mod.blues(.x))) |>
  # remove raw data
  dplyr::select(-data) |>
  # unnest blues
  unnest(blues) |>
  # clean and rename variables
  clean_names() |>
  rename(prd_val=predicted_value) |>
  # remove status column
  dplyr::select(-status) |>
  # add year and location columns by env
  left_join(dat1l |>
              group_by(env) |>
              summarise(year=unique(year),
                        loc=unique(loc)),
            by='env') |>
  glimpse()

#### Test weight imputation
# In 5 trials test weight was evaluated in only one replication. Because of that it was excluded from the single trials models. To be able to add these observations to the genetic trend analyses we used the test weight from the single replication and the average standard error from all trials before 2014. We did not use data from 2014 and later because the measurement method was different and all trials with information from only one replication were prior to 2014.

# Trials with test weight data from only one replication
dat1l|>
  group_by(env,var) |>
  mutate(no_blk=length(unique(blk))) |> 
  filter(no_blk==1) |>
  ungroup() |>
  group_by(year, loc) |>
  summarise() |>
  ungroup() |>
  droplevels() |>
  glimpse()

# Average standard error for test weight standard error from 2001 to 2013 for imputation
tw.std_error <- out_blues |>
  ungroup() |>
  filter(var=='test_weight'&as.numeric(as.character(year))<2014) |>
  summarise(se=sqrt(mean((std_error^2)*3))) |>
  unlist(use.names = F)

# Create a long format data set with blues and standard errors
blues <- # long format df with blues and standard errors
  out_blues |> # blues output
  # add trials with test weight data from only one block and mean standard error
  bind_rows(dat1l|>
              group_by(env,var) |>
              mutate(no_blk=length(unique(blk))) |> 
              filter(no_blk==1) |>
              mutate(prd_val=val,
                     std_error=tw.std_error) |>
              dplyr::select(env,var,gen,prd_val,std_error,year,loc)) |>
  # group by gen to create new columns and organize the data
  group_by(gen) |>
  mutate(gen=as.factor(gen), # change to factor
         # create gidyr (first year the genotype entered the trial)
         gidyr=min(as.numeric(as.character(year))), 
         # first year 2001 we consider IL-98 lines as first entry
         gidyr=ifelse(gidyr!=2001,gidyr,ifelse(grepl('^98',gen),2001,NA)), 
         # add NA to gidyr of gen we will exclude from the analysis
         ## NA if not starts with number (not IL) 
         gidyr=ifelse(!grepl('^\\d',gen),NA,gidyr), 
         #IL gen that starts with US
         gidyr=ifelse(grepl('^US',gen),
                      min(as.numeric(as.character(year))),gidyr), 
         gidyr=ifelse(grepl('B-B',gen),NA,gidyr), # (not IL)
         gidyr=ifelse(grepl('^96',gen),NA,gidyr), # (first entered before 2001)
         gidyr=ifelse(grepl('^97',gen),NA,gidyr), # (first entered before 2001)
         ## gidyr for Kaskaskia = 1993 based on the line name (90-7514)
         gidyr=ifelse(gen=='Kaskaskia',1993,gidyr), 
         gidyr=as.factor(gidyr), # change to factor
         year_n=as.numeric(as.character(year)), # numeric year
         gidyr_n=as.numeric(as.character(gidyr)), # numeric gidyr
         wt=(1/(std_error^2)))|> # weights to include in the models
  ungroup() |> # ungroup
  dplyr::select(env,year,year_n,loc,gen,#check,check_f,
                gidyr,gidyr_n,
                var,prd_val,std_error,wt) |> # select variables
  arrange(year,loc,gen,var) |> # arrange data
  droplevels.data.frame()

## Genetic correlation among traits ----
asreml.options(trace=T) # Report convergence monitoring

### All years ----

#### Data
# Include the NA's to the blues data set for the genetic correlation model
blues.na <- blues |>
  dplyr::select(-c(std_error,wt)) |> # filter out to make it easier to pivot wider
  pivot_wider(values_from = prd_val, names_from = var) |> #pivot wider (add NA)
  pivot_longer(names_to='var', values_to='prd_val',
               cols=c(grain_yield:heading_time), 
               values_drop_na = F) |> # pivot longer without dropping NA's
  left_join(blues) |> # add back std_error and weights
  mutate_if(is.character,as.factor) |>
  arrange(var,env,gen) |>
  glimpse()

#### Model
# Multi-trait model for the genetic correlation among traits
mt_mod <- asreml(fixed = prd_val~var,
                 random = ~us(var):gen + env:var + gen:env:var,
                 family=asr_gaussian(dispersion = 1),
                 residual= ~units:var, 
                 na.action = na.method(y='include', x='include'), 
                 asmv=var, data=blues.na, weights=wt,
                 workspace="8gb")

# Extract unique variables
var <- as.vector(unique(blues.na$var)) 
# Extract variance components and convert to matrix
mat_varcomp <- 
  as.data.frame(summary(mt_mod)$varcomp) |> # Get variance components
  dplyr::select(component) |>  # Select the column containing the variance components
  rownames_to_column('varcomp') |>  # Convert row names to a separate column
  filter(str_detect(varcomp, paste(var, collapse = '|'))) |>  # Filter rows that contain variables
  mutate(varcomp=str_replace_all(varcomp,'var:gen!var_','')) |>  # Replace specific patterns
  separate(varcomp, c('var1','var2'),sep = ':') |>  # Split the 'varcomp' column into two separate columns
  spread(var2,component) |>  # Spread the 'component' values into separate columns
  column_to_rownames('var1') |>  # Use 'var1' column as row names
  as.matrix() |>  # Convert the resulting data frame to a matrix
  t()  # Transpose the matrix

mat_varcomp

# Function to calculate genetic correlations from the matrix of variance components
calc_gen_corr <- function(mat) {
  # Get the number of rows in the matrix
  n <- nrow(mat)  
  # Create a square matrix to store the genetic correlations, initialized with 1's
  genetic_corr_mat <- matrix(1, nrow = n, ncol = n,  
                             dimnames = list(rownames(mat), colnames(mat)))  # Set the row and column names
  
  for (i in 1:(n - 1)) {  # Loop over the rows of the matrix (excluding the last row)
    for (j in (i + 1):n) {  # Loop over the columns of the matrix starting from i+1 to n
      
      genetic_corr <- mat[i, j] / (sqrt(mat[i, i]) * sqrt(mat[j, j]))  # Calculate the genetic correlation
      # Assign the calculated genetic correlation to the corresponding cell in the upper triangle of the matrix
      genetic_corr_mat[i, j] <- genetic_corr
      # Assign the same value to the corresponding cell in the lower triangle of the matrix
      genetic_corr_mat[j, i] <- genetic_corr
    }
  }
  return(genetic_corr_mat)  # Return the matrix containing the genetic correlations
}

# Calculate genetic correlations
genetic_corr <- calc_gen_corr(mat_varcomp)  # Call the function with the input matrix
genetic_corr

#### Plot all ----
colnames(genetic_corr) <- c('Grain yield','Heading time',
                            'Plant height','Test weight')
rownames(genetic_corr) <- c('Grain yield','Heading time',
                            'Plant height','Test weight')

# Define a color palette transitioning from orange to blue
my_palette <- colorRampPalette(c("#13294B", "white","#FF552E"))(n = 100)

# Plot the correlation matrix with the specified color palette
png("figures/gen_corr.png", width = 6, height = 4, units = "in", res = 320)
corrplot(genetic_corr, method = "color", type = "upper", addCoef.col = TRUE, 
         addgrid = FALSE, outline = "black",  
         tl.cex = 1, tl.srt = 0, tl.col = "black", 
         cl.cex = 1, cl.align.text = "l",
         number.cex = 1,
         diag = FALSE, mar = c(0.2, 0.2, 0.2, 0.2),
         col = my_palette)
dev.off()

### First 3yrs ----

#### Data
blues_1st3yrs <- blues.na |>
  filter(year_n <2004) |>
  glimpse()

#### Model 
mt_mod_1st3yr <- asreml(fixed = prd_val~var,
                        random = ~us(var):gen + env:var + gen:env:var,
                        family=asr_gaussian(dispersion = 1),
                        residual= ~units:var, 
                        na.action = na.method(y='include', x='include'), 
                        asmv=var, data=blues_1st3yrs, weights=wt,
                        workspace="2gb")

# Extract unique variables
var <- as.vector(unique(blues.na$var)) 
# Extract variance components and convert to matrix
mat_varcomp_1st3yr <- 
  as.data.frame(summary(mt_mod_1st3yr)$varcomp) |> # Get variance components
  dplyr::select(component) |>  # Select the column containing the variance components
  rownames_to_column('varcomp') |>  # Convert row names to a separate column
  filter(str_detect(varcomp, paste(var, collapse = '|'))) |>  # Filter rows that contain variables
  mutate(varcomp=str_replace_all(varcomp,'var:gen!var_','')) |>  # Replace specific patterns
  separate(varcomp, c('var1','var2'),sep = ':') |>  # Split the 'varcomp' column into two separate columns
  spread(var2,component) |>  # Spread the 'component' values into separate columns
  column_to_rownames('var1') |>  # Use 'var1' column as row names
  as.matrix() |>  # Convert the resulting data frame to a matrix
  t()  # Transpose the matrix

mat_varcomp_1st3yr

# Calculate genetic correlations
genetic_corr_1st3yr <- calc_gen_corr(mat_varcomp_1st3yr)  # Call the function with the input matrix
genetic_corr_1st3yr

#### Plot 1st 3yrs----
png("figures/gen_corr_1st3yr.png", width = 6, height = 4, units = "in", res = 320)
corrplot(genetic_corr_1st3yr, method = "color", type = "upper", addCoef.col = TRUE, 
         addgrid = FALSE, outline = "black",  
         tl.cex = 1, tl.srt = 0, tl.col = "black", 
         cl.cex = 1, cl.align.text = "l",
         number.cex = 1,
         diag = FALSE, mar = c(0.2, 0.2, 0.2, 0.2),
         col = my_palette)
dev.off()

### Last 3yrs

# Data
blues_lst3yrs <- blues.na |>
  filter(year_n >2018) |>
  glimpse()

# Model
mt_mod_lst3yr <- asreml(fixed = prd_val~var,
                        random = ~us(var):gen + env:var + gen:env:var,
                        family=asr_gaussian(dispersion = 1),
                        residual= ~units:var, 
                        na.action = na.method(y='include', x='include'), 
                        asmv=var, data=blues_lst3yrs, weights=wt,
                        workspace="2gb")
# Extract unique variables
var <- as.vector(unique(blues.na$var)) 
# Extract variance components and convert to matrix
mat_varcomp_lst3yr <- 
  as.data.frame(summary(mt_mod_lst3yr)$varcomp) |> # Get variance components
  dplyr::select(component) |>  # Select the column containing the variance components
  rownames_to_column('varcomp') |>  # Convert row names to a separate column
  filter(str_detect(varcomp, paste(var, collapse = '|'))) |>  # Filter rows that contain variables
  mutate(varcomp=str_replace_all(varcomp,'var:gen!var_','')) |>  # Replace specific patterns
  separate(varcomp, c('var1','var2'),sep = ':') |>  # Split the 'varcomp' column into two separate columns
  spread(var2,component) |>  # Spread the 'component' values into separate columns
  column_to_rownames('var1') |>  # Use 'var1' column as row names
  as.matrix() |>  # Convert the resulting data frame to a matrix
  t()  # Transpose the matrix

mat_varcomp_lst3yr

# Calculate genetic correlations
genetic_corr_lst3yr <- calc_gen_corr(mat_varcomp_lst3yr)  # Call the function with the input matrix
genetic_corr_lst3yr

#### Plot lst 3yrs----
png("figures/gen_corr_lst3yr.png", width = 6, height = 4, units = "in", res = 320)
corrplot(genetic_corr_lst3yr, method = "color", type = "upper", addCoef.col = TRUE, 
         addgrid = FALSE, outline = "black",  
         tl.cex = 1, tl.srt = 0, tl.col = "black", 
         cl.cex = 1, cl.align.text = "l",
         number.cex = 1,
         diag = FALSE, mar = c(0.2, 0.2, 0.2, 0.2),
         col = my_palette)
dev.off()

## Trend analysis ----

### Data ----
dat_trend_21yr <- 
  blues |> # blues from single trial analysis
  mutate(check=ifelse(gen=='Kaskaskia',T,F), # Assign check & entry genotypes
         check_f=as.factor(check)) |> # Check as factor
  filter(!is.na(gidyr)|check==T) |> # Remove genotypes without gidyr
  filter(check==T|year_n==gidyr_n) |> # filter checks & entries to keep entry data only for the first year it was tested
  glimpse()

### Models ----
# Toy data used to create and test the function
dat<- dat_trend_21yr|>filter(var=='grain_yield') |>
  glimpse()

# Function to run the genetic, environmental, and phenotypic trends
mod.trends <- 
  function(dat){
    var <- unique(as.character(dat$variable)) # variable being analyzed
    if('heading_time' %in% var) { # models for heading time
      
      # Genetic trend
      mod_gen <- asreml(fixed=prd_val ~gidyr_n + check_f,
                        random=~year + at(year,'2021'):gen + at(year,'2021'):loc,
                        weights=wt, family=asr_gaussian(dispersion=1),
                        na.action=na.method('omit'), data=dat)
      mod_gen <- update(mod_gen)
      ## Point estimates of gidyr
      mod_gen_pt <- asreml(fixed=prd_val ~gidyr + check_f,
                           random=~year+at(year,'2021'):gen+at(year,'2021'):loc,
                           weights=wt, family=asr_gaussian(dispersion=1),
                           na.action=na.method('omit'), data=dat)
      mod_gen_pt <- update(mod_gen_pt)
      
      # Environmental trend - check only
      mod_env <- asreml(fixed=prd_val ~year_n,
                        random=~at(year,'2021'):gen+at(year,'2021'):loc,
                        weights=wt, family=asr_gaussian(dispersion=1),
                        na.action=na.method('omit'),data=dat|>filter(check==T))
      mod_env <- update(mod_env)
      
      # Phenotypic trend - entry only
      mod_pheno <- asreml(fixed=prd_val ~gidyr_n,
                          random=~at(year,'2021'):gen+at(year,'2021'):loc,
                          weights=wt, family=asr_gaussian(dispersion=1),
                          na.action=na.method('omit'),data=dat|>filter(check==F))
      mod_pheno <- update(mod_pheno)
      
    }else{ # models for grain yield, test weight, and plant height
      
      # Genetic trend
      mod_gen <- asreml(fixed=prd_val~ gidyr_n+ check_f + loc + check_f:loc,
                        random=~ year + gen:year + gen:loc,
                        weights=wt, family=asr_gaussian(dispersion=1),
                        na.action=na.method('omit'), data=dat)
      mod_gen <- update(mod_gen)
      ## Point estimates of gidyr
      mod_gen_pt <- asreml(fixed=prd_val~ gidyr+ check_f + loc + check_f:loc,
                           random=~ year + gen:year + gen:loc,
                           weights=wt, family=asr_gaussian(dispersion=1),
                           na.action=na.method('omit'), data=dat)
      mod_gen_pt <- update(mod_gen_pt)
      
      # Environmental trend - check only
      mod_env <- asreml(fixed=prd_val~ year_n, 
                        random=~ loc:year + gen:loc + gen:year,
                        weights=wt, family=asr_gaussian(dispersion=1),
                        na.action=na.method('omit'),
                        data=dat|>filter(check==T))
      mod_env <- update(mod_env)
      
      # Phenotypic trend - entry only
      mod_pheno <- asreml(fixed=prd_val~ gidyr_n + loc, 
                          random=~ year + gen:loc + gen:year,
                          weights=wt, family=asr_gaussian(dispersion=1),
                          na.action=na.method('omit'),
                          data=dat|>filter(check==F))
      mod_pheno <- update(mod_pheno)
    }
    
    # Genetic trend figure data
    ## Predict genetic trend line
    gen_line <- predict(mod_gen,  classify='gidyr_n:check_f',
                        levels=list(check_f="FALSE",
                                    gidyr_n=c(unique(dat$year_n))))$pvals |>
      as.data.frame() |> clean_names() |> 
      rename(pval_line=predicted_value, se_line=std_error) |>
      mutate(gidyr=factor(gidyr_n),
             ci1=pval_line+1.96*se_line, ci2=pval_line-1.96*se_line)
    ## Predict genetic trend points
    gen_pts <- predict(mod_gen_pt, classify='gidyr:check_f',
                       levels=list(check_f="FALSE",
                                   gidyr=factor(c(unique(dat$year_n)))))$pvals |>
      as.data.frame() |>clean_names() |> 
      rename(pval_point=predicted_value,se_point=std_error)
    ## Join trend line and point predictions
    dat_gen_trend <- left_join(gen_line, gen_pts) |>
      nest(dat_gen_trend=everything())
    
    # Wald test of fixed effects
    ## Function to convert asreml wald to df
    wald_df <- function(mod,trend){
      aov_df <- wald(mod) |> as.data.frame() |> rownames_to_column('term') |> 
        mutate(trend=trend) |> relocate(trend)
      return(aov_df)
    }
    ## Bind all trends to the same df
    wald <- bind_rows(wald_df(mod_gen,'gen'), 
                      wald_df(mod_env,'env'), 
                      wald_df(mod_pheno,'pheno'))
    
    # Fixed effects estimates
    ## Function to convert asreml summary to df
    fixeff_df <- function(mod,trend){
      fixeff_df <- summary(mod, coef=TRUE)$coef.fixed |> as.data.frame() |> 
        rownames_to_column('term') |> mutate(trend=trend) |> relocate(trend)
      return(fixeff_df)
    }
    ## Bind all trends to the same df
    coef_fixeff <- bind_rows(fixeff_df(mod_gen,'gen'), 
                             fixeff_df(mod_env,'env'), 
                             fixeff_df(mod_pheno,'pheno'))
    ## Raw fixed effects
    raw_fixeff <- bind_rows(list(wald=wald,coef_fixeff=coef_fixeff),.id='id') |>
      clean_names()|>nest(raw_fixeff=everything())
    ## Summarized fixed effects
    fixeff <- full_join(wald,coef_fixeff,by=c('trend', 'term')) |> 
      clean_names() |> filter(!grepl('Intercept|residual',term)) |> 
      select(-c(sum_of_sq,z_ratio)) |> nest(fixeff=everything())
    
    # Variance components
    ## Function to convert varcomp to df
    varcomp_df <- function(mod,trend){
      varcomp_df <- summary(mod)$varcomp |> as.data.frame() |> 
        rownames_to_column('term') |> mutate(trend=trend) |> relocate(trend)
      return(varcomp_df)
    }
    ## Bind rows
    varcomp <- bind_rows(varcomp_df(mod_gen,'gen'), 
                         varcomp_df(mod_env,'env'), 
                         varcomp_df(mod_pheno,'pheno')) |>
      clean_names() |> nest(varcomp=everything())
    
    return(data.frame(dat_gen_trend=dat_gen_trend,
                      fixeff=fixeff,
                      varcomp=varcomp))
  }

#### Run model
# Run the model for all variables
out_trend <- dat_trend_21yr |>
  mutate(variable=as.character(var)) |>
  group_by(var) |>
  nest() |>
  mutate(trends=map(data,
                    ~mod.trends(.x))) |>
  dplyr::select(-data) |>
  unnest(trends) |>
  ungroup()

### Figures ----
# Data
gt_dat <- out_trend |>
  dplyr::select(var,dat_gen_trend) |>
  unnest() |>
  glimpse()

trends <- out_trend |>
  dplyr::select(var, fixeff) |>
  unnest(fixeff) |>
  filter(term%in%c('year_n','gidyr_n','at(check_f, FALSE):gidyr_n')) |>
  mutate(slope=round(solution,3),
         se=round(std_error,3),
         sign=ifelse(pr_chisq<0.001,'***',
                     ifelse(pr_chisq<0.01,'**',
                            ifelse(pr_chisq<0.05,'*','')))) |>
  glimpse()

#### Genetic trends ----

##### Grain yield ----
slope_gy <- trends |>
  filter(var=="grain_yield"&trend=="gen") |>
  mutate(label=paste(slope, '±', se, 'kg ha<sup>-1</sup> yr<sup>-1</sup>', sign)) |>
  glimpse()

gt_dat |> filter(var=="grain_yield") |>
  ggplot(aes(gidyr_n)) +
  geom_line(aes(y=ci1),color='#13294B', linetype = 2, size=0.5, alpha=0.5) +
  geom_line(aes(y=ci2),color='#13294B', linetype = 2, size=0.5, alpha=0.5) +
  geom_ribbon(aes(ymin = ci2, ymax = ci1), fill = '#13294B', alpha = 0.1) +
  geom_point(aes(y=pval_point), color='#13294B', fill='#13294B', size=2) +
  geom_line(aes(x=gidyr_n,y=pval_line),color='#FF552E', size=1) +
  scale_x_continuous(name= "Year", breaks = seq(2001,2021, by = 1)) +
  scale_y_continuous(name = expression("Grain yield (kg ha"^-1*")"),
                     limits = c(4200,5800), expand = c(0,0), breaks = seq(3000, 7000, by = 250)) +
  geom_richtext(data=slope_gy,aes(x=2001,y=5625,label=label), hjust=0, label.color = NA, size =4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = 'none', panel.grid = element_blank())
ggsave("figures/grain_yield.png", width = 6, height = 4, units = "in", dpi=320)

##### Test weight ----
slope_tw <- trends |>
  filter(var=="test_weight"&trend=="gen") |>
  mutate(label=paste(slope, '±', se,'g L<sup>-1</sup> yr<sup>-1</sup>', sign)) |>
  glimpse()

gt_dat |> filter(var=="test_weight") |>
  ggplot(aes(gidyr_n)) +
  geom_line(aes(y=ci1),color='#13294B', linetype = 2, size=0.5, alpha=0.5) +
  geom_line(aes(y=ci2),color='#13294B', linetype = 2, size=0.5, alpha=0.5) +
  geom_ribbon(aes(ymin = ci2, ymax = ci1), fill = '#13294B', alpha = 0.1) +
  geom_point(aes(y=pval_point), color='#13294B', fill='#13294B', size=2) +
  geom_line(aes(x=gidyr_n,y=pval_line),color='#FF552E', size=1) +
  scale_x_continuous(name= "Year", breaks = seq(2001,2021, by = 1)) +
  scale_y_continuous(name = expression("Test weight (g L"^-1*")"),
                     limits = c(718,768), expand = c(0,0), breaks = seq(720, 770, by = 5)) +
  geom_richtext(data=slope_tw,aes(x=2001,y=762.5,label=label), hjust=0, label.color = NA, size =4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = 'none', panel.grid = element_blank())
ggsave("figures/test_weight.png", width = 6, height = 4, units = "in", dpi=320)

##### Heading time ----
slope_ht <- trends |>
  filter(var=="heading_time"&trend=="gen") |>
  mutate(label=paste(slope, '±', se,'d yr<sup>-1</sup>', sign)) |>
  glimpse()

gt_dat |> filter(var=="heading_time") |>
  ggplot(aes(gidyr_n)) +
  geom_line(aes(y=ci1),color='#13294B', linetype = 2, size=0.5, alpha=0.5) +
  geom_line(aes(y=ci2),color='#13294B', linetype = 2, size=0.5, alpha=0.5) +
  geom_ribbon(aes(ymin = ci2, ymax = ci1), fill = '#13294B', alpha = 0.1) +
  geom_point(aes(y=pval_point), color='#13294B', fill='#13294B', size=2) +
  geom_line(aes(x=gidyr_n,y=pval_line),color='#FF552E', size=1) +
  scale_x_continuous(name= "Year", breaks = seq(2001,2021, by = 1)) +
  scale_y_continuous(name = expression("Heading time (JD)"),
                     limits = c(129,139), expand = c(0,0), breaks = seq(130, 138, by = 1)) +
  geom_richtext(data=slope_ht,aes(x=2001,y=138,label=label), hjust=0, label.color = NA, size =4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = 'none', panel.grid = element_blank())
ggsave("figures/heading_time.png", width = 6, height = 4, units = "in", dpi=320)

### Table ----
trends |>
  mutate(value=paste(slope, '±', se, sign)) |>
  dplyr::select(trend, var, value) |>
  pivot_wider(names_from = trend, values_from = value) |>
  mutate(var=c( 'Grain yield', 'Plant height', 'Test weight', 'Heading time')) |> 
  write_csv("tables/table2.csv")

# Save ----

save.image("data/wheat_genetic_trends_R.RData")
#load("data/wheat_genetic_trends_R.RData")