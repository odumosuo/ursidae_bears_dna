####Introductions####
#The topic of interest for this paper is the Ursidae. The Ursidae family is more commonly referred to as bears (Darin, 2015). The higher taxonomic classifications are Animalia, Chordata, Mammalia and Carnivora for the kingdom, phylum, class and order respectively (Darin, 2015). Of the 13 species that have been recorded, only eight of them remain extant. Though most acknowledge the rarity of the giant panda, with an estimated population around 1000 (Swaisgood et al., 2016), most are unaware that several of the other species within this family are just as rare. The International Union for Conservation of Nature (IUCN) lists seven species as vulnerable: Giant Panda (Ailuropoda melanoleuca), Spectacled Bear (Tremarctos ornatus), Sun bear (Helarctos malayanus), Polar bear (Ursus maritimus), Sloth bear (Melursus ursinus), and Asiatic Black Bear (Ursus thibetanus)( IUCN red list of threatened species., n.d.). The IUCN states a vulnerable species as “A species considered to be facing a high risk of extinction in the wild” (IUCN red list of threatened species., n.d.). Because of the rarity of species in this family, much of the data exploration and inquiries in this script are based upon geographical distributions of Ursidae as well as the species richness.
#For data exploration, a subset of variables that might yield interesting information was created. Outliers were checked for errors in data entry. The expected values of latitude should be between -90 and 90, whereas, the expected values of longitude should be between -180 and 180. Outliers in the elevation and depth data were also examined. The Ursidae are usually found in the habitats of deserts, forests, grasslands, ice floes, and mountains (Darin, 2015). All elevation data should be above zero and depth data should be virtually non-existent. Additional investigations included checking the number of country data and the number of coordinate data, exploring the distribution of the Ursidae family, species and BINs amongst the countries and exploring the distribution of species and BINs amongst each other 
#For this paper it’s hypothesized that the distribution of data by country is different for the Ursus (genus with largest specimen record in the BOLD data) and for the Ailuropoda (genus with the second largest specimen record in the BOLD data). That is, there is a relationship between the genus variable and the country variable.  This hypothesis is predicated on the fact that the Ailuropoda genus includes a single species, the Giant Panda (Ailuropoda melanoleuca) usually found China (Swaisgood et al., 2016) On the other hand, the Ursus genus includes four extant species found in a variety of continents (Darin, 2015). A chi-square test was used to analyze this claim. The paper inspects the accuracy of the genus locations as expected from literature and the adequacy of the sampling. 

####Setting up####
### load package(S)
library(tidyverse)
library(vegan)

###Get present Working Directory. Previously set working directory with setwd(). Excluded runnable function to prevent error.
getwd()


###Get bold file from online
ursidae_bold <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Ursidae&format=tsv")


###save bold file to current working directory
write_tsv(ursidae_bold, "ursidae_bold_data.txt")
View(ursidae_bold)


###Explore variables in the dataframe 
names(ursidae_bold)


###Create a subset of interested variables from ursidae_bold to form into new tibble called ursidae_bold2 
ursidae_bold2 <- ursidae_bold[ , c("processid", "institution_storing", "bin_uri", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name", "identification_method", "tissue_type", "sampling_protocol", "lifestage", "sex", "reproduction", "habitat", "extrainfo", "notes", "lat", "lon", "elev", "depth", "country", "province_state", "region", "sector", "exactsite", "markercode", "run_dates")]
names(ursidae_bold2)
View(ursidae_bold2)





####Made Functions####
#create afunction to make a contingency table of the distributions of two categorical variables. x = first variable position in ursidae_bold2 .y = second variable position in ursidae_bold2. Exclude na's
make_table_two_variables <- function(x, y) 
  {
  ursidae_bold2[ , c(x, y)] %>%
    filter(!is.na(x)) %>%
    filter(!is.na(y)) %>%
    table()
  }





####Data Exploration####
###get summary and classes of data
class(ursidae_bold2)
summary(ursidae_bold2)


###check which markercodes were used?
unique(na.omit(ursidae_bold2$markercode))


###check genus and number of genus. Exclude na's
unique(na.omit(ursidae_bold2$genus_name))
length(na.omit(unique(ursidae_bold2$genus_name)))
###check species and number of species. Exclude na's
unique(na.omit(ursidae_bold2$species_name))
length(na.omit(unique(ursidae_bold2$species_name)))
###check BINs number of BINs. Exclude na's
unique(na.omit(ursidae_bold2$bin_uri))
length(na.omit(unique(ursidae_bold2$bin_uri)))
###check number of records with coordinate data vs number of records with country data.Exclude na's  
length(na.omit(ursidae_bold2$lon)) 
length(na.omit(ursidae_bold2$country))


###Check for outliers 
#latitude needs to be -90 <= x <= 90.longitude needs to be -180 <= to <= 180. Get summary of data
summary(ursidae_bold2$lat)
summary(ursidae_bold2$lon)
#elevation needs to be above 0.Exclude na's
sum(ursidae_bold2$lat < 0, na.rm = TRUE)
#all data for depth should be 0 or na.Exclude na's
sum(ursidae_bold2$lat < 0, na.rm = TRUE)



  

####Figure out how data is distributed among countries. 
countries <- table(ursidae_bold2$country)
countries
class(countries)
#Check how many countries are represented.
length(countries)



### map of the genus locations on world map
#make data frame of longitude and latitude from Ursidae_bold2. Exclude records with na values for either lon or la.
names(ursidae_bold2)
genus_map_table <- ursidae_bold2[ , c("lon","lat", "genus_name")] %>%
  filter(!is.na(lon)) %>%
  filter(!is.na(lat)) %>%
  filter(!is.na(genus_name))
#make world map. Get world data, Map world data, then map Ursidae_bold2 data. 
world <- map_data("world") 
ggplot() +
  geom_map(data = world, map = world, aes(x = long, y = lat , map_id = region), fill = "#CCCCCC", colour = "black", size = 0.3) + 
  geom_point(data = genus_map_table, aes(x = lon, y = lat, color = genus_name, size = 0.7), alpha = 0.8) +
  labs(title ="World map of Ursidae distribution", x = "Longitude", y = "Latitude", colour = "Genus")


###Explore species and BIN distribution among themselves  
#create a table that lists unique BINs for each species.0 Frequency values were excluded
names(ursidae_bold2)
bin_vs_species_contingencytable <- make_table_two_variables(3, 9)
df_bin_vs_species_contingencytable <- data.frame(bin_vs_species_contingencytable) %>%
  filter(Freq != 0)
#create a plot of BINs and species. With species on x axis to see occurrences of BINs within each species 
ggplot(df_bin_vs_species_contingencytable, aes(x = species_name, y = Freq , fill = bin_uri)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title ="Stacked bar graph of BINs and species with species on x axis", x = "Species name", y = "Frequency", fill = "BIN URI")
  
#create a plot of BINs and species. with BINs on x axis to see occurrences of species within each BIN 
ggplot(df_bin_vs_species_contingencytable, aes(x = bin_uri, y = Freq , fill = species_name)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title ="Stacked bar graph of BINs and species with BINs on x axis", x = "BIN URI", y = "Frequency", fill = "Species name")




###Figure out how species and BINs are distributed among countries
#make a table and dataframe of the distributions of species among countries. Used made function make_table_two_variables Exclude 0 frequency of occurence
species_vs_country_contingencytable <- make_table_two_variables(9,23)
df_species_vs_country_contingencytable <- data.frame(species_vs_country_contingencytable) %>%
  filter(Freq != 0)
#make a table and dataframe of the distributions of BINs among countries. Used made function make_table_two_variables. Exclude 0 frequency of occurence 
bins_vs_country_contingencytable <- make_table_two_variables(3,23)
df_bins_vs_country_contingencytable <- data.frame(bins_vs_country_contingencytable) %>%
  filter(Freq != 0)
#create a plot of species vs countries to visualize the countries where the species can be found
names(df_species_vs_country_contingencytable)
ggplot(df_species_vs_country_contingencytable, aes(x = species_name, y = Freq , fill = country)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  coord_flip() +
  labs(title ="Stacked bar graph of country and species", x = "Species name", y = "Frequency", fill = "Country")
  
#create a plot of BINs vs countries to visualize the countries where the BINs can be found
names(df_bins_vs_country_contingencytable)
ggplot(df_bins_vs_country_contingencytable, aes(x = bin_uri, y = Freq , fill = country)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  coord_flip() +
  labs(title ="Stacked bar graph of country and BINs", x = "BIN URI", y = "Frequency", fill = "Country")
  


  
  
  
  
  
  
  






####Analysis to Address Questions####
###Chi square test of the distribution of countries in Ursus vs Ailuropoda genus
#create contingency table with Ursus and Ailuropoda along with the countries the they are found in. Exclude records with na values for either country or genus_name. 
ursus_and_ailuropoda_vs_country_contingencytable <- ursidae_bold2[ , c("country", "genus_name")] %>%
  filter(!is.na(country)) %>%
  filter(!is.na(genus_name)) %>%
  filter(genus_name %in% c("Ursus", "Ailuropoda")) %>%
  table()
#perform chi square analysis
chisq.test(ursus_and_ailuropoda_vs_country_contingencytable)


  
###How are genus spread among countries? Create a table of countries with associated genus. Create a figure with countries vs associated genus to visualize which countries each genus was found in. Used made function make_table_two_variables
#make a table and dataframe of the distributions of genus among countries. Used made function make_table_two_variables. Exclude 0 frequency of occurence
names(ursidae_bold2)
genus_vs_country_contingencytable <- make_table_two_variables(8,23)
df_genus_vs_country_contingencytable <- data.frame(genus_vs_country_contingencytable) %>%
  filter(Freq != 0)
#create a plot of genus vs countries to visualize the countries where the genus can be found
ggplot(df_genus_vs_country_contingencytable, aes(x = genus_name, y = Freq , fill = country)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title ="Stacked bar graph of country and genus with genus on x axis", x = "Genus", y = "Frequency", fill = "Country")



###How well sampled is my group? Create a rarefraction plot to view curve
#create a data frame with the counts of species.Exclude na's
df_countby_species <- ursidae_bold2 %>%
  filter(!is.na(species_name)) %>%
  group_by(species_name) %>%
  count(species_name)
#reshape data to increase columns and decrease rows
df_countby_species_spread <- pivot_wider(data = df_countby_species, names_from  = species_name, values_from = n)
#plot rarefraction curve
rarecurve(df_countby_species_spread, main = "Rarefraction curve of species",  xlab = "Number of Individuals", ylab = "Species Richness") 



###What is the  species richness? 
diversity_index_for_species <- estimateR(df_countby_species_spread)




  

####Plots to submit####
### map of the genus locations on world map
#make data frame of lon and lat. Exclude records with na values for either lon or lan.
#subset lon, lat and genus name data to make world map. Exclude NA's
names(ursidae_bold2)
genus_map_table <- ursidae_bold2[ , c("lon","lat", "genus_name")] %>%
  filter(!is.na(lon)) %>%
  filter(!is.na(lat)) %>%
  filter(!is.na(genus_name))
#make world map of Ursidae distribution to view distribution in the world.Get world data, Map world data, then map Ursidae_bold2 data.
world <- map_data("world") 
ggplot() +
  geom_map(data = world, map = world, aes(x = long, y = lat , map_id = region), fill = "#CCCCCC", colour = "black", size = 0.3) + 
  geom_point(data = genus_map_table, aes(x = lon, y = lat, color = genus_name, size = 0.7), alpha = 0.8) +
  labs(title ="World map of Ursidae distribution", x = "Longitude", y = "Latitude", colour = "Genus")


###A figure of ursus and ailuropoda vs country frequency.(x axis countries, why axis frquency and 2 symbols for ursus and ailuropoda). Figure to view distribution asked in chi-square analysis.  
#Make ursus, ailuropoda contingency table as data frame. Exclude 0 frequency of occurence
df_ursus_and_ailuropoda_vs_country_contingencytable <- data.frame(ursus_and_ailuropoda_vs_country_contingencytable) %>%
  filter(Freq != 0)
#plot ursus and ailuropoda against countries
ggplot(df_ursus_and_ailuropoda_vs_country_contingencytable, aes(x = country , y = Freq , fill = genus_name)) +
  geom_bar(position = "stack", stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title ="Bar graph showing distribution of Ursus and Ailuropoda among countries", x = "country", y = "Frequency", fill = "Genus") +
  coord_flip()

 
 
###Figure showing species richness rarefraction curve.
rarecurve(df_countby_species_spread, main = "Rarefraction curve of species",  xlab = "Number of Individuals", ylab = "Species Richness") 


  
  




####Results and Discussion####
#Ff the seven genus, 13 species, and 16 BINs, there were 38 records with coordinate data and 415 records with information about the country. The coordinate data is noticeably underepresented compared to the country data. This could be accounted for by the fact that a lot of the species in the family are at high risk of extinction (IUCN red list of threatened species., n.d.); researchers are careful to exclude coordinate data to avoid poaching or invasion. Russia had the most amount of bears (90), whereas, Columbia, Greenland, and Japan had the least amount of bears (one each). The species count is less than the BIN count. Before exploration it was believed that this was a result of the presence of cryptic species not identified with the traditional taxonomic grouping but identifiable with the BIN algorithm.  That is, the expectation was that unique species would contain multiple BINs.  This was found. However, with further data exploration, it was also found that multiple species were found in unique BINs. With the incorporation of each successive data point, the reliability of the BIN algorithm improves. As such, increasing the volume of sampling would inevitably lead to a more powerful tool for investigating biodiversity and to understanding this dissonance. From the bar graphs (included in data exploration and not in plots), it is clear that the species Ursus arctos was found in the most countries, as well as the BIN BOLD:AAC3976. A caveat to this is that both Ursus arctos and BIN BOLD:AAC3976 were the most heavily sampled. Unsurprisingly, Ursus arctos is the predominant species in the BIN BOLD:AAC3976 and vis-vera. From figure 1 We can see that only two genus (Tremarctos and Ursus) are represented on the map. This is largely due to the lack of coordinate data for reasons previously explained.
#The chi-square test showed that there was a significant relationship between the genus variable (Ursus and Ailuropoda) the country variable, X2 (22, N = 400) = 297.55, p < 2.2e-16. A graph of the distribution of the two genus among the countries can be seen in figure 2. To explore whether the different genus were found in their expected locations, a bar plot with the genus on x axis and their respective countries stacked on each genus was created (in data exploration not plots). Most genus were found in their respective habitat. However, the extinct genus of Arctotherium had no country information and the giant panda (Ailuropoda melanoleuca) was found in a non native country, the United States. The giant panda is native to China (Swaisgood et al., 2016). It is likely the specimen sampled was being kept in captivity or kept for conservation purposes in the United States. Sampling completeness and species richness was assessed with a rarefaction curve. As seen in figure 3. the curve is reaching an asymptotic value around 13. Complemented with the Chao index and a value of 13.000, we come to the conclusion that our data set is adequately sampled with a species richness of 13. For future studies, it would be interesting to investigate the area with the most density of the Ursidae family where conservation can be concentrated.




####Acknowledgments####
#I would like to acknowledge the following people for their contribution to the completion of this project. These people might not have helped me directly but at some point, we engaged in discussion pertaining to this project: Jessica Labarcena, Nishita Sharif, Daniel Amoako, Amjad Osman, Jesse, and Omar Khan. 





####References####
#Chi-Square Test in R | Explore the Examples and Essential concepts! (n.d.). Data fair training. Retrieved October 3, 2022, from https://data-flair.training/blogs/chi-square-test-in-r/
  
#Creating Functions. (n.d). Swcarpentry. Retrieved October 3, 2022, from https://swcarpentry.github.io/r-novice-inflammation/02-func-R/
  
#Darin M. C. (2015). Chapter 50 – Ursidae. In R. E. Miller & M. E. Fowler (Eds.), Fowler's Zoo and Wild Animal Medicine, Volume 8 (pp. 498-508) W.B. Saunders. https://doi.org/10.1016/B978-1-4557-7397-8.00050-5 

#Filter multiple values on a string column in dplyr. (2014). Stackoverflow. Retrieved October 3, 2022, from https://stackoverflow.com/questions/25647470/filter-multiple-values-on-a-string-column-in-dplyr

#Grouped, stacked and percent stacked barplot in ggplot2 (2018). The R graph gallery. Retrieved October 3, 2022, from https://r-graph-gallery.com/48-grouped-barplot-with-ggplot2

#How To Make World Map with ggplot2 in R? (2020). Datavizpyr. Retrieved October 3, 2022, from https://datavizpyr.com/how-to-make-world-map-with-ggplot2-in-r/#Simple_World_Map_in_R

#IUCN red list of threatened species. (n.d.). IUCN redlist. Retrieved October 3, 2022, from https://www.iucnredlist.org/search?taxonomies=101398&searchType=species 

#IUCN red list of threatened species. (n.d.). IUCN redlist. Retrieved October 3, 2022, from https://www.iucnredlist.org/

#Rotating and spacing axis labels in ggplot2 (2020). Stackoverflow. Retrieved October 3, 2022, from https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2

#Swaisgood, R., Wang, D., & Wei, F. (2016). Ailuropoda melanoleuca. The IUCN Red List of Threatened Species, 2016, e-T712A102080907.

